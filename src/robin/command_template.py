"""Run a user's shell command template on a watched file, safely.

Templates such as ``samtools index {file}`` are the operator's own shell
commands, so they run through the shell on purpose. The file name comes from
the watched folder, though, so it must never become shell code. It is never
pasted into the command text: ``{file}`` becomes a reference to the
``ROBIN_FILE`` environment variable, and the shell never runs what a variable
holds.

Limit: this protects the shell Robin starts. A template that hands text to a
second interpreter (``sh -c "... {file}"``, ``bash -c``, ``eval``,
``ssh host "... {file}"``, ``xargs sh -c``) gives it the expanded file name,
and that interpreter parses it again. In such templates use
``"$ROBIN_FILE"`` inside the inner command instead of ``{file}``; it is
exported to the command's environment.
"""

import os
import subprocess

FILE_ENV_VAR = "ROBIN_FILE"
_PLACEHOLDER = "{file}"


def shell_command_for_file(template: str) -> str:
    """Return the template with each ``{file}`` turned into a ``$ROBIN_FILE`` reference.

    The reference is quoted to suit where the placeholder sits, so the path
    stays one word: ``"${ROBIN_FILE}"`` outside quotes, ``${ROBIN_FILE}``
    inside double quotes, and ``'"${ROBIN_FILE}"'`` inside single quotes
    (close the quote, insert, reopen).

    "Where it sits" is the innermost context. ``$(...)``, a backtick command
    and a ``(...)`` group each start a new unquoted context, even inside
    double quotes, so in ``"$(dirname "{file}")"`` the inner quotes are the
    ones that count. A single quote-state misread that and left the path
    unquoted, so it split at spaces. Not followed: ``case`` patterns inside
    ``$(...)``, here-documents and comments. None of this affects safety:
    the file name is never shell text, whatever the quoting.
    """
    ref = "${" + FILE_ENV_VAR + "}"
    out = []
    # Innermost context last: "top", "$(", "(" and "`" are unquoted;
    # "'" and '"' are quoted.
    stack = ["top"]
    i = 0
    while i < len(template):
        ctx = stack[-1]
        if template.startswith(_PLACEHOLDER, i):
            if ctx == '"':
                out.append(ref)
            elif ctx == "'":
                out.append(f"'\"{ref}\"'")
            else:
                out.append(f'"{ref}"')
            i += len(_PLACEHOLDER)
            continue
        char = template[i]
        out.append(char)
        if ctx == "'":
            if char == "'":
                stack.pop()
            i += 1
            continue
        if char == "\\" and i + 1 < len(template):
            out.append(template[i + 1])  # escaped character: never opens or closes anything
            i += 2
            continue
        if template.startswith("$(", i):
            out.append("(")
            stack.append("$(")
            i += 2
            continue
        if char == "`":
            if ctx == "`":
                stack.pop()
            else:
                stack.append("`")
        elif ctx == '"':
            if char == '"':
                stack.pop()
        elif char in "'\"":
            stack.append(char)
        elif char == "(":
            stack.append("(")
        elif char == ")" and ctx in ("$(", "("):
            stack.pop()
        i += 1
    return "".join(out)


def run_file_command(template: str, file_path: str, **kwargs) -> subprocess.CompletedProcess:
    """Run ``template`` through the shell with ``{file}`` standing for ``file_path``.

    Extra keyword arguments go to ``subprocess.run``.
    """
    env = {**os.environ, FILE_ENV_VAR: file_path}
    return subprocess.run(shell_command_for_file(template), shell=True, env=env, **kwargs)
