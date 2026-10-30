"""Run a user's shell command template on a watched file, safely.

Templates such as ``samtools index {file}`` are the operator's own shell
commands, so they run through the shell on purpose. The file name comes from
the watched folder, though, so it must never become shell code. It is never
pasted into the command text: ``{file}`` becomes a reference to the
``ROBIN_FILE`` environment variable, and the shell never runs what a variable
holds.
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
    """
    ref = "${" + FILE_ENV_VAR + "}"
    out = []
    quote = None  # None, "'" or '"'
    i = 0
    while i < len(template):
        if template.startswith(_PLACEHOLDER, i):
            if quote is None:
                out.append(f'"{ref}"')
            elif quote == '"':
                out.append(ref)
            else:
                out.append(f"'\"{ref}\"'")
            i += len(_PLACEHOLDER)
            continue
        char = template[i]
        out.append(char)
        if char == "\\" and quote != "'" and i + 1 < len(template):
            out.append(template[i + 1])  # escaped character: never opens or closes a quote
            i += 2
            continue
        if quote is None and char in "'\"":
            quote = char
        elif char == quote:
            quote = None
        i += 1
    return "".join(out)


def run_file_command(template: str, file_path: str, **kwargs) -> subprocess.CompletedProcess:
    """Run ``template`` through the shell with ``{file}`` standing for ``file_path``.

    Extra keyword arguments go to ``subprocess.run``.
    """
    env = {**os.environ, FILE_ENV_VAR: file_path}
    return subprocess.run(shell_command_for_file(template), shell=True, env=env, **kwargs)
