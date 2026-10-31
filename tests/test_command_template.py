"""SEC-9: a watched file's name reaches the user's command as data, never as code."""

import pytest

from robin.command_template import run_file_command

PLAIN = "/d/a b.bam"
# Every character the shell treats specially, plus a command that leaves a trace.
HOSTILE = "x$(touch MARKER);y'z\"w`touch MARKER2`\\v*.bam"


def _printed_args(template, file_path, cwd):
    """Run the template with printf so each argument the shell built prints on its own line."""
    result = run_file_command(template, file_path, capture_output=True, text=True, cwd=cwd)
    assert result.returncode == 0, result.stderr
    return result.stdout.split("\n")[:-1]


@pytest.mark.parametrize(
    "template, expected",
    [
        ("printf '%s\\n' {file}", [PLAIN]),
        ("printf '%s\\n' \"{file}\"", [PLAIN]),
        ("printf '%s\\n' '{file}'", [PLAIN]),
        ("printf '%s\\n' {file}.bai", [PLAIN + ".bai"]),
        ("printf '%s\\n' \"{file}.bai\"", [PLAIN + ".bai"]),
        ("printf '%s\\n' 'in={file}'", ["in=" + PLAIN]),
        ("printf '%s\\n' {file} \"{file}\"", [PLAIN, PLAIN]),
        ("printf '%s\\n' \\\"{file}\\\"", ['"' + PLAIN + '"']),
        # Guards for the context tracking (fix-series review, finding 5)
        ("printf '%s\\n' \"$((1+2))\" \"{file}\"", ["3", PLAIN]),
        ("( printf '%s\\n' \"{file}\" )", [PLAIN]),
        ("d=`dirname \"{file}\"`; printf '%s\\n' \"$d\"", ["/d"]),
    ],
)
def test_placeholder_becomes_the_exact_path_in_every_quoting(tmp_path, template, expected):
    assert _printed_args(template, PLAIN, tmp_path) == expected


@pytest.mark.parametrize(
    "template",
    ["printf '%s\\n' {file}", "printf '%s\\n' \"{file}\"", "printf '%s\\n' 'pre {file} post'"],
)
def test_filename_with_shell_syntax_does_not_run(tmp_path, template):
    path = str(tmp_path / HOSTILE)
    out = _printed_args(template, path, tmp_path)
    assert path in out[0]
    assert not (tmp_path / "MARKER").exists()
    assert not (tmp_path / "MARKER2").exists()


def test_filename_does_not_run_inside_command_substitution(tmp_path):
    # Inside "$(...)" the shell starts a new quoting context. Nothing in the
    # name runs there, and {file} is quoted for that context, so it stays whole.
    path = str(tmp_path / HOSTILE)
    _printed_args("printf '%s\\n' \"$(printf '%s ' {file})\"", path, tmp_path)
    assert not (tmp_path / "MARKER").exists()
    assert not (tmp_path / "MARKER2").exists()


def test_template_without_placeholder_runs_unchanged(tmp_path):
    assert _printed_args("printf '%s\\n' {other} 'a b'", PLAIN, tmp_path) == ["{other}", "a b"]


# Inside "$(...)" the shell starts a new quoting context, so a quoted
# placeholder there must stay one word (fix-series review, finding 5).


def test_a_quoted_placeholder_inside_command_substitution_stays_one_word(tmp_path):
    path = str(tmp_path / "a b" / "s.bam")
    template = "printf '%s\\n' \"$(dirname \"{file}\")/output.bai\""
    assert _printed_args(template, path, tmp_path) == [str(tmp_path / "a b" / "output.bai")]


def test_nested_command_substitutions_keep_the_path_whole(tmp_path):
    # Three levels: with two, an odd count of '"' before {file} made the old
    # one-level tracker land on the right answer by luck.
    path = str(tmp_path / "a b" / "c d" / "s.bam")
    template = "printf '%s\\n' \"$(basename \"$(dirname \"$(dirname \"{file}\")\")\")\""
    assert _printed_args(template, path, tmp_path) == ["a b"]


def test_a_shell_variable_named_file_is_left_alone(tmp_path):
    # "${file}" is the shell's own parameter expansion, not the placeholder.
    # It used to become "$"${ROBIN_FILE}"", i.e. the PID followed by text.
    template = "file=x; printf '%s\\n' \"${file}\" ${file} {file}"
    assert _printed_args(template, PLAIN, tmp_path) == ["x", "x", PLAIN]


def test_a_hostile_name_inside_command_substitution_comes_back_exact(tmp_path):
    path = str(tmp_path / HOSTILE)
    template = "printf '%s\\n' \"$(basename \"{file}\")\""
    assert _printed_args(template, path, tmp_path) == [HOSTILE]
    assert not (tmp_path / "MARKER").exists()
    assert not (tmp_path / "MARKER2").exists()


# The two places Robin runs a user's command on a watched file.


def test_watcher_command_does_not_run_the_filename(tmp_path, monkeypatch, capsys):
    from robin.watcher import FileChangeHandler

    monkeypatch.chdir(tmp_path)
    handler = FileChangeHandler(["*"], [], "printf '%s\\n' {file}")
    path = str(tmp_path / HOSTILE)
    handler._run_command(path)
    assert capsys.readouterr().out.strip() == path
    assert not (tmp_path / "MARKER").exists()
    assert not (tmp_path / "MARKER2").exists()


def test_workflow_command_handler_does_not_run_the_filename(tmp_path, monkeypatch):
    from robin.workflow_simple import Job, WorkflowContext, command_handler

    monkeypatch.chdir(tmp_path)
    path = str(tmp_path / HOSTILE)
    job = Job(job_id=1, job_type="index", context=WorkflowContext(filepath=path),
              origin="fast", workflow=["index"])
    command_handler(job, "printf '%s\\n' {file}")
    result = job.context.results["index"]
    assert result["returncode"] == 0
    assert result["stdout"].strip() == path
    assert not (tmp_path / "MARKER").exists()
    assert not (tmp_path / "MARKER2").exists()
