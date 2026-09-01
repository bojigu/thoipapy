"""Tests for utils.Command, which runs every external tool the pipeline depends on.

Command used to run a single string through /bin/sh. The paths interpolated into that string are
built from the protein name, and the protein name is whatever the user put in their FASTA header,
so every character of it reached the shell. These tests pin the fix: arguments are passed to the
program, not to a shell, and nothing in them is interpreted.
"""

import sys

import pytest
from thoipapy.utils import Command


def test_a_string_command_is_refused():
    """The old calling convention must fail loudly rather than becoming a missing-file error."""
    with pytest.raises(TypeError, match="list of arguments"):
        Command("echo hello")


def test_shell_metacharacters_in_an_argument_are_not_interpreted(tmp_path):
    """A path carrying a payload reaches the program as one argument, and nothing runs.

    This is the live failure: a user submitted a UniProt header containing pipes, thoipapy built
    "cd-hit -i /path/sp|Q8WWF3|... -o ..." and the shell split it. The same mechanism runs
    anything else.
    """
    marker = tmp_path / "EXECUTED"
    hostile = f"x; touch {marker}; echo `touch {marker}` $(touch {marker}) | touch {marker}"
    written = tmp_path / "argv.txt"

    command = Command(
        [sys.executable, "-c", "import sys; open(sys.argv[2],'w').write(sys.argv[1])", hostile, str(written)]
    )
    command.run(timeout=60)

    assert command.succeeded(), command.stderr
    assert not marker.exists(), "the injected command ran"
    assert written.read_text() == hostile, "the argument did not arrive intact as one argv element"


def test_stdout_path_writes_the_programs_output_to_a_file(tmp_path):
    """Replaces the shell "> file" redirects."""
    out = tmp_path / "captured.txt"
    command = Command([sys.executable, "-c", "print('hello')"], stdout_path=out)
    command.run(timeout=60)

    assert command.succeeded(), command.stderr
    assert out.read_text().strip() == "hello"


def test_stdin_bytes_reach_the_program(tmp_path):
    """Replaces the shell pipelines."""
    out = tmp_path / "echoed.txt"
    command = Command(
        [sys.executable, "-c", "import sys; sys.stdout.write(sys.stdin.read())"],
        stdin_bytes=b"ACDEF\nGHIKL\n",
        stdout_path=out,
    )
    command.run(timeout=60)

    assert command.succeeded(), command.stderr
    assert out.read_bytes() == b"ACDEF\nGHIKL\n"


def test_cwd_is_where_the_program_runs(tmp_path):
    """rate4site writes r4s.res into the working directory whatever -o says."""
    workdir = tmp_path / "work"
    workdir.mkdir()
    command = Command([sys.executable, "-c", "open('written_here.txt','w').write('x')"], cwd=workdir)
    command.run(timeout=60)

    assert command.succeeded(), command.stderr
    assert (workdir / "written_here.txt").is_file()


def test_a_nonzero_exit_is_not_success():
    command = Command([sys.executable, "-c", "raise SystemExit(3)"])
    command.run(timeout=60)

    assert command.returncode == 3
    assert not command.succeeded()


def test_a_missing_program_is_reported_as_a_failure_not_a_traceback():
    """Without a shell there is no 127 to inherit, so Command has to supply the failure itself.

    A tool absent from the container image is the most common way this fails in production, and it
    has to arrive at the caller's succeeded() check like any other failure.
    """
    command = Command(["thoipapy-no-such-program", "-x"])
    command.run(timeout=60)

    assert not command.succeeded()
    assert command.returncode == 127
    assert "FileNotFoundError" in command.stderr


def test_a_timeout_is_not_success():
    command = Command([sys.executable, "-c", "import time; time.sleep(30)"])
    command.run(timeout=1)

    assert command.timed_out
    assert not command.succeeded()


def test_command_string_is_copy_pasteable():
    command = Command(["cd-hit", "-i", "/tmp/a b|c.fas"])
    assert command.command_string == "cd-hit -i '/tmp/a b|c.fas'"
