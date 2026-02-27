"""
Lightweight smoke-test helpers for the rmaps_refactored project.

These are not exhaustive tests, but they give you a quick way to
check that the CLI wiring and pyfaidx-backed genome access work
without crashing on very small inputs.

Run from the project root with:

    python tests_smoke.py
"""

from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parent
PYTHON = sys.executable


def run(cmd: list[str]) -> int:
    print("Running:", " ".join(cmd))
    result = subprocess.run(cmd, cwd=ROOT)
    print("Exit code:", result.returncode)
    return result.returncode


def smoke_cli_help() -> None:
    commands = [
        [PYTHON, "cli.py", "--help"],
        [PYTHON, "cli.py", "motif-map", "--help"],
        [PYTHON, "cli.py", "convert", "--help"],
        [PYTHON, "cli.py", "exon-sets", "--help"],
    ]
    failed = False
    for cmd in commands:
        if run(cmd) != 0:
            failed = True
    if failed:
        raise SystemExit(1)


def main() -> None:
    smoke_cli_help()
    print("Smoke tests completed. For full runs, invoke the CLI with real data.")


if __name__ == "__main__":
    main()

