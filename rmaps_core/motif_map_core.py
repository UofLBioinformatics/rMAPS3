from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict

import sys
import subprocess
import shutil

from rmaps_core.input_utils import maybe_prepare_rmats_input
from rmaps_core.path_utils import build_subprocess_env, repo_root, resolve_user_path
from rmaps_core.stat_utils import normalize_stat_method


PYTHON = sys.executable
REPO_ROOT = repo_root()


@dataclass(frozen=True)
class EventSpec:
    """
    Minimal shared description for each event type.

    Right now this encodes:
    - which motifMap script to run
    - which miso2rMATS converter to use (for the standalone convert CLI)
    """

    name: str
    script: str
    miso_converter: str


EVENT_SPECS: Dict[str, EventSpec] = {
    "se": EventSpec(
        name="SE",
        script="legacy/motifMapSE_MP.py",
        miso_converter="bin/miso2rMATS.SE.pl",
    ),
    "a3ss": EventSpec(
        name="A3SS",
        script="legacy/motifMapA3SS_MP.py",
        miso_converter="bin/miso2rMATS.A3SS.pl",
    ),
    "a5ss": EventSpec(
        name="A5SS",
        script="legacy/motifMapA5SS_MP.py",
        miso_converter="bin/miso2rMATS.A5SS.pl",
    ),
    "ri": EventSpec(
        name="RI",
        script="legacy/motifMapRI_MP.py",
        miso_converter="bin/miso2rMATS.RI.pl",
    ),
    "mxe": EventSpec(
        name="MXE",
        script="legacy/motifMapMXE_MP.py",
        miso_converter="bin/miso2rMATS.MXE.pl",
    ),
}


def event_script(event: str) -> Path:
    """
    Return the path to the motifMap script for the given event (se, a3ss, a5ss, ri, mxe).
    """
    key = event.lower()
    if key not in EVENT_SPECS:
        raise ValueError(f"Unsupported event type: {event}")
    return REPO_ROOT / EVENT_SPECS[key].script


def miso_converter_script(event: str) -> Path:
    """
    Return the path to the Perl miso2rMATS converter for the given event.
    """
    key = event.lower()
    if key not in EVENT_SPECS:
        raise ValueError(f"Unsupported event type: {event}")
    return REPO_ROOT / EVENT_SPECS[key].miso_converter


def run_subprocess(
    cmd: list[str],
    env_overrides: dict[str, str] | None = None,
    base_cwd: Path | None = None,
) -> int:
    """
    Shared helper for launching child processes from this project.
    """
    env = build_subprocess_env(REPO_ROOT, env_overrides)
    result = subprocess.run(cmd, cwd=base_cwd or Path.cwd(), env=env)
    return result.returncode


def run_motif_map(
    event: str,
    known_motifs: Path,
    motifs: str,
    fasta_root: Path,
    genome: str,
    output: Path,
    rmats: str,
    miso: str,
    up: str,
    down: str,
    background: str,
    label: str,
    intron: int,
    exon: int,
    window: int,
    step: int,
    sig_fdr: float,
    sig_delta_psi: float,
    separate: bool,
    stat_method: str = "fisher",
    stat_permutations: int | None = None,
    stat_seed: int | None = None,
    keep_temp: bool = False,
    base_cwd: Path | None = None,
) -> int:
    """
    Build and run the legacy motifMap* script for a given event type.

    This keeps all event-specific wiring in one place so the Typer CLI
    can call a single entrypoint per event.
    """
    script_path = event_script(event)
    stat_method = normalize_stat_method(stat_method)
    base_cwd = base_cwd or Path.cwd()
    known_motifs = resolve_user_path(known_motifs, base_cwd)
    motifs = resolve_user_path(motifs, base_cwd)
    fasta_root = resolve_user_path(fasta_root, base_cwd)
    output = Path(resolve_user_path(output, base_cwd))
    rmats = resolve_user_path(rmats, base_cwd)
    miso = resolve_user_path(miso, base_cwd)
    up = resolve_user_path(up, base_cwd)
    down = resolve_user_path(down, base_cwd)
    background = resolve_user_path(background, base_cwd)
    rmats = maybe_prepare_rmats_input(rmats, output)
    cmd: list[str] = [
        PYTHON,
        str(script_path),
        "-k",
        known_motifs,
        "-m",
        motifs,
        "--fasta-root",
        fasta_root,
        "-g",
        genome,
        "-o",
        str(output),
        "-r",
        rmats,
        "-mi",
        miso,
        "-u",
        up,
        "-d",
        down,
        "-b",
        background,
        "--label",
        label,
        "--intron",
        str(intron),
        "--exon",
        str(exon),
        "--window",
        str(window),
        "--step",
        str(step),
        "--sigFDR",
        str(sig_fdr),
        "--sigDeltaPSI",
        str(sig_delta_psi),
    ]
    if separate:
        cmd.append("--separate")

    env_overrides = {"RMAPS_STAT_METHOD": stat_method}
    if stat_permutations is not None:
        env_overrides["RMAPS_STAT_PERMUTATIONS"] = str(stat_permutations)
    if stat_seed is not None:
        env_overrides["RMAPS_STAT_SEED"] = str(stat_seed)

    code = run_subprocess(cmd, env_overrides, base_cwd=base_cwd)

    # Keep temp data on failure for debugging; clean on success unless requested.
    if code == 0 and not keep_temp:
        shutil.rmtree(output / "temp", ignore_errors=True)

    return code

