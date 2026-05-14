"""
Core dispatcher for CLIP-seq RNA map analysis.

This module provides the dispatcher logic for running CLIP-seq analysis
across all 5 alternative splicing event types (SE, A3SS, A5SS, RI, MXE).
"""

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
class ClipEventSpec:
    """
    Specification for each CLIP-seq event type.

    Attributes:
        name: Full event name (e.g., "SE", "A5SS")
        script: Path to the legacy clip processing script
        description: Human-readable description of the event type
    """

    name: str
    script: str
    description: str


CLIP_EVENT_SPECS: Dict[str, ClipEventSpec] = {
    "se": ClipEventSpec(
        name="SE",
        script="legacy/clipSeqSE.py",
        description="Skipped Exon (Exon Skipping)",
    ),
    "a3ss": ClipEventSpec(
        name="A3SS",
        script="legacy/clipSeqA3SS.py",
        description="Alternative 3' Splice Site",
    ),
    "a5ss": ClipEventSpec(
        name="A5SS",
        script="legacy/clipSeqA5SS.py",
        description="Alternative 5' Splice Site",
    ),
    "ri": ClipEventSpec(
        name="RI",
        script="legacy/clipSeqRI.py",
        description="Retained Intron",
    ),
    "mxe": ClipEventSpec(
        name="MXE",
        script="legacy/clipSeqMXE.py",
        description="Mutually Exclusive Exons",
    ),
}


def clip_event_script(event: str) -> Path:
    """
    Return the path to the CLIP-seq processing script for the given event.

    Args:
        event: Event type (se, a3ss, a5ss, ri, mxe)

    Returns:
        Path to the clip processing script

    Raises:
        ValueError: If event type is not supported
    """
    key = event.lower()
    if key not in CLIP_EVENT_SPECS:
        raise ValueError(f"Unsupported event type: {event}")
    return REPO_ROOT / CLIP_EVENT_SPECS[key].script


def run_subprocess(
    cmd: list[str],
    env_overrides: dict[str, str] | None = None,
    base_cwd: Path | None = None,
) -> int:
    """
    Shared helper for launching child processes from this project.

    Args:
        cmd: Command list to execute

    Returns:
        Return code from the subprocess
    """
    env = build_subprocess_env(REPO_ROOT, env_overrides)
    result = subprocess.run(cmd, cwd=base_cwd or Path.cwd(), env=env)
    return result.returncode


def run_clip_map(
    event: str,
    peak: Path,
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
    Run CLIP-seq RNA map analysis for the specified event type.

    Args:
        event: Event type (se, a3ss, a5ss, ri, mxe)
        peak: Path to CLIP-seq peak file
        output: Output directory path
        rmats: rMATS file path or "NA"
        miso: MISO file path or "NA"
        up: Upregulated exons coordinate file or "NA"
        down: Downregulated exons coordinate file or "NA"
        background: Background exons coordinate file or "NA"
        label: RNA-binding protein label
        intron: Intron length to examine (default: 250)
        exon: Exon length to examine (default: 50)
        window: Window size for analysis (default: 10)
        step: Step size for sliding window (default: 1)
        sig_fdr: FDR cutoff for significant events
        sig_delta_psi: Delta PSI cutoff for significant events
        separate: Whether to create separate p-value plots

    Returns:
        Return code from the subprocess
    """
    script = clip_event_script(event)
    stat_method = normalize_stat_method(stat_method)
    base_cwd = base_cwd or Path.cwd()
    peak = resolve_user_path(peak, base_cwd)
    output = Path(resolve_user_path(output, base_cwd))
    rmats = resolve_user_path(rmats, base_cwd)
    miso = resolve_user_path(miso, base_cwd)
    up = resolve_user_path(up, base_cwd)
    down = resolve_user_path(down, base_cwd)
    background = resolve_user_path(background, base_cwd)
    rmats = maybe_prepare_rmats_input(rmats, output)

    cmd = [
        PYTHON,
        str(script),
        "--peak", peak,
        "--output", str(output),
        "--rMATS", rmats,
        "--miso", miso,
        "--up", up,
        "--down", down,
        "--background", background,
        "--label", label,
        "--intron", str(intron),
        "--exon", str(exon),
        "--window", str(window),
        "--step", str(step),
        "--sigFDR", str(sig_fdr),
        "--sigDeltaPSI", str(sig_delta_psi),
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


def get_event_description(event: str) -> str:
    """
    Get the human-readable description for an event type.

    Args:
        event: Event type (se, a3ss, a5ss, ri, mxe)

    Returns:
        Description string

    Raises:
        ValueError: If event type is not supported
    """
    key = event.lower()
    if key not in CLIP_EVENT_SPECS:
        raise ValueError(f"Unsupported event type: {event}")
    return CLIP_EVENT_SPECS[key].description


def list_supported_events() -> list[str]:
    """
    Get list of all supported event types.

    Returns:
        List of event type keys (lowercase)
    """
    return list(CLIP_EVENT_SPECS.keys())
