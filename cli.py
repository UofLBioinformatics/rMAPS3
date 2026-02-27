import sys
import subprocess
from pathlib import Path
from typing import Optional

import typer

from motif_map_core import event_script, miso_converter_script, run_subprocess, run_motif_map


app = typer.Typer(help="RNA motif maps and exon set utilities")

motif_app = typer.Typer(help="Motif map generation for splicing event types")
convert_app = typer.Typer(help="Converters (e.g., MISO to rMATS)")
exon_app = typer.Typer(help="Build up/down/background exon sets from rMATS+cuffdiff")

app.add_typer(motif_app, name="motif-map")
app.add_typer(convert_app, name="convert")
app.add_typer(exon_app, name="exon-sets")

REPO_ROOT = Path(__file__).resolve().parent
PYTHON = sys.executable


def run_cmd(cmd: list[str]) -> None:
    """
    Run a subprocess in the repository root and exit with its return code.
    """
    code = run_subprocess(cmd)
    raise typer.Exit(code=code)


#
# motif-map subcommands
#


@motif_app.command("se")
def motif_map_se(
    known_motifs: Path = typer.Option(
        ...,
        "--known-motifs",
        "--knownMotifs",
        help="Known motif table (MotifName, regularExpression, ...).",
    ),
    motifs: str = typer.Option(
        "NA",
        "--motifs",
        "--motif",
        help="Optional motif list file; use 'NA' for none.",
    ),
    fasta_root: Path = typer.Option(
        ...,
        "--fasta-root",
        "--fastaRoot",
        help="FASTA root directory (genomedata root).",
    ),
    genome: str = typer.Option(
        ...,
        "--genome",
        help="Genome build (e.g. hg19, hg38, mm10, dm3).",
    ),
    output: Path = typer.Option(
        ...,
        "--output",
        "-o",
        help="Output directory.",
    ),
    rmats: str = typer.Option(
        "NA",
        "--rMATS",
        "--rmats",
        help="rMATS SE event file (or 'NA').",
    ),
    miso: str = typer.Option(
        "NA",
        "--miso",
        help="MISO SE event file (or 'NA').",
    ),
    up: str = typer.Option(
        "NA",
        "--up",
        help="Upregulated exon coordinate file (or 'NA').",
    ),
    down: str = typer.Option(
        "NA",
        "--down",
        "--dn",
        help="Downregulated exon coordinate file (or 'NA').",
    ),
    background: str = typer.Option(
        "NA",
        "--background",
        "--bg",
        help="Background exon coordinate file (or 'NA').",
    ),
    label: str = typer.Option(
        "RBP",
        "--label",
        help="Motif label, e.g. 'TG-rich'.",
    ),
    intron: int = typer.Option(
        250,
        "--intron",
        help="Intron length to examine around junctions.",
    ),
    exon: int = typer.Option(
        50,
        "--exon",
        help="Exon length to examine at each end.",
    ),
    window: int = typer.Option(
        50,
        "--window",
        help="Window size for motif scanning.",
    ),
    step: int = typer.Option(
        1,
        "--step",
        help="Step size for sliding window.",
    ),
    sig_fdr: float = typer.Option(
        0.05,
        "--sig-fdr",
        "--sigFDR",
        help="FDR cutoff for significant events.",
    ),
    sig_delta_psi: float = typer.Option(
        0.05,
        "--sig-delta-psi",
        "--sigDeltaPSI",
        help="Delta PSI cutoff for significant events.",
    ),
    separate: bool = typer.Option(
        False,
        "--separate",
        help="Draw separate p-value panels under motif maps.",
    ),
) -> None:
    """
    Generate motif maps for SE events.
    """
    code = run_motif_map(
        "se",
        known_motifs=known_motifs,
        motifs=motifs,
        fasta_root=fasta_root,
        genome=genome,
        output=output,
        rmats=rmats,
        miso=miso,
        up=up,
        down=down,
        background=background,
        label=label,
        intron=intron,
        exon=exon,
        window=window,
        step=step,
        sig_fdr=sig_fdr,
        sig_delta_psi=sig_delta_psi,
        separate=separate,
    )
    raise typer.Exit(code=code)


@motif_app.command("a3ss")
def motif_map_a3ss(
    known_motifs: Path = typer.Option(
        ...,
        "--known-motifs",
        "--knownMotifs",
        help="Known motif table.",
    ),
    motifs: str = typer.Option(
        "NA",
        "--motifs",
        "--motif",
        help="Optional motif list file; use 'NA' for none.",
    ),
    fasta_root: Path = typer.Option(
        ...,
        "--fasta-root",
        "--fastaRoot",
        help="FASTA root directory.",
    ),
    genome: str = typer.Option(..., "--genome"),
    output: Path = typer.Option(..., "--output", "-o"),
    rmats: str = typer.Option("NA", "--rMATS", "--rmats"),
    miso: str = typer.Option("NA", "--miso"),
    up: str = typer.Option("NA", "--up"),
    down: str = typer.Option("NA", "--down", "--dn"),
    background: str = typer.Option("NA", "--background", "--bg"),
    label: str = typer.Option("RBP", "--label"),
    intron: int = typer.Option(250, "--intron"),
    exon: int = typer.Option(50, "--exon"),
    window: int = typer.Option(50, "--window"),
    step: int = typer.Option(1, "--step"),
    sig_fdr: float = typer.Option(0.05, "--sig-fdr", "--sigFDR"),
    sig_delta_psi: float = typer.Option(0.05, "--sig-delta-psi", "--sigDeltaPSI"),
    separate: bool = typer.Option(False, "--separate"),
) -> None:
    """
    Generate motif maps for A3SS events.
    """
    code = run_motif_map(
        "a3ss",
        known_motifs=known_motifs,
        motifs=motifs,
        fasta_root=fasta_root,
        genome=genome,
        output=output,
        rmats=rmats,
        miso=miso,
        up=up,
        down=down,
        background=background,
        label=label,
        intron=intron,
        exon=exon,
        window=window,
        step=step,
        sig_fdr=sig_fdr,
        sig_delta_psi=sig_delta_psi,
        separate=separate,
    )
    raise typer.Exit(code=code)


@motif_app.command("a5ss")
def motif_map_a5ss(
    known_motifs: Path = typer.Option(..., "--known-motifs", "--knownMotifs"),
    motifs: str = typer.Option("NA", "--motifs", "--motif"),
    fasta_root: Path = typer.Option(..., "--fasta-root", "--fastaRoot"),
    genome: str = typer.Option(..., "--genome"),
    output: Path = typer.Option(..., "--output", "-o"),
    rmats: str = typer.Option("NA", "--rMATS", "--rmats"),
    miso: str = typer.Option("NA", "--miso"),
    up: str = typer.Option("NA", "--up"),
    down: str = typer.Option("NA", "--down", "--dn"),
    background: str = typer.Option("NA", "--background", "--bg"),
    label: str = typer.Option("RBP", "--label"),
    intron: int = typer.Option(250, "--intron"),
    exon: int = typer.Option(50, "--exon"),
    window: int = typer.Option(50, "--window"),
    step: int = typer.Option(1, "--step"),
    sig_fdr: float = typer.Option(0.05, "--sig-fdr", "--sigFDR"),
    sig_delta_psi: float = typer.Option(0.05, "--sig-delta-psi", "--sigDeltaPSI"),
    separate: bool = typer.Option(False, "--separate"),
) -> None:
    """
    Generate motif maps for A5SS events.
    """
    code = run_motif_map(
        "a5ss",
        known_motifs=known_motifs,
        motifs=motifs,
        fasta_root=fasta_root,
        genome=genome,
        output=output,
        rmats=rmats,
        miso=miso,
        up=up,
        down=down,
        background=background,
        label=label,
        intron=intron,
        exon=exon,
        window=window,
        step=step,
        sig_fdr=sig_fdr,
        sig_delta_psi=sig_delta_psi,
        separate=separate,
    )
    raise typer.Exit(code=code)


@motif_app.command("ri")
def motif_map_ri(
    known_motifs: Path = typer.Option(..., "--known-motifs", "--knownMotifs"),
    motifs: str = typer.Option("NA", "--motifs", "--motif"),
    fasta_root: Path = typer.Option(..., "--fasta-root", "--fastaRoot"),
    genome: str = typer.Option(..., "--genome"),
    output: Path = typer.Option(..., "--output", "-o"),
    rmats: str = typer.Option("NA", "--rMATS", "--rmats"),
    miso: str = typer.Option("NA", "--miso"),
    up: str = typer.Option("NA", "--up"),
    down: str = typer.Option("NA", "--down", "--dn"),
    background: str = typer.Option("NA", "--background", "--bg"),
    label: str = typer.Option("RBP", "--label"),
    intron: int = typer.Option(250, "--intron"),
    exon: int = typer.Option(50, "--exon"),
    window: int = typer.Option(50, "--window"),
    step: int = typer.Option(1, "--step"),
    sig_fdr: float = typer.Option(0.05, "--sig-fdr", "--sigFDR"),
    sig_delta_psi: float = typer.Option(0.05, "--sig-delta-psi", "--sigDeltaPSI"),
    separate: bool = typer.Option(False, "--separate"),
) -> None:
    """
    Generate motif maps for RI events.
    """
    code = run_motif_map(
        "ri",
        known_motifs=known_motifs,
        motifs=motifs,
        fasta_root=fasta_root,
        genome=genome,
        output=output,
        rmats=rmats,
        miso=miso,
        up=up,
        down=down,
        background=background,
        label=label,
        intron=intron,
        exon=exon,
        window=window,
        step=step,
        sig_fdr=sig_fdr,
        sig_delta_psi=sig_delta_psi,
        separate=separate,
    )
    raise typer.Exit(code=code)


@motif_app.command("mxe")
def motif_map_mxe(
    known_motifs: Path = typer.Option(..., "--known-motifs", "--knownMotifs"),
    motifs: str = typer.Option("NA", "--motifs", "--motif"),
    fasta_root: Path = typer.Option(..., "--fasta-root", "--fastaRoot"),
    genome: str = typer.Option(..., "--genome"),
    output: Path = typer.Option(..., "--output", "-o"),
    rmats: str = typer.Option("NA", "--rMATS", "--rmats"),
    miso: str = typer.Option("NA", "--miso"),
    up: str = typer.Option("NA", "--up"),
    down: str = typer.Option("NA", "--down", "--dn"),
    background: str = typer.Option("NA", "--background", "--bg"),
    label: str = typer.Option("RBP", "--label"),
    intron: int = typer.Option(250, "--intron"),
    exon: int = typer.Option(50, "--exon"),
    window: int = typer.Option(50, "--window"),
    step: int = typer.Option(1, "--step"),
    sig_fdr: float = typer.Option(0.05, "--sig-fdr", "--sigFDR"),
    sig_delta_psi: float = typer.Option(0.05, "--sig-delta-psi", "--sigDeltaPSI"),
    separate: bool = typer.Option(False, "--separate"),
) -> None:
    """
    Generate motif maps for MXE events.
    """
    code = run_motif_map(
        "mxe",
        known_motifs=known_motifs,
        motifs=motifs,
        fasta_root=fasta_root,
        genome=genome,
        output=output,
        rmats=rmats,
        miso=miso,
        up=up,
        down=down,
        background=background,
        label=label,
        intron=intron,
        exon=exon,
        window=window,
        step=step,
        sig_fdr=sig_fdr,
        sig_delta_psi=sig_delta_psi,
        separate=separate,
    )
    raise typer.Exit(code=code)


#
# convert subcommands
#


@convert_app.command("miso")
def convert_miso(
    event: str = typer.Option(
        ...,
        "--event",
        case_sensitive=False,
        help="Event type: se, a3ss, a5ss, ri, mxe.",
    ),
    miso: Path = typer.Option(..., "--in", help="MISO input file."),
    output: Path = typer.Option(..., "--out", help="Output rMATS-like file."),
    bayes_low: int = typer.Option(1, "--bayes-low"),
    bayes_high: int = typer.Option(100, "--bayes-high"),
) -> None:
    """
    Run the appropriate Perl miso2rMATS converter directly.
    """
    event_l = event.lower()
    try:
        script = miso_converter_script(event_l)
    except ValueError as e:
        typer.echo(str(e), err=True)
        raise typer.Exit(code=1)
    cmd = [
        "perl",
        str(script),
        str(bayes_low),
        str(bayes_high),
        str(miso),
        str(output),
    ]
    run_cmd(cmd)


#
# exon-sets subcommands
#


@exon_app.command("se")
def exon_sets_se(
    input_file: Path = typer.Option(
        ...,
        "--input",
        help="Combined rMATS+cuffdiff file (common.txt).",
    ),
    sample1: str = typer.Option(..., "--sample1"),
    sample2: str = typer.Option(..., "--sample2"),
    out_dir: Path = typer.Option(..., "--out"),
) -> None:
    """
    Build up/down/background SE exon sets using getExonSets.py.
    """
    script = REPO_ROOT / "bin" / "getExonSets.py"
    cmd = [
        PYTHON,
        str(script),
        str(input_file),
        sample1,
        sample2,
        str(out_dir),
    ]
    run_cmd(cmd)


if __name__ == "__main__":
    app()

