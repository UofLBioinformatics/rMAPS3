# rmaps_refactored

Refactored version of the original RBP motif-mapping pipeline.

Key changes:

- **Python 3 only**: All core scripts have been updated to run on Python 3.
- **pyfaidx genome access layer**: Genome access now goes through `pyfaidx.Fasta` via `genome_access.py`, using a `genomedata/`-style layout (e.g. `genomedata/hg19/hg19.fa`).
- **Unified CLI**: The preferred entry point is the Typer-based CLI defined in `cli.py`.

## Installation

Create and activate a virtualenv, then install:

```bash
pip install -r requirements.txt
```

You will also need indexed FASTA files laid out like:

```text
genomedata/
  hg19/hg19.fa
  hg38/hg38.fa
  mm10/mm10.fa
  ...
```

## Usage

From the project root:

```bash
python cli.py --help
python cli.py motif-map se --help
python cli.py convert miso --help
python cli.py exon-sets se --help
```

Motif-map commands expect a `--fasta-root` / `--fastaRoot` argument that now points to your `genomedata/` FASTA root (kept for compatibility with existing scripts and logs).

The legacy `motifMap*.py` scripts are still present but should be considered deprecated in favor of the new CLI.

## Testing

Detailed test documentation is in [`tests/README.md`](tests/README.md), including:
- smoke checks (`python tests/smoke_cli.py`)
- per-event run scripts
- full end-to-end matrix (`bash tests/run_all_events.sh`)
- required test input files and expected outputs

