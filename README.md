# rmaps_refactored

Refactored version of the original RBP motif-mapping pipeline.

Key changes:

- **Python 3 only**: All core scripts have been updated to run on Python 3.
- **pyfaidx instead of pygr**: Genome access now goes through `pyfaidx.Fasta` via `genome_access.py`, using a `genomedata/`-style layout (e.g. `genomedata/hg19/hg19.fa`).
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

Motif-map commands expect a `--pygr-station` / `--pygrStation` argument that now points to your `genomedata/` FASTA root (kept for compatibility with existing scripts and logs).

The legacy `motifMap*.py` scripts are still present but should be considered deprecated in favor of the new CLI.

