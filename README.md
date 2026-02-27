# rMAPS Refactored

Production-oriented Python 3 refactor of the original rMAPS motif-mapping pipeline.

## What This Project Provides

- Python 3 support across legacy motif-map event engines.
- `pyfaidx`-based genome access (no `pygr` runtime dependency).
- Unified CLI (`cli.py`) for:
  - motif-map generation (`se`, `a3ss`, `a5ss`, `ri`, `mxe`)
  - MISO-to-rMATS conversion
  - exon-set generation
- Cross-platform-safe multiprocessing flow in event scripts.
- Output export fallbacks:
  - PDF: native PyX output, then PNG-to-PDF fallback if needed
  - PNG: native Ghostscript path, then PDF-to-PNG fallback via `pypdfium2`

## Requirements

- Python 3.10+ recommended
- OS: Linux/macOS/Windows
- Perl available on `PATH` for MISO conversion (`miso2rMATS*.pl`)
- Optional Ghostscript for native PNG export (not required with Python fallback)

Install dependencies:

```bash
pip install -r requirements.txt
```

## Genome Data Layout

This repo expects a FASTA root directory passed as `--fasta-root` (or `--fastaRoot`) with:

```text
genomedata/
  hg19/hg19.fa
  hg38/hg38.fa
  mm10/mm10.fa
  dm3/dm3.fa
  ...
```

Important:
- Build name and FASTA name must match (`<build>/<build>.fa`).
- If you run with `--genome hg19`, the loader uses `genomedata/hg19/hg19.fa`.

## Project Structure

- `cli.py`: primary user entrypoint
- `motif_map_core.py`: centralized routing/execution for event engines and converters
- `genome_access.py`: shared `pyfaidx` sequence fetch layer
- `drawutils.py`: plotting helpers and export fallback logic
- `legacy/`: event-specific motif engines (`motifMap*_MP.py`)
- `bin/`: converter/util scripts
- `tests/`: smoke and integration test runners

## CLI Usage

Top-level help:

```bash
python cli.py --help
```

Motif-map help:

```bash
python cli.py motif-map --help
python cli.py motif-map se --help
```

Converters:

```bash
python cli.py convert --help
python cli.py convert miso --help
```

Exon sets:

```bash
python cli.py exon-sets --help
python cli.py exon-sets se --help
```

### Example: SE Motif Map from rMATS

```bash
python cli.py motif-map se \
  --known-motifs data/knownMotifs.human.mouse.txt \
  --motifs data/ESRP.like.motif.txt \
  --fasta-root genomedata \
  --genome hg19 \
  --output temp/example_se \
  --rMATS testData/SE.MATS.ReadsOnTargetAndJunctionCounts.txt \
  --miso NA \
  --up NA --down NA --background NA
```

### Example: SE Motif Map from MISO

```bash
python cli.py motif-map se \
  --known-motifs data/knownMotifs.human.mouse.txt \
  --motifs data/ESRP.like.motif.txt \
  --fasta-root genomedata \
  --genome hg19 \
  --output temp/example_se_miso \
  --rMATS NA \
  --miso testData/ESRP.OE.miso_bf \
  --up NA --down NA --background NA
```

### Example: A3SS Event

```bash
python cli.py motif-map a3ss \
  --known-motifs data/testMotifs.txt \
  --motifs NA \
  --fasta-root genomedata \
  --genome mm10 \
  --output temp/example_a3ss \
  --rMATS temp/A3SS.MATS.ReadsOnTargetAndJunctionCounts.txt \
  --miso NA \
  --up NA --down NA --background NA
```

## Outputs

Each run writes to the given `--output` directory:

- `exon/`: `up.coord.txt`, `dn.coord.txt`, `bg.coord.txt`
- `fasta/`: per-region FASTA windows
- `maps/`: motif map outputs (`.pdf`, `.png` where available)
- `temp/`: intermediate files
- `log.motifMap.txt`: run log

## Testing

Detailed test documentation and commands are in:
- [`tests/README.md`](tests/README.md)

Quick smoke:

```bash
python tests/smoke_cli.py
```

Full event matrix:

```bash
bash tests/run_all_events.sh
```

Windows PowerShell SE legacy matrix:

```powershell
powershell -ExecutionPolicy Bypass -File tests/run_legacy_se_matrix.ps1
```

## Operational Notes

- Use the CLI (`cli.py`) as the public interface.
- Legacy scripts under `legacy/` are internal engines; direct invocation is not recommended.
- `--fasta-root` is the canonical argument; `--fastaRoot` is accepted for compatibility.
- `__pycache__` and genome data contents are ignored in git (folder skeleton retained with `.gitkeep`).

## Troubleshooting

- `FastaNotFoundError`:
  - Verify `--genome` and `--fasta-root` match expected layout (`<build>/<build>.fa`).
- MISO conversion fails:
  - Ensure Perl is installed and available in `PATH`.
  - Verify `bin/miso2rMATS.*.pl` scripts exist.
- PNG files missing:
  - Native Ghostscript export may be unavailable.
  - Ensure `pypdfium2` is installed for PNG fallback from PDF.
- Pylance warnings for dynamic globals in legacy files:
  - Some globals are initialized at runtime by `setup_runtime()`; this is expected in legacy engines.

## Migration Status

- Fully migrated user surface: CLI/core/genome access/testing.
- Legacy event logic is retained for behavioral parity and gradual modernization.

