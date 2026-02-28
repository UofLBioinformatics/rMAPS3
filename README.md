# rMAPS Refactored

A Python 3 refactor of the original rMAPS motif-mapping pipeline.

## What This Project Provides

- Python 3 support across legacy motif-map event engines.
- `pyfaidx`-based genome access (no `pygr` runtime dependency).
- Unified CLI (`cli.py`) for:
  - motif-map generation (`se`, `a3ss`, `a5ss`, `ri`, `mxe`)
  - MISO-to-rMATS conversion
  - exon-set generation

## Requirements

- Python 3.10+ recommended
- OS: Linux/macOS/Windows
- Perl available on `PATH` for MISO conversion (`miso2rMATS*.pl`)
- Flask for Web UI (`run_web.py`)
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
- `motif_map_core.py`: centralized routing/execution for event workflows and converters
- `genome_access.py`: shared `pyfaidx` sequence fetch layer
- `drawutils.py`: plotting/export helpers
- `legacy/`: event-specific motif scripts (`motifMap*_MP.py`)
- `bin/`: converter/util scripts
- `tests/`: smoke and integration test runners

## CLI Usage

CLI reference is documented in:

- [`CLI_USAGE.md`](CLI_USAGE.md)

## Web UI (Local)

Run the local web server:

```bash
python run_web.py
```

Open:

- `http://127.0.0.1:5000`

Full Web UI documentation is in:

- [`webui/README.md`](webui/README.md)

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

Note: shell test runners under `tests/` require a Bash-compatible environment (Linux/macOS, WSL, or Git Bash on Windows).

## Documentation

- Installation: [`INSTALL.md`](INSTALL.md)
- Quick start: [`QUICKSTART.md`](QUICKSTART.md)
- CLI reference: [`CLI_USAGE.md`](CLI_USAGE.md)
- Web UI usage and API: [`webui/README.md`](webui/README.md)
- Test scripts and matrix: [`tests/README.md`](tests/README.md)

## Operational Notes

- Use the CLI (`cli.py`) as the public interface.
- For local browser-based runs, use `run_web.py`.
- Scripts under `legacy/` are internal implementation details; use the CLI for production runs.
- `--fasta-root` is the canonical argument; `--fastaRoot` is accepted for compatibility.

## Troubleshooting

- `FastaNotFoundError`:
  - Verify `--genome` and `--fasta-root` match expected layout (`<build>/<build>.fa`).
- MISO conversion fails:
  - Ensure Perl is installed and available in `PATH`.
  - Verify `bin/miso2rMATS.*.pl` scripts exist.
- PNG files missing:
  - Native Ghostscript export may be unavailable.
  - Install optional rendering dependencies from `requirements.txt`.
