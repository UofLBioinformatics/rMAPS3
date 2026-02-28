# Testing Guide

All maintained test-run shell scripts live in this `tests/` folder.

## Test data and prerequisites

- Run test scripts from repository root.
- Bash-compatible shell required for `.sh` runners (Linux/macOS, WSL, or Git Bash on Windows).
- Install dependencies from `requirements.txt`.
- Ensure Perl is available on `PATH` for MISO conversion tests (`miso2rMATS.*.pl`).
- Ensure FASTA files exist under your selected `--fasta-root` with this layout:
  - `genomedata/hg19/hg19.fa`
  - `genomedata/hg38/hg38.fa`
  - `genomedata/mm10/mm10.fa`
  - `genomedata/dm3/dm3.fa`
- Event-specific rMATS test inputs used by scripts are expected at:
  - `temp/SE.MATS.ReadsOnTargetAndJunctionCounts.txt`
  - `temp/A3SS.MATS.ReadsOnTargetAndJunctionCounts.txt`
  - `temp/A5SS.MATS.ReadsOnTargetAndJunctionCounts.txt`
  - `temp/RI.MATS.ReadsOnTargetAndJunctionCounts.txt`
  - `temp/MXE.MATS.ReadsOnTargetAndJunctionCounts.txt`
- MISO test inputs used by scripts are expected at:
  - `testData/ESRP.OE.miso_bf`
  - `testData/example.dm3.miso_bf`

## Script overview

- `tests/run_all_events.sh`
  - Master entrypoint.
  - Runs all event rMATS/user-input test matrices and the SE MISO matrix.

- `tests/run_event_pair.sh <event> <genome> <rmats_file> <label>`
  - Shared event test helper called by per-event wrappers.
  - Runs four checks for one event:
    - rMATS with optional motifs (`data/ESRP.like.motif.txt`)
    - user-input replay using generated `exon/up.coord.txt`, `dn.coord.txt`, `bg.coord.txt` from previous run
    - rMATS with no optional motifs (`--motifs NA`)
    - user-input replay for no-optional branch
  - Uses known motifs from `data/testMotifs.txt`.
  - Writes outputs to `temp/tests_<event>_<timestamp>/...`.

- Per-event wrappers:
  - `tests/run_se.sh`
  - `tests/run_a3ss.sh`
  - `tests/run_a5ss.sh`
  - `tests/run_ri.sh`
  - `tests/run_mxe.sh`
  - Each wrapper invokes `run_event_pair.sh` with event-specific genome + rMATS file.

- `tests/run_se_miso.sh`
  - Dedicated SE MISO matrix.
  - Runs four checks:
    - hg19 with optional motifs
    - hg19 with no optional motifs
    - dm3 with optional motifs
    - dm3 with no optional motifs
  - Uses known motifs from `data/knownMotifs.human.mouse.txt`.
  - Writes outputs to `temp/tests_se_miso_<timestamp>/...`.

- `tests/smoke_cli.py`
  - Lightweight smoke test.
  - Checks CLI entrypoints/help commands without running long motif-map jobs.

## Common commands

Run everything:

```bash
bash tests/run_all_events.sh
```

Run one event matrix:

```bash
bash tests/run_se.sh
bash tests/run_a3ss.sh
bash tests/run_a5ss.sh
bash tests/run_ri.sh
bash tests/run_mxe.sh
```

Run SE MISO matrix only:

```bash
bash tests/run_se_miso.sh
```

Run helper directly for a custom event/input:

```bash
bash tests/run_event_pair.sh se hg38 ./temp/SE.MATS.ReadsOnTargetAndJunctionCounts.txt motif
```

Run lightweight smoke checks:

```bash
python tests/smoke_cli.py
```

## Expected outputs

- Each script creates a timestamped folder under `temp/`.
- Typical outputs include:
  - `exon/*.coord.txt`
  - `fasta/*.fasta`
  - `maps/*.pdf`
  - `log.motifMap.txt`
- Non-fatal PNG export warnings can appear if Ghostscript is not installed; PDF generation is the primary expected artifact.
