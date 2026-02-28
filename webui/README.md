# Web UI Guide

This folder contains the local Flask web interface for running motif-map analyses.

## Start

From repository root:

```bash
python run_web.py
```

Open:

- `http://127.0.0.1:5000`

## Features

- Run motif-map jobs for: `SE`, `A3SS`, `A5SS`, `RI`, `MXE`
- Input modes:
  - rMATS upload
  - MISO upload
  - coordinate file upload (`up`, `down`, `bg`)
- Live status polling
- Live log tail during job execution
- One-click local test run
- Load existing job by ID
- Recent jobs dropdown
- Browser history support (`back`/`forward`) for job views

## Environment Variables

- `RMAPS_FASTA_ROOT`:
  FASTA root directory (default: `genomedata/`)
- `RMAPS_RESULTS_DIR`:
  output directory for web jobs (default: `results/`)
- `RMAPS_QUICKTEST_DIR`:
  one-click test data directory (default: `testData/`)

## Key Endpoints

- `GET /`:
  Web UI
- `GET /api/genomes`:
  Supported genomes + local availability
- `GET /api/events`:
  Supported event types
- `POST /api/submit`:
  Submit user job
- `POST /api/quick-test/run`:
  Submit one-click test job
- `GET /api/status/<job_id>`:
  Job status
- `GET /api/logs/<job_id>`:
  Filtered live logs
- `GET /api/results/<job_id>`:
  Top summary rows (with legacy fallback parsing)
- `GET /api/jobs`:
  Recent jobs list

## Job Output Layout

Each web job writes under:

```text
results/<job_id>/           or results/quicktest_<job_id>/
  exon/
  fasta/
  maps/
  temp/
  log.motifMap.txt
  pVal.up.vs.bg.RNAmap.txt
  pVal.dn.vs.bg.RNAmap.txt
```

## Troubleshooting

- Recent jobs dropdown empty:
  - Verify `RMAPS_RESULTS_DIR` points to the correct results directory.
  - Check `GET /api/jobs` response.
- Job load says not found:
  - Ensure the results folder exists in `RMAPS_RESULTS_DIR`.
  - Use the recent jobs dropdown value directly.
- Live logs not updating:
  - Confirm `log.motifMap.txt` is being written in the job folder.
  - Hard refresh browser after JS updates (`Ctrl+F5`).
- MISO run fails:
  - Ensure Perl is available in the server process `PATH`.

