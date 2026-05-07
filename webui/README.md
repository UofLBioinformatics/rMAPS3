# Web UI Guide

This folder contains the local Flask web interface for rMAPS 3 motif-map and CLIP-seq RNA-map analyses.

## Start

From repository root:

```bash
python run_web.py
```

Open:

- `http://127.0.0.1:5000`

## Features

- Top toggle between:
  - Motif Map
  - CLIP-seq Map
- Run jobs for: `SE`, `A3SS`, `A5SS`, `RI`, `MXE`
- Input modes:
  - rMATS upload
  - MISO upload
  - coordinate file upload (`up`, `down`, `bg`)
- CLIP-seq mode peak upload (`.bed`) for `clip-map`
- P-value method selector (`fisher`, `mannwhitney_greater`, `brunnermunzel_greater`, `permutation_one_sided`)
- Optional permutation controls (`stat_permutations`, `stat_seed`) for permutation runs
- Live status polling
- Live log tail during job execution
- One-click local test run
- Load existing job by ID
- Recent jobs dropdown
- Browser history support (`back`/`forward`) for job views

## One-Click Test Genome Behavior

- Motif quick test uses the genome selected in the UI genome dropdown.
- CLIP quick test does not require genome FASTA.
- If no genome is provided by the client, the server default is `hg19`.
- Required FASTA layout for the selected genome is:
  - `genomedata/<build>/<build>.fa`
  - Example for default: `genomedata/hg19/hg19.fa`

## Environment Variables

- `RMAPS_FASTA_ROOT`:
  FASTA root directory (default: `genomedata/`)
- `RMAPS_RESULTS_DIR`:
  output directory for web jobs (default: `results/`)
- `RMAPS_QUICKTEST_DIR`:
  one-click test data directory (default: `data/test/`)
- `RMAPS_MAX_UPLOAD_MB`:
  maximum upload size in MB (default: `200`)

## P-value Method

- The Analysis Parameters panel includes a **P-value Method** selector.
- Methods available: `fisher` (default), `mannwhitney_greater`, `brunnermunzel_greater`, `permutation_one_sided`.
- Optional fields `stat_permutations` and `stat_seed` are used by `permutation_one_sided` and ignored by other methods.
- Selected method is submitted as `stat_method` in `/api/submit` and used by both motif-map and clip-map runs.

Method guidance:
- Use `fisher` for backward-compatible and fast baseline runs.
- Use `mannwhitney_greater` or `brunnermunzel_greater` for rank/distribution-based alternatives.
- Use `permutation_one_sided` for empirical confirmation when runtime is less constrained.

Performance note:
- `permutation_one_sided` can take significantly longer than other methods, especially with high `stat_permutations`.

## Key Endpoints

- `GET /`:
  Web UI
- `GET /api/genomes`:
  Supported genomes + local availability
- `GET /api/events?analysis_type=motif|clip`:
  Supported event types for selected analysis mode
- `POST /api/submit`:
  Submit user job (includes `stat_method` field)
- `POST /api/quick-test/run`:
  Submit one-click test job
- `GET /api/status/<job_id>`:
  Job status
- `GET /api/logs/<job_id>`:
  Filtered live logs
- `GET /api/results/<job_id>`:
  Motif summary rows or CLIP output file list
- `GET /api/jobs`:
  Recent jobs list

## Job Output Layout

Each web job writes under:

```text
results/<job_id>/           or results/quicktest_<job_id>/
  motif jobs:
    exon/
    fasta/
    maps/
    temp/
    log.motifMap.txt
    pVal.up.vs.bg.RNAmap.txt
    pVal.dn.vs.bg.RNAmap.txt

  clip jobs:
    exon/
    temp/
    log.CLIPSeq3.0.1.txt
    *.RNAmap.txt
    *.pdf / *.eps
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
