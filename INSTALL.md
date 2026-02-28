# Installation

## 1. Prerequisites

Required:
- Python 3.10+ (recommended)
- `pip`

Optional, depending on workflow:
- `perl` for MISO conversion
- `flask` for Web UI (`run_web.py`)
- TeX (MiKTeX/TeX Live) for native PyX text rendering
- Ghostscript for native PyX PNG export

## 2. Install Dependencies

From repository root:

```bash
pip install -r requirements.txt
```

## 3. Verify CLI

```bash
python cli.py --help
python cli.py motif-map --help
python cli.py convert --help
python cli.py exon-sets --help
python tests/smoke_cli.py
```

## 4. Genome FASTA Setup

Use this layout:

```text
genomedata/
  hg19/hg19.fa
  hg38/hg38.fa
  mm10/mm10.fa
  dm3/dm3.fa
  ...
```

Important:
- folder name must match genome build
- FASTA filename must match build (`<build>.fa`)

## 5. Verify Web UI (Optional)

```bash
python run_web.py
```

Open:
- `http://127.0.0.1:5000`

## 6. MISO Support (Perl)

MISO conversion calls bundled Perl converters in `bin/`.

If `perl` is not in `PATH`, MISO jobs fail.

PowerShell example for portable Perl:

```powershell
$env:PATH = "C:\myperl\bin;C:\myperl\perl\bin;C:\myperl\c\bin;" + $env:PATH
```

## 7. Common Issues

### FASTA not found

Confirm `--fasta-root` and `--genome` point to `genomedata/<genome>/<genome>.fa`.

### MISO conversion fails

Confirm Perl is available in the shell/server environment.

### PNG export warnings

Ghostscript may be missing. PDF output is still produced.

### Port 5000 in use

Stop the process using port `5000`, then restart `run_web.py`.

## 8. Next Docs

- CLI reference: `CLI_USAGE.md`
- Quick start: `QUICKSTART.md`
- Web UI details: `webui/README.md`
