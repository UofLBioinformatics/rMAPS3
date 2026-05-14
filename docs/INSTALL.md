# Installation

## 1. Python Environment (Recommended)

Using a virtual environment is recommended but not required.

```bash
python -m venv .venv
```

Activate:

- Windows PowerShell:

```powershell
.\.venv\Scripts\Activate.ps1
```

- Linux/macOS:

```bash
source .venv/bin/activate
```

## 2. Install Dependencies

```bash
python -m pip install -r requirements.txt
```

Optional system dependencies (not installed via `pip`):

- MiKTeX or TeX Live: improves PyX text rendering in motif-map outputs.
  If TeX fonts are unavailable, rMAPS3 falls back to simplified Pillow-based
  motif-map PDF/PNG rendering.
- Ghostscript: enables native PNG export path

## 3. Prepare Genome FASTA Data

Layout must be:

```text
genomedata/
  <build>/<build>.fa
```

Examples:

```text
genomedata/hg19/hg19.fa
genomedata/hg38/hg38.fa
genomedata/mm10/mm10.fa
```

### Download genomes using fetch scripts

**Windows (PowerShell):**

```powershell
# Download specific genomes
.\scripts\fetch_genomes.ps1 -Genomes hg19,hg38

# Or download all available genomes
.\scripts\fetch_genomes.ps1 -Genomes all

# Preview what would download without actually downloading
.\scripts\fetch_genomes.ps1 -Genomes hg19 -DryRun

# Force re-download if already exists
.\scripts\fetch_genomes.ps1 -Genomes hg19 -Force
```

**Linux/macOS:**

```bash
# Download specific genomes
./scripts/fetch_genomes.sh --genomes hg19,hg38

# Or download all available genomes
./scripts/fetch_genomes.sh --genomes all

# Preview what would download
./scripts/fetch_genomes.sh --genomes hg19 --dry-run

# Force re-download if already exists
./scripts/fetch_genomes.sh --genomes hg19 --force
```

The script automatically:
- Downloads from UCSC/Ensembl reference sources
- Verifies SHA256 checksums
- Decompresses `.fa.gz` → `.fa`
- Creates `<build>/<build>.fa` structure

Available genomes: `hg19`, `hg38`, `mm10`, `dm3`, `dm6`, `rn6`, `rn7`, `galGal5`, `galGal6`, `danRer10`, `danRer11`, `ce11`, `xenLae2`, `xenTro7`, `xenTro9`, `bosTau9`, `susScr11`, `araTha1`, `oSa7`

For complete list with URLs and verification hashes, see `scripts/genomes.manifest.tsv`

## 4. Verify Installation

Quick CLI smoke check:

```bash
python tests/smoke_cli.py
```

Recommended pre-deployment validation:

```bash
python tests/test_clip.py
python tests/test_motif.py --fasta-root genomedata --genome hg19
```

Optional command checks:

```bash
python cli.py --help
python cli.py motif-map --help
python cli.py convert --help
python cli.py exon-sets --help
```

## 5. Shared HPC / Read-Only Installs

rMAPS3 can be installed once in a shared location and run from a separate working directory, such as a scratch or job directory.

Recommended pattern:

```bash
cd /scratch/$USER/rmaps-job-001
python /path/to/rMAPS3/cli.py motif-map se \
  --known-motifs /path/to/rMAPS3/data/knownMotifs.human.mouse.txt \
  --motifs NA \
  --fasta-root /shared/genomedata \
  --genome hg38 \
  --output ./results/motif_se \
  --rMATS ./inputs/SE.MATS.JC.txt \
  --miso NA --up NA --down NA --background NA
```

Path behavior:
- rMAPS3 auto-detects its installed repository/package location for internal scripts.
- Relative user input and output paths are resolved from the directory where you run the command.
- Absolute paths are used as-is.
- Optional file arguments still accept `NA`.
- Child processes run from the user's invocation directory, not the installed repository directory.

This avoids writing job outputs into a shared or read-only rMAPS3 installation.

## 6. Optional: Start Local Web UI

Run local web UI:

```bash
python run_web.py
```

Then open `http://127.0.0.1:5000`.

## 7. Next Docs

- CLI examples and command reference: [`CLI_USAGE.md`](CLI_USAGE.md)
- Web UI details: [`../webui/README.md`](../webui/README.md)
- Test scripts and matrix: [`../tests/README.md`](../tests/README.md)
