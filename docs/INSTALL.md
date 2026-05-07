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

- MiKTeX or TeX Live: improves PyX text rendering in motif-map outputs
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

## 5. Optional: Start Local Web UI

Run local web UI:

```bash
python run_web.py
```

Then open `http://127.0.0.1:5000`.

## 6. Next Docs

- CLI examples and command reference: [`CLI_USAGE.md`](CLI_USAGE.md)
- Web UI details: [`../webui/README.md`](../webui/README.md)
- Test scripts and matrix: [`../tests/README.md`](../tests/README.md)
