# CLI Usage

This project's command-line interface is `cli.py`.

## Help

```bash
python cli.py --help
python cli.py motif-map --help
python cli.py convert --help
python cli.py exon-sets --help
```

## Motif Map Command Shape

```bash
python cli.py motif-map {se|a3ss|a5ss|ri|mxe} \
  --genome GENOME \
  --output OUTPUT_DIR \
  --known-motifs KNOWN_MOTIFS \
  --fasta-root FASTA_ROOT \
  [--motifs CUSTOM_MOTIFS_OR_NA] \
  [--label LABEL] [--window N] [--step N] [--intron N] [--exon N] \
  [--sig-fdr F] [--sig-delta-psi F] \
  (--rMATS RMATS_FILE | --miso MISO_FILE | --up UP --down DOWN --background BG)
```

Input modes are mutually exclusive:
- rMATS mode: `--rMATS`
- MISO mode: `--miso`
- coordinate mode: `--up --down --background`

## Required Arguments

- event subcommand: `se`, `a3ss`, `a5ss`, `ri`, or `mxe`
- `--genome`
- `--output`
- `--known-motifs`
- `--fasta-root`
- exactly one input mode

## Common Options

- `--motifs` (default `NA`)
- `--label` (default `RBP`)
- `--window` (default `50`)
- `--step` (default `1`)
- `--intron` (default `250`)
- `--exon` (default `50`)
- `--sig-fdr` / `--sigFDR` (default `0.05`)
- `--sig-delta-psi` / `--sigDeltaPSI` (default `0.05`)
- `--separate`

## Examples

### 1) SE from rMATS

```bash
python cli.py motif-map se \
  --known-motifs data/knownMotifs.human.mouse.txt \
  --motifs data/ESRP.like.motif.txt \
  --fasta-root genomedata \
  --genome hg19 \
  --output temp/se_rmats \
  --rMATS testData/SE.MATS.ReadsOnTargetAndJunctionCounts.txt \
  --miso NA --up NA --down NA --background NA
```

### 2) SE from MISO

```bash
python cli.py motif-map se \
  --known-motifs data/knownMotifs.human.mouse.txt \
  --motifs data/ESRP.like.motif.txt \
  --fasta-root genomedata \
  --genome hg19 \
  --output temp/se_miso \
  --rMATS NA \
  --miso testData/ESRP.OE.miso_bf \
  --up NA --down NA --background NA
```

### 3) A3SS from coordinate files

```bash
python cli.py motif-map a3ss \
  --known-motifs data/knownMotifs.human.mouse.txt \
  --motifs NA \
  --fasta-root genomedata \
  --genome hg19 \
  --output temp/a3ss_coords \
  --rMATS NA --miso NA \
  --up path/to/up.txt \
  --down path/to/down.txt \
  --background path/to/bg.txt
```

## Convert Commands

```bash
python cli.py convert miso --event se --in input.miso_bf --out output.rmats.txt
```

## Exon Set Commands

```bash
python cli.py exon-sets se \
  --input path/to/common.txt \
  --sample1 SAMPLE_1 \
  --sample2 SAMPLE_2 \
  --out temp/exon_sets
```

## Output Files

For each motif-map run (`--output out_dir`), typical outputs include:
- `out_dir/pVal.up.vs.bg.RNAmap.txt`
- `out_dir/pVal.dn.vs.bg.RNAmap.txt`
- `out_dir/log.motifMap.txt`
- `out_dir/exon/{up,dn,bg}.coord.txt`
- `out_dir/fasta/*.fasta`
- `out_dir/maps/*`
- `out_dir/temp/*.txt`
