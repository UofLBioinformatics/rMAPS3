# Quickstart

## 1) Run an SE analysis from rMATS

```bash
python cli.py motif-map se \
  --known-motifs data/knownMotifs.human.mouse.txt \
  --motifs data/ESRP.like.motif.txt \
  --fasta-root genomedata \
  --genome hg19 \
  --output temp/quickstart_se_rmats \
  --rMATS testData/SE.MATS.ReadsOnTargetAndJunctionCounts.txt \
  --miso NA --up NA --down NA --background NA
```

## 2) Run an SE analysis from MISO

```bash
python cli.py motif-map se \
  --known-motifs data/knownMotifs.human.mouse.txt \
  --motifs data/ESRP.like.motif.txt \
  --fasta-root genomedata \
  --genome hg19 \
  --output temp/quickstart_se_miso \
  --rMATS NA \
  --miso testData/ESRP.OE.miso_bf \
  --up NA --down NA --background NA
```

Note: MISO mode requires Perl in `PATH`.

## 3) Start local web UI

```bash
python run_web.py
```

Open:
- `http://127.0.0.1:5000`

## 4) One-click web test

In the web page:
- choose event/genome
- click `Run One-Click Test`

Default quick-test data root is `testData/` (or `RMAPS_QUICKTEST_DIR`).

## 5) Check outputs

Inspect the output directory shown by CLI or web UI.

Typical files:
- `log.motifMap.txt`
- `pVal.up.vs.bg.RNAmap.txt`
- `pVal.dn.vs.bg.RNAmap.txt`
- `maps/`
- `fasta/`
- `exon/`

## More

- Installation: `INSTALL.md`
- Full CLI reference: `CLI_USAGE.md`
- Web UI details: `webui/README.md`
