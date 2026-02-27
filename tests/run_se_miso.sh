#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

TS="$(date +%Y%m%d_%H%M%S)"
BASE="./temp/tests_se_miso_${TS}"
mkdir -p "$BASE"

echo "[1/4] se MISO hg19 (with optional motifs)"
python cli.py motif-map se \
  --genome hg19 \
  --fasta-root ./genomedata \
  --output "${BASE}/hg19_with_optional" \
  --rMATS NA \
  --miso ./testData/ESRP.OE.miso_bf \
  --up NA --down NA --background NA \
  --intron 250 --exon 50 --window 50 --step 1 \
  --label "hg19.miso ESRP-like" \
  --motifs ./data/ESRP.like.motif.txt \
  --known-motifs ./data/knownMotifs.human.mouse.txt

echo "[2/4] se MISO hg19 (no optional motifs)"
python cli.py motif-map se \
  --genome hg19 \
  --fasta-root ./genomedata \
  --output "${BASE}/hg19_no_optional" \
  --rMATS NA \
  --miso ./testData/ESRP.OE.miso_bf \
  --up NA --down NA --background NA \
  --intron 250 --exon 50 --window 50 --step 1 \
  --label "hg19.miso ESRP-like" \
  --motifs NA \
  --known-motifs ./data/knownMotifs.human.mouse.txt

echo "[3/4] se MISO dm3 (with optional motifs)"
python cli.py motif-map se \
  --genome dm3 \
  --fasta-root ./genomedata \
  --output "${BASE}/dm3_with_optional" \
  --rMATS NA \
  --miso ./testData/example.dm3.miso_bf \
  --up NA --down NA --background NA \
  --intron 250 --exon 50 --window 50 --step 1 \
  --label "dm3.miso ESRP-like" \
  --motifs ./data/ESRP.like.motif.txt \
  --known-motifs ./data/knownMotifs.human.mouse.txt

echo "[4/4] se MISO dm3 (no optional motifs)"
python cli.py motif-map se \
  --genome dm3 \
  --fasta-root ./genomedata \
  --output "${BASE}/dm3_no_optional" \
  --rMATS NA \
  --miso ./testData/example.dm3.miso_bf \
  --up NA --down NA --background NA \
  --intron 250 --exon 50 --window 50 --step 1 \
  --label "dm3.miso ESRP-like" \
  --motifs NA \
  --known-motifs ./data/knownMotifs.human.mouse.txt

echo "Done: ${BASE}"
