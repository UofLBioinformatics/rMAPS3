#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 4 ]]; then
  echo "Usage: $0 <event> <genome> <rmats_file> <label>"
  exit 1
fi

EVENT="$1"
GENOME="$2"
RMATS_FILE="$3"
LABEL="$4"

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

TS="$(date +%Y%m%d_%H%M%S)"
BASE="./temp/tests_${EVENT}_${TS}"
mkdir -p "$BASE"

echo "[1/4] ${EVENT} rMATS run (with optional motifs)"
python cli.py motif-map "$EVENT" \
  --genome "$GENOME" \
  --fasta-root ./genomedata \
  --output "${BASE}/result_${EVENT}_rMATS_with_optional" \
  --rMATS "$RMATS_FILE" \
  --miso NA \
  --up NA --down NA --background NA \
  --intron 250 --exon 50 --window 50 --step 1 \
  --label "$LABEL" \
  --motifs ./data/ESRP.like.motif.txt \
  --known-motifs ./data/testMotifs.txt

echo "[2/4] ${EVENT} user-input replay run (with optional motifs)"
python cli.py motif-map "$EVENT" \
  --genome "$GENOME" \
  --fasta-root ./genomedata \
  --output "${BASE}/result_${EVENT}_userInput_with_optional" \
  --rMATS NA \
  --miso NA \
  --up "${BASE}/result_${EVENT}_rMATS_with_optional/exon/up.coord.txt" \
  --down "${BASE}/result_${EVENT}_rMATS_with_optional/exon/dn.coord.txt" \
  --background "${BASE}/result_${EVENT}_rMATS_with_optional/exon/bg.coord.txt" \
  --intron 250 --exon 50 --window 50 --step 1 \
  --label "$LABEL" \
  --motifs ./data/ESRP.like.motif.txt \
  --known-motifs ./data/testMotifs.txt

echo "[3/4] ${EVENT} rMATS run (no optional motifs)"
python cli.py motif-map "$EVENT" \
  --genome "$GENOME" \
  --fasta-root ./genomedata \
  --output "${BASE}/result_${EVENT}_rMATS_no_optional" \
  --rMATS "$RMATS_FILE" \
  --miso NA \
  --up NA --down NA --background NA \
  --intron 250 --exon 50 --window 50 --step 1 \
  --label "$LABEL" \
  --motifs NA \
  --known-motifs ./data/testMotifs.txt

echo "[4/4] ${EVENT} user-input replay run (no optional motifs)"
python cli.py motif-map "$EVENT" \
  --genome "$GENOME" \
  --fasta-root ./genomedata \
  --output "${BASE}/result_${EVENT}_userInput_no_optional" \
  --rMATS NA \
  --miso NA \
  --up "${BASE}/result_${EVENT}_rMATS_no_optional/exon/up.coord.txt" \
  --down "${BASE}/result_${EVENT}_rMATS_no_optional/exon/dn.coord.txt" \
  --background "${BASE}/result_${EVENT}_rMATS_no_optional/exon/bg.coord.txt" \
  --intron 250 --exon 50 --window 50 --step 1 \
  --label "$LABEL" \
  --motifs NA \
  --known-motifs ./data/testMotifs.txt

echo "Done: ${BASE}"
