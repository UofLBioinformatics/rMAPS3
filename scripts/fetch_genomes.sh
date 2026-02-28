#!/usr/bin/env bash
set -euo pipefail

MANIFEST="scripts/genomes.manifest.tsv"
FASTA_ROOT="${RMAPS_FASTA_ROOT:-genomedata}"
GENOMES="all"
FORCE=0
KEEP_ARCHIVE=0
ALLOW_UNVERIFIED=0
DRY_RUN=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest) MANIFEST="$2"; shift 2 ;;
    --fasta-root) FASTA_ROOT="$2"; shift 2 ;;
    --genomes) GENOMES="$2"; shift 2 ;;
    --force) FORCE=1; shift ;;
    --keep-archive) KEEP_ARCHIVE=1; shift ;;
    --allow-unverified) ALLOW_UNVERIFIED=1; shift ;;
    --dry-run) DRY_RUN=1; shift ;;
    *) echo "Unknown arg: $1" >&2; exit 1 ;;
  esac
done

if [[ ! -f "$MANIFEST" ]]; then
  echo "Manifest not found: $MANIFEST" >&2
  exit 1
fi

mkdir -p "$FASTA_ROOT"

contains_genome() {
  local list="$1"
  local item="$2"
  [[ ",$list," == *",$item,"* ]]
}

verify_sha256() {
  local file="$1"
  local expected="$2"
  local got=""
  if command -v sha256sum >/dev/null 2>&1; then
    got="$(sha256sum "$file" | awk '{print $1}')"
  elif command -v shasum >/dev/null 2>&1; then
    got="$(shasum -a 256 "$file" | awk '{print $1}')"
  else
    echo "No sha256 tool found (sha256sum/shasum)." >&2
    exit 1
  fi
  [[ "${got,,}" == "${expected,,}" ]]
}

while IFS=$'\t' read -r genome url sha256 _; do
  [[ -z "${genome:-}" ]] && continue
  [[ "${genome:0:1}" == "#" ]] && continue

  if [[ "$GENOMES" != "all" ]] && ! contains_genome "$GENOMES" "$genome"; then
    continue
  fi

  if [[ -z "${url:-}" ]]; then
    echo "Missing URL for $genome in manifest" >&2
    exit 1
  fi
  if [[ -z "${sha256:-}" && "$ALLOW_UNVERIFIED" -ne 1 ]]; then
    echo "Missing sha256 for $genome. Add checksum or use --allow-unverified." >&2
    exit 1
  fi

  target_dir="$FASTA_ROOT/$genome"
  mkdir -p "$target_dir"
  leaf="$(basename "$url")"
  archive="$target_dir/$leaf"
  final_fa="$target_dir/$genome.fa"

  echo "==> $genome"
  echo "    URL: $url"
  echo "    Target: $final_fa"
  if [[ "$DRY_RUN" -eq 1 ]]; then
    continue
  fi

  if [[ -f "$final_fa" && "$FORCE" -ne 1 ]]; then
    echo "    Skip: FASTA exists"
  else
    if [[ "$FORCE" -eq 1 || ! -f "$archive" ]]; then
      echo "    Downloading..."
      if command -v curl >/dev/null 2>&1; then
        curl -L --retry 3 -o "$archive" "$url"
      elif command -v wget >/dev/null 2>&1; then
        wget -O "$archive" "$url"
      else
        echo "Need curl or wget" >&2
        exit 1
      fi
    else
      echo "    Reusing existing archive"
    fi

    if [[ -n "${sha256:-}" ]]; then
      if ! verify_sha256 "$archive" "$sha256"; then
        echo "SHA256 mismatch for $genome" >&2
        exit 1
      fi
      echo "    Checksum OK"
    else
      echo "    Warning: no checksum for $genome (allowed by --allow-unverified)"
    fi

    if [[ "$archive" == *.gz ]]; then
      echo "    Decompressing..."
      gzip -dc "$archive" > "$final_fa"
    else
      cp -f "$archive" "$final_fa"
    fi
  fi

  if [[ ! -f "$final_fa.fai" ]]; then
    if command -v python >/dev/null 2>&1; then
      echo "    Indexing with pyfaidx..."
      python -c "from pyfaidx import Fasta; Fasta(r'$final_fa'); print('indexed')"
    else
      echo "    Warning: python not found; skipping index"
    fi
  fi

  if [[ "$KEEP_ARCHIVE" -ne 1 && "$archive" != "$final_fa" ]]; then
    rm -f "$archive"
  fi
done < "$MANIFEST"

echo "Done."
