#!/usr/bin/env bash
set -euo pipefail
"$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/run_event_pair.sh" se hg38 ./temp/SE.MATS.ReadsOnTargetAndJunctionCounts.txt motif
