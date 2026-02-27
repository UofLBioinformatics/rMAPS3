#!/usr/bin/env bash
set -euo pipefail
"$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/run_event_pair.sh" ri hg38 ./temp/RI.MATS.ReadsOnTargetAndJunctionCounts.txt motif
