#!/usr/bin/env bash
set -euo pipefail
"$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/run_event_pair.sh" a5ss hg38 ./temp/A5SS.MATS.ReadsOnTargetAndJunctionCounts.txt motif
