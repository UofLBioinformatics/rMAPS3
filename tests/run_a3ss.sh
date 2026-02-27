#!/usr/bin/env bash
set -euo pipefail
"$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/run_event_pair.sh" a3ss mm10 ./temp/A3SS.MATS.ReadsOnTargetAndJunctionCounts.txt motif
