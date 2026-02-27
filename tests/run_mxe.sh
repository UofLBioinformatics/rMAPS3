#!/usr/bin/env bash
set -euo pipefail
"$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/run_event_pair.sh" mxe mm10 ./temp/MXE.MATS.ReadsOnTargetAndJunctionCounts.txt motif
