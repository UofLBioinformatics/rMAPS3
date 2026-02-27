#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
"$ROOT/run_se.sh"
"$ROOT/run_a3ss.sh"
"$ROOT/run_a5ss.sh"
"$ROOT/run_ri.sh"
"$ROOT/run_mxe.sh"
"$ROOT/run_se_miso.sh"
echo "All event tests completed"
