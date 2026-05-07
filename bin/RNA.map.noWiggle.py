#!/usr/bin/env python3
"""Compatibility dispatcher for event-specific CLIP RNA-map renderers."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


EVENT_RENDERERS = {
    "SE": "RNA.map.noWiggle.SE.py",
    "A3SS": "RNA.map.noWiggle.A3SS.py",
    "A5SS": "RNA.map.noWiggle.A5SS.py",
    "RI": "RNA.map.noWiggle.RI.py",
    "MXE": "RNA.map.noWiggle.MXE.py",
}


def main() -> int:
    event_type = sys.argv[14].upper() if len(sys.argv) > 14 else "SE"
    renderer_name = EVENT_RENDERERS.get(event_type)
    if renderer_name is None:
        supported = ", ".join(sorted(EVENT_RENDERERS))
        print(f"Unsupported event type for RNA map renderer: {event_type}")
        print(f"Supported event types: {supported}")
        return 2

    renderer = Path(__file__).resolve().with_name(renderer_name)
    cmd = [sys.executable, str(renderer), *sys.argv[1:]]
    return subprocess.run(cmd).returncode


if __name__ == "__main__":
    raise SystemExit(main())
