from __future__ import annotations

import os
from pathlib import Path


SENTINEL_VALUES = {"NA"}


def repo_root() -> Path:
    """Return the rMAPS3 repository/install root."""
    return Path(__file__).resolve().parents[1]


def resolve_user_path(path: str | Path, base_cwd: Path) -> str:
    """Resolve a user-supplied path relative to the original invocation cwd."""
    if isinstance(path, str) and path in SENTINEL_VALUES:
        return path

    user_path = Path(path).expanduser()
    if user_path.is_absolute():
        return str(user_path)
    return str((base_cwd / user_path).resolve())


def build_subprocess_env(
    root: Path,
    env_overrides: dict[str, str] | None = None,
) -> dict[str, str]:
    """Build a subprocess environment that can import rMAPS modules."""
    env = dict(os.environ)
    if env_overrides:
        env.update(env_overrides)

    root_str = str(root)
    existing = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = root_str if not existing else f"{root_str}{os.pathsep}{existing}"
    return env
