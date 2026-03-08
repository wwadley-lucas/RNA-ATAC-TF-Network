"""
Shared utilities for RNA-ATAC-TF-Network pipeline scripts.
"""

from pathlib import Path

import yaml


def load_config(config_path: str) -> dict:
    """Load YAML configuration file."""
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def resolve_path(path: str, base_dir: Path) -> Path:
    """Resolve relative paths from config file location."""
    p = Path(path)
    if p.is_absolute():
        return p
    return (base_dir / p).resolve()
