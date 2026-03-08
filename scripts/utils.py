"""
Shared utilities for RNA-ATAC-TF-Network pipeline scripts.
"""

from pathlib import Path

import yaml


_REQUIRED_CONFIG_KEYS = {'data', 'contrasts', 'scoring', 'output'}


def load_config(config_path: str) -> dict:
    """Load and validate YAML configuration file.

    Raises:
        ValueError: If required top-level keys are missing from the config.
    """
    with open(config_path, 'r') as f:
        cfg = yaml.safe_load(f)

    if cfg is None:
        raise ValueError(f"Config file is empty: {config_path}")

    missing = _REQUIRED_CONFIG_KEYS - set(cfg.keys())
    if missing:
        raise ValueError(
            f"Config file {config_path} missing required keys: {sorted(missing)}. "
            f"Required keys: {sorted(_REQUIRED_CONFIG_KEYS)}"
        )

    return cfg


def resolve_path(path: str, base_dir: Path) -> Path:
    """Resolve relative paths from config file location."""
    p = Path(path)
    if p.is_absolute():
        return p
    return (base_dir / p).resolve()
