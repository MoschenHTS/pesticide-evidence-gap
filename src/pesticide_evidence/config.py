"""Configuration loading and path resolution."""

import os
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]


def load_config():
    """Parameters, seeds and path layout from config/config.yaml."""
    with open(ROOT / "config" / "config.yaml") as fh:
        return yaml.safe_load(fh)


CFG = load_config()


DATA_ROOT = Path(os.environ.get("PESTICIDE_DATA_DIR", ROOT)).expanduser().resolve()
DATA_LAYERS = {"raw", "interim", "processed"}


def path(layer, *parts):
    """Absolute path inside a data layer (raw, interim, processed, under PESTICIDE_DATA_DIR) or results."""
    base = DATA_ROOT if layer in DATA_LAYERS else ROOT
    return base.joinpath(CFG["paths"][layer], *parts)


def treaty_ingredients():
    """Ingredients listed under the Rotterdam or Stockholm conventions."""
    import csv

    with open(ROOT / "config" / "treaty_pesticides.csv") as fh:
        return {row["ingredient"] for row in csv.DictReader(fh)}
