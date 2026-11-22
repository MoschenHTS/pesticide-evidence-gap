"""Reader for the global atlas of pesticide use (supplementary data file 2, Table S1)."""

import pandas as pd

from pesticide_evidence.config import path
from pesticide_evidence.names import to_canonical

COLUMNS = {
    "Country/region name": "country",
    "Chemical active ingredient name": "name",
    "CAS RN": "cas",
    "Annual applied mass [kg]": "kg",
    "Crop name": "crop",
    "Crop group": "cg",
    "Primary target class": "target",
    "Data origin": "origin",
}


def read_atlas():
    """All atlas rows with short column names."""
    return pd.read_excel(path("raw", "atlas", "sdata2.xlsx"), sheet_name="Table S1", header=1).rename(columns=COLUMNS)


def reported_atlas():
    """Rows reported by the country (not modelled), with the canonical ingredient name in `ai`."""
    x = read_atlas()
    x = x[x.origin == "Reported"].copy()
    x["ai"] = x.name.map(lambda s: to_canonical(str(s)))
    return x.dropna(subset=["ai"])
