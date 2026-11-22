"""Regulatory status tables (EU Pesticides Database substances, treaty listings)."""

import json
import re

import pandas as pd

from .config import path
from .names import to_canonical

STATUS_PRIORITY = {"Approved": 0, "Pending": 1, "Not approved": 2}


def year_of(text):
    """Four-digit year at the end of a date text, or None."""
    return int(text[-4:]) if text and re.search(r"\d{4}$", text) else None


def load_eu(index=None):
    """One row per canonical ingredient (status, approval year, expiry year).

    When an ingredient has several EU entries the most favourable status is kept
    (Approved before Pending before Not approved). `index` optionally restricts the rows.
    """
    rows = []
    with open(path("raw", "eu", "eu_active_substances.json")) as fh:
        for line in fh:
            if not line.strip():
                continue
            rec = json.loads(line)
            ai = to_canonical(rec["substance_name"])
            if ai:
                rows.append(dict(ai=ai, status=rec["substance_status"], appr=year_of(rec.get("approval_date")), exp=year_of(rec.get("expiry_date"))))
    table = pd.DataFrame(rows).sort_values("status", key=lambda s: s.map(STATUS_PRIORITY)).drop_duplicates("ai").set_index("ai")
    return table if index is None else table[table.index.isin(index)]
