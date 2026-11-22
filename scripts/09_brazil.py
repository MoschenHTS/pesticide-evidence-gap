"""Stage 09. Brazilian subsample: the gap model with Brazilian applied mass in place of Southern use, and the five best-selling pesticides.

Reads  data/processed/ingredient_table.parquet and data/processed/use_atlas_country_ingredient.csv.
Writes results/brazil.json (consumed by 11_write_ledger.py) and results/tables/Table_2.csv.
"""

import json

import numpy as np
import pandas as pd

from pesticide_evidence.config import CFG, path
from pesticide_evidence.stats import AUTHORIZED, GROUP_KEYS, GROUPS, analysis_sample, design, ppml, round5

# the five best-selling pesticides in Brazil named as not approved in the EU by Bombardi (2023)
BEST_SELLERS = ["mancozeb", "atrazine", "acephate", "chlorothalonil", "chlorpyrifos"]


def brazil_sample(min_t):
    """Analysis sample with Brazilian applied mass above a threshold and indicators for absence from each register."""
    table = pd.read_parquet(path("processed", "ingredient_table.parquet"))
    sample = analysis_sample(table, min_south_t=CFG["sample"]["min_south_use_t"])
    use = pd.read_csv(path("processed", "use_atlas_country_ingredient.csv"))
    brazil = use[(use.country == "Brazil") & (use.origin != "Imputed")].groupby("ai").tonnes.sum()
    sample["brazil_t"] = sample.index.map(brazil).astype(float).fillna(0)
    sample["not_EU"] = (~sample.eu_ok.astype(bool)).astype(float)
    sample["not_US"] = (~sample.us_ok.astype(bool)).astype(float)
    sample["not_CA"] = (~sample.ca_ok.astype(bool)).astype(float)
    return sample, sample[sample.brazil_t >= min_t].copy()


def fit_brazil(data, terms):
    """Gap model with Brazilian applied mass in place of Southern use; returns ratios with intervals and the use elasticity."""
    d = data.copy()
    d["South"] = d["brazil_t"]
    X, names = design(d, terms=terms)
    b, s = ppml(d.en_zh_tox.values.astype(float), X, np.arange(len(d)))
    out = {}
    for k in terms:
        i = names.index(k)
        out[k] = [round5(np.exp(b[i])), round5(np.exp(b[i] - 1.96 * s[i])), round5(np.exp(b[i] + 1.96 * s[i]))]
    return out, round5(b[names.index("ln_south")])


def run(data):
    """Group ratios, single-register ratios and the joint register model on the Brazilian sample."""
    d = data.copy()
    for key, group in zip(GROUP_KEYS, GROUPS):
        d[key] = (d.grp4 == group).astype(float)
    groups, elasticity = fit_brazil(d, GROUP_KEYS)
    registers = {k: fit_brazil(d, [k])[0][k] for k in ("not_EU", "not_US", "not_CA")}
    registers["joint"] = fit_brazil(d, ["not_EU", "not_US", "not_CA"])[0]
    return {
        "n": len(d),
        "by_group": d.grp4.value_counts().to_dict(),
        "groups": groups,
        "elasticity_brazil_use": elasticity,
        "registers": registers,
        "not_eu_n": int(d.not_EU.sum()),
        "not_eu_us_registered_n": int(((d.not_EU == 1) & d.us_ok.astype(bool)).sum()),
    }


def best_sellers(sample):
    """Observed and expected papers for the best-selling pesticides, with the expectation from the model fitted to authorized pesticides."""
    X, _ = design(sample, terms=[])
    authorized = (sample.grp4 == AUTHORIZED).values
    b, _ = ppml(sample.en_zh_tox.values.astype(float)[authorized], X[authorized], np.arange(authorized.sum()))
    expected = pd.Series(np.exp(X @ b), index=sample.index)
    return {
        a: {
            "eu": bool(sample.eu_ok[a]),
            "us": bool(sample.us_ok[a]),
            "ca": bool(sample.ca_ok[a]),
            "brazil_t": round5(sample.brazil_t[a]),
            "observed": float(sample.en_zh_tox[a]),
            "expected_authorized_fit": round5(expected[a]),
            "ratio": round5(sample.en_zh_tox[a] / expected[a]),
        }
        for a in BEST_SELLERS
    }


def table_2(rows):
    """Table 2 of the paper: register status and observed against expected papers for the best sellers."""
    yes = {True: "Yes", False: "No"}
    return pd.DataFrame(
        [
            {
                "Pesticide": a.capitalize(),
                "EU": yes[v["eu"]],
                "US": yes[v["us"]],
                "Canada": yes[v["ca"]],
                "Brazilian use (t)": round(v["brazil_t"]),
                "Papers": int(v["observed"]),
                "Expected papers": round(v["expected_authorized_fit"]),
                "Ratio": round(v["ratio"], 2),
            }
            for a, v in rows.items()
        ]
    )


if __name__ == "__main__":
    min_t = CFG["sample"]["min_south_use_t"]
    sample, main = brazil_sample(min_t)
    out = {"min_brazil_use_t": min_t, "primary": run(main), "best_sellers": best_sellers(sample)}
    out["sensitivity"] = {}
    for label, threshold in (("at_least_10_t", 10), ("any_reported_use", 0.001)):
        out["sensitivity"][label] = run(sample[sample.brazil_t >= threshold])
    out["sensitivity"]["at_least_10_t"]["threshold_t"] = 10
    out["sensitivity"]["any_reported_use"]["threshold_t"] = 0
    json.dump(out, open(path("results", "brazil.json"), "w"), indent=1)
    tables = path("results", "tables")
    tables.mkdir(parents=True, exist_ok=True)
    frame = table_2(out["best_sellers"])
    frame.to_csv(tables / "Table_2.csv", index=False)
    (tables / "Table_2.md").write_text(frame.to_markdown(index=False) + "\n")
    print(json.dumps(out["primary"], indent=1))
