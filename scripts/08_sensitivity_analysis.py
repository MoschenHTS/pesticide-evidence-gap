"""Stage 08. Sensitivity analyses of the primary result.

Refits the primary model under alternative samples and specifications (outcome definitions, use threshold, market age, exclusions, China share),
leave-one-out over the absent group and a within-target-class label permutation.

Reads data/processed/ingredient_table.parquet; writes results/sensitivity.json.
"""

import json

import numpy as np
import pandas as pd

from pesticide_evidence.config import CFG, path
from pesticide_evidence.stats import ABSENT, analysis_sample, design, fit



def corpus_counts(sample, flag):
    """Selected-corpus paper counts per ingredient, from papers flagged by the named column."""
    corpus = pd.read_parquet(path("interim", "corpus", "core_corpus.parquet"), columns=["pmid", "compounds", flag])
    corpus = corpus[corpus[flag] & corpus.compounds.notna() & (corpus.compounds != "")]
    pairs = corpus[["pmid"]].assign(ai=corpus.compounds.str.split(";")).explode("ai")
    pairs["ai"] = pairs.ai.str.strip().str.lower().map({a.lower(): a for a in sample.index})
    return sample.index.map(pairs.dropna(subset=["ai"]).groupby("ai").pmid.nunique()).fillna(0).astype(float).values


def sensitivity(sample, table):
    """Primary-model ratios under alternative samples, outcome definitions, thresholds, exclusions and China share, plus leave-one-out and label permutation."""
    atlas = pd.read_csv(path("processed", "use_atlas_country_ingredient.csv"))
    atlas = atlas[atlas.origin == "Reported"]

    def prepare(s):
        """Sample copy with each ingredient's China share of Southern use."""
        s = s.copy()
        s["cn"] = s.index.map(atlas[atlas.country == "China"].groupby("ai").tonnes.sum()).fillna(0).values / s.South.values
        return s

    def run(s, **kw):
        """Fit the primary model on a sample and return its size and ratios."""
        out, _, _ = fit(s, "en_zh_tox", **kw)
        return {"n": len(s), **{k: list(v) for k, v in out.items()}}

    s = prepare(sample)
    res = {"baseline": run(s)}
    res["PubMed toxicity only"] = {"n": len(s), **{k: list(v) for k, v in fit(s, "pm_tox")[0].items()}}
    res["all topics"] = {"n": len(s), **{k: list(v) for k, v in fit(s, "en_zh_all")[0].items()}}
    res[">=200 t Southern use"] = run(prepare(analysis_sample(table, 200)))
    res["excluding reinstated (malathion, fenvalerate, chlorpyrifos)"] = run(s[~s.index.isin(["malathion", "fenvalerate", "chlorpyrifos"])])
    res["dropping five largest Southern uses"] = run(s.drop(s.sort_values("South", ascending=False).index[:5]))
    res["excluding fumigants"] = run(s[~s.index.isin(["metam", "1,3-dichloropropene", "chloropicrin", "dazomet"])])
    res["first PubMed record 1985 or later (age<=40)"] = run(s[s.age <= 40])
    for label, flag in (("first-filter corpus counts", "core"), ("selected corpus counts", "selected")):
        counted = s.assign(corpus_papers=corpus_counts(s, flag))
        res[label] = {"n": len(counted), **{k: list(v) for k, v in fit(counted, "corpus_papers")[0].items()}}
    res["market age of at least 20 years"] = run(s[s.age >= 20])
    res["market age of at least 30 years"] = run(s[s.age >= 30])
    res["+ China share of Southern use"] = run(s, extra=["cn"])
    drops = [fit(s.drop(a), "en_zh_tox")[0]["absent"][0] for a in s.index[s.grp4 == ABSENT]]
    res["leave-one-out (absent), range of ratio"] = {"min": min(drops), "max": max(drops)}
    rng = np.random.default_rng(CFG["seeds"]["within_class_permutation"])
    observed = res["baseline"]["absent"][0]
    permuted = []
    for _ in range(CFG["within_class_permutations"]):
        p = s.copy()
        p["grp4"] = p.groupby("target").grp4.transform(lambda x: rng.permutation(x.values))
        permuted.append(fit(p, "en_zh_tox")[0]["absent"][0])
    res["permutation within target class (400)"] = {"reached_observed": int(sum(v <= observed for v in permuted)), "median_permuted": float(np.median(permuted))}
    return {"sensitivity": res}


def model_specification(sample):
    """Coefficients and variance inflation factors of the primary model, and the model with the squared ln Southern use."""
    s = sample.assign(ln_south_sq=np.log(sample.South.values) ** 2)
    out, coef, se = fit(s, "en_zh_tox")
    X, names = design(s)
    x = X[:, 1:]
    vif = {}
    for j, name in enumerate(names[1:]):
        others = np.delete(x, j, axis=1)
        beta = np.linalg.lstsq(np.column_stack([np.ones(len(x)), others]), x[:, j], rcond=None)[0]
        resid = x[:, j] - np.column_stack([np.ones(len(x)), others]) @ beta
        vif[name] = float(1 / max(1e-12, resid.var() / x[:, j].var()))
    coefficients = {
        k: {
            "estimate": float(coef[k]),
            "se": float(se[k]),
            "ratio": float(np.exp(coef[k])),
            "lo": float(np.exp(coef[k] - 1.96 * se[k])),
            "hi": float(np.exp(coef[k] + 1.96 * se[k])),
        }
        for k in ("ln_south", "ln_north", "ln_age", "withdrawn", "lapsed", "absent")
    }
    forms = {}
    for label, kwargs, data in (
        ("primary specification", {}, s),
        ("without crop-group shares", {}, s[[c for c in s.columns if not c.startswith("crop_")]]),
        ("without Northern use", {"drop": ("ln_north",)}, s),
        ("without market age", {"drop": ("ln_age",)}, s),
        ("without target class", {"drop": tuple(n for n in names if n.startswith("target_"))}, s),
        ("with squared ln Southern use", {"extra": ("ln_south_sq",)}, s),
    ):
        res = fit(data, "en_zh_tox", **kwargs)[0]
        forms[label] = {k: list(v) for k, v in res.items()}
    return {
        "model_specification": {
            "n": len(s),
            "coefficients": coefficients,
            "vif_max": {"name": max(vif, key=vif.get), "value": max(vif.values())},
            "vif_median": float(np.median(list(vif.values()))),
            "n_covariates": len(names),
            "alternative_forms": forms,
        }
    }


if __name__ == "__main__":
    table = pd.read_parquet(path("processed", "ingredient_table.parquet"))
    sample = analysis_sample(table, CFG["sample"]["min_south_use_t"])
    out = sensitivity(sample, table)
    out.update(model_specification(sample))
    json.dump(out, open(path("results", "sensitivity.json"), "w"), indent=1, default=str)
    print("wrote results/sensitivity.json")
