"""Stage 07. Event studies of research after Northern market exit.

Design 1, stacked matched cohorts: an ingredient's first Northern exit year (EU non-renewal, end of all US registrations,
end of all Canadian registrations) in 2010-2020; three nearest controls (EU-approved, no US or Canadian exit) matched on pre-period
paper volume and year of first PubMed record; window -5..+5 years; Poisson pseudo-ML with unit-by-cohort and cohort-by-year effects,
errors clustered by ingredient; separate outcomes for Northern- and Southern-based authors.
Calibration: 200 placebo runs (random never-exited ingredients, random exit years).

Design 2, EU review timing: ingredients not renewed at their scheduled expiry (treated) against old ingredients renewed within one year
of the same date (controls); calibration by 300 label permutations.

Writes results/event_studies.json.
"""

import json

import numpy as np
import pandas as pd

from pesticide_evidence import geography
from pesticide_evidence.config import CFG, path, treaty_ingredients
from pesticide_evidence.names import to_canonical
from pesticide_evidence.registers import load_eu, year_of
from pesticide_evidence.stats import pct_change, ppml_cov

ES = CFG["event_study"]
YEARS = list(range(CFG["sample"]["years"][0], CFG["sample"]["years"][1] + 1))
TREATY = treaty_ingredients()


class Inputs:
    """Panels, registers and exit events shared by the event-study designs."""

    def __init__(self):
        """Load the paper panels, registers, use table and author-country panels, and derive each ingredient's first Northern exit and the control pool."""
        self.panel = pd.read_parquet(path("interim", "corpus", "panel_ingredient_year.parquet"))
        self.panel.columns = [int(c) for c in self.panel.columns]
        self.first_record = json.load(open(path("interim", "corpus", "first_record.json")))
        self.eu = load_eu(self.panel.index)
        us = pd.read_csv(path("processed", "us_registrations.csv"), index_col=0)
        ca = pd.read_csv(path("processed", "canada_registrations.csv"), index_col=0)
        self.south_use = pd.read_parquet(path("processed", "ingredient_table.parquet")).South
        self.events = self._exits(us, ca)
        self.first_exit = {a: min(d.values()) for a, d in self.events.items()}

        def ok_us(a):
            """True when the ingredient has a current US product or no US record."""
            return a not in us.index or us.loc[a, "n_current"] > 0

        def ok_ca(a):
            """True when the ingredient has a current Canadian product or no Canadian record."""
            return a not in ca.index or ca.loc[a, "n_current"] > 0

        approved = self.eu[self.eu.status == "Approved"].index
        self.controls = [a for a in approved if a not in self.events and ok_us(a) and ok_ca(a) and a in self.panel.index]
        self.author_panels = self._author_panels()

    def _exits(self, us, ca):
        """Exit year per register for ingredients that left the EU, US or Canadian register between 2008 and 2020."""
        lo, hi = 2008, 2020
        events = {}
        for a, r in self.eu[(self.eu.status == "Not approved") & self.eu.appr.notna() & self.eu.exp.between(lo, hi)].iterrows():
            events.setdefault(a, {})["EU"] = int(r.exp)
        for a, r in us[(us.n_current == 0) & us.last_cancel.between(lo, hi)].iterrows():
            events.setdefault(a, {})["US"] = int(r.last_cancel)
        for a, r in ca[(ca.n_current == 0) & ca.last_exp.between(lo, hi)].iterrows():
            events.setdefault(a, {})["CA"] = int(r.last_exp)
        return events

    def _author_panels(self):
        """Papers per ingredient and year for Northern-authored and Southern-authored papers, by the country of the main affiliation."""
        corpus = pd.read_parquet(path("interim", "corpus", "core_corpus.parquet"))
        aff = pd.read_parquet(path("interim", "corpus", "affiliations_countries.parquet"))
        d = corpus[corpus.selected][["pmid", "year", "compounds"]].merge(aff[["pmid", "c_main"]], on="pmid")
        d["grp"] = np.where(d.c_main.isna(), "unk", np.where(d.c_main.isin(geography.NORTH_AUTHORS), "North", "South"))

        def panel(mask):
            """Papers per ingredient and year for the rows selected by the mask."""
            x = d[mask]
            e = x.assign(_a=x.compounds.replace("", np.nan).str.split(";")).explode("_a").dropna(subset=["_a"])
            return e.groupby(["_a", "year"]).pmid.nunique().unstack(fill_value=0).reindex(columns=YEARS, fill_value=0)

        return {"all": self.panel, "North": panel(d.grp == "North"), "South": panel(d.grp == "South")}

    def exits_in_window(self):
        """First exit year per ingredient, for exits inside the event window set in the configuration."""
        lo, hi = ES["exit_years"]
        return {a: c for a, c in self.first_exit.items() if lo <= c <= hi}


def pre_papers(panel, a, c):
    """Papers on ingredient a in the five years before year c."""
    return panel.loc[a, [y for y in YEARS if c - 5 <= y <= c - 1]].sum()


def build_stack(inp, cohorts, pool, outcome, window=ES["window"], n_match=ES["n_matched_controls"], min_pre=ES["min_pre_papers"]):
    """Stacked cohort panel: treated ingredients and their nearest controls matched on pre-exit paper volume and year of first record."""
    panel, first = inp.panel, inp.first_record
    rows = []
    for c in sorted(set(cohorts.values())):
        treated = [a for a, y in cohorts.items() if y == c and a in panel.index]

        treated = [a for a in treated if pre_papers(panel, a, c) >= min_pre]
        candidates = [a for a in pool if pre_papers(panel, a, c) >= min_pre and a in first]
        if not treated or not candidates:
            continue
        c_pre = np.log1p(np.array([pre_papers(panel, a, c) / 5 for a in candidates]))
        c_first = np.array([first[a] for a in candidates], float)
        sd_pre, sd_first = c_pre.std() or 1, c_first.std() or 1
        chosen = []
        for a in treated:
            dist = np.abs(c_pre - np.log1p(pre_papers(panel, a, c) / 5)) / sd_pre + np.abs(c_first - first.get(a, 1980)) / sd_first
            chosen += list(np.array(candidates)[np.argsort(dist)[:n_match]])
        for a in treated + chosen:
            for y in range(c - window, c + window + 1):
                if YEARS[0] <= y <= YEARS[-1]:
                    rows.append((a, c, y, y - c, a in treated, float(outcome.loc[a, y]) if a in outcome.index else 0.0))
    L = pd.DataFrame(rows, columns=["ai", "coh", "year", "k", "treated", "n"])
    L["treated"] = L.treated.astype(bool)
    L["unit"] = L.ai + "|" + L.coh.astype(str) + "|" + L.treated.astype(int).astype(str)
    L["cy"] = L.coh.astype(str) + "|" + L.year.astype(str)
    return L


def joint_effect(LN, LS):
    """North and South post effects in one model (clustered by ingredient) and their difference."""
    L = pd.concat([LN.assign(g="N"), LS.assign(g="S")], ignore_index=True)
    U = pd.get_dummies(L.unit + "|" + L.g, drop_first=True).astype(float).values
    CY = pd.get_dummies(L.cy + "|" + L.g, drop_first=True).astype(float).values
    pn = ((L.treated) & (L.k >= 1) & (L.g == "N")).astype(float).values
    ps = ((L.treated) & (L.k >= 1) & (L.g == "S")).astype(float).values
    b, V = ppml_cov(L.n.values, np.column_stack([np.ones(len(L)), pn, ps, U, CY]), L.ai.values, iters=60)
    d = b[1] - b[2]
    return dict(bN=b[1], seN=np.sqrt(V[1, 1]), bS=b[2], seS=np.sqrt(V[2, 2]), d=d, sed=np.sqrt(max(V[1, 1] + V[2, 2] - 2 * V[1, 2], 0)))


def stacked_calibration(inp):
    """Joint North and South effect of exits for two ingredient sets, with placebo and permutation calibration."""
    exits = inp.exits_in_window()

    def used(a):
        """True when the ingredient is still used in the South above the configured threshold."""
        return inp.south_use.get(a, 0) >= ES["still_used_t"]

    sets = {
        "any exit, still used >=200 t": {a: c for a, c in exits.items() if used(a)},
        "not treaty-listed, any use": {a: c for a, c in exits.items() if a not in TREATY},
    }
    eligible = [a for a in inp.controls if a in inp.first_record]
    years = np.array(list(exits.values()))
    rng = np.random.default_rng(CFG["seeds"]["stacked_placebo"])
    out = {}
    for label, selection in sets.items():
        LN = build_stack(inp, selection, inp.controls, inp.author_panels["North"])
        LS = build_stack(inp, selection, inp.controls, inp.author_panels["South"])
        r = joint_effect(LN, LS)
        rec = dict(
            set=label,
            events=len(selection),
            nT=int(LN[LN.treated].ai.nunique()),
            North=pct_change(r["bN"]),
            North_ci=[pct_change(r["bN"] - 1.96 * r["seN"]), pct_change(r["bN"] + 1.96 * r["seN"])],
            South=pct_change(r["bS"]),
            South_ci=[pct_change(r["bS"] - 1.96 * r["seS"]), pct_change(r["bS"] + 1.96 * r["seS"])],
            diff_log=r["d"],
            diff_se=r["sed"],
            diff_z=r["d"] / r["sed"],
        )
        placebo = []
        for _ in range(ES["placebo_runs"]):
            pseudo = rng.choice(eligible, size=len(selection), replace=False)
            cohorts = {a: int(rng.choice(years)) for a in pseudo}
            pool = [a for a in eligible if a not in cohorts]
            try:
                pn = build_stack(inp, cohorts, pool, inp.author_panels["North"])
                ps = build_stack(inp, cohorts, pool, inp.author_panels["South"])
                if pn[pn.treated].ai.nunique() < 8:
                    continue
                q = joint_effect(pn, ps)
                placebo.append((q["bN"], q["bS"], q["d"]))
            except Exception:
                pass
        placebo = np.array(placebo)
        rec["placebo_n"] = len(placebo)
        for j, name in enumerate(["North", "South", "diff"]):
            x = placebo[:, j]
            rec["placebo_" + name] = dict(mean=float(x.mean()), sd=float(x.std()), p2_5=float(np.percentile(x, 2.5)), p97_5=float(np.percentile(x, 97.5)))
        rec["p_North_le"] = float(np.mean(placebo[:, 0] <= np.log1p(rec["North"] / 100)))
        rec["p_South_le"] = float(np.mean(placebo[:, 1] <= np.log1p(rec["South"] / 100)))
        rec["p_diff_le"] = float(np.mean(placebo[:, 2] <= rec["diff_log"]))
        pN, pS = rec["placebo_North"], rec["placebo_South"]
        rec["North_cal_ci"] = [pct_change(np.log1p(rec["North"] / 100) - pN["p97_5"]), pct_change(np.log1p(rec["North"] / 100) - pN["p2_5"])]
        rec["South_cal_ci"] = [pct_change(np.log1p(rec["South"] / 100) - pS["p97_5"]), pct_change(np.log1p(rec["South"] / 100) - pS["p2_5"])]
        out[label] = rec
    return out


def balance(inp):
    """Means and standardized differences of treated and control units in the stacked design."""
    table = pd.read_parquet(path("processed", "ingredient_table.parquet"))
    exits = inp.exits_in_window()
    selection = {a: c for a, c in exits.items() if a not in TREATY}
    stack = build_stack(inp, selection, inp.controls, inp.author_panels["North"])
    units = stack.drop_duplicates(["ai", "coh", "treated"])[["ai", "coh", "treated"]]

    def papers(a, c, lo, hi):
        """Papers on ingredient a from lo to hi years after the exit year c."""
        return float(inp.panel.loc[a, [y for y in YEARS if c + lo <= y <= c + hi]].sum())

    rows = []
    for u in units.itertuples():
        t = table.loc[u.ai] if u.ai in table.index else None
        rows.append(
            {
                "treated": bool(u.treated),
                "papers_per_year_before": papers(u.ai, u.coh, -5, -1) / 5,
                "first_record_year": float(inp.first_record[u.ai]),
                "change_in_papers_before": np.log1p(papers(u.ai, u.coh, -2, -1) / 2) - np.log1p(papers(u.ai, u.coh, -5, -4) / 2),
                "ln_south_use": float(np.log(t.South)) if t is not None and t.South > 0 else np.nan,
                "ln_north_use": float(np.log(t.North + 1)) if t is not None else np.nan,
                "insecticide": float(t is not None and t.target == "Insecticides"),
                "herbicide": float(t is not None and t.target == "Herbicides"),
                "fungicide": float(t is not None and t.target == "Fungicides"),
            }
        )
    d = pd.DataFrame(rows)
    out = {}
    for col in d.columns.drop("treated"):
        x, y = d.loc[d.treated, col].dropna(), d.loc[~d.treated, col].dropna()
        pooled = np.sqrt((x.var() + y.var()) / 2)
        out[col] = {"treated": float(x.mean()), "control": float(y.mean()), "std_diff": float((x.mean() - y.mean()) / pooled) if pooled > 0 else 0.0}
    out["n_treated_units"] = int(d.treated.sum())
    out["n_control_units"] = int((~d.treated).sum())
    return out


EVENT_BINS = {"-5 to -4": (-5, -4), "-3 to -2": (-3, -2), "0": (0, 0), "1 to 2": (1, 2), "3 to 5": (3, 5)}


def event_time(inp):
    """Event-time effects by period around the exit, for Northern-authored and Southern-authored papers."""
    exits = inp.exits_in_window()
    selection = {a: c for a, c in exits.items() if a not in TREATY}
    out = {}
    for who in ("North", "South"):
        L = build_stack(inp, selection, inp.controls, inp.author_panels[who])
        U = pd.get_dummies(L.unit, drop_first=True).astype(float).values
        CY = pd.get_dummies(L.cy, drop_first=True).astype(float).values
        bins = np.column_stack([((L.treated) & L.k.between(lo, hi)).astype(float).values for lo, hi in EVENT_BINS.values()])
        b, V = ppml_cov(L.n.values, np.column_stack([np.ones(len(L)), bins, U, CY]), L.ai.values, iters=60)
        se = np.sqrt(np.diag(V)[1 : 1 + len(EVENT_BINS)])
        out[who] = {name: dict(pct=pct_change(b[1 + i]), lo=pct_change(b[1 + i] - 1.96 * se[i]), hi=pct_change(b[1 + i] + 1.96 * se[i])) for i, name in enumerate(EVENT_BINS)}
    return out


def review_cohorts(inp):
    """Non-renewed ingredients (treated) and old ingredients renewed at their scheduled date (controls), keyed to the event year."""
    rows = []
    with open(path("raw", "eu", "eu_active_substances.json")) as fh:
        for line in fh:
            if line.strip():
                rec = json.loads(line)

                rows.append(
                    dict(ai=to_canonical(rec["substance_name"]), status=rec["substance_status"], appr=year_of(rec.get("approval_date")), exp=year_of(rec.get("expiry_date")))
                )
    raw = pd.DataFrame(rows).dropna(subset=["ai"])
    renewed = raw[(raw.status == "Approved") & raw.appr.between(2011, 2021)].drop_duplicates("ai").set_index("ai")
    not_renewed = raw[(raw.status == "Not approved") & raw.exp.between(2011, 2021) & raw.appr.notna()].drop_duplicates("ai").set_index("ai")

    def is_old(a):
        """True when the first PubMed record is from 1998 or earlier."""
        return inp.first_record.get(a, 2000) <= 1998

    treated = {a: int(r.exp) for a, r in not_renewed.iterrows() if a in inp.panel.index}
    controls = {a: int(r.appr) for a, r in renewed.iterrows() if a in inp.panel.index and is_old(a)}
    return treated, controls


def build_review(inp, treated, controls, outcome, min_pre=3, window=5):
    """Cohort panel for the review-timing design, matched on pre-period paper volume."""
    panel = inp.panel

    rows = []
    for c in sorted(set(treated.values())):
        T = [a for a, y in treated.items() if y == c and pre_papers(panel, a, c) >= min_pre]
        C = [a for a, y in controls.items() if abs(y - c) <= 1 and pre_papers(panel, a, c) >= min_pre]
        if not T or not C:
            continue
        for a in T + C:
            for y in range(c - window, c + window + 1):
                if YEARS[0] <= y <= YEARS[-1]:
                    rows.append((a, c, y, y - c, a in T, float(outcome.loc[a, y]) if a in outcome.index else 0.0))
    L = pd.DataFrame(rows, columns=["ai", "coh", "year", "k", "treated", "n"])
    L["unit"] = L.ai + "|" + L.coh.astype(str)
    L["cy"] = L.coh.astype(str) + "|" + L.year.astype(str)
    return L


def review_timing(inp):
    """North and South effects of non-renewal against renewal at the scheduled review date."""
    treated, controls = review_cohorts(inp)

    def run(t, c):
        """Fit the joint North and South effect for one treated and control set; None when a side is empty."""
        LN = build_review(inp, t, c, inp.author_panels["North"])
        LS = build_review(inp, t, c, inp.author_panels["South"])
        if LN.empty or LN.treated.sum() == 0 or (~LN.treated).sum() == 0:
            return None
        r = joint_effect(LN, LS)
        r["nT"], r["nC"] = int(LN[LN.treated].ai.nunique()), int(LN[~LN.treated].ai.nunique())
        return r

    real = run(treated, controls)
    reference = {**treated, **{a: y for a, y in controls.items() if a not in treated}}
    ids = list(reference)
    rng = np.random.default_rng(CFG["seeds"]["review_permutation"])
    perm = []
    for _ in range(ES["review_permutations"]):
        chosen = set(rng.choice(ids, size=len(treated), replace=False))
        try:
            r = run({a: reference[a] for a in ids if a in chosen}, {a: reference[a] for a in ids if a not in chosen})
            if r and r["nT"] >= 8 and r["nC"] >= 4:
                perm.append((r["bN"], r["bS"], r["d"]))
        except Exception:
            pass
    perm = np.array(perm)
    return dict(
        real={k: float(v) for k, v in real.items()},
        North_pct=pct_change(real["bN"]),
        South_pct=pct_change(real["bS"]),
        perm_n=len(perm),
        perm_sd_N=float(perm[:, 0].std()),
        perm_sd_S=float(perm[:, 1].std()),
        p_N_le=float((perm[:, 0] <= real["bN"]).mean()),
        p_S_le=float((perm[:, 1] <= real["bS"]).mean()),
        p_diff_le=float((perm[:, 2] <= real["d"]).mean()),
    )


if __name__ == "__main__":
    inputs = Inputs()
    result = {"stacked": stacked_calibration(inputs), "review_timing": review_timing(inputs), "event_time": event_time(inputs), "balance": balance(inputs)}
    out = path("results", "event_studies.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    json.dump(result, open(out, "w"), indent=1)
    print(f"wrote {out}")
