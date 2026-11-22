"""Stage 06. Cross-sectional analysis: every reported number that does not come from the event studies.

Reads  data/processed/ingredient_table.parquet and the cached counts (data/interim) plus the registers and use tables.
Writes results/numbers.json (consumed by 11_write_ledger.py).

Sections: sample and headline ratios, five-year windows, registers, funders, human-health evidence, group shares, selection flow,
named ingredients, missing papers, country shares, atlas imputation, applied mass and crops, coverage against FAOSTAT, register cross-check.
"""

import json
import zipfile

import numpy as np
import pandas as pd

from pesticide_evidence import geography
from pesticide_evidence.atlas import reported_atlas
from pesticide_evidence.config import CFG, path
from pesticide_evidence.funders import classify_funder, classify_funder_by_rule
from pesticide_evidence.selection import STEPS
from pesticide_evidence.stats import ABSENT, GROUP_KEYS, analysis_sample, design, fit, ppml, ppml_cov, round5




def triple(t):
    """Estimate and interval limits rounded to five decimals."""
    return [round5(t[0]), round5(t[1]), round5(t[2])]


def ratios(data, outcome="en_zh_tox", **kwargs):
    """Primary-model ratios for the three non-authorized groups, with the Southern-use elasticity."""
    out, coefs, _ = fit(data, outcome, **kwargs)
    return {k: triple(v) for k, v in out.items()}, round5(coefs["ln_south"])


def load_counts(table, sample):
    """Attach the cached literature counts that are not part of the ingredient table."""
    cas = json.load(open(path("interim", "pubmed", "cas_registry_counts.json")))
    sample["rn"] = [float((cas.get(a) or {}).get("rn", np.nan)) for a in sample.index]
    return sample


def headline(sample):
    """Sample sizes and the group ratios under the primary outcome, PubMed only, all topics and the CAS registry count."""
    out = {"sample": {"n": len(sample), **sample.grp4.value_counts().to_dict()}, "headline": {}}
    for label, outcome in (
        ("combined_tox_pubmed_plus_zh", "en_zh_tox"),
        ("pubmed_tox_only", "pm_tox"),
        ("all_topics_combined", "en_zh_all"),
        ("cas_registry_pubmed", "rn"),
    ):
        r, e = ratios(sample, outcome)
        out["headline"][label] = {**r, "elasticity_south_use": e}
    out["raw_median_papers"] = sample.groupby("grp4").en_zh_tox.median().to_dict()
    return out


def registers(sample):
    """Ratios for ingredients missing from each register, separately, in the none-of-three contrast and in the joint model."""
    s = sample.copy()
    s["not_EU"], s["not_US"], s["not_CA"] = (~s.eu_ok).astype(float), (~s.us_ok).astype(float), (~s.ca_ok).astype(float)
    s["none3"] = s.not_EU * s.not_US * s.not_CA
    out = {}
    for label, terms in (
        ("not_EU", ["not_EU"]),
        ("not_US", ["not_US"]),
        ("not_CA", ["not_CA"]),
        ("none_of_three", ["none3"]),
        ("joint_EU_US_CA", ["not_EU", "not_US", "not_CA"]),
    ):
        r, _ = ratios(s, terms=terms)
        out[label] = r
    rng = np.random.default_rng(CFG["seeds"]["register_bootstrap"])
    boot = []
    for _ in range(CFG["bootstrap_replicates"]):
        draw = s.iloc[rng.integers(0, len(s), len(s))]
        try:
            r, _ = ratios(draw, terms=["not_EU", "not_US"])
            boot.append(np.log(r["not_US"][0] / r["not_EU"][0]))
        except Exception:
            pass
    boot = np.array(boot)
    out["US_over_EU_joint"] = {
        "median": round5(np.exp(np.median(boot))),
        "ci": [round5(np.exp(np.percentile(boot, 2.5))), round5(np.exp(np.percentile(boot, 97.5)))],
        "share_US_lt_EU": round5(np.mean(boot < 0)),
        "reps": len(boot),
    }
    return {"registers": out}


def windows(sample):
    """Group ratios for toxicity papers in each five-year window."""
    counts = json.load(open(path("interim", "pubmed", "window_counts.json")))
    s = sample[sample.index.isin(counts)].copy()
    out = {"n": len(s)}
    for j, (lo, hi) in enumerate(CFG["windows"]):
        s["w"] = [counts[a][j] for a in s.index]
        r, _ = ratios(s, "w")
        out[f"{lo}-{str(hi)[-2:]}"] = {"absent": r["absent"], "lapsed": r["lapsed"], "withdrawn": r["withdrawn"]}
    return {"windows": out}


def funder_table(sample, classifier):
    """Per-ingredient counts of toxicity works by funder class, from the OpenAlex funder records."""
    counts = json.load(open(path("interim", "openalex", "funder_counts.json")))
    rows = []
    for ai, rec in counts.items():
        tally = {"N": 0, "S": 0, "I": 0, "O": 0}
        for _, name, n in rec["funders"]:
            tally[classifier(name)] += n
        rows.append(dict(ai=ai, n=rec["n"], **tally))
    return sample.join(pd.DataFrame(rows).set_index("ai"), how="inner")


def funder_gap(s):
    """Joint model of Northern-funded and Southern-funded works, with the difference between the two absent-group effects."""
    X, names = design(s)
    k = names.index("absent")
    X_stack = np.block([[X, np.zeros_like(X)], [np.zeros_like(X), X]])
    y = np.concatenate([s.N.values.astype(float), s.S.values.astype(float)])
    init = np.zeros(X_stack.shape[1])
    init[0], init[X.shape[1]] = np.log(max(s.N.mean(), 0.1)), np.log(max(s.S.mean(), 0.1))
    b, cov = ppml_cov(y, X_stack, np.concatenate([s.index.values, s.index.values]), init=init)
    d = b[k] - b[X.shape[1] + k]
    sd = np.sqrt(cov[k, k] + cov[X.shape[1] + k, X.shape[1] + k] - 2 * cov[k, X.shape[1] + k])
    return {
        "absent_ratio_northern_funded": round(float(np.exp(b[k])), 4),
        "absent_ratio_southern_funded": round(float(np.exp(b[X.shape[1] + k])), 4),
        "ratio_of_ratios_N_over_S": round(float(np.exp(d)), 3),
        "ci": [round(float(np.exp(d - 1.96 * sd)), 3), round(float(np.exp(d + 1.96 * sd)), 3)],
        "z": round(float(d / sd), 2),
        "n_ingredients": len(s),
    }


def funders(sample):
    """Absent-group ratios by funder class, and the share of works with a Northern or Southern funder in each group."""
    s = funder_table(sample, classify_funder)
    out = {"n_ingredients": len(s), "groups": s.grp4.value_counts().to_dict()}
    for label, col in (("all", "n"), ("northern", "N"), ("southern", "S"), ("industry", "I"), ("unclassified", "O")):
        out[label] = ratios(s, col)[0]["absent"]
    shares = s.groupby("grp4")[["n", "N", "S"]].sum()
    out["share_works_northern_funder"] = {g: round5(shares.loc[g, "N"] / shares.loc[g, "n"]) for g in shares.index}
    out["share_works_southern_funder"] = {g: round5(shares.loc[g, "S"] / shares.loc[g, "n"]) for g in shares.index}
    diff = funder_gap(s)
    rule = funder_table(sample, classify_funder_by_rule)
    out["rule_only"] = {
        "northern": ratios(rule, "N")[0]["absent"],
        "southern": ratios(rule, "S")[0]["absent"],
        "unclassified": ratios(rule, "O")[0]["absent"],
        "diff": funder_gap(rule),
    }
    funded = json.load(open(path("interim", "openalex", "funded_counts.json")))
    if all(a in funded for a in s.index):
        s = s.assign(funded=[funded.get(a, 0) for a in s.index])
        coverage = s.groupby("grp4")[["funded", "n"]].sum()
        out["share_works_with_listed_funder"] = {g: round5(coverage.loc[g, "funded"] / coverage.loc[g, "n"]) for g in coverage.index}
        out["absent_ratio_works_with_listed_funder"] = ratios(s, "funded")[0]["absent"]
        out["northern_funder_share_of_funded_works"] = {g: round5(shares.loc[g, "N"] / coverage.loc[g, "funded"]) for g in shares.index}
        out["southern_funder_share_of_funded_works"] = {g: round5(shares.loc[g, "S"] / coverage.loc[g, "funded"]) for g in shares.index}
    return {"funders_openalex": out, "funders_diff": diff}


def human_evidence(sample):
    """Absent-group ratios for human, epidemiology, poisoning and rodent studies, and the interaction with WHO hazard class."""
    counts = json.load(open(path("interim", "pubmed", "human_evidence_counts.json")))
    s = sample[sample.index.isin(counts)].copy()
    for key in ("human", "epi", "poison", "mamm"):
        s[key] = [counts[a][key] for a in s.index]
    who = pd.read_csv(path("processed", "who_hazard_class.csv")).set_index("ai").who_class
    s["who"] = who.reindex(s.index)
    s["haz"] = s.who.isin(["Ia", "Ib", "II"]).astype(float)
    s["abs_x_haz"] = (s.grp4 == ABSENT).astype(float) * s.haz
    out = {"sample_n": len(s), "total_papers": {k: int(s[k].sum()) for k in ("human", "epi", "poison", "mamm")}, "zero_paper_share_by_group": {}}
    for key in ("human", "epi", "poison", "mamm"):
        r, _, _ = fit(s, key)
        out[key] = {"n": len(s), **{k: [round(v[0], 4), round(v[1], 4), round(v[2], 4)] for k, v in r.items()}}
        out["zero_paper_share_by_group"][key] = s.groupby("grp4")[key].apply(lambda x: round(float((x == 0).mean()), 3)).to_dict()

    def interaction(data):
        """Primary model with an absent-by-hazard interaction."""
        full, c2, s2 = fit(data, "en_zh_tox", terms=GROUP_KEYS + ["abs_x_haz"], extra=["haz"])
        res = {"n": len(data), **{k: [round(v[0], 4), round(v[1], 4), round(v[2], 4)] for k, v in full.items() if k in GROUP_KEYS}}
        res["hazardous_main"] = [
            round(float(np.exp(c2["haz"])), 3),
            round(float(np.exp(c2["haz"] - 1.96 * s2["haz"])), 3),
            round(float(np.exp(c2["haz"] + 1.96 * s2["haz"])), 3),
        ]
        res["absent_x_hazardous"] = [round(v, 3) for v in full["abs_x_haz"]]
        return res

    out["tox_hazard_interaction"] = interaction(s)
    out["n_hazardous_by_group"] = s.groupby("grp4").haz.sum().astype(int).to_dict()
    thin = s[(s.haz == 1) & (s.human <= 20)]
    out["mass_hazardous_thin_human"] = {
        "share_of_sample_South_mass": round(float(thin.South.sum() / s.South.sum()), 4),
        "by_group_kt": {g: round(float(thin[thin.grp4 == g].South.sum() / 1000), 1) for g in s.grp4.unique()},
    }
    return {"human": out}


BOMBARDI_TOP5 = ["mancozeb", "atrazine", "acephate", "chlorothalonil", "chlorpyrifos"]


def named_ingredients(sample):
    """Observed and expected papers for the named ingredients, with the expectation from the model fitted to authorized pesticides."""
    X, _ = design(sample, terms=[])
    authorized = (sample.grp4 == "authorized").values
    b, _ = ppml(sample.en_zh_tox.values.astype(float)[authorized], X[authorized], np.arange(authorized.sum()))
    expected = pd.Series(np.exp(X @ b), index=sample.index)
    rows = {
        a: {
            "group": sample.grp4[a],
            "eu": bool(sample.eu_ok[a]),
            "us": bool(sample.us_ok[a]),
            "ca": bool(sample.ca_ok[a]),
            "observed": float(sample.en_zh_tox[a]),
            "expected_authorized_fit": round5(expected[a]),
            "ratio": round5(sample.en_zh_tox[a] / expected[a]),
        }
        for a in BOMBARDI_TOP5
        if a in sample.index
    }
    return {"bombardi_top5": rows}


def group_shares(sample):
    """Share of Southern mass and of papers in each group, and the median market age."""
    mass = sample.groupby("grp4").South.sum() / sample.South.sum()
    papers = sample.groupby("grp4").en_zh_tox.sum() / sample.en_zh_tox.sum()
    return {
        "group_shares": {
            g: {"mass": round5(mass[g]), "papers": round5(papers[g]), "papers_per_mass": round5(papers[g] / mass[g]), "median_age": float(sample.age[sample.grp4 == g].median())}
            for g in mass.index
        }
    }


def selection_flow():
    """Counts of PubMed queries, records, rule exclusions and selected papers at each step."""
    corpus = pd.read_parquet(path("interim", "corpus", "core_corpus.parquet"), columns=["pmid", "core", "excl", "selected", "selection_step"])
    ruled_out = corpus.excl != ""
    removed = corpus.selection_step.value_counts()
    return {
        "selection_flow": {
            "queries": len(list(path("interim", "pubmed", "queries").glob("*.json"))),
            "records": len(corpus),
            "excluded_by_rule": int(ruled_out.sum()),
            "excluded_other": int((~corpus.core & ~ruled_out).sum()),
            "rule_retained": int(corpus.core.sum()),
            "removed_by_step": {step: int(removed.get(step, 0)) for step in STEPS},
            "selected": int(corpus.selected.sum()),
        }
    }


def missing_papers(sample):
    """Papers expected for the absent group had it been authorized, against the papers observed."""
    X, names = design(sample)
    b, _ = ppml(sample.en_zh_tox.values.astype(float), X, np.arange(len(sample)))
    absent = (sample.grp4 == ABSENT).values
    X0 = X.copy()
    X0[:, names.index("absent")] = 0
    expected = pd.Series(np.exp(X0 @ b), index=sample.index)[absent]
    observed = sample.en_zh_tox[absent]
    by_ingredient = {a: {"south_t": round5(sample.South[a]), "observed": float(observed[a]), "expected": round5(expected[a])} for a in sample.index[absent]}
    return {
        "missing_papers": {
            "n_absent": int(absent.sum()),
            "observed": float(observed.sum()),
            "expected": round5(expected.sum()),
            "missing": round5(expected.sum() - observed.sum()),
            "by_ingredient": by_ingredient,
        }
    }


def country_absent_share(table):
    """Share of each Southern country's reported mass that is in pesticides absent from the registers."""
    atlas = pd.read_csv(path("processed", "use_atlas_country_ingredient.csv"))
    atlas = atlas[(atlas.origin == "Reported") & atlas.ai.isin(table.index)].copy()
    atlas["absent"] = atlas.ai.map(table.grp4) == ABSENT
    south = atlas[~atlas.country.isin(geography.NORTH_USE)]
    by = south.groupby("country").apply(lambda d: pd.Series({"t": d.tonnes.sum(), "absent_t": d.tonnes[d.absent].sum()}), include_groups=False)
    by["share"] = by.absent_t / by.t
    iso = geography.atlas_with_regions(south[["country"]].drop_duplicates()).set_index("country").iso3
    by["iso3"] = iso
    top = by.sort_values("absent_t", ascending=False).head(8)
    return {
        "country_absent_share": {
            "countries": int(len(by)),
            "by_country": {c: {"iso3": r.iso3, "t": round5(r.t), "absent_t": round5(r.absent_t), "share": round5(r.share)} for c, r in by.iterrows()},
            "top_by_absent_mass": {c: round5(r.share) for c, r in top.iterrows()},
            "share_in_top_three_of_absent_mass": round5(top.absent_t.head(3).sum() / by.absent_t.sum()),
        }
    }


def atlas_imputation_check(table, sample):
    """Reported and imputed countries and mass in the atlas, and the group ratios when imputed use is added to reported use."""
    countries = pd.read_excel(path("raw", "atlas", "sdata1.xlsx"), sheet_name="Table S12")
    countries.columns = ["n", "country", "iso3", "origin", "kg", "continent", "remarks"]
    countries = countries.dropna(subset=["iso3"])
    wb = geography.world_bank()
    countries["region"] = countries.iso3.map(lambda i: wb.get(i, (None,))[0])
    ssa = countries[countries.region == "Sub-Saharan Africa"]
    atlas = pd.read_csv(path("processed", "use_atlas_country_ingredient.csv"))
    atlas["south"] = ~atlas.country.isin(geography.NORTH_USE)
    alt = sample.copy()
    for label, rows in (("reported", atlas[atlas.origin == "Reported"]), ("reported_and_imputed", atlas)):
        south = rows[rows.south].groupby("ai").tonnes.sum().reindex(alt.index).fillna(0)
        north = rows[~rows.south].groupby("ai").tonnes.sum().reindex(alt.index).fillna(0)
        alt["South_" + label], alt["North_" + label] = south, north
    check = alt.assign(South=alt.South_reported_and_imputed, North=alt.North_reported_and_imputed)
    r, _ = ratios(check)
    imputed = atlas[(atlas.origin == "Imputed") & atlas.south & atlas.ai.isin(sample.index)]
    return {
        "atlas_imputation": {
            "countries_reported": int((countries.origin == "Reported").sum()),
            "countries_imputed": int((countries.origin == "Imputed").sum()),
            "ssa_countries_reported": int((ssa.origin == "Reported").sum()),
            "ssa_countries_imputed": int((ssa.origin == "Imputed").sum()),
            "ssa_imputed_kt_all_ingredients": round5(ssa[ssa.origin == "Imputed"].kg.sum() / 1e6),
            "ssa_reported_kt_all_ingredients": round5(ssa[ssa.origin == "Reported"].kg.sum() / 1e6),
            "imputed_share_of_atlas_mass": round5(countries[countries.origin == "Imputed"].kg.sum() / countries.kg.sum()),
            "sample_imputed_south_kt_by_group": {g: round5(v / 1e3) for g, v in imputed.groupby(imputed.ai.map(sample.grp4)).tonnes.sum().items()},
            "ratios_with_imputed_use": r,
        }
    }


def mass_and_crops(table):
    """Atlas mass, Southern share, China share of absent mass, and per-crop absent shares and papers per kilotonne."""
    atlas = pd.read_csv(path("processed", "use_atlas_country_ingredient.csv"))
    atlas = atlas[(atlas.origin == "Reported") & atlas.ai.isin(table.index)].copy()
    atlas["grp4"] = atlas.ai.map(table.grp4)
    atlas["south"] = ~atlas.country.isin(geography.NORTH_USE)
    total = atlas.tonnes.sum()
    absent_south = atlas[atlas.south & (atlas.grp4 == ABSENT)]
    mass = {
        "atlas_kt_446_ingredients": round5(total / 1e3),
        "south_share": round5(atlas[atlas.south].tonnes.sum() / total),
        "china_share_of_absent_south_mass": round5(absent_south[absent_south.country == "China"].tonnes.sum() / absent_south.tonnes.sum()),
    }
    x = reported_atlas()
    x = x[x.ai.isin(table.index)].copy()
    x["grp4"] = x.ai.map(table.grp4)
    x["t"] = x.kg / 1000
    papers = table.pm_tox.fillna(0) + table.oa_zh_tox.fillna(0)
    x["ppk"] = x.ai.map(papers / (x.groupby("ai").t.sum() / 1e3))
    crops = {}
    for crop in ["Rice", "Soybean", "Maize", "Wheat", "Cotton", "Citrus", "Sugarcane", "Potato", "Grape", "Apple", "Rapeseed"]:
        d = x[x.crop.str.lower().str.startswith(crop.lower())]
        if len(d) == 0:
            continue
        a, au = d[d.grp4 == ABSENT], d[d.grp4 == "authorized"]
        crops[crop] = {
            "kt": round5(d.t.sum() / 1e3),
            "absent_share": round5(a.t.sum() / d.t.sum()),
            "ppk_absent": round5((a.ppk * a.t).sum() / a.t.sum()) if a.t.sum() > 0 else None,
            "ppk_authorized": round5((au.ppk * au.t).sum() / au.t.sum()),
        }
    return {"mass": mass, "crops": crops}


def coverage():
    """Atlas coverage against FAOSTAT and GloPUT by region."""
    atlas = pd.read_csv(path("processed", "use_atlas_country_ingredient.csv"))
    atlas = geography.atlas_with_regions(atlas[atlas.origin == "Reported"])
    atlas["wbreg"] = atlas.region.astype(str).str.strip()
    world_bank = geography.world_bank()
    wb_names = {v[2].upper(): v for v in world_bank.values()}
    alias = {
        "UNITED REPUBLIC OF TANZANIA": "TANZANIA",
        "CÔTE D'IVOIRE": "COTE D'IVOIRE",
        "EGYPT": "EGYPT, ARAB REP.",
        "TÜRKIYE": "TURKIYE",
        "REPUBLIC OF MOLDOVA": "MOLDOVA",
        "VIET NAM": "VIETNAM",
        "BOLIVIA (PLURINATIONAL STATE OF)": "BOLIVIA",
        "CHINA, MAINLAND": "CHINA",
        "REPUBLIC OF KOREA": "KOREA, REP.",
        "IRAN (ISLAMIC REPUBLIC OF)": "IRAN, ISLAMIC REP.",
        "VENEZUELA (BOLIVARIAN REPUBLIC OF)": "VENEZUELA, RB",
        "UNITED STATES OF AMERICA": "UNITED STATES",
        "CZECHIA": "CZECH REPUBLIC",
        "SLOVAKIA": "SLOVAK REPUBLIC",
        "LAO PEOPLE'S DEMOCRATIC REPUBLIC": "LAO PDR",
        "DEMOCRATIC REPUBLIC OF THE CONGO": "CONGO, DEM. REP.",
        "CONGO": "CONGO, REP.",
        "GAMBIA": "GAMBIA, THE",
    }
    with zipfile.ZipFile(path("raw", "fao", "fao_pesticides_use.zip")) as z:
        fao = pd.read_csv(z.open("Inputs_Pesticides_Use_E_All_Data_(Normalized).csv"), encoding="latin-1")
    fao = fao[(fao.Item == "Pesticides (total)") & (fao.Element == "Agricultural Use") & fao.Year.between(2016, 2020)]
    mean_use = fao.groupby("Area").Value.mean().to_frame("t")
    region = {}
    for area in mean_use.index:
        key = area.upper().replace(" (THE)", "").strip()
        hit = wb_names.get(key) or wb_names.get(alias.get(key, "").upper())
        region[area] = hit[0] if hit else None
    mean_use["region"] = pd.Series(region)
    fao_kt = mean_use.dropna().groupby("region").t.sum() / 1000
    glopu = pd.read_csv(path("raw", "glopu", "GloPUT_1995-2019.csv"))
    glopu["region"] = glopu.ISO3.map(lambda iso: world_bank.get(iso, (None,))[0])
    glopu_kt = glopu.assign(t=glopu[["2016", "2017", "2018", "2019"]].mean(axis=1)).groupby("region").t.sum() / 1000
    atlas_kt = atlas.groupby("wbreg").tonnes.sum() / 1000
    tab = pd.DataFrame({"fao_kt": fao_kt, "glopu_kt": glopu_kt, "atlas_kt": atlas_kt}).fillna(0)
    tab["atlas_over_fao"] = tab.atlas_kt / tab.fao_kt.replace(0, np.nan)
    tab["atlas_over_glopu"] = tab.atlas_kt / tab.glopu_kt.replace(0, np.nan)
    africa = atlas[atlas.wbreg == "Sub-Saharan Africa"]
    return {
        "coverage": {
            "fao_total_kt_region_sum": round(float(fao_kt.sum())),
            "by_region": {k: {kk: (round(float(vv), 2) if vv == vv else None) for kk, vv in r.items()} for k, r in tab.to_dict("index").items()},
            "africa": {
                "atlas_countries": int(africa.country.nunique()),
                "wb_countries": 48,
                "atlas_kt": round(float(africa.tonnes.sum() / 1000), 1),
                "share_atlas_mass": round(float(africa.tonnes.sum() / atlas.tonnes.sum()), 4),
                "fao_share_world": round(float(fao_kt.get("Sub-Saharan Africa", 0) / fao_kt.sum()), 4),
                "glopu_ssa_kt": round(float(glopu_kt.get("Sub-Saharan Africa", 0)), 1),
                "nigeria_glopu_kt": round(float(glopu.loc[glopu.ISO3 == "NGA", ["2016", "2017", "2018", "2019"]].mean(axis=1).iloc[0] / 1000), 1),
                "nigeria_fao_kt": round(float(mean_use.loc["Nigeria", "t"] / 1000), 1),
            },
        }
    }


def register_crosscheck(sample):
    """Agreement of the coded register status with the independently compiled table."""
    other = pd.read_csv(path("interim", "crosscheck", "independent_register_table.csv"))
    other["key"] = other.name_norm.astype(str).str.lower().str.strip()
    m = sample.join(other.drop_duplicates("key").set_index("key")[["eu_status", "us_n_active", "ca_n_current"]], how="left")
    m = m[m.eu_status.notna()].copy()
    m["o_eu"], m["o_us"], m["o_ca"] = m.eu_status.isin(["Approved", "Pending"]), m.us_n_active.fillna(0) > 0, m.ca_n_current.fillna(0) > 0
    m["o_any"] = m.o_eu | m.o_us | m.o_ca
    m["our_any"] = m.eu_ok | m.us_ok | m.ca_ok
    absent = m[m.grp4 == ABSENT]
    return {
        "register_crosscheck": {
            "sample_n": int(len(sample)),
            "matched_n": int(len(m)),
            "agreement": {
                r: round(float((m[a] == m[b]).mean()), 4)
                for r, a, b in (("EU", "eu_ok", "o_eu"), ("US", "us_ok", "o_us"), ("Canada", "ca_ok", "o_ca"), ("any register", "our_any", "o_any"))
            },
            "absent_matched_n": int(len(absent)),
            "absent_registered_in_independent_table": int(absent.o_any.sum()),
            "disagreements_any_register": sorted(m.index[m.our_any != m.o_any].tolist()),
        }
    }


if __name__ == "__main__":
    table = pd.read_parquet(path("processed", "ingredient_table.parquet"))
    sample = load_counts(table, analysis_sample(table, CFG["sample"]["min_south_use_t"]))
    numbers = {}
    for section in (
        headline(sample),
        registers(sample),
        windows(sample),
        funders(sample),
        human_evidence(sample),
        group_shares(sample),
        selection_flow(),
        named_ingredients(sample),
        missing_papers(sample),
        country_absent_share(table),
        atlas_imputation_check(table, sample),
        mass_and_crops(table),
        coverage(),
        register_crosscheck(sample),
    ):
        numbers.update(section)
    out = path("results", "numbers.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    json.dump(numbers, open(out, "w"), indent=1, default=str)
    print(f"wrote {out} with {len(numbers)} sections")
