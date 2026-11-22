"""Stage 05. Build the ingredient-level analysis table and the corpus panels from the cached inputs.

Inputs  (data/raw, data/interim, data/processed): use atlas, EU/US/Canada registers, EU act scan, name table,
        PubMed and OpenAlex counts, selected corpus, first-record years, PubMed records.
Outputs data/interim/corpus/panel_ingredient_year.parquet      selected papers per ingredient and year
        data/interim/corpus/affiliations_raw.parquet           first- and last-author affiliation text of the selected papers
        data/interim/corpus/affiliations_countries.parquet     author country per paper
        data/processed/ingredient_table.parquet                one row per ingredient (the analysis table)
"""

import json

import numpy as np
import pandas as pd

from pesticide_evidence import geography
from pesticide_evidence.atlas import reported_atlas
from pesticide_evidence.config import CFG, ROOT, path, treaty_ingredients
from pesticide_evidence.pubmed import load_records
from pesticide_evidence.registers import load_eu

YEARS = list(range(CFG["sample"]["years"][0], CFG["sample"]["years"][1] + 1))


def build_panel():
    """Selected toxicology papers per canonical ingredient and year."""
    corpus = pd.read_parquet(path("interim", "corpus", "core_corpus.parquet"))
    core = corpus[corpus.selected]
    long = core.assign(_a=core.compounds.replace("", np.nan).str.split(";")).explode("_a").dropna(subset=["_a"])
    panel = long.groupby(["_a", "year"]).pmid.nunique().unstack(fill_value=0).reindex(columns=YEARS, fill_value=0)
    panel.to_parquet(path("interim", "corpus", "panel_ingredient_year.parquet"))
    return panel


def build_affiliations_raw():
    """First- and last-author affiliation text and author count of every selected paper."""
    corpus = pd.read_parquet(path("interim", "corpus", "core_corpus.parquet"))
    records = load_records(list(corpus.pmid[corpus.selected]))
    rows = []
    for pmid, rec in records.items():
        authors = rec["authors"] or []
        rows.append((pmid, authors[-1]["aff"] if authors else "", authors[0]["aff"] if authors else "", len(authors)))
    raw = pd.DataFrame(rows, columns=["pmid", "aff_last", "aff_first", "n"])
    raw.to_parquet(path("interim", "corpus", "affiliations_raw.parquet"), index=False)
    return raw


def build_affiliation_countries():
    """Author country of the first author, the last author and the main author (last author when there is more than one)."""
    raw = pd.read_parquet(path("interim", "corpus", "affiliations_raw.parquet"))
    years = pd.read_parquet(path("interim", "corpus", "core_corpus.parquet"))[["pmid", "year"]]
    out = raw.merge(years, on="pmid")
    out["c_first"] = out.aff_first.map(geography.parse)
    out["c_last"] = out.aff_last.map(geography.parse)
    out["c_main"] = np.where(out.c_last.notna() & (out.n > 1), out.c_last, out.c_first)
    out["c_main"] = out.c_main.where(out.c_main.notna(), out.c_first)
    out.to_parquet(path("interim", "corpus", "affiliations_countries.parquet"))
    return out


def use_by_region():
    """Reported applied mass in tonnes per ingredient for Southern and Northern countries."""
    atlas = pd.read_csv(path("processed", "use_atlas_country_ingredient.csv"))
    atlas = atlas[atlas.origin == "Reported"].copy()
    atlas["north"] = atlas.country.isin(geography.NORTH_USE)
    use = atlas.pivot_table(index="ai", columns="north", values="tonnes", aggfunc="sum", fill_value=0)
    use.columns = ["South", "North"]
    return use


def atlas_attributes():
    """Crop-group shares of Southern applied mass, and the modal target class per ingredient."""
    x = reported_atlas()
    south = x[~x.country.isin(geography.NORTH_USE)]
    shares = south.pivot_table(index="ai", columns="cg", values="kg", aggfunc="sum", fill_value=0)
    shares = shares.div(shares.sum(axis=1), axis=0)
    shares.columns = ["crop_" + c for c in shares.columns]
    target = x.groupby("ai").target.agg(lambda s: s.value_counts().index[0])
    return shares, target


def group_labels(table, us_reg, eu_scan, treaty):
    """Northern-register groups: authorized, withdrawn for hazard or treaty reasons, lapsed for other reasons, absent."""
    scan = eu_scan[eu_scan.nonren].groupby("ai")[["market", "hazard"]].sum()
    table["eu_hazard"] = (scan.hazard > scan.market).reindex(table.index).fillna(False)
    table["eu_market"] = (scan.market >= scan.hazard).reindex(table.index).fillna(False) & table.index.isin(scan.index)
    table["us_reg"] = (us_reg.share_regulatory >= 0.5).reindex(table.index).fillna(False) & (us_reg.n_current.reindex(table.index) == 0)
    table["treaty"] = table.index.isin(treaty)
    table["grp"] = np.where(
        table.north_any, "allowed", np.where(table.treaty | table.eu_hazard | table.us_reg, "not allowed: hazard/treaty basis", "not allowed: other/unknown")
    )
    return table


def build_table(panel):
    """One row per ingredient: use, register status and dates, authorization group, literature counts and covariates."""
    first_record = json.load(open(path("interim", "corpus", "first_record.json")))
    use = use_by_region()
    eu = load_eu()
    us = pd.read_csv(path("processed", "us_registrations.csv"), index_col=0)
    ca = pd.read_csv(path("processed", "canada_registrations.csv"), index_col=0)
    table = use.copy()
    table["papers"] = panel.sum(axis=1).reindex(table.index).fillna(0)
    table["eu"] = eu.status.reindex(table.index)
    table["us"] = us.n_current.reindex(table.index)
    table["ca"] = ca.n_current.reindex(table.index)
    table["eu_ok"] = table.eu.isin(["Approved", "Pending"])
    table["us_ok"] = table.us.fillna(0) > 0
    table["ca_ok"] = table.ca.fillna(0) > 0
    table["known"] = table.eu.notna()
    table["north_any"] = table.eu_ok | table.us_ok | table.ca_ok
    table["age"] = [max(2025 - first_record[a], 1) if a in first_record else 20 for a in table.index]
    shares, target = atlas_attributes()
    table = table.join(shares).join(target)
    name_counts = pd.DataFrame(json.load(open(path("interim", "pubmed", "name_counts.json")))).T.rename(columns={"all": "pm_all", "tox": "pm_tox"})
    openalex = json.load(open(path("interim", "openalex", "chinese_counts.json")))
    openalex = pd.DataFrame({a: {k: v for k, v in rec.items() if not k.endswith("_lang")} for a, rec in openalex.items()}).T.rename(columns=lambda c: "oa_" + c)
    table = table.join(name_counts).join(openalex)
    table["cas"] = pd.read_csv(ROOT / "config" / "ingredient_names.csv", index_col=0).cas.reindex(table.index)
    table = group_labels(table, us, pd.read_csv(path("processed", "eu_act_text_scan.csv")), treaty_ingredients())
    table["en_zh_tox"] = table.pm_tox + table.oa_zh_tox
    table["en_zh_all"] = table.pm_all + table.oa_zh_all
    table["us_ever"] = table.index.isin(us.index)
    table["ca_ever"] = table.index.isin(ca.index)
    table["eu_dated"] = table.index.map(lambda a: a in eu.index and pd.notna(eu.loc[a, "appr"]))
    table["grp4"] = np.where(
        table.grp == "allowed",
        "authorized",
        np.where(
            table.grp.str.contains("hazard"),
            "withdrawn: hazard/treaty",
            np.where(~table.us_ever & ~table.ca_ever & ~table.eu_dated, "absent from Northern registers", "lapsed: other reason"),
        ),
    )
    return table


if __name__ == "__main__":
    panel = build_panel()
    build_affiliations_raw()
    build_affiliation_countries()
    table = build_table(panel)
    table.to_parquet(path("processed", "ingredient_table.parquet"))
    print(f"ingredients: {len(table)}; groups: {table.grp4.value_counts().to_dict()}")
