"""Stage 04. Literature and funder counts from PubMed and OpenAlex.

Every step caches to data/interim and resumes from the cache. Run one step with `--only <step>`; limit the ingredients with `--limit N` for a test.

Steps and outputs (data/interim):
    first_record     corpus/first_record.json            year of the first PubMed record per ingredient
    name_counts      pubmed/name_counts.json             PubMed papers mentioning the ingredient, with and without toxicity terms, 2005-2025
    cas_counts       pubmed/cas_registry_counts.json     PubMed papers indexed with the ingredient's CAS registry number
    window_counts    pubmed/window_counts.json           toxicity papers in the five-year windows
    human_counts     pubmed/human_evidence_counts.json   human, epidemiology, poisoning and rodent toxicity papers
    chinese_counts   openalex/chinese_counts.json        OpenAlex Chinese-language papers, with and without toxicity terms
    funder_counts    openalex/funder_counts.json         OpenAlex toxicity papers by funder
    funded_counts    openalex/funded_counts.json         OpenAlex toxicity papers with at least one listed funder

Keys are read from the environment: OPENALEX_API_KEY and OPENALEX_MAILTO.
"""

import argparse
import json
import os
import time

import pandas as pd
import requests

from pesticide_evidence.config import CFG, ROOT, path
from pesticide_evidence.names import AMBIG, canonical, class_of_name

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
TOX_PUBMED = '(toxic*[tiab] OR toxicity[tiab] OR ecotox*[tiab] OR genotox*[tiab] OR cytotox*[tiab] OR "oxidative stress"[tiab] OR mortality[tiab] OR LC50[tiab])'
TOX_OPENALEX = '(toxic OR toxicity OR ecotoxicity OR genotoxicity OR cytotoxicity OR "oxidative stress" OR mortality OR LC50)'
TOX_ZH = "(毒性 OR 毒理 OR 急性毒性 OR 毒害)"
EPIDEMIOLOGY = '(epidemiolog*[tiab] OR cohort*[tiab] OR "case-control"[tiab] OR "cross-sectional"[tiab] OR biomonitor*[tiab] OR urinary[tiab])'
POISONING = '(poisoning[tiab] OR intoxication[tiab] OR overdose[tiab] OR "self-poisoning"[tiab] OR suicide[tiab])'
YEARS = CFG["sample"]["years"]


def load_json(file):
    """Cached JSON object, or an empty dict when the cache file does not exist."""
    return json.load(open(file)) if file.exists() else {}


def save_json(obj, file):
    """Write a JSON object to a cache file, creating its folder."""
    file.parent.mkdir(parents=True, exist_ok=True)
    json.dump(obj, open(file, "w"))


def pubmed_count(term, mindate=YEARS[0], maxdate=YEARS[1], pause=0.34):
    """Number of PubMed records matching a query in a publication-date window, retried with back-off; None after repeated failures."""
    for attempt in range(6):
        try:
            r = requests.get(EUTILS, params=dict(db="pubmed", term=term, retmode="json", retmax=0, datetype="pdat", mindate=mindate, maxdate=maxdate), timeout=60)
            time.sleep(pause)
            return int(r.json()["esearchresult"]["count"])
        except (requests.RequestException, KeyError, ValueError):
            time.sleep(2**attempt)
    return None


def openalex_get(params, endpoint="works"):
    """GET request to the OpenAlex API with the contact address and key from the environment, retried with back-off."""
    params = {**params, "mailto": os.environ.get(CFG["openalex"]["mailto_env"], "")}
    key = os.environ.get(CFG["openalex"]["key_env"])
    if key:
        params["api_key"] = key
    for attempt in range(6):
        r = requests.get(f"https://api.openalex.org/{endpoint}", params=params, timeout=60)
        if r.status_code == 200:
            return r.json()
        time.sleep(2**attempt * 5 if r.status_code == 429 else 2**attempt)
    return None


def ingredient_names():
    """Table of curated ingredient names (CAS number and Chinese name), indexed by ingredient."""
    return pd.read_csv(ROOT / "config" / "ingredient_names.csv", index_col=0)


def query_ingredients(limit=None):
    """Ingredients with at least the minimum Southern applied mass that have a curated name entry."""
    from pesticide_evidence import geography

    atlas = pd.read_csv(path("processed", "use_atlas_country_ingredient.csv"))
    atlas = atlas[(atlas.origin == "Reported") & ~atlas.country.isin(geography.NORTH_USE)]
    south = atlas.groupby("ai").tonnes.sum()
    names = ingredient_names()
    ids = [a for a in south[south >= CFG["sample"]["min_south_use_t"]].sort_values(ascending=False).index if a in names.index]
    return ids[:limit] if limit else ids


def first_record(limit=None):
    """Year of the first PubMed record naming each ingredient, found by bisection over publication years."""
    file = path("interim", "corpus", "first_record.json")
    result = load_json(file)
    generic = " OR ".join(f"{t}*[tiab]" for t in ("pesticide", "herbicide", "insecticide", "fungicide", "acaricide", "nematicide", "rodenticide", "agrochemical"))
    ingredients = sorted({canonical(n) for n in class_of_name()})
    todo = [a for a in ingredients if a not in result and a != "copper"][:limit]
    for ai in todo:
        term = f'("{ai}"[tiab] AND ({generic}))' if ai in AMBIG or len(ai) < 5 else f'"{ai}"[tiab]'
        if pubmed_count(term, "1900/01/01", "2025/12/31") == 0:
            result[ai] = None
            continue
        lo, hi = 1900, 2025
        while lo < hi:
            mid = (lo + hi) // 2
            if pubmed_count(term, "1900/01/01", f"{mid}/12/31") > 0:
                hi = mid
            else:
                lo = mid + 1
        result[ai] = lo
        save_json(result, file)
    save_json(result, file)


def name_counts(limit=None):
    """PubMed records naming the ingredient in title or abstract, with and without a toxicity term."""
    file = path("interim", "pubmed", "name_counts.json")
    result = load_json(file)
    for ai in list(ingredient_names().index)[:limit]:
        if ai not in result:
            result[ai] = {"all": pubmed_count(f'"{ai}"[tiab]'), "tox": pubmed_count(f'"{ai}"[tiab] AND {TOX_PUBMED}')}
            save_json(result, file)


def cas_counts(limit=None):
    """PubMed records indexed with the ingredient's CAS registry number, with and without toxicity subheadings."""
    file = path("interim", "pubmed", "cas_registry_counts.json")
    result = load_json(file)
    names = ingredient_names()
    for ai in query_ingredients(limit):
        cas = names.loc[ai, "cas"]
        if ai in result or not isinstance(cas, str):
            continue
        result[ai] = {"rn": pubmed_count(f'"{cas}"[rn]'), "rn_tox": pubmed_count(f'"{cas}"[rn] AND (toxicity[sh] OR poisoning[sh] OR "adverse effects"[sh])')}
        save_json(result, file)


def window_counts(limit=None):
    """PubMed toxicity records per five-year window."""
    file = path("interim", "pubmed", "window_counts.json")
    result = load_json(file)
    for ai in query_ingredients(limit):
        if ai not in result:
            result[ai] = [pubmed_count(f'"{ai}"[tiab] AND {TOX_PUBMED}', lo, hi) for lo, hi in CFG["windows"]]
            save_json(result, file)


def human_counts(limit=None):
    """PubMed records of human, epidemiology, poisoning and rodent toxicity evidence for each ingredient."""
    file = path("interim", "pubmed", "human_evidence_counts.json")
    result = load_json(file)
    for ai in query_ingredients(limit):
        if ai in result:
            continue
        q = f'"{ai}"[tiab]'
        result[ai] = {
            "human": pubmed_count(f"{q} AND humans[mh]"),
            "epi": pubmed_count(f"{q} AND humans[mh] AND {EPIDEMIOLOGY}"),
            "poison": pubmed_count(f"{q} AND {POISONING}"),
            "mamm": pubmed_count(f"{q} AND (rats[mh] OR mice[mh]) AND {TOX_PUBMED}"),
        }
        save_json(result, file)


def chinese_counts(limit=None):
    """OpenAlex works in Chinese that name the ingredient, with and without a toxicity term."""
    file = path("interim", "openalex", "chinese_counts.json")
    result = load_json(file)
    names = ingredient_names()
    window = f"publication_year:{YEARS[0]}-{YEARS[1]}"
    for ai in query_ingredients(limit):
        if ai in result:
            continue
        record = {}
        zh = names.loc[ai, "zh_name"]
        if isinstance(zh, str):
            for key, query in (("zh_all", f'"{zh}"'), ("zh_tox", f'("{zh}" AND {TOX_ZH})')):
                data = openalex_get({"filter": f"title_and_abstract.search:{query},{window}", "per-page": 1})
                if data:
                    record[key] = data["meta"]["count"]
        result[ai] = record
        save_json(result, file)


def funder_counts(limit=None):
    """OpenAlex toxicity works for each ingredient, grouped by funder."""
    file = path("interim", "openalex", "funder_counts.json")
    result = load_json(file)
    for ai in query_ingredients(limit):
        if ai in result:
            continue
        data = openalex_get({"filter": f'title_and_abstract.search:("{ai}" AND {TOX_OPENALEX}),publication_year:{YEARS[0]}-{YEARS[1]}', "group_by": "funders.id"})
        if data is None:
            break
        result[ai] = {"n": data["meta"]["count"], "funders": [(g["key"], g["key_display_name"], g["count"]) for g in data["group_by"]]}
        save_json(result, file)


def funded_counts(limit=None):
    """OpenAlex toxicity works for each ingredient that list at least one funder."""
    file = path("interim", "openalex", "funded_counts.json")
    result = load_json(file)
    for ai in query_ingredients(limit):
        if ai in result:
            continue
        data = openalex_get({"filter": f'title_and_abstract.search:("{ai}" AND {TOX_OPENALEX}),publication_year:{YEARS[0]}-{YEARS[1]},funders.id:!null', "per-page": 1})
        if data is None:
            break
        result[ai] = data["meta"]["count"]
        save_json(result, file)


STEPS = {
    "first_record": first_record,
    "name_counts": name_counts,
    "cas_counts": cas_counts,
    "window_counts": window_counts,
    "human_counts": human_counts,
    "chinese_counts": chinese_counts,
    "funder_counts": funder_counts,
    "funded_counts": funded_counts,
}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--only", choices=sorted(STEPS))
    parser.add_argument("--limit", type=int)
    args = parser.parse_args()
    for name, step in STEPS.items():
        if args.only in (None, name):
            step(args.limit)
            print(f"{name}: done")
