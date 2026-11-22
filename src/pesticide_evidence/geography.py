"""Country groupings: Global North for use and authorship, World Bank regions, and author country from affiliation text."""

import json
import re

import pandas as pd

from .config import CFG, path

EU27 = {
    "AUT",
    "BEL",
    "BGR",
    "HRV",
    "CYP",
    "CZE",
    "DNK",
    "EST",
    "FIN",
    "FRA",
    "DEU",
    "GRC",
    "HUN",
    "IRL",
    "ITA",
    "LVA",
    "LTU",
    "LUX",
    "MLT",
    "NLD",
    "POL",
    "PRT",
    "ROU",
    "SVK",
    "SVN",
    "ESP",
    "SWE",
}
NORTH_ISO = EU27 | {"GBR", "NOR", "CHE", "ISL", "USA", "CAN", "JPN", "KOR", "AUS", "NZL", "ISR"}
NORTH_USE = set(CFG["north_use_countries"])
NORTH_AUTHORS = set(CFG["north_author_countries"])

COUNTRY_ALIASES = {
    "usa": "United States",
    "united states of america": "United States",
    "u.s.a": "United States",
    "the united states": "United States",
    "u.s.": "United States",
    "uk": "United Kingdom",
    "england": "United Kingdom",
    "scotland": "United Kingdom",
    "wales": "United Kingdom",
    "northern ireland": "United Kingdom",
    "china": "China",
    "p.r. china": "China",
    "pr china": "China",
    "people's republic of china": "China",
    "p. r. china": "China",
    "prc": "China",
    "hong kong": "Hong Kong SAR, China",
    "macau": "Macao SAR, China",
    "taiwan": "Taiwan, China",
    "taiwan, china": "Taiwan, China",
    "republic of china": "Taiwan, China",
    "republic of korea": "Korea, Rep.",
    "south korea": "Korea, Rep.",
    "korea": "Korea, Rep.",
    "korea (republic of)": "Korea, Rep.",
    "iran": "Iran, Islamic Rep.",
    "egypt": "Egypt, Arab Rep.",
    "turkey": "Turkiye",
    "türkiye": "Turkiye",
    "turkiye": "Turkiye",
    "russia": "Russian Federation",
    "russian federation": "Russian Federation",
    "vietnam": "Viet Nam",
    "czech republic": "Czechia",
    "czechia": "Czechia",
    "slovak republic": "Slovak Republic",
    "slovakia": "Slovak Republic",
    "méxico": "Mexico",
    "mexico": "Mexico",
    "tunisie": "Tunisia",
    "brasil": "Brazil",
    "venezuela": "Venezuela, RB",
    "syria": "Syrian Arab Republic",
    "laos": "Lao PDR",
    "the netherlands": "Netherlands",
    "holland": "Netherlands",
    "deutschland": "Germany",
    "españa": "Spain",
    "italia": "Italy",
    "république du cameroun": "Cameroon",
    "cote d'ivoire": "Cote d'Ivoire",
    "ivory coast": "Cote d'Ivoire",
    "palestine": "West Bank and Gaza",
    "kyrgyzstan": "Kyrgyz Republic",
    "north macedonia": "North Macedonia",
    "macedonia": "North Macedonia",
    "bolivia": "Bolivia",
    "tanzania": "Tanzania",
    "gambia": "Gambia, The",
    "yemen": "Yemen, Rep.",
    "congo": "Congo, Rep.",
    "democratic republic of the congo": "Congo, Dem. Rep.",
    "singapore": "Singapore",
    "kosovo": "Kosovo",
    "lebanon": "Lebanon",
}
US_STATE_CODES = set(
    "al ak az ar ca co ct de fl ga hi id il in ia ks ky la me md ma mi mn ms mo mt ne nv nh nj nm ny nc nd oh ok or pa ri sc sd tn tx ut vt va wa wv wi wy dc".split()
)
US_STATE_NAMES = {
    "alabama",
    "alaska",
    "arizona",
    "arkansas",
    "california",
    "colorado",
    "connecticut",
    "delaware",
    "florida",
    "georgia",
    "hawaii",
    "idaho",
    "illinois",
    "indiana",
    "iowa",
    "kansas",
    "kentucky",
    "louisiana",
    "maine",
    "maryland",
    "massachusetts",
    "michigan",
    "minnesota",
    "mississippi",
    "missouri",
    "montana",
    "nebraska",
    "nevada",
    "new hampshire",
    "new jersey",
    "new mexico",
    "new york",
    "north carolina",
    "north dakota",
    "ohio",
    "oklahoma",
    "oregon",
    "pennsylvania",
    "rhode island",
    "south carolina",
    "south dakota",
    "tennessee",
    "texas",
    "utah",
    "vermont",
    "virginia",
    "washington",
    "west virginia",
    "wisconsin",
    "wyoming",
}


def world_bank():
    """ISO3 -> (region, income level, name) for countries (aggregates excluded)."""
    rows = json.load(open(path("raw", "worldbank", "wb_countries.json")))[1]
    return {z["id"]: (z["region"]["value"].strip(), z["incomeLevel"]["value"], z["name"]) for z in rows if z["region"]["value"] != "Aggregates"}


def atlas_with_regions(atlas):
    """Add ISO3 and World Bank region to the country-by-ingredient atlas table."""
    names = pd.read_excel(path("raw", "atlas", "sdata1.xlsx"), sheet_name="Table S4")
    iso = dict(zip(names["Country/region, harmonized"], names["ISO 3166-1 alpha-3"]))
    wb = world_bank()
    atlas = atlas.copy()
    atlas["iso3"] = atlas.country.map(iso)
    atlas["region"] = atlas.iso3.map(lambda i: wb.get(i, (None,))[0])
    return atlas


def _country_names():
    """Lower-case name -> World Bank country name, for full names, names before the first comma, and the aliases."""
    rows = [x for x in json.load(open(path("raw", "worldbank", "wb_countries.json")))[1] if x["region"]["value"] != "Aggregates"]
    names = {x["name"].lower(): x["name"] for x in rows}
    for x in rows:
        names[re.sub(r",.*", "", x["name"]).lower()] = x["name"]
    names.update(COUNTRY_ALIASES)
    return names


COUNTRY_NAMES = _country_names()
NAMES_LONGEST_FIRST = sorted(COUNTRY_NAMES, key=len, reverse=True)


def parse(affiliation):
    """Country of the first affiliation in a PubMed affiliation string, or None."""
    if not affiliation:
        return None
    text = affiliation.split(" | ")[0].strip().rstrip(".")
    text = re.sub(r"(electronic address|email|e-mail)\s*:?.*$", "", text, flags=re.I)
    text = re.sub(r"\S+@\S+", "", text).strip().rstrip(".;, ")
    parts = [t.strip().lower() for t in re.split(r"[,;]", text) if t.strip()]
    if not parts:
        return None
    last = re.sub(r"^[\d\s]+", "", re.sub(r"[\d\.\-]+$", "", parts[-1]).strip())
    if last in COUNTRY_NAMES:
        return COUNTRY_NAMES[last]
    for part in reversed(parts[-3:]):
        part = re.sub(r"[\d]+", "", part).strip(" .")
        if part in COUNTRY_NAMES:
            return COUNTRY_NAMES[part]
    if parts[-1].strip() in US_STATE_CODES or parts[-1].strip(" .") in US_STATE_NAMES or (len(parts) > 1 and parts[-2].strip() in US_STATE_NAMES):
        return "United States"
    tail = text.lower()[-60:]
    for name in NAMES_LONGEST_FIRST:
        if len(name) >= 5 and re.search(r"\b" + re.escape(name) + r"\b", tail):
            return COUNTRY_NAMES[name]
    return None
