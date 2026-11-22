"""Stage 03. Regulatory status and use tables from the raw downloads.

Writes to data/processed:
    us_registrations.csv             EPA product register (PPIS): products, current products, first approval, last cancellation and reason class per ingredient
    canada_registrations.csv         Health Canada PMRA product register: products, current products, last expiry per ingredient
    who_hazard_class.csv             WHO hazard class (2019 classification) per ingredient
    use_atlas_country_ingredient.csv applied mass per country and ingredient from the global atlas of pesticide use
    eu_act_text_scan.csv             reason stated in the EU non-renewal acts (market or hazard)
"""

import json
import re
import time

import numpy as np
import pandas as pd
import requests

from pesticide_evidence.atlas import read_atlas
from pesticide_evidence.config import path
from pesticide_evidence.names import CORE_GROUPS, canonical, class_of_name, to_canonical

REGULATORY_REASONS = {"10", "07", "21", "22", "23", "26", "27", "43", "44", "45", "47", "AS", "CS", "CR", "AR", "24", "35"}
MARKET_REASONS = {"01", "02", "05", "09", "32", "33", "34", "48", "49", "50", "11", "41", "42", "08", "03"}
MARKET_TEXT = re.compile(
    r"(did not (?:submit|provide)|not submitted|no longer (?:support|wish)|withdr[ae]w (?:its|their|the) (?:support|application|interest)|withdrawn (?:its|their) support|not (?:been )?supported|does not (?:wish|intend)|no (?:notifier|applicant)|did not (?:wish|apply)|declared that it (?:does|did) not|informed the commission that (?:it|they) (?:no longer|did not|do not))",
    re.I,
)
HAZARD_TEXT = re.compile(
    r"(concern|high risk|unacceptable|not (?:been )?(?:demonstrated|established|shown)|could not be (?:excluded|concluded)|risk assessment (?:could not|was not)|critical area|toxic|endocrine|genotox|carcinogen|persistent|bioaccumul|data gap|not satisf|failed to meet|do not satisfy)",
    re.I,
)


def eu_substances():
    """Records of the EU Pesticides Database active-substance file."""
    with open(path("raw", "eu", "eu_active_substances.json")) as fh:
        return [json.loads(line) for line in fh if line.strip()]


def us_registrations():
    """Per ingredient: US product counts, current products, first approval, last cancellation, and the share of recent cancellations that are regulatory or market-driven."""
    raw = path("raw", "epa")
    products = pd.read_fwf(
        raw / "product.txt",
        colspecs=[(0, 11), (11, 13), (13, 14), (14, 22), (22, 30), (30, 32), (32, 102)],
        names=["reg", "form", "tox", "appr", "can", "ct", "name"],
        dtype=str,
        header=None,
        encoding="latin-1",
    )
    formula = pd.read_fwf(raw / "formula.txt", colspecs=[(0, 11), (11, 17), (17, 24)], names=["reg", "pc", "pct"], dtype=str, header=None, encoding="latin-1")
    chem = pd.read_fwf(
        raw / "chemname.txt", colspecs=[(0, 6), (6, 26), (26, 296), (296, 297)], names=["pc", "type", "name", "disp"], dtype=str, header=None, encoding="latin-1"
    )
    chem["name"] = chem.name.str.strip()
    best = chem[chem.type.str.strip().isin(["ISO", "ANSI", "Other Common", "BSI"])].copy()
    best["rank"] = best.type.str.strip().map({"ISO": 0, "ANSI": 1, "BSI": 2, "Other Common": 3})
    best = best.sort_values(["pc", "rank"])
    pc_to_ai = {}
    for pc, names in best.groupby("pc").name.apply(list).items():
        for name in names:
            ai = to_canonical(name)
            if ai:
                pc_to_ai[pc] = ai
                break
    formula["ai"] = formula.pc.map(pc_to_ai)
    formula = formula.dropna(subset=["ai"])
    products["appr_y"] = pd.to_numeric(products.appr.str[:4], errors="coerce")
    products["can_y"] = pd.to_numeric(products.can.str[:4], errors="coerce")
    products["current"] = products.can.fillna("0").str.strip().isin(["00000000", "0", ""]) | (products.can_y.fillna(0) == 0)
    products.loc[products.current, "can_y"] = np.nan
    merged = formula.merge(products[["reg", "appr_y", "can_y", "current", "ct"]], on="reg")
    merged["reason"] = np.where(merged.ct.isin(REGULATORY_REASONS), "regulatory", np.where(merged.ct.isin(MARKET_REASONS), "market", "other"))
    rows = []
    for ai, g in merged.groupby("ai"):
        ended = g.loc[~g.current]
        latest = ended.sort_values("can_y").tail(max(1, int(0.2 * len(ended)))) if len(ended) else ended
        rows.append(
            dict(
                ai=ai,
                n_products=len(g),
                n_current=int(g.current.sum()),
                first_appr=g.appr_y.min(),
                last_cancel=g.can_y.max(),
                share_regulatory=(latest.reason == "regulatory").mean() if len(latest) else np.nan,
                share_market=(latest.reason == "market").mean() if len(latest) else np.nan,
            )
        )
    table = pd.DataFrame(rows).set_index("ai")
    table.to_csv(path("processed", "us_registrations.csv"))
    return table


def canada_registrations():
    """Per ingredient: Canadian product counts, current products, first registration and last expiry year."""
    try:
        products = pd.read_csv(path("raw", "pmra", "product.csv"), low_memory=False)
    except UnicodeDecodeError:
        products = pd.read_csv(path("raw", "pmra", "product.csv"), encoding="latin-1", low_memory=False)
    products["ai_text"] = products["Active ingredients - English"].astype(str).str.upper()
    products["current"] = products["Current / Historical"].eq("Current")
    products["exp"] = pd.to_datetime(products["Expiry date"], errors="coerce").dt.year
    products["first"] = pd.to_datetime(products["Date first registered"], errors="coerce").dt.year
    focal = sorted({canonical(n) for n, (group, _) in class_of_name().items() if group in CORE_GROUPS})
    rows = []
    for ai in focal:
        if len(ai) < 5 or ai in ("copper", "metam", "nicotine"):
            continue
        pattern = re.compile(r"(?<![A-Z])" + re.escape(ai.upper()) + r"(?![A-Z])")
        hit = products[products.ai_text.map(lambda s: bool(pattern.search(s)))]
        if len(hit):
            rows.append(dict(ai=ai, n_products=len(hit), n_current=int(hit.current.sum()), last_exp=hit.exp.max(), first_reg=hit["first"].min()))
    table = pd.DataFrame(rows).set_index("ai")
    table.to_csv(path("processed", "canada_registrations.csv"))
    return table


def who_hazard_class():
    """WHO hazard class of each ingredient, parsed from the 2019 classification text."""
    lines = open(path("raw", "who", "who_2019.txt")).read().splitlines()
    start = next(i for i, line in enumerate(lines) if "INDEX. Classification" in line or "COMMON NAME" in line and "CLASS" in line)
    token = re.compile(r"\b(Ia|la|Ib|lb|II|III|U|O|FM)\s+(\d{1,3}(?:, \d{1,3})?)(?=\s|$)")
    rows = []
    for line in lines[start:]:
        pos = 0
        for m in token.finditer(line):
            name = line[pos : m.start()].strip(" .")
            pos = m.end()
            if name and len(name) > 2 and not name.startswith(("INDEX", "COMMON")):
                rows.append((name, {"la": "Ia", "lb": "Ib"}.get(m.group(1), m.group(1))))
    who = pd.DataFrame(rows, columns=["name", "who_class"])
    who["ai"] = who.name.map(to_canonical)
    who = who.dropna(subset=["ai"]).drop_duplicates("ai")
    who.to_csv(path("processed", "who_hazard_class.csv"), index=False)
    return who


def cas_to_canonical():
    """Map from CAS number to canonical ingredient name for numbers that identify one name in the EU file."""
    by_cas = {}
    for rec in eu_substances():
        cas = (rec.get("as_cas_number") or "").strip()
        ai = to_canonical(rec["substance_name"])
        if cas and ai:
            by_cas.setdefault(cas, set()).add(ai)
    return {c: next(iter(v)) for c, v in by_cas.items() if len(v) == 1}


def use_atlas():
    """Applied mass in tonnes by country, ingredient and data origin from the use atlas."""
    atlas = read_atlas()
    atlas["kg"] = pd.to_numeric(atlas.kg, errors="coerce")
    cas_map = cas_to_canonical()
    names = atlas[["name", "cas"]].drop_duplicates()
    names["ai"] = [cas_map.get(str(c).strip()) or to_canonical(n) for n, c in zip(names.name, names.cas)]
    atlas = atlas.merge(names, on=["name", "cas"], how="left")
    table = atlas.dropna(subset=["ai"]).groupby(["country", "ai", "origin"]).kg.sum().div(1000).rename("tonnes").reset_index()
    table.to_csv(path("processed", "use_atlas_country_ingredient.csv"), index=False)
    return table


def celex_ids(label):
    """CELEX identifiers of the EU regulations cited in a register legislation label."""
    ids = []
    for m in re.finditer(r"(?:Reg\.?|Regulation)\s*\(E[UC]\)\s*(?:No\.?\s*)?(\d{1,4})/(\d{2,4})", label, re.I):
        a, b = m.group(1), m.group(2)
        if len(a) == 4:
            year, number = a, b
        elif len(b) == 4:
            year, number = b, a
        else:
            year, number = ("19" if int(b) > 50 else "20") + b, a
        ids.append(f"3{year}R{int(number):04d}")
    for m in re.finditer(r"(\d{2,4})/(\d{1,4})/(?:EC|EU|EEC)", label):
        year, number = m.group(1), m.group(2)
        if len(year) == 2:
            year = ("19" if int(year) > 50 else "20") + year
        ids.append(f"3{year}D{int(number):04d}")
    return ids


def act_text(celex):
    """Text of an EU act from the cache, or from the Publications Office when not cached."""
    cache = path("raw", "eu", "acts")
    cache.mkdir(parents=True, exist_ok=True)
    file = cache / f"{celex}.txt"
    if file.exists():
        return file.read_text()
    try:
        r = requests.get(f"http://publications.europa.eu/resource/celex/{celex}", headers={"Accept": "application/xhtml+xml", "Accept-Language": "eng"}, timeout=40)
        text = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", r.text)) if r.status_code == 200 else ""
    except requests.RequestException:
        text = ""
    file.write_text(text)
    time.sleep(0.3)
    return text


def eu_act_scan():
    """For each not-approved EU substance, counts of market-related and hazard-related terms in the text of its non-renewal act."""
    rows = []
    for rec in eu_substances():
        if rec["substance_status"] != "Not approved":
            continue
        ai = to_canonical(rec["substance_name"])
        label = re.sub("<[^>]+>", " ", str(rec.get("legislations_actives") or "")) + " " + re.sub("<[^>]+>", " ", str(rec.get("remark") or ""))
        ids = celex_ids(label)
        if not ai or not ids:
            continue
        for celex in dict.fromkeys(ids):
            text = act_text(celex)
            if not text:
                continue
            title = text[:400]
            rows.append(
                dict(
                    ai=ai,
                    celex=celex,
                    title=title[:200],
                    market=len(MARKET_TEXT.findall(text)),
                    hazard=len(HAZARD_TEXT.findall(text)),
                    nonren=bool(re.search(r"non-?renewal|not renew|non-?inclusion|not to renew|not approved|not be included", title, re.I)),
                    length=len(text),
                )
            )
    table = pd.DataFrame(rows)
    table.to_csv(path("processed", "eu_act_text_scan.csv"), index=False)
    return table


if __name__ == "__main__":
    for step in (us_registrations, canada_registrations, who_hazard_class, use_atlas, eu_act_scan):
        result = step()
        print(f"{step.__name__}: {len(result)} rows")
