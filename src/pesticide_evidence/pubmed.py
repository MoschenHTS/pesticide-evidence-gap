"""PubMed retrieval helpers: query builders, throttled E-utilities calls, complete ID retrieval and a SQLite record cache."""

import json
import re
import sqlite3
import threading
import time
import datetime as dt
import xml.etree.ElementTree as ET

import requests

from .config import path

EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
S = requests.Session()
PAUSE = 0.36  # seconds between requests, below the E-utilities limit of 3 per second without an API key
MAXPAGE = 9000  # a date range returning more records than this is split in two before fetching


def tiab(*terms):
    """PubMed title-or-abstract query for the terms, quoting those that contain spaces or punctuation."""
    return "(" + " OR ".join(f'"{t}"[tiab]' if re.search(r"[ \-/,.()]", t) else f"{t}[tiab]" for t in terms) + ")"


def AND(*xs):
    """Parenthesized AND of query parts."""
    return "(" + " AND ".join(xs) + ")"


def OR(*xs):
    """Parenthesized OR of query parts."""
    return "(" + " OR ".join(xs) + ")"


_lock = threading.Lock()
_last = [0.0]


def _throttle():
    """Wait so that requests stay within the request rate."""
    with _lock:
        wait = _last[0] + PAUSE - time.time()
        if wait > 0:
            time.sleep(wait)
        _last[0] = time.time()


def _post(endpoint, data, tries=7, timeout=180, **kw):
    """POST to an E-utilities endpoint with retries and back-off."""
    last = None
    for a in range(tries):
        try:
            _throttle()
            r = S.post(EU + endpoint, data=data, timeout=timeout, **kw)
            if r.status_code == 429:
                raise RuntimeError("429")
            r.raise_for_status()
            return r
        except Exception as e:
            last = e
            time.sleep(min(60, 2**a))
    raise RuntimeError(f"POST {endpoint} failed after {tries} tries: {last}")


def esearch_count(term, d0=None, d1=None):
    """Number of PubMed records matching a query, optionally within publication dates."""
    p = dict(db="pubmed", term=term, retmode="json", retmax=0)
    if d0:
        p.update(datetype="pdat", mindate=d0, maxdate=d1)
    js = _post("esearch.fcgi", p).json()["esearchresult"]
    if "count" not in js:
        raise RuntimeError("esearch error: " + str(js)[:300])
    return int(js["count"])


def _esearch_ids_once(term, d0, d1, n):
    """Identifiers returned by one search request for a date range."""
    p = dict(db="pubmed", term=term, retmode="json", retmax=n, datetype="pdat", mindate=d0, maxdate=d1)
    js = _post("esearch.fcgi", p).json()["esearchresult"]
    return js.get("idlist", [])


def esearch_all(term, y0=2005, y1=2025):
    """All PMIDs for `term` with publication date in [y0, y1]. Returns (ids, count_total, n_slices); date windows are split until every slice fits one page."""
    ids, slices = [], 0
    whole = (f"{y0}/01/01", f"{y1}/12/31")
    stack = [whole]
    full = None
    while stack:
        d0, d1 = stack.pop()
        n = esearch_count(term, d0, d1)
        if (d0, d1) == whole:
            full = n
        if n == 0:
            continue
        if n > MAXPAGE:
            a = dt.date.fromisoformat(d0.replace("/", "-"))
            b = dt.date.fromisoformat(d1.replace("/", "-"))
            if a == b:
                raise RuntimeError(f"single day {d0} has {n} hits for {term[:80]}")
            if (b - a).days > 366:
                stack.extend((f"{y}/01/01", f"{y}/12/31") for y in range(a.year, b.year + 1))
            else:
                m = a + (b - a) / 2
                stack.append((d0, m.strftime("%Y/%m/%d")))
                stack.append(((m + dt.timedelta(days=1)).strftime("%Y/%m/%d"), d1))
            continue
        got = _esearch_ids_once(term, d0, d1, n)
        assert len(got) == n, f"truncation: slice {d0}-{d1} count {n} but {len(got)} IDs"
        ids.extend(got)
        slices += 1
    uniq = list(dict.fromkeys(ids))
    assert len(uniq) == full, f"slice mismatch: unique IDs {len(uniq)} vs whole-range count {full}"
    return uniq, full, slices


def _txt(e):
    """Text content of an XML element, or an empty string."""
    return "".join(e.itertext()).strip() if e is not None else ""


def parse_article(a):
    """Record fields (identifiers, year, journal, title, abstract, MeSH, substances, keywords, authors, language) from one PubMed article element."""
    pm = a.findtext(".//MedlineCitation/PMID")
    pd = a.find(".//Article/Journal/JournalIssue/PubDate")
    y = None
    if pd is not None:
        y = pd.findtext("Year") or (pd.findtext("MedlineDate") or "")[:4]
    ad = a.find(".//Article/ArticleDate")
    if ad is not None and ad.findtext("Year"):
        y = ad.findtext("Year")
    auth = []
    for x in a.findall(".//AuthorList/Author"):
        ln = x.findtext("LastName")
        ini = x.findtext("Initials")
        if ln:
            auth.append({"n": f"{ln} {ini or ''}".strip(), "aff": " | ".join(_txt(z) for z in x.findall(".//AffiliationInfo/Affiliation"))})
    ids = {x.get("IdType"): _txt(x) for x in a.findall(".//PubmedData/ArticleIdList/ArticleId")}
    return dict(
        pmid=pm,
        year=int(y) if y and y[:4].isdigit() else None,
        journal=a.findtext(".//Journal/ISOAbbreviation") or a.findtext(".//Journal/Title"),
        title=_txt(a.find(".//ArticleTitle")),
        abstract=" ".join(_txt(p) for p in a.findall(".//Abstract/AbstractText")),
        mesh=[_txt(m) for m in a.findall(".//MeshHeading/DescriptorName")],
        chem=[_txt(c) for c in a.findall(".//ChemicalList/Chemical/NameOfSubstance")],
        kw=[_txt(k) for k in a.findall(".//KeywordList/Keyword")],
        ptype=[_txt(p) for p in a.findall(".//PublicationType")],
        doi=ids.get("doi"),
        pmc=ids.get("pmc"),
        authors=auth,
        lang=a.findtext(".//Language"),
    )


def _db(db_path=None):
    """Connection to the SQLite record cache, created when absent."""
    con = sqlite3.connect(str(db_path or path("raw", "pubmed", "records.sqlite")))
    con.execute("CREATE TABLE IF NOT EXISTS records (pmid TEXT PRIMARY KEY, js TEXT)")
    con.execute("CREATE TABLE IF NOT EXISTS missing (pmid TEXT PRIMARY KEY)")
    return con


def _efetch_batch(b):
    """Fetch and parse one batch of records; returns the batch and the records found."""
    r = _post("efetch.fcgi", {"db": "pubmed", "id": ",".join(b), "retmode": "xml"}, timeout=240)
    root = ET.fromstring(r.content)
    got = {}
    for a in root.findall(".//PubmedArticle"):
        rec = parse_article(a)
        got[rec["pmid"]] = rec
    return b, got


def fetch_records(pmids, db_path=None, batch=200, workers=3, log=print):
    """Fetch PubMed records into the SQLite cache (resumable, parallel); PMIDs that PubMed does not return are stored in the `missing` table."""
    import concurrent.futures as cf

    con = _db(db_path)
    have = {r[0] for r in con.execute("SELECT pmid FROM records")} | {r[0] for r in con.execute("SELECT pmid FROM missing")}
    todo = [p for p in dict.fromkeys(pmids) if p not in have]
    log(f"fetch: {len(pmids)} requested, {len(todo)} to fetch")
    t0 = time.time()
    chunks = [todo[i : i + batch] for i in range(0, len(todo), batch)]
    done = 0
    with cf.ThreadPoolExecutor(workers) as ex:
        for k, fut in enumerate(cf.as_completed([ex.submit(_efetch_batch, c) for c in chunks])):
            b, got = fut.result()
            for p, rec in got.items():
                con.execute("INSERT OR REPLACE INTO records VALUES (?,?)", (p, json.dumps(rec)))
            for p in set(b) - set(got):
                con.execute("INSERT OR REPLACE INTO missing VALUES (?)", (p,))
            con.commit()
            done += len(b)
            if k % 25 == 0:
                log(f"  {done}/{len(todo)}  {time.time() - t0:.0f}s")
    n_rec = con.execute("SELECT COUNT(*) FROM records").fetchone()[0]
    n_mis = con.execute("SELECT COUNT(*) FROM missing").fetchone()[0]
    con.close()
    return n_rec, n_mis


def load_records(pmids=None, db_path=None):
    """Cached records as a dictionary keyed by PMID, for all records or the given identifiers."""
    con = _db(db_path)
    if pmids is None:
        rows = con.execute("SELECT pmid, js FROM records").fetchall()
    else:
        pm = list(pmids)
        rows = []
        for i in range(0, len(pm), 900):
            q = ",".join("?" * len(pm[i : i + 900]))
            rows += con.execute(f"SELECT pmid, js FROM records WHERE pmid IN ({q})", pm[i : i + 900]).fetchall()
    con.close()
    return {p: json.loads(j) for p, j in rows}
