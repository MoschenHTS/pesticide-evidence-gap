"""Stage 01. PubMed retrieval, 2005-2025.

Runs the frozen query set (a generic pesticide block, MeSH supplementary-concept names in batches, and ingredient names by class in batches), each combined
with the animal-or-cell block and the toxicity block. Every query is split by date until all slices fit one result page and the retrieved identifiers
are checked against the PubMed count. The union of identifiers is fetched into the SQLite record cache.

Writes data/interim/pubmed/queries/*.json, data/interim/pubmed/pmid_sources.json and data/raw/pubmed/records.sqlite.
"""

import concurrent.futures as cf
import hashlib
import json
import re

from pesticide_evidence.config import CFG, ROOT, path
from pesticide_evidence.names import focal_classes
from pesticide_evidence.pubmed import esearch_all, fetch_records
from pesticide_evidence.text_rules import ANY_MODEL, FR, GEN_Q, TOX, batches, names_for_class, tiab

Y0, Y1 = CFG["sample"]["years"]
OUT = path("interim", "pubmed", "queries")
CORE = f"AND {ANY_MODEL} AND {TOX}"


def jobs():
    """List of (label, query) pairs: the generic block, MeSH concept batches and ingredient-name batches by class."""
    queries = [("generic", f"({GEN_Q} {CORE}){FR}")]
    concepts = sorted({x["label"] for x in json.load(open(ROOT / "config" / "mesh_pesticide_concepts.json"))})
    concepts = [s for s in concepts if len(s) < 45 and "(" not in s and re.search("[A-Za-z]{3}", s)]
    for i, batch in enumerate(batches(concepts, 12)):
        queries.append((f"scr_{i:02d}", "((" + " OR ".join(f'"{x}"[tiab] OR "{x}"[nm]' for x in batch) + f") {CORE}){FR}"))
    for group, cls in focal_classes():
        for i, batch in enumerate(batches(names_for_class(group, cls), 25)):
            queries.append((f"cls_{re.sub('[^A-Za-z]+', '_', cls)}_{i}", f"(({tiab(*batch)}) {CORE}){FR}"))
    return queries


def run(label, term):
    """Run one query, split by date, and cache its identifiers; returns the label and a status text."""
    file = OUT / f"{label}.json"
    digest = hashlib.md5(term.encode()).hexdigest()
    if file.exists() and json.loads(file.read_text()).get("hash") == digest:
        return label, "cached"
    ids, n, slices = esearch_all(term, Y0, Y1)
    file.write_text(json.dumps({"label": label, "hash": digest, "count": n, "n_ids": len(ids), "slices": slices, "term": term, "ids": ids}))
    return label, f"{n} ids, {slices} slices"


def union_of_queries():
    """Write each retrieved PMID with the labels of the queries that returned it."""
    sources = {}
    for file in OUT.glob("*.json"):
        record = json.loads(file.read_text())
        assert record["count"] == record["n_ids"] == len(record["ids"]), file
        for pmid in record["ids"]:
            sources.setdefault(pmid, []).append(record["label"])
    path("interim", "pubmed", "pmid_sources.json").write_text(json.dumps(sources))
    return sources


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    queries = jobs()
    print(len(queries), "queries", flush=True)
    with cf.ThreadPoolExecutor(4) as pool:
        for future in cf.as_completed([pool.submit(run, *q) for q in queries]):
            print(*future.result(), flush=True)
    sources = union_of_queries()
    print("union", len(sources))
    print(fetch_records(list(sources), log=lambda m: print(m, flush=True)))
