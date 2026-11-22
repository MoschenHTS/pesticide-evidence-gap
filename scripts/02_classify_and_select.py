"""Stage 02. Rule-based classification and selection of the toxicology corpus.

Steps (run one with `--step`):
    classify    rule-based classification of every retrieved record: ingredients, classes, model groups, outcome, method and design flags, exclusions
    select      selection rules applied to the records that pass the first rule filter; marks the selected corpus and the rule that removed each other record

Writes data/interim/corpus/core_corpus.parquet.
"""

import argparse
import collections
import json
import re

import pandas as pd

from pesticide_evidence import text_rules as K
from pesticide_evidence.config import ROOT, path
from pesticide_evidence.names import AMBIG, CASE, CLASS_WORDS_AMBIGUOUS, CLASS_WORDS_SAFE, NAMES, canonical, class_of_name
from pesticide_evidence.pubmed import load_records
from pesticide_evidence.selection import removal_step

CORPUS = path("interim", "corpus", "core_corpus.parquet")
C_OF = class_of_name()
CLASS_GROUP = {c: g for (g, c) in NAMES}
# biocide and fumigant names are recorded as classes but not counted as pesticide ingredient hits
BIOCIDE_GROUPS = ["Biocide", "Fumigant"]
PA_MESH = {
    "Pesticides",
    "Herbicides",
    "Insecticides",
    "Fungicides, Industrial",
    "Acaricides",
    "Rodenticides",
    "Molluscacides",
    "Nematocides",
    "Insect Repellents",
    "Pesticide Residues",
}
MESH_MODEL = {
    "Zebrafish": "Zebrafish",
    "Fishes": "Other fish",
    "Rats": "Rodents",
    "Mice": "Rodents",
    "Rats, Wistar": "Rodents",
    "Rats, Sprague-Dawley": "Rodents",
    "Mice, Inbred C57BL": "Rodents",
    "Xenopus laevis": "Amphibians",
    "Amphibians": "Amphibians",
    "Anura": "Amphibians",
    "Larva": None,
    "Bees": "Bees and pollinators",
    "Daphnia": "Aquatic invertebrates",
    "Crustacea": "Aquatic invertebrates",
    "Birds": "Birds",
    "Chickens": "Birds",
    "Quail": "Birds",
    "Oligochaeta": "Soil and other invertebrates",
    "Caenorhabditis elegans": "Soil and other invertebrates",
    "Drosophila melanogaster": "Soil and other invertebrates",
    "Drosophila": "Soil and other invertebrates",
    "Cell Line": "Cells",
    "Cells, Cultured": "Cells",
    "Cell Line, Tumor": "Cells",
    "Hep G2 Cells": "Cells",
    "Rabbits": "Other mammals",
    "Swine": "Other mammals",
    "Cattle": "Other mammals",
    "Sheep": "Other mammals",
    "Dogs": "Other mammals",
}
OTHER_FISH = re.compile(r"carp|trout|medaka|tilapia|salmon|catfish|goldfish|minnow|killifish|stickleback|sturgeon|oryzias|bass|guppy|loach|croaker|flounder|bream")


def _alternation(words, case_sensitive):
    """Regular-expression alternation of the words, longest first, with word-boundary and plural handling."""
    parts = []
    for w in sorted(set(words), key=len, reverse=True):
        escaped = re.escape(w)
        if w.endswith("-"):
            parts.append(escaped)
        elif len(w) <= 5:
            parts.append(escaped + r"(?:s|es)?(?![A-Za-z])")
        else:
            parts.append(escaped)
    return re.compile(r"(?<![A-Za-z0-9])(" + "|".join(parts) + ")", 0 if case_sensitive else re.I)


RX_CASE = _alternation([n for n in C_OF if n in CASE], True)
NOCASE_NAMES = [n for n in C_OF if n not in CASE]
RX_NOCASE = _alternation(NOCASE_NAMES, False)
LOWER_TO_NAME = {n.lower(): n for n in NOCASE_NAMES}
CLASS_WORD = {}
for table, ambiguous in ((CLASS_WORDS_SAFE, False), (CLASS_WORDS_AMBIGUOUS, True)):
    for cls, words in table.items():
        for w in words:
            CLASS_WORD[w.lower()] = (cls, ambiguous)
RX_CLASS_WORD = _alternation([w for w in CLASS_WORD if w != "syp-"] + ["SYP-"], False)
MESH_CONCEPTS = {x["label"].lower() for x in json.load(open(ROOT / "config" / "mesh_pesticide_concepts.json"))}
MODEL_RX = {k: re.compile(v, re.I) for k, v in K.MODEL_RE.items()}
FLAG_RX = {
    prefix: {k: re.compile(p, re.I) for k, p in table.items()} for prefix, table in (("o_", K.OUTCOME_RE), ("m_", K.METHOD_RE), ("d_", K.DESIGN_RE), ("x_", K.EXTRA_RE))
}
HUMAN = re.compile(K.HUMAN_RE, re.I)
TARGET = re.compile(K.TARGET_PEST_RE, re.I)


def exclusions(text, title, models, named, generic, target_pest, year, retracted, mesh):
    """Names of the exclusion rules that a record meets (delivery, synthesis, tool use, non-research and others)."""
    found = []
    pesticide = named or generic
    if pesticide and K.DELIV.search(text):
        found.append("delivery_detection")
    if pesticide and K.SYNTH.search(text):
        found.append("synthesis")
    if pesticide and K.PD_TOOL_T.search(text) and K.PD_TOOL_C.search(text) and not K.PEST_CTX.search(title):
        found.append("pd_tool")
    if (
        pesticide
        and K.EFFIC.search(text)
        and not K.BENEFICIAL.search(text)
        and (not models or target_pest and not (set(models) - {"Soil and other invertebrates", "Aquatic invertebrates"}))
    ):
        found.append("efficacy_resistance")
    if pesticide and K.RESID.search(title) and not K.STRONG_TOX.search(text):
        found.append("residues_analytical")
    core = bool(pesticide and (models or "Animals" in mesh) and not found and year and 2005 <= year <= 2025 and not retracted)
    return ";".join(found), core


def classify_record(r):
    """Rule-based flags for one PubMed record: pesticide context, named ingredients and classes, model groups, outcomes, methods and design."""
    title, abstract = r["title"] or "", r["abstract"] or ""
    text = " ".join([title, abstract, " ".join(r["kw"]), " ".join(r["chem"])])
    lower = text.lower()
    generic = (
        bool(K.PEST_CTX.search(text))
        or any(c in PA_MESH for c in r["mesh"])
        or any(c in PA_MESH for c in r["chem"])
        or any(c.lower() in MESH_CONCEPTS for c in r["chem"])
    )
    ingredient_hits, class_hits, biocide = collections.Counter(), set(), set()
    for rx in (RX_CASE, RX_NOCASE):
        for m in rx.finditer(text):
            token = m.group(1)
            if token in C_OF:
                name = token
            else:
                low = token.lower()
                name = LOWER_TO_NAME.get(low) or LOWER_TO_NAME.get(re.sub(r"(s|es)$", "", low))
            if name is None or name not in C_OF:
                continue
            if name in AMBIG and not generic:
                continue
            group, cls = C_OF[name]
            if group in BIOCIDE_GROUPS:
                biocide.add(cls)
                continue
            ingredient_hits[canonical(name)] += 1
            class_hits.add(cls)
    word_classes = set()
    for m in RX_CLASS_WORD.finditer(text):
        word = m.group(1).lower()
        cls, ambiguous = CLASS_WORD.get(word, (None, False))
        if cls is None and word.startswith("syp-"):
            cls, ambiguous = "Other fungicides", False
        if cls and (not ambiguous or generic):
            word_classes.add(cls)
    classes = sorted(class_hits | word_classes)
    groups = sorted({CLASS_GROUP[c] for c in classes})
    named = bool(classes)
    models = {k for k, rx in MODEL_RX.items() if rx.search(text)}
    for term in r["mesh"]:
        mapped = MESH_MODEL.get(term)
        if mapped:
            models.add(mapped)
    if "Other fish" in models and "Zebrafish" in models and not OTHER_FISH.search(lower):
        models.discard("Other fish")
    counts = {k: len(MODEL_RX[k].findall(text)) for k in models if k in MODEL_RX}
    primary = max(sorted(counts), key=counts.get) if counts else None
    out = dict(
        pmid=r["pmid"],
        year=r["year"],
        journal=r["journal"],
        title=title,
        doi=r["doi"],
        n_auth=len(r["authors"]),
        generic=generic,
        named=named,
        classes=";".join(classes),
        groups=";".join(groups),
        compounds=";".join(sorted(ingredient_hits)),
        n_ai=len(ingredient_hits),
        top_ai=(max(sorted(ingredient_hits), key=ingredient_hits.get) if ingredient_hits else None),
        class_word_only=bool(classes) and not ingredient_hits,
        biocide_cls=";".join(sorted(biocide)),
        models=";".join(sorted(models)),
        primary_model=primary,
        model_any=bool(models),
        animals_mesh="Animals" in r["mesh"],
        human_ctx=bool(HUMAN.search(text)),
        human_mesh="Humans" in r["mesh"],
        target_pest=bool(TARGET.search(text)),
        review=any("Review" in t for t in r["ptype"]),
        retracted=any("Retract" in t for t in r["ptype"]),
        lang=r["lang"],
        has_abstract=bool(abstract),
    )
    for prefix, table in FLAG_RX.items():
        for k, rx in table.items():
            out[prefix + k] = bool(rx.search(text))
    out["excl"], out["core"] = exclusions(text, title, models, named, generic, out["target_pest"], r["year"], out["retracted"], r["mesh"])
    out["tier"] = "named" if named else ("generic" if generic else "none")
    return out


def classify_corpus(out=CORPUS):
    """Classify every cached record, attach the queries that retrieved it, and write the corpus table."""
    records = load_records()
    df = pd.DataFrame([classify_record(r) for r in records.values()])
    sources = json.load(open(path("interim", "pubmed", "pmid_sources.json")))
    df["n_sources"] = df.pmid.map(lambda p: len(sources.get(p, [])))
    df["src_generic"] = df.pmid.map(lambda p: "generic" in sources.get(p, []))
    df["src_scr"] = df.pmid.map(lambda p: any(s.startswith("scr_") for s in sources.get(p, [])))
    df["src_names"] = df.pmid.map(lambda p: any(s.startswith("cls_") for s in sources.get(p, [])))
    df.to_parquet(out, index=False)
    return df


def select_corpus(out=CORPUS):
    """Apply the selection rules to the records that pass the first filter and mark the selected corpus and the rule that removed each other record."""
    df = pd.read_parquet(out)
    for col in ("selected", "selection_step"):
        if col in df.columns:
            df = df.drop(columns=col)
    core = df[df.core]
    records = load_records(list(core.pmid))
    steps = {
        row.pmid: removal_step(records[row.pmid]["title"] or "", records[row.pmid]["abstract"] or "", records[row.pmid]["ptype"], row.compounds or "")
        for row in core.itertuples()
    }
    df["selection_step"] = df.pmid.map(steps)
    df["selected"] = df.core & df.selection_step.isna()
    df.to_parquet(out, index=False)
    return df


STEPS = {"classify": classify_corpus, "select": select_corpus}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--step", choices=sorted(STEPS), required=True)
    args = parser.parse_args()
    STEPS[args.step]()
    print(f"{args.step}: done")
