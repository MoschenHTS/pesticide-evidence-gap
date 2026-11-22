"""Selection rules applied to the records that pass the first rule filter."""

import re

from pesticide_evidence import text_rules as K

TARGET = re.compile(K.TARGET_PEST_RE, re.I)
PARKINSON = re.compile(r"parkinson", re.I)
TOOL_COMPOUNDS = re.compile(r"rotenone|paraquat|mptp", re.I)
DRUG_ONLY = {"warfarin", "ivermectin", "niclosamide", "moxidectin", "nicotine"}
EXCLUDED_TYPES = {"Editorial", "Letter", "Comment", "News", "Biography", "Interview", "Published Erratum", "Newspaper Article"}
MIN_MENTIONS = 2  # a record must name the pesticide at least this many times in title and abstract

STEPS = ("no_ingredient", "drug_only", "single_mention", "publication_type", "tool_use", "efficacy", "pest_in_title", "residue")


def removal_step(title, abstract, publication_types, compounds):
    """Name of the first rule that removes the record, or None when the record is selected."""
    text = title + " " + abstract
    names = {x.strip().lower() for x in compounds.split(";") if x.strip()}
    non_target = bool(K.BENEFICIAL.search(text))
    checks = (
        ("no_ingredient", not names),
        ("drug_only", names <= DRUG_ONLY),
        ("single_mention", sum(text.lower().count(name) for name in names) < MIN_MENTIONS),
        ("publication_type", bool(EXCLUDED_TYPES & set(publication_types))),
        ("tool_use", bool(PARKINSON.search(text) and TOOL_COMPOUNDS.search(text))),
        ("efficacy", bool(K.EFFIC.search(text) and TARGET.search(text) and not non_target)),
        ("pest_in_title", bool(TARGET.search(title) and not non_target)),
        ("residue", bool(K.RESID.search(text) and not K.STRONG_TOX.search(text))),
    )
    return next((name for name, hit in checks if hit), None)
