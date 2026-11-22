# Data dictionary

## data/processed/ingredient_table.parquet

One row per canonical active ingredient (446 ingredients with reported use in the atlas). Built by `scripts/05_build_ingredient_table.py`.

| Column | Meaning |
|---|---|
| `South`, `North` | Applied mass (t) in Global South and Global North countries, atlas reported rows, 2016-2020 |
| `papers` | Selected toxicology papers per ingredient, 2005-2025 |
| `eu`, `us`, `ca` | EU approval status; number of current US and Canadian products |
| `eu_ok`, `us_ok`, `ca_ok` | Authorized in the EU (approved or pending), registered in the US, registered in Canada |
| `known` | EU status available |
| `north_any` | Authorized in at least one register |
| `age` | Years since the first PubMed record (minimum 1; 20 when unknown) |
| `crop_*` | Share of Southern applied mass by crop group |
| `target` | Modal primary target class in the atlas |
| `pm_all`, `pm_tox` | PubMed records mentioning the ingredient in title or abstract, all topics and with toxicity terms, 2005-2025 |
| `oa_zh_all`, `oa_zh_tox` | OpenAlex works mentioning the Chinese name, all topics and with toxicity terms |
| `cas` | CAS registry number |
| `eu_hazard`, `eu_market` | EU non-renewal act states a hazard or a market reason |
| `us_reg` | US cancellation carried a regulatory reason code and no current registration |
| `treaty` | Listed under the Rotterdam or Stockholm conventions (`config/treaty_pesticides.csv`) |
| `grp` | Authorized, or not authorized with a hazard or treaty basis, or not authorized for another reason |
| `en_zh_tox`, `en_zh_all` | PubMed plus Chinese-language OpenAlex counts: primary outcome (toxicity) and all topics |
| `us_ever`, `ca_ever`, `eu_dated` | Ever registered in the US or Canada; EU entry with an approval date |
| `grp4` | Analysis group: authorized, withdrawn (hazard or treaty), lapsed (other reason), absent from Northern registers |

## data/interim

| File | Content |
|---|---|
| `corpus/core_corpus.parquet` | One row per retrieved paper: ingredients, classes, model groups, outcome, method and design flags, exclusions, first-filter flag (`core`), selected-corpus flag (`selected`), first rule that removed the record (`selection_step`) |
| `corpus/panel_ingredient_year.parquet` | Selected papers per ingredient and year |
| `corpus/affiliations_raw.parquet`, `affiliations_countries.parquet` | First- and last-author affiliation text; parsed author country |
| `corpus/first_record.json` | Year of the first PubMed record per ingredient |
| `pubmed/name_counts.json`, `cas_registry_counts.json`, `window_counts.json`, `human_evidence_counts.json` | PubMed counts by name, CAS registry number, five-year window, evidence type |
| `openalex/chinese_counts.json`, `funder_counts.json` | OpenAlex Chinese-language counts and counts by funder |
| `crosscheck/independent_register_table.csv` | Independently compiled register table used for the cross-check |

## results

| File | Content |
|---|---|
| `numbers.json` | All cross-sectional numbers |
| `event_studies.json` | Event-study estimates with placebo and permutation calibration |
| `NUMBERS_LEDGER.md` | The numbers quoted in the manuscript |
