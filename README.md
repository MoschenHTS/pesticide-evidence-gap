# Pesticide evidence gap — bibliometric pipeline for toxicity research by register status

Code, configuration and result summaries supporting the paper *"Pesticide science studies what Northern registers authorize, not what the South sprays."*

This repository covers the full analysis: PubMed and OpenAlex retrieval, rule-based classification and selection of toxicology papers, the ingredient-level table that joins use, register status and literature counts, the cross-sectional gap model, the event studies, the sensitivity analyses, the Brazilian subsample, and the figures and tables. A pesticide here is an active ingredient (the chemical substance, not the commercial product). Third-party data are not redistributed.

## Layout

```
config/            config.yaml (parameters, seeds, country sets), treaty_pesticides.csv, ingredient_names.csv,
                   funder_classes_checked.csv, mesh_pesticide_concepts.json
src/pesticide_evidence/
                   config.py  stats.py  atlas.py  names.py  geography.py  registers.py  funders.py  pubmed.py  text_rules.py  selection.py
scripts/           pipeline stages, run in numeric order (below)
docs/              data_sources.md (inputs, sources, terms), data_dictionary.md (columns and files)
results/           final numbers and tables (below); figures are not included
data/              not in the repository; see "Primary data (not included)"
```

## Pipeline stages (`scripts/`)

```
01_retrieve_pubmed            PubMed retrieval, 2005-2025
02_classify_and_select        rule-based classification, then selection of the toxicology corpus (--step classify | select)
03_regulatory_and_use_data    EU, US and Canadian registers, treaty lists, use atlas
04_acquire_counts             PubMed and OpenAlex counts per ingredient
05_build_ingredient_table     ingredient-level analysis table and corpus panels
06_cross_sectional_analysis   gap model and every cross-sectional number
07_event_studies              stacked event study and EU review timing, with placebo calibration
08_sensitivity_analysis       alternative samples, counts, periods and model forms
09_brazil                     gap model on the Brazilian subsample and the five best-selling pesticides
10_make_figures               figures, Table 1 and the supplementary tables
11_write_ledger               collects the numbers of stages 06-09 into one list
```

Scripts are meant to be run from the repository root:

```bash
make acquire     # stages 01-04, cached and resumable
make analysis    # stages 05-11
python scripts/04_acquire_counts.py --only name_counts --limit 5      # test a single step
```

## Setup

```bash
mamba env create -f environment.yml
mamba activate pesticide-evidence
pip install -e .
cp .env.example .env      # OPENALEX_API_KEY and OPENALEX_MAILTO, needed only for stage 04
export PESTICIDE_DATA_DIR=/path/to/folder/containing/data
```

## Results (`results/`)

Final, human-readable numbers and tables, enough to check every reported number without re-running anything. Figures are not included; the published figures are in the paper, and stage 10 regenerates them.

| File | Contents |
|---|---|
| `NUMBERS_LEDGER.md` | Every number quoted in the paper, with the file it comes from |
| `numbers.json` | Cross-sectional numbers (stage 06) |
| `event_studies.json` | Event-study estimates with placebo and permutation calibration (stage 07) |
| `sensitivity.json` | Sensitivity analyses (stage 08, Table S1) |
| `brazil.json` | Brazilian subsample (stage 09) |
| `tables/Table_1.csv`, `Table_2.csv` | Sample by authorization group; the five best-selling pesticides in Brazil |
| `tables/Table_S7.csv` | CAS number, group, register status and basis of removal for each of the 295 pesticides |
| `tables/ingredient_sample.csv` | The 295-ingredient analysis sample |
| `supplement/Supplementary_Tables.xlsx` | Supplementary tables, one per sheet |

## Primary data (not included)

The third-party files are not redistributed. `docs/data_sources.md` gives each source, its access date, its terms and the file names that the scripts expect under `data/raw/`. PubMed and OpenAlex are queried through their public interfaces and add records over time, so a later run can return counts that differ from the ones reported; the reported counts are fixed in `results/numbers.json`. The checksums of the files used are written by `make manifest` to `data/MANIFEST.sha256`.

Scripts read the data from the folder named by `PESTICIDE_DATA_DIR`, which contains `data/` and can sit outside the repository.

## Third-party data

The inputs listed in `docs/data_sources.md` keep their original terms rather than this repository's license.

## License

Code and original data in this repository are released under the [MIT License](LICENSE), except for the third-party inputs above.

## Citation

If you use this code, please cite the paper (full citation to be added on publication). For citing this repository specifically, see [CITATION.cff](CITATION.cff).
