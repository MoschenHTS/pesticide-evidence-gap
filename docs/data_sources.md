# Data sources

| Data | Source | Accessed | Licence or terms |
|---|---|---|---|
| Applied mass per country, ingredient and crop, 2016-2020 (99 countries, 1,077 ingredients) | Global atlas of pesticide use, Research Square preprint rs-10665000, supplementary data files 1 and 2 | 2026-10-25 | See the preprint; the raw files are not redistributed here |
| EU approval status and dates | EU Pesticides Database, active substances (API v3) | 2026-10-14 | European Commission reuse policy |
| Text of EU non-renewal acts | EUR-Lex / Publications Office, Cellar (CELEX) | 2026-10-14 | EU reuse policy |
| US registrations | EPA Pesticide Product Information System (FIFRA section 3), www3.epa.gov/pesticides/PPISdata | 2026-10-25 | US government work |
| Canadian registrations | Health Canada PMRA product register (current and historical) | 2026-10-25 | Open Government Licence - Canada |
| WHO hazard class | WHO Recommended Classification of Pesticides by Hazard, 2019 edition | 2026-10-14 | WHO terms |
| Countries, regions, income groups | World Bank API | 2026-10-25 | CC BY 4.0 |
| National and regional pesticide use totals | FAOSTAT, Pesticides Use (1990-2022) | 2026-10-13 | CC BY-NC-SA 3.0 IGO |
| Treaty-listed pesticides | Rotterdam Convention Annex III; Stockholm Convention Annexes A and B | 2026-10-25 | Public |
| Total pesticide use by country, 1995-2019, corrected for FAO reporting gaps (GloPUT) | Shattuck et al. 2023, Open Science Framework project dyu38, `GloPUT_1995-2019.csv` | 2026-10-26 | See the OSF project |
| Literature records and counts | PubMed E-utilities | 2026-10-13 and 2026-10-25/26 | NLM terms |
| Literature counts by language; funder counts | OpenAlex API | 2026-10-13 and 2026-10-25/26 | CC0 |
| MeSH supplementary concept names for pesticides | NLM MeSH (`config/mesh_pesticide_concepts.json`) | 2026-10-11 | NLM terms |
| Chinese ingredient names | Wikidata labels matched by CAS number, completed by hand (`config/ingredient_names.csv`) | 2026-10-13 | CC0 |
| Independent register table used for the cross-check | Compiled 2026-10-02 from the EU Pesticides Database API v3, EPA PPIS and PMRA files | 2026-10-02 | Derived table |

## Country boundaries (Figure 1)
Natural Earth 1:110m Admin 0 countries, public domain, `https://naciscdn.org/naturalearth/110m/cultural/ne_110m_admin_0_countries.zip`, stored in `data/raw/geography/`.

## Files expected in `data/raw/`
No third-party file is redistributed in this repository. Download each source from the address in the table above and save it under the name below, inside the folder named by `PESTICIDE_DATA_DIR`.

| Path under `data/raw/` | Content | Created by |
|---|---|---|
| `atlas/sdata1.xlsx`, `atlas/sdata2.xlsx` | Supplementary data files 1 and 2 of the global atlas of pesticide use | download |
| `eu/eu_active_substances.json` | EU Pesticides Database, active substances (API v3) | download |
| `eu/acts/<CELEX>.txt` | Text of the EU non-renewal acts | stage 03 |
| `epa/product.txt`, `epa/chemname.txt`, `epa/formula.txt` | EPA PPIS tables | download |
| `pmra/product.csv` | Health Canada PMRA product register | download |
| `who/who_2019.txt` | WHO Recommended Classification of Pesticides by Hazard, 2019, as plain text | download |
| `fao/fao_pesticides_use.zip` | FAOSTAT Pesticides Use bulk file | download |
| `glopu/GloPUT_1995-2019.csv` | GloPUT total use by country | download |
| `worldbank/wb_countries.json` | World Bank country API response | download |
| `geography/ne_110m_admin_0_countries.zip` | Natural Earth country boundaries | download |
| `pubmed/records.sqlite` | PubMed records retrieved by the queries in `data/interim/pubmed/queries/` | stage 01 |

Stages 01 and 04 query PubMed and OpenAlex, which add records over time, so a later run returns counts that can differ from those reported in the paper. The counts used in the paper are in `results/numbers.json` and `results/NUMBERS_LEDGER.md`.

## CAS registry numbers in `config/ingredient_names.csv`
The number is the one under which PubMed indexes the named active ingredient, checked against PubChem. For twelve ingredients whose ISO name denotes an ester (metsulfuron-methyl, fenoxaprop-ethyl, quizalofop-ethyl, bensulfuron-methyl, chlorimuron-ethyl, thifensulfuron-methyl, diclofop-methyl, fluazifop-butyl, tribenuron-methyl, clodinafop-propargyl, mesosulfuron-methyl, thiencarbazone-methyl) the number of that ester is kept; for the others it is the parent compound. The PubMed window is 2005 to 2025.
