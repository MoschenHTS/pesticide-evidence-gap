PYTHON ?= python
export PYTHONPATH := src

.PHONY: all analysis acquire manifest

all: analysis

analysis:
	$(PYTHON) scripts/05_build_ingredient_table.py
	$(PYTHON) scripts/06_cross_sectional_analysis.py
	$(PYTHON) scripts/07_event_studies.py
	$(PYTHON) scripts/08_sensitivity_analysis.py
	$(PYTHON) scripts/09_brazil.py
	$(PYTHON) scripts/10_make_figures.py
	$(PYTHON) scripts/11_write_ledger.py

acquire:
	$(PYTHON) scripts/01_retrieve_pubmed.py
	$(PYTHON) scripts/02_classify_and_select.py --step classify
	$(PYTHON) scripts/02_classify_and_select.py --step select
	$(PYTHON) scripts/03_regulatory_and_use_data.py
	$(PYTHON) scripts/04_acquire_counts.py

manifest:
	cd $${PESTICIDE_DATA_DIR:-.} && find data -type f ! -name '.DS_Store' -exec shasum -a 256 {} \; | sort -k2 > data/MANIFEST.sha256
