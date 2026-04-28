PYTHON ?= python3
PIP := $(PYTHON) -m pip

.PHONY: setup pipeline dashboard

setup:
	$(PIP) install --upgrade pip
	$(PIP) install -r requirements.txt

pipeline:
	$(PYTHON) load_data.py
	$(PYTHON) run_pipeline.py

dashboard:
	$(PYTHON) -m streamlit run app.py --server.address 0.0.0.0 --server.port 8501 --server.headless true --browser.gatherUsageStats false
