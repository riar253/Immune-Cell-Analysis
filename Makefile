PYTHON = python3
VENV   = venv
BIN    = $(VENV)/bin

.DEFAULT_GOAL := all

.PHONY: all setup pipeline dashboard clean

all: setup pipeline dashboard

setup: ## Create a virtual environment and install dependencies
	$(PYTHON) -m venv $(VENV)
	$(BIN)/pip install --upgrade pip
	$(BIN)/pip install -r requirements.txt

pipeline: ## Initialize the DB, load data, generate tables + plots
	$(BIN)/python load_data.py
	$(BIN)/python analysis.py

dashboard: ## Start the local server for the interactive dashboard
	$(BIN)/streamlit run app.py

clean: ## Remove generated artifacts and the virtual environment
	rm -rf $(VENV) cell_count.db
	rm -rf output/
