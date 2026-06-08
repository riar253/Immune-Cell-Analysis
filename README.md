# Immune Cell Analysis — Loblaw Bio

Analysis of immune cell population frequencies from a clinical trial of the drug
candidate **miraclib**. The project loads per-sample cell counts into a normalized
SQLite database, computes relative frequencies, statistically compares treatment
responders against non-responders, profiles a baseline patient subset, and serves
the results through an interactive Streamlit dashboard.

## Quick start

The project is driven by a `Makefile`. From the repository root:

```bash
make setup       # create a virtual environment and install dependencies
make pipeline    # build the database and run all analyses (Parts 1–4)
make dashboard   # launch the interactive dashboard
make clean       # remove the venv and generated outputs
```

`make setup` creates a `venv/` and installs `requirements.txt`. Every target
invokes the virtual environment's interpreter directly, so you do **not** need to
activate the venv yourself. (Running `make` with no target is equivalent to
`make setup pipeline dashboard`.)

### Reproducing the outputs

`make pipeline` runs `load_data.py` then `analysis.py`. It prints the summary
tables to the terminal and writes all generated artifacts to `output/`:

| File | Contents |
|------|----------|
| `cell_frequencies.csv` | Part 2 — relative frequency of each population per sample |
| `responder_boxplot.png` | Part 3 — boxplot of frequencies, responders vs non-responders |
| `responder_stats.csv` | Part 3 — per-population statistical comparison |
| `baseline_samples.csv` | Part 4 — the identified baseline subset |
| `baseline_per_project.csv` | Part 4 — sample counts per project |
| `baseline_per_response.csv` | Part 4 — subject counts by response |
| `baseline_per_sex.csv` | Part 4 — subject counts by sex |

### Running the dashboard

`make dashboard` starts a local Streamlit server (default `http://localhost:8501`).
Run `make pipeline` first so the database exists.

**In GitHub Codespaces:** Streamlit's port is forwarded automatically — open the
URL from the "Ports" tab (or click the toast notification) to view the dashboard
in your browser.

<!-- Hosted dashboard (optional): if deployed to Streamlit Community Cloud, the
public link goes here, e.g. https://<your-app>.streamlit.app -->

## Database schema

The flat `cell-count.csv` mixes three levels of information that repeat
redundantly on every row: patient facts, sample facts, and the cell-count
measurements. The schema normalizes these into three tables.

```
subjects (1) ──< samples (1) ──< cell_counts
```

**`subjects`** — one row per patient.

| Column | Type | Notes |
|--------|------|-------|
| `subject_id` | TEXT PK | |
| `project` | TEXT | study/project the patient belongs to |
| `condition` | TEXT | indication, e.g. melanoma |
| `age` | INTEGER | |
| `sex` | TEXT | M / F |
| `treatment` | TEXT | e.g. miraclib |
| `response` | TEXT | yes / no / NULL |

**`samples`** — one row per biological sample.

| Column | Type | Notes |
|--------|------|-------|
| `sample_id` | TEXT PK | |
| `subject_id` | TEXT NOT NULL, FK → subjects | |
| `sample_type` | TEXT | e.g. PBMC |
| `time_from_treatment_start` | INTEGER | timepoint; baseline = 0 |

**`cell_counts`** — **long format**, one row per (sample, population).

| Column | Type | Notes |
|--------|------|-------|
| `sample_id` | TEXT NOT NULL, FK → samples | |
| `population` | TEXT NOT NULL | b_cell, cd8_t_cell, cd4_t_cell, nk_cell, monocyte |
| `count` | INTEGER NOT NULL | |
| | | composite PK = (`sample_id`, `population`) |

### Design rationale

- **Normalization (subjects / samples).** A patient contributes several samples
  across timepoints. Storing demographics, treatment, and response once per
  patient (rather than repeated on every sample row) removes redundancy and the
  risk of inconsistent duplicates.

- **Long-format cell counts.** Counts are stored as rows (`population`, `count`)
  rather than as five fixed columns. This was chosen deliberately because:
  - the analytical questions are all *per population* (relative frequency,
    per-population boxplots, per-population statistics), and the long shape
    matches that grain — Part 2's required output is a single straightforward
    aggregation query rather than an unpivot;
  - adding a new cell population later means inserting rows, not an
    `ALTER TABLE` plus query rewrites;
  - it generalizes to other per-sample measurements without schema change.

- **Constraints and indexes.** Foreign keys enforce referential integrity
  (`PRAGMA foreign_keys = ON`); the composite primary key on `cell_counts`
  prevents duplicate measurements. Indexes back the columns that filters and
  joins hit (`samples.subject_id`, `samples.sample_type`,
  `samples.time_from_treatment_start`, `cell_counts.population`). At the current
  data size these indexes are not yet performance-critical — they express the
  access patterns the design anticipates at scale.

- **Assumption.** `treatment` and `response` are modeled at the subject level
  because they are constant across each subject's samples in this dataset. If a
  trial recorded per-timepoint treatment, those fields would move to `samples`.

### Scaling considerations

The design is intended to hold up at hundreds of projects, thousands of samples,
and varied analytics:

- **Promote `project` to its own table.** It is currently a column on `subjects`.
  With hundreds of projects carrying their own attributes (sponsor, site, assay
  platform, dates), a `projects` table with `subjects.project_id` as a foreign key
  keeps project metadata normalized and queryable.
- **Long format scales naturally.** New cell populations or new assay readouts
  become new rows, never schema migrations — important when many projects measure
  overlapping-but-different panels.
- **Indexing / partitioning.** The existing indexes target the common filter and
  join columns. As volume grows, composite indexes tuned to frequent query shapes
  (e.g. `samples(sample_type, time_from_treatment_start)`) and, on a larger engine
  (Postgres), table partitioning by project would keep queries fast.
- **Migration path.** SQLite is ideal for a self-contained, reproducible
  submission. The schema is standard relational SQL and ports directly to
  Postgres/MySQL if the data outgrows a single file or needs concurrent writers.

## Code structure

```
├── load_data.py     # Part 1: build schema + load CSV into cell_count.db
├── analysis.py      # Parts 2–4: analysis functions + CLI pipeline (main)
├── app.py           # Streamlit dashboard (imports analysis.py)
├── Makefile         # setup / pipeline / dashboard / clean targets
├── requirements.txt
├── cell-count.csv   # input data
└── output/          # generated tables and plots
```

The code separates three concerns:

- **`load_data.py`** owns data ingestion only: it defines the schema, reshapes
  the CSV into the three tables (the wide-to-long `melt` for cell counts), and
  loads them. It is idempotent — it drops and recreates the tables on each run, so
  the pipeline can be re-run safely.

- **`analysis.py`** holds the analytical logic as small, focused functions
  (`cell_population_frequencies`, `responder_comparison_data`,
  `compare_responders`, `baseline_subset`, `baseline_breakdowns`, plotting). Each
  returns a tidy DataFrame or figure. Its `main()` runs the full pipeline and
  writes the `output/` artifacts.

- **`app.py`** is presentation only. It imports the functions from `analysis.py`
  and renders them in the dashboard, with query results cached via
  `st.cache_data`.

The key design choice is that **the analysis logic lives in one place**
(`analysis.py`) and is shared by both the command-line pipeline and the
dashboard. The two never compute results differently, and the dashboard can never
drift from the graded pipeline output.

## Methodology notes (Part 3)

Responders and non-responders are compared with the **Mann–Whitney U test** (a
non-parametric test that makes no normality assumption — appropriate for
relative-frequency data), and p-values are corrected for testing five populations
using the **Benjamini–Hochberg** false-discovery-rate procedure. Significance is
assessed on the adjusted p-values at α = 0.05.

In the provided data, no population shows a significant difference after
correction (cd4_t_cell is the strongest raw signal but does not survive FDR
adjustment). One caveat worth stating: each subject contributes multiple PBMC
samples across timepoints, so the samples are not fully independent observations;
the effective sample size is smaller than the raw count, and a repeated-measures
or per-timepoint model would be the rigorous next step.
