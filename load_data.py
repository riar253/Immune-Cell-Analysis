"""
Initializes the SQLite database and loads cell-count.csv into a normalized schema.
Creates `cell_count.db` in the repository root with three tables:
  subjects     - one row per patient (demographics, treatment, response)
  samples      - one row per biological sample (timepoint, sample type)
  cell_counts  - long format: one row per (sample, cell population) measurement
"""

import sqlite3
import pandas as pd

DB_PATH = "cell_count.db"
CSV_PATH = "cell-count.csv"

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

# Create the relational schema, dropping any existing tables first
def create_schema(conn):
    cursor = conn.cursor()
    cursor.execute("PRAGMA foreign_keys = ON;")

    cursor.executescript(
        """
        DROP TABLE IF EXISTS cell_counts;
        DROP TABLE IF EXISTS samples;
        DROP TABLE IF EXISTS subjects;

        CREATE TABLE subjects (
            subject_id  TEXT PRIMARY KEY,
            project     TEXT,
            condition   TEXT,
            age         INTEGER,
            sex         TEXT,
            treatment   TEXT,
            response    TEXT
        );

        CREATE TABLE samples (
            sample_id                  TEXT PRIMARY KEY,
            subject_id                 TEXT NOT NULL,
            sample_type                TEXT,
            time_from_treatment_start  INTEGER,
            FOREIGN KEY (subject_id) REFERENCES subjects (subject_id)
        );

        CREATE TABLE cell_counts (
            sample_id   TEXT NOT NULL,
            population  TEXT NOT NULL,
            count       INTEGER NOT NULL,
            PRIMARY KEY (sample_id, population),
            FOREIGN KEY (sample_id) REFERENCES samples (sample_id)
        );

        CREATE INDEX idx_samples_subject     ON samples (subject_id);
        CREATE INDEX idx_samples_type        ON samples (sample_type);
        CREATE INDEX idx_samples_timepoint   ON samples (time_from_treatment_start);
        CREATE INDEX idx_cell_counts_pop     ON cell_counts (population);
        """
    )
    conn.commit()

# Split the flat CSV dataframe into the three normalized frames
def build_tables(df):
    subjects = (
        df[["subject", "project", "condition", "age", "sex", "treatment", "response"]]
        .drop_duplicates(subset="subject")
        .rename(columns={"subject": "subject_id"})
    )

    samples = df[
        ["sample", "subject", "sample_type", "time_from_treatment_start"]
    ].rename(columns={"sample": "sample_id", "subject": "subject_id"})

    cell_counts = df.melt(
        id_vars="sample",
        value_vars=POPULATIONS,
        var_name="population",
        value_name="count",
    ).rename(columns={"sample": "sample_id"})

    return subjects, samples, cell_counts


def main():
    df = pd.read_csv(CSV_PATH)

    with sqlite3.connect(DB_PATH) as conn:
        create_schema(conn)
        subjects, samples, cell_counts = build_tables(df)

        subjects.to_sql("subjects", conn, if_exists="append", index=False)
        samples.to_sql("samples", conn, if_exists="append", index=False)
        cell_counts.to_sql("cell_counts", conn, if_exists="append", index=False)

        print(f"Loaded {len(subjects)} subjects")
        print(f"Loaded {len(samples)} samples")
        print(f"Loaded {len(cell_counts)} cell-count rows into {DB_PATH}")


if __name__ == "__main__":
    main()
