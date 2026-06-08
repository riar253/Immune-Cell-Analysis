"""
Analysis pipeline for the immune-cell dataset. Reads from
cell_count.db and writes result tables/plots into the output/ directory.
"""

import os
import sqlite3

import matplotlib

matplotlib.use("Agg")  # non-interactive backend so plots render headless in the pipeline

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu

DB_PATH = "cell_count.db"
OUTPUT_DIR = "output"

# Significance threshold applied to the FDR-adjusted p-values in Part 3.
ALPHA = 0.05

# Fixed display order for the five cell populations (keeps plots consistent).
POPULATION_ORDER = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

# Open a connection to the project database
def get_connection():
    return sqlite3.connect(DB_PATH)


"""
Returns one row per (sample, population) with the sample's total cell
count and that population's percentage of the total. Percentage is kept at
full precision here; callers round for display.
"""
def cell_population_frequencies(conn):
    
    query = """
        SELECT
            sample_id AS sample,
            SUM(count) OVER (PARTITION BY sample_id) AS total_count,
            population,
            count,
            100.0 * count / SUM(count) OVER (PARTITION BY sample_id) AS percentage
        FROM cell_counts
        ORDER BY sample_id, population;
    """
    return pd.read_sql_query(query, conn)


# Part 3 cohort: per-sample population frequencies for melanoma patients on
# miraclib, PBMC samples only, restricted to clear responders/non-responders.
# Returns tidy rows (sample, subject_id, response, population, percentage) that
# both the stats and the boxplot consume, and that the dashboard can reuse.
def responder_comparison_data(conn):
    query = """
        SELECT
            cc.sample_id AS sample,
            s.subject_id,
            sub.response,
            cc.population,
            100.0 * cc.count / SUM(cc.count) OVER (PARTITION BY cc.sample_id) AS percentage
        FROM cell_counts cc
        JOIN samples s    ON s.sample_id = cc.sample_id
        JOIN subjects sub ON sub.subject_id = s.subject_id
        WHERE sub.condition = 'melanoma'
          AND sub.treatment = 'miraclib'
          AND s.sample_type = 'PBMC'
          AND sub.response IN ('yes', 'no')
        ORDER BY cc.population, cc.sample_id;
    """
    return pd.read_sql_query(query, conn)


# Benjamini-Hochberg FDR adjustment for a list of p-values. Returns adjusted
# p-values aligned to the input order, clipped to [0, 1].
def benjamini_hochberg(pvals):
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = p.argsort()
    ranked = p[order] * n / (np.arange(n) + 1)
    # Enforce monotonicity from the largest p-value downward.
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty(n)
    adjusted[order] = np.clip(ranked, 0, 1)
    return adjusted


# Part 3 stats: per population, Mann-Whitney U comparing responder vs
# non-responder relative frequencies, with BH-adjusted p-values. Sorted by raw
# p-value so the strongest signals surface first.
def compare_responders(df):
    rows = []
    for population, group in df.groupby("population"):
        responders = group.loc[group["response"] == "yes", "percentage"]
        non_responders = group.loc[group["response"] == "no", "percentage"]
        u_stat, p_value = mannwhitneyu(
            responders, non_responders, alternative="two-sided"
        )
        rows.append(
            {
                "population": population,
                "n_responder": len(responders),
                "n_non_responder": len(non_responders),
                "median_responder": responders.median(),
                "median_non_responder": non_responders.median(),
                "u_statistic": u_stat,
                "p_value": p_value,
            }
        )

    stats = pd.DataFrame(rows)
    stats["p_adjusted"] = benjamini_hochberg(stats["p_value"])
    stats["significant"] = stats["p_adjusted"] < ALPHA
    return stats.sort_values("p_value").reset_index(drop=True)


# Part 3 plot: side-by-side boxplots of relative frequency per population,
# split by response. Returns the figure so callers can save it (pipeline) or
# render it (dashboard).
def plot_responder_boxplots(df):
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.boxplot(
        data=df,
        x="population",
        y="percentage",
        hue="response",
        order=POPULATION_ORDER,
        ax=ax,
    )
    ax.set_title(
        "Relative frequency by cell population: responders vs non-responders\n"
        "(melanoma, miraclib, PBMC samples)"
    )
    ax.set_xlabel("Cell population")
    ax.set_ylabel("Relative frequency (%)")
    ax.legend(title="Responder")
    fig.tight_layout()
    return fig


# Part 4 subset: baseline (time_from_treatment_start = 0) melanoma PBMC samples
# from miraclib-treated patients. One row per sample, carrying the subject-level
# fields the breakdowns need. Reused by the dashboard.
def baseline_subset(conn):
    query = """
        SELECT
            s.sample_id,
            s.subject_id,
            sub.project,
            sub.response,
            sub.sex
        FROM samples s
        JOIN subjects sub ON sub.subject_id = s.subject_id
        WHERE sub.condition = 'melanoma'
          AND sub.treatment = 'miraclib'
          AND s.sample_type = 'PBMC'
          AND s.time_from_treatment_start = 0
        ORDER BY s.sample_id;
    """
    return pd.read_sql_query(query, conn)


# Part 4 breakdowns of the baseline subset. Returns three small tables:
#   samples_per_project  - sample counts (one row per sample)
#   subjects_per_response - distinct-subject counts by response
#   subjects_per_sex      - distinct-subject counts by sex
# Response/sex count subjects (deduplicated) since those are subject attributes.
def baseline_breakdowns(df):
    subjects = df.drop_duplicates(subset="subject_id")

    samples_per_project = (
        df.groupby("project").size().reset_index(name="n_samples")
    )
    subjects_per_response = (
        subjects.groupby("response").size().reset_index(name="n_subjects")
    )
    subjects_per_sex = (
        subjects.groupby("sex").size().reset_index(name="n_subjects")
    )
    return samples_per_project, subjects_per_response, subjects_per_sex


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    with get_connection() as conn:
        # Cell population relative frequencies
        freq = cell_population_frequencies(conn)
        freq_path = os.path.join(OUTPUT_DIR, "cell_frequencies.csv")
        freq.round({"percentage": 2}).to_csv(freq_path, index=False)

        print("\n=== Part 2: Cell population relative frequencies ===")
        
        # Render a pipe-separated grid table. to_markdown has no row limit, so
        # show a head + ellipsis + tail preview; the full table is in the CSV.
        preview = freq.round({"percentage": 2})
        head, tail = preview.head(5), preview.tail(5)
        ellipsis = pd.DataFrame([["..."] * preview.shape[1]], columns=preview.columns)
        preview = pd.concat([head, ellipsis, tail], ignore_index=True)
        print(preview.to_markdown(index=False))
        print(f"\n({len(freq)} rows total — full table written to {freq_path})")

        # Part 3: responders vs non-responders (melanoma, miraclib, PBMC).
        comparison = responder_comparison_data(conn)

        boxplot_path = os.path.join(OUTPUT_DIR, "responder_boxplot.png")
        fig = plot_responder_boxplots(comparison)
        fig.savefig(boxplot_path, dpi=150)
        plt.close(fig)

        stats = compare_responders(comparison)
        stats_path = os.path.join(OUTPUT_DIR, "responder_stats.csv")
        stats.to_csv(stats_path, index=False)

        n_samples = comparison["sample"].nunique()
        print("\n=== Part 3: Responders vs non-responders (melanoma, miraclib, PBMC) ===")
        print(f"Cohort: {n_samples} samples")
        print(stats.round({
            "median_responder": 2, "median_non_responder": 2,
            "p_value": 4, "p_adjusted": 4,
        }).to_markdown(index=False))

        significant = stats.loc[stats["significant"], "population"].tolist()
        if significant:
            print(f"\nSignificant populations (BH-adjusted p < {ALPHA}): "
                  f"{', '.join(significant)}")
        else:
            print(f"\nNo populations significant at BH-adjusted p < {ALPHA}.")
        print(f"Boxplot written to {boxplot_path}")
        print(f"Stats table written to {stats_path}")

        # Part 4: baseline melanoma/miraclib/PBMC subset breakdowns.
        baseline = baseline_subset(conn)
        baseline_path = os.path.join(OUTPUT_DIR, "baseline_samples.csv")
        baseline.to_csv(baseline_path, index=False)

        per_project, per_response, per_sex = baseline_breakdowns(baseline)
        per_project.to_csv(os.path.join(OUTPUT_DIR, "baseline_per_project.csv"), index=False)
        per_response.to_csv(os.path.join(OUTPUT_DIR, "baseline_per_response.csv"), index=False)
        per_sex.to_csv(os.path.join(OUTPUT_DIR, "baseline_per_sex.csv"), index=False)

        n_subjects = baseline["subject_id"].nunique()
        print("\n=== Part 4: Baseline subset (melanoma, miraclib, PBMC, time = 0) ===")
        print(f"{len(baseline)} samples from {n_subjects} subjects")
        print("\nSamples per project:")
        print(per_project.to_markdown(index=False))
        print("\nSubjects by response:")
        print(per_response.to_markdown(index=False))
        print("\nSubjects by sex:")
        print(per_sex.to_markdown(index=False))
        print(f"\nSubset and breakdowns written to {OUTPUT_DIR}/baseline_*.csv")


if __name__ == "__main__":
    main()
