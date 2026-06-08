"""Interactive Streamlit dashboard for the immune-cell analysis (Parts 2-4).

Launch with `make dashboard` (or `streamlit run app.py`). Reads cell_count.db,
so run `python load_data.py` (or `make pipeline`) first to build it.

All analysis logic is imported from analysis.py so the dashboard and the
command-line pipeline stay in sync.
"""

import os

import streamlit as st

import analysis

st.set_page_config(page_title="Loblaw Bio — Immune Cell Analysis", layout="wide")


# Cache the queries so re-renders (tab switches, widget changes) don't re-hit
# the database. Each helper opens and closes its own connection.
@st.cache_data
def load_frequencies():
    with analysis.get_connection() as conn:
        return analysis.cell_population_frequencies(conn)


@st.cache_data
def load_responder_data():
    with analysis.get_connection() as conn:
        return analysis.responder_comparison_data(conn)


@st.cache_data
def load_baseline():
    with analysis.get_connection() as conn:
        return analysis.baseline_subset(conn)


def part2_tab():
    st.header("Part 2 — Cell population relative frequencies")
    st.caption(
        "For each sample, the relative frequency of each of the five immune "
        "cell populations as a percentage of that sample's total cell count."
    )
    freq = load_frequencies()

    samples = sorted(freq["sample"].unique())
    selected = st.multiselect(
        "Filter by sample (leave empty to show all):", samples, default=[]
    )
    table = freq[freq["sample"].isin(selected)] if selected else freq

    st.write(f"Showing **{table['sample'].nunique():,}** samples "
             f"({len(table):,} rows).")
    st.dataframe(
        table.round({"percentage": 2}), use_container_width=True, hide_index=True
    )
    st.download_button(
        "Download as CSV",
        table.round({"percentage": 2}).to_csv(index=False),
        file_name="cell_frequencies.csv",
        mime="text/csv",
    )


def part3_tab():
    st.header("Part 3 — Responders vs non-responders")
    st.caption(
        "Melanoma patients treated with miraclib, PBMC samples only. "
        "Comparison of cell population relative frequencies between responders "
        "(response = yes) and non-responders (response = no)."
    )
    df = load_responder_data()

    n_samples = df["sample"].nunique()
    n_resp = df.loc[df["response"] == "yes", "sample"].nunique()
    n_non = df.loc[df["response"] == "no", "sample"].nunique()
    c1, c2, c3 = st.columns(3)
    c1.metric("Samples in cohort", f"{n_samples:,}")
    c2.metric("Responder samples", f"{n_resp:,}")
    c3.metric("Non-responder samples", f"{n_non:,}")

    st.subheader("Boxplot")
    fig = analysis.plot_responder_boxplots(df)
    st.pyplot(fig)

    st.subheader("Statistical comparison")
    st.caption(
        "Mann–Whitney U test per population with Benjamini–Hochberg FDR "
        f"correction. Significance at adjusted p < {analysis.ALPHA}."
    )
    stats = analysis.compare_responders(df)
    st.dataframe(
        stats.round({
            "median_responder": 2, "median_non_responder": 2,
            "p_value": 4, "p_adjusted": 4,
        }),
        use_container_width=True,
        hide_index=True,
    )

    significant = stats.loc[stats["significant"], "population"].tolist()
    if significant:
        st.success(f"Significant populations (adjusted p < {analysis.ALPHA}): "
                   f"{', '.join(significant)}")
    else:
        st.info(
            f"No populations show a significant difference at adjusted "
            f"p < {analysis.ALPHA} after correcting for multiple comparisons."
        )


def part4_tab():
    st.header("Part 4 — Baseline subset analysis")
    st.caption(
        "Melanoma PBMC samples at baseline (time_from_treatment_start = 0) "
        "from patients treated with miraclib."
    )
    baseline = load_baseline()
    per_project, per_response, per_sex = analysis.baseline_breakdowns(baseline)

    c1, c2 = st.columns(2)
    c1.metric("Baseline samples", f"{len(baseline):,}")
    c2.metric("Distinct subjects", f"{baseline['subject_id'].nunique():,}")

    col1, col2, col3 = st.columns(3)
    with col1:
        st.subheader("Samples per project")
        st.dataframe(per_project, hide_index=True, use_container_width=True)
        st.bar_chart(per_project, x="project", y="n_samples")
    with col2:
        st.subheader("Subjects by response")
        st.dataframe(per_response, hide_index=True, use_container_width=True)
        st.bar_chart(per_response, x="response", y="n_subjects")
    with col3:
        st.subheader("Subjects by sex")
        st.dataframe(per_sex, hide_index=True, use_container_width=True)
        st.bar_chart(per_sex, x="sex", y="n_subjects")


def main():
    st.title("Loblaw Bio — Immune Cell Analysis")

    if not os.path.exists(analysis.DB_PATH):
        st.error(
            f"Database `{analysis.DB_PATH}` not found. "
            "Run `python load_data.py` (or `make pipeline`) first."
        )
        st.stop()

    tab2, tab3, tab4 = st.tabs([
        "Part 2 · Frequencies",
        "Part 3 · Responders",
        "Part 4 · Baseline subset",
    ])
    with tab2:
        part2_tab()
    with tab3:
        part3_tab()
    with tab4:
        part4_tab()


if __name__ == "__main__":
    main()
