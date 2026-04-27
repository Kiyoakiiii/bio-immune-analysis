"""Interactive Streamlit dashboard for the Loblaw Bio immune analysis."""

from __future__ import annotations

import sqlite3
from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st


PROJECT_ROOT = Path(__file__).resolve().parent
DB_PATH = PROJECT_ROOT / "cell_counts.db"
OUTPUT_DIR = PROJECT_ROOT / "outputs"


st.set_page_config(
    page_title="Loblaw Bio Immune Analysis",
    layout="wide",
)


@st.cache_data
def read_csv(path: str) -> pd.DataFrame:
    return pd.read_csv(path)


@st.cache_data
def read_text(path: str) -> str:
    return Path(path).read_text(encoding="utf-8").strip()


@st.cache_data
def response_frequency_data() -> pd.DataFrame:
    with sqlite3.connect(DB_PATH) as conn:
        query = """
            WITH sample_totals AS (
                SELECT sample_id, SUM(count) AS total_count
                FROM cell_counts
                GROUP BY sample_id
            )
            SELECT
                p.project_code AS project,
                subj.subject_code AS subject,
                subj.response AS response,
                subj.gender AS gender,
                s.sample_code AS sample,
                s.sample_type AS sample_type,
                s.time_from_treatment_start AS time_from_treatment_start,
                cp.population_name AS population,
                cc.count AS count,
                (100.0 * cc.count / st.total_count) AS percentage
            FROM samples s
            JOIN subjects subj ON subj.subject_id = s.subject_id
            JOIN projects p ON p.project_id = subj.project_id
            JOIN sample_totals st ON st.sample_id = s.sample_id
            JOIN cell_counts cc ON cc.sample_id = s.sample_id
            JOIN cell_populations cp ON cp.population_id = cc.population_id
            WHERE subj.indication = 'melanoma'
              AND subj.treatment = 'miraclib'
              AND s.sample_type = 'PBMC'
              AND subj.response IN ('yes', 'no')
        """
        return pd.read_sql_query(query, conn)


def require_pipeline_outputs() -> bool:
    required_paths = [
        DB_PATH,
        OUTPUT_DIR / "frequency_summary.csv",
        OUTPUT_DIR / "miraclib_response_stats.csv",
        OUTPUT_DIR / "baseline_miraclib_melanoma_pbmc_samples.csv",
        OUTPUT_DIR / "baseline_subset_counts.csv",
        OUTPUT_DIR / "extra_question_answer.txt",
    ]
    missing = [path.name for path in required_paths if not path.exists()]
    if missing:
        st.error("Run `make pipeline` before opening the dashboard.")
        st.write("Missing files:", ", ".join(missing))
        return False
    return True


def overview_page() -> None:
    frequency = read_csv(str(OUTPUT_DIR / "frequency_summary.csv"))
    st.subheader("Cell Population Frequency By Sample")

    left, middle, right = st.columns(3)
    left.metric("Samples", f"{frequency['sample'].nunique():,}")
    middle.metric("Cell populations", f"{frequency['population'].nunique():,}")
    right.metric("Summary rows", f"{len(frequency):,}")

    populations = sorted(frequency["population"].unique())
    selected_population = st.multiselect(
        "Population filter",
        populations,
        default=populations,
    )
    filtered = frequency[frequency["population"].isin(selected_population)]

    sample_options = sorted(frequency["sample"].unique())
    sample_index = 0
    selected_sample = st.selectbox("Inspect one sample", sample_options, index=sample_index)
    sample_data = frequency[frequency["sample"] == selected_sample]

    fig = px.bar(
        sample_data,
        x="population",
        y="percentage",
        color="population",
        text="percentage",
        title=f"Relative frequency for {selected_sample}",
    )
    fig.update_traces(texttemplate="%{text:.2f}%", textposition="outside")
    fig.update_layout(showlegend=False, yaxis_title="Relative frequency (%)")
    st.plotly_chart(fig, use_container_width=True)

    st.dataframe(filtered, use_container_width=True, hide_index=True)


def responder_page() -> None:
    stats = read_csv(str(OUTPUT_DIR / "miraclib_response_stats.csv"))
    response_data = response_frequency_data()

    st.subheader("Miraclib Response Comparison")
    st.caption("Melanoma PBMC samples, all time points, responder vs non-responder.")

    fig = px.box(
        response_data,
        x="population",
        y="percentage",
        color="response",
        points="all",
        hover_data=["sample", "subject", "project", "time_from_treatment_start"],
        category_orders={"response": ["yes", "no"]},
        labels={"percentage": "Relative frequency (%)", "response": "Response"},
    )
    st.plotly_chart(fig, use_container_width=True)

    significant = stats[stats["significant_fdr_0_05"] == True]
    st.metric("Significant populations at FDR < 0.05", len(significant))
    st.dataframe(stats, use_container_width=True, hide_index=True)


def baseline_page() -> None:
    samples = read_csv(str(OUTPUT_DIR / "baseline_miraclib_melanoma_pbmc_samples.csv"))
    counts = read_csv(str(OUTPUT_DIR / "baseline_subset_counts.csv"))
    extra_answer = read_text(str(OUTPUT_DIR / "extra_question_answer.txt"))

    st.subheader("Baseline Melanoma PBMC Miraclib Subset")

    left, middle, right = st.columns(3)
    left.metric("Baseline samples", f"{len(samples):,}")
    middle.metric("Subjects", f"{samples['subject'].nunique():,}")
    right.metric("Male melanoma responder B-cell average at time 0", extra_answer)

    fig = px.bar(
        counts,
        x="value",
        y="count",
        color="metric",
        facet_col="metric",
        text="count",
        labels={"value": "", "count": "Samples"},
    )
    fig.update_xaxes(matches=None)
    st.plotly_chart(fig, use_container_width=True)

    st.dataframe(counts, use_container_width=True, hide_index=True)
    st.dataframe(samples, use_container_width=True, hide_index=True)


def schema_page() -> None:
    st.subheader("Database Schema")
    st.write(
        """
        The SQLite model separates projects, subjects, samples, cell populations, and
        counts. This keeps repeated metadata out of the measurement table and makes
        population-level analytics easier to extend as new cell types or projects are
        added.
        """
    )
    schema = pd.DataFrame(
        [
            {"table": "projects", "purpose": "One row per clinical project."},
            {"table": "subjects", "purpose": "Subject-level metadata and response."},
            {"table": "samples", "purpose": "Sample identifiers, type, and treatment time."},
            {"table": "cell_populations", "purpose": "Controlled vocabulary for immune populations."},
            {"table": "cell_counts", "purpose": "Long-format measurement table."},
        ]
    )
    st.dataframe(schema, use_container_width=True, hide_index=True)


def main() -> None:
    st.title("Loblaw Bio Immune Cell Analysis")
    if not require_pipeline_outputs():
        return

    page = st.sidebar.radio(
        "View",
        [
            "Overview",
            "Responder comparison",
            "Baseline subset",
            "Database schema",
        ],
    )

    if page == "Overview":
        overview_page()
    elif page == "Responder comparison":
        responder_page()
    elif page == "Baseline subset":
        baseline_page()
    else:
        schema_page()


if __name__ == "__main__":
    main()
