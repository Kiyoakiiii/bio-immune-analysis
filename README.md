# Loblaw Bio Immune Cell Analysis

This repository analyzes immune cell count data from a clinical trial and provides a reproducible SQLite pipeline plus an interactive Streamlit dashboard.

Dashboard URL after starting the server:

<http://localhost:8501>

## Quick Start

```bash
make setup
make pipeline
make dashboard
```

`make pipeline` creates `cell_counts.db` in the repository root and writes analysis outputs to `outputs/`.

## Project Structure

- `load_data.py` initializes the SQLite database and loads `data/cell-count.csv`.
- `run_pipeline.py` generates the frequency summary, statistical comparison, baseline subset tables, static plot, and extra question answer.
- `app.py` starts the Streamlit dashboard.
- `outputs/` contains generated CSV, PNG, and TXT outputs from the pipeline.

## Database Schema

The database uses a normalized relational design:

- `projects`: one row per project code.
- `subjects`: subject-level metadata, including indication, age, gender, treatment, and response.
- `samples`: sample-level metadata, including sample type and time from treatment start.
- `cell_populations`: controlled vocabulary for immune cell populations.
- `cell_counts`: long-format counts keyed by sample and population.

This schema avoids repeating subject and project metadata on every measurement row. It also scales cleanly if the trial grows to hundreds of projects, thousands of samples, or additional immune cell populations. New cell types can be added as rows in `cell_populations`, and downstream analyses can query `cell_counts` without changing the table shape. Indexes on treatment, response, indication, gender, sample type, time, and population support the filters used by the dashboard and analysis scripts.

The source CSV uses `condition`, `sex`, and `sample`; these are mapped to `indication`, `gender`, and `sample_code` in the database. Empty response values are stored as NULL.

## Analysis

Part 2 computes the relative frequency of each immune cell population in each sample:

```text
sample,total_count,population,count,percentage
```

Part 3 compares melanoma PBMC samples from miraclib-treated responders and non-responders across all available time points. The analysis uses a two-sided Mann-Whitney U test for each immune population, followed by Benjamini-Hochberg FDR correction at `q < 0.05`.

Part 4 filters to melanoma PBMC baseline samples from miraclib-treated subjects and summarizes counts by project, response, and gender.

The extra question answer is written to `outputs/extra_question_answer.txt`.

## Results Summary

The baseline melanoma PBMC miraclib subset contains 656 samples. Counts by project are `prj1=384` and `prj3=272`; counts by response are `yes=331` and `no=325`; counts by gender are `M=344` and `F=312`.

For the responder versus non-responder comparison, no immune cell population is significant after Benjamini-Hochberg FDR correction at `q < 0.05`. The lowest nominal p-value is for `cd4_t_cell` (`p=0.013344`, `q=0.066721`), which is suggestive but does not meet the adjusted significance threshold.

Considering melanoma males, responders at `time_from_treatment_start=0` have an average B-cell count of `10206.15`.

## Generated Outputs

- `outputs/frequency_summary.csv`
- `outputs/miraclib_response_stats.csv`
- `outputs/miraclib_response_boxplot.png`
- `outputs/baseline_miraclib_melanoma_pbmc_samples.csv`
- `outputs/baseline_subset_counts.csv`
- `outputs/extra_question_answer.txt`

## Reproducibility Notes

The pipeline is deterministic. Re-running `make pipeline` recreates `cell_counts.db` and regenerates all outputs from `data/cell-count.csv`.
