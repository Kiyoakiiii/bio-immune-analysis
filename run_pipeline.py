"""Run all analysis steps for the Loblaw Bio immune cell project."""

from __future__ import annotations

import csv
import math
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median


PROJECT_ROOT = Path(__file__).resolve().parent
DB_PATH = PROJECT_ROOT / "cell_counts.db"
OUTPUT_DIR = PROJECT_ROOT / "outputs"

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]


def connect() -> sqlite3.Connection:
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    return conn


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def frequency_rows(conn: sqlite3.Connection) -> list[dict[str, object]]:
    query = """
        WITH sample_totals AS (
            SELECT sample_id, SUM(count) AS total_count
            FROM cell_counts
            GROUP BY sample_id
        )
        SELECT
            s.sample_code AS sample,
            st.total_count AS total_count,
            cp.population_name AS population,
            cc.count AS count,
            (100.0 * cc.count / st.total_count) AS percentage
        FROM samples s
        JOIN sample_totals st ON st.sample_id = s.sample_id
        JOIN cell_counts cc ON cc.sample_id = s.sample_id
        JOIN cell_populations cp ON cp.population_id = cc.population_id
        ORDER BY s.sample_code, cp.population_name
    """
    return [
        {
            "sample": row["sample"],
            "total_count": row["total_count"],
            "population": row["population"],
            "count": row["count"],
            "percentage": round(row["percentage"], 6),
        }
        for row in conn.execute(query)
    ]


def response_frequency_rows(conn: sqlite3.Connection) -> list[sqlite3.Row]:
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
            s.sample_code AS sample,
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
        ORDER BY cp.population_name, subj.response, s.sample_code
    """
    return list(conn.execute(query))


def average_ranks(values: list[float]) -> list[float]:
    indexed = sorted(enumerate(values), key=lambda item: item[1])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(indexed):
        j = i + 1
        while j < len(indexed) and indexed[j][1] == indexed[i][1]:
            j += 1
        average_rank = (i + 1 + j) / 2.0
        for k in range(i, j):
            ranks[indexed[k][0]] = average_rank
        i = j
    return ranks


def mann_whitney_u_pvalue(group_a: list[float], group_b: list[float]) -> tuple[float, float]:
    """Two-sided Mann-Whitney U test using a tie-corrected normal approximation."""

    n_a = len(group_a)
    n_b = len(group_b)
    if n_a == 0 or n_b == 0:
        return math.nan, math.nan

    combined = group_a + group_b
    ranks = average_ranks(combined)
    rank_sum_a = sum(ranks[:n_a])
    u_a = rank_sum_a - (n_a * (n_a + 1) / 2.0)
    u_b = n_a * n_b - u_a
    u_stat = min(u_a, u_b)

    n_total = n_a + n_b
    tie_counts = Counter(combined)
    tie_adjustment = sum(count**3 - count for count in tie_counts.values())
    variance = (n_a * n_b / 12.0) * (
        (n_total + 1) - tie_adjustment / (n_total * (n_total - 1))
    )
    if variance <= 0:
        return u_stat, 1.0

    mean_u = n_a * n_b / 2.0
    continuity = 0.5 if abs(u_stat - mean_u) >= 0.5 else 0.0
    z_score = (abs(u_stat - mean_u) - continuity) / math.sqrt(variance)
    p_value = math.erfc(z_score / math.sqrt(2.0))
    return u_stat, min(max(p_value, 0.0), 1.0)


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    valid = [(index, p_value) for index, p_value in enumerate(p_values) if not math.isnan(p_value)]
    adjusted = [math.nan] * len(p_values)
    if not valid:
        return adjusted

    valid.sort(key=lambda item: item[1])
    m = len(valid)
    running_min = 1.0
    for rank, (index, p_value) in reversed(list(enumerate(valid, start=1))):
        running_min = min(running_min, p_value * m / rank)
        adjusted[index] = min(running_min, 1.0)
    return adjusted


def response_statistics(rows: list[sqlite3.Row]) -> list[dict[str, object]]:
    by_population: dict[str, dict[str, list[float]]] = defaultdict(lambda: {"yes": [], "no": []})
    for row in rows:
        by_population[row["population"]][row["response"]].append(float(row["percentage"]))

    stats_rows: list[dict[str, object]] = []
    p_values: list[float] = []
    for population in POPULATIONS:
        responder = by_population[population]["yes"]
        non_responder = by_population[population]["no"]
        u_stat, p_value = mann_whitney_u_pvalue(responder, non_responder)
        p_values.append(p_value)
        stats_rows.append(
            {
                "population": population,
                "n_responder": len(responder),
                "n_non_responder": len(non_responder),
                "median_responder_percentage": round(median(responder), 6),
                "median_non_responder_percentage": round(median(non_responder), 6),
                "mann_whitney_u": round(u_stat, 6),
                "p_value": p_value,
            }
        )

    q_values = benjamini_hochberg(p_values)
    for row, q_value in zip(stats_rows, q_values):
        row["q_value_fdr"] = q_value
        row["significant_fdr_0_05"] = bool(q_value < 0.05)
        row["p_value"] = round(row["p_value"], 12)
        row["q_value_fdr"] = round(row["q_value_fdr"], 12)

    return stats_rows


def create_response_boxplot(rows: list[sqlite3.Row], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    grouped: dict[str, dict[str, list[float]]] = defaultdict(lambda: {"yes": [], "no": []})
    for row in rows:
        grouped[row["population"]][row["response"]].append(float(row["percentage"]))

    positions: list[float] = []
    values: list[list[float]] = []
    colors: list[str] = []
    tick_positions: list[float] = []
    tick_labels: list[str] = []

    for index, population in enumerate(POPULATIONS):
        base = index * 3.0
        responder_position = base + 1.0
        non_responder_position = base + 1.8
        positions.extend([responder_position, non_responder_position])
        values.extend([grouped[population]["yes"], grouped[population]["no"]])
        colors.extend(["#2f80ed", "#f2994a"])
        tick_positions.append((responder_position + non_responder_position) / 2.0)
        tick_labels.append(population)

    fig, ax = plt.subplots(figsize=(11, 6))
    boxplot = ax.boxplot(values, positions=positions, widths=0.55, patch_artist=True)
    for patch, color in zip(boxplot["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.68)
    for median_line in boxplot["medians"]:
        median_line.set_color("#1f2933")
        median_line.set_linewidth(1.5)

    ax.set_title("Melanoma PBMC miraclib samples: responders vs non-responders")
    ax.set_ylabel("Cell population relative frequency (%)")
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, rotation=15, ha="right")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(
        handles=[
            Patch(facecolor="#2f80ed", alpha=0.68, label="Responder"),
            Patch(facecolor="#f2994a", alpha=0.68, label="Non-responder"),
        ],
        loc="upper right",
    )
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def baseline_subset_rows(conn: sqlite3.Connection) -> list[dict[str, object]]:
    query = """
        SELECT
            p.project_code AS project,
            subj.subject_code AS subject,
            subj.indication AS indication,
            subj.age AS age,
            subj.gender AS gender,
            subj.treatment AS treatment,
            subj.response AS response,
            s.sample_code AS sample,
            s.sample_type AS sample_type,
            s.time_from_treatment_start AS time_from_treatment_start,
            SUM(CASE WHEN cp.population_name = 'b_cell' THEN cc.count END) AS b_cell,
            SUM(CASE WHEN cp.population_name = 'cd8_t_cell' THEN cc.count END) AS cd8_t_cell,
            SUM(CASE WHEN cp.population_name = 'cd4_t_cell' THEN cc.count END) AS cd4_t_cell,
            SUM(CASE WHEN cp.population_name = 'nk_cell' THEN cc.count END) AS nk_cell,
            SUM(CASE WHEN cp.population_name = 'monocyte' THEN cc.count END) AS monocyte
        FROM samples s
        JOIN subjects subj ON subj.subject_id = s.subject_id
        JOIN projects p ON p.project_id = subj.project_id
        JOIN cell_counts cc ON cc.sample_id = s.sample_id
        JOIN cell_populations cp ON cp.population_id = cc.population_id
        WHERE subj.indication = 'melanoma'
          AND subj.treatment = 'miraclib'
          AND s.sample_type = 'PBMC'
          AND s.time_from_treatment_start = 0
        GROUP BY
            p.project_code,
            subj.subject_code,
            subj.indication,
            subj.age,
            subj.gender,
            subj.treatment,
            subj.response,
            s.sample_code,
            s.sample_type,
            s.time_from_treatment_start
        ORDER BY p.project_code, subj.subject_code, s.sample_code
    """
    return [dict(row) for row in conn.execute(query)]


def baseline_counts(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for metric, source_column in [
        ("project", "project"),
        ("response", "response"),
        ("gender", "gender"),
    ]:
        counts = Counter(str(row[source_column]) for row in rows)
        for value, count in sorted(counts.items()):
            output.append({"metric": metric, "value": value, "count": count})
    return output


def validate_outputs(
    conn: sqlite3.Connection,
    frequency: list[dict[str, object]],
    baseline_rows: list[dict[str, object]],
) -> None:
    sample_count = conn.execute("SELECT COUNT(*) FROM samples").fetchone()[0]
    population_count = conn.execute("SELECT COUNT(*) FROM cell_populations").fetchone()[0]
    if sample_count != 10500:
        raise AssertionError(f"Expected 10,500 samples; found {sample_count:,}.")
    if population_count != 5:
        raise AssertionError(f"Expected 5 cell populations; found {population_count}.")
    if len(frequency) != sample_count * population_count:
        raise AssertionError("Frequency summary row count does not match samples x populations.")

    totals: dict[str, float] = defaultdict(float)
    for row in frequency:
        totals[str(row["sample"])] += float(row["percentage"])
    worst_delta = max(abs(total - 100.0) for total in totals.values())
    if worst_delta > 0.0001:
        raise AssertionError(f"Frequency percentages do not sum to 100%; max delta {worst_delta}.")

    expected_counts = {
        ("project", "prj1"): 384,
        ("project", "prj3"): 272,
        ("response", "yes"): 331,
        ("response", "no"): 325,
        ("gender", "M"): 344,
        ("gender", "F"): 312,
    }
    observed = {
        (row["metric"], row["value"]): row["count"] for row in baseline_counts(baseline_rows)
    }
    for key, expected in expected_counts.items():
        if observed.get(key) != expected:
            raise AssertionError(f"Expected {key}={expected}; found {observed.get(key)}.")
    if len(baseline_rows) != 656:
        raise AssertionError(f"Expected 656 baseline subset rows; found {len(baseline_rows)}.")


def main() -> None:
    if not DB_PATH.exists():
        raise FileNotFoundError("cell_counts.db is missing. Run python load_data.py first.")

    OUTPUT_DIR.mkdir(exist_ok=True)
    with connect() as conn:
        frequency = frequency_rows(conn)
        write_csv(
            OUTPUT_DIR / "frequency_summary.csv",
            ["sample", "total_count", "population", "count", "percentage"],
            frequency,
        )

        response_rows = response_frequency_rows(conn)
        stats_rows = response_statistics(response_rows)
        write_csv(
            OUTPUT_DIR / "miraclib_response_stats.csv",
            [
                "population",
                "n_responder",
                "n_non_responder",
                "median_responder_percentage",
                "median_non_responder_percentage",
                "mann_whitney_u",
                "p_value",
                "q_value_fdr",
                "significant_fdr_0_05",
            ],
            stats_rows,
        )
        create_response_boxplot(
            response_rows,
            OUTPUT_DIR / "miraclib_response_boxplot.png",
        )

        baseline_rows = baseline_subset_rows(conn)
        write_csv(
            OUTPUT_DIR / "baseline_miraclib_melanoma_pbmc_samples.csv",
            [
                "project",
                "subject",
                "indication",
                "age",
                "gender",
                "treatment",
                "response",
                "sample",
                "sample_type",
                "time_from_treatment_start",
                "b_cell",
                "cd8_t_cell",
                "cd4_t_cell",
                "nk_cell",
                "monocyte",
            ],
            baseline_rows,
        )

        baseline_count_rows = baseline_counts(baseline_rows)
        write_csv(
            OUTPUT_DIR / "baseline_subset_counts.csv",
            ["metric", "value", "count"],
            baseline_count_rows,
        )

        validate_outputs(conn, frequency, baseline_rows)

    print("Pipeline complete. Outputs written to outputs/.")


if __name__ == "__main__":
    main()
