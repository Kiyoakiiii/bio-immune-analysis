"""Load immune cell count data into a normalized SQLite database."""

from __future__ import annotations

import csv
import sqlite3
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent
DATA_PATH = PROJECT_ROOT / "data" / "cell-count.csv"
DB_PATH = PROJECT_ROOT / "cell_counts.db"

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

REQUIRED_COLUMNS = {
    "project",
    "subject",
    "condition",
    "age",
    "sex",
    "treatment",
    "response",
    "sample",
    "sample_type",
    "time_from_treatment_start",
    *POPULATIONS,
}


def connect() -> sqlite3.Connection:
    conn = sqlite3.connect(DB_PATH)
    conn.execute("PRAGMA foreign_keys = ON")
    return conn


def initialize_database(conn: sqlite3.Connection) -> None:
    conn.executescript(
        """
        CREATE TABLE projects (
            project_id INTEGER PRIMARY KEY AUTOINCREMENT,
            project_code TEXT NOT NULL UNIQUE
        );

        CREATE TABLE subjects (
            subject_id INTEGER PRIMARY KEY AUTOINCREMENT,
            project_id INTEGER NOT NULL,
            subject_code TEXT NOT NULL,
            indication TEXT NOT NULL,
            age INTEGER,
            gender TEXT NOT NULL,
            treatment TEXT NOT NULL,
            response TEXT,
            FOREIGN KEY (project_id) REFERENCES projects(project_id),
            UNIQUE (project_id, subject_code)
        );

        CREATE TABLE samples (
            sample_id INTEGER PRIMARY KEY AUTOINCREMENT,
            subject_id INTEGER NOT NULL,
            sample_code TEXT NOT NULL UNIQUE,
            sample_type TEXT NOT NULL,
            time_from_treatment_start INTEGER NOT NULL,
            FOREIGN KEY (subject_id) REFERENCES subjects(subject_id)
        );

        CREATE TABLE cell_populations (
            population_id INTEGER PRIMARY KEY AUTOINCREMENT,
            population_name TEXT NOT NULL UNIQUE
        );

        CREATE TABLE cell_counts (
            sample_id INTEGER NOT NULL,
            population_id INTEGER NOT NULL,
            count INTEGER NOT NULL CHECK (count >= 0),
            PRIMARY KEY (sample_id, population_id),
            FOREIGN KEY (sample_id) REFERENCES samples(sample_id),
            FOREIGN KEY (population_id) REFERENCES cell_populations(population_id)
        );

        CREATE INDEX idx_subjects_indication_treatment_response
            ON subjects(indication, treatment, response);
        CREATE INDEX idx_subjects_gender
            ON subjects(gender);
        CREATE INDEX idx_samples_type_time
            ON samples(sample_type, time_from_treatment_start);
        CREATE INDEX idx_cell_counts_population
            ON cell_counts(population_id);
        """
    )

    conn.executemany(
        "INSERT INTO cell_populations (population_name) VALUES (?)",
        [(population,) for population in POPULATIONS],
    )


def get_or_create_project(conn: sqlite3.Connection, project_code: str) -> int:
    conn.execute(
        "INSERT OR IGNORE INTO projects (project_code) VALUES (?)",
        (project_code,),
    )
    row = conn.execute(
        "SELECT project_id FROM projects WHERE project_code = ?",
        (project_code,),
    ).fetchone()
    if row is None:
        raise RuntimeError(f"Could not load project {project_code!r}")
    return int(row[0])


def get_or_create_subject(
    conn: sqlite3.Connection,
    *,
    project_id: int,
    subject_code: str,
    indication: str,
    age: int | None,
    gender: str,
    treatment: str,
    response: str | None,
) -> int:
    conn.execute(
        """
        INSERT OR IGNORE INTO subjects (
            project_id, subject_code, indication, age, gender, treatment, response
        )
        VALUES (?, ?, ?, ?, ?, ?, ?)
        """,
        (project_id, subject_code, indication, age, gender, treatment, response),
    )
    row = conn.execute(
        """
        SELECT subject_id
        FROM subjects
        WHERE project_id = ? AND subject_code = ?
        """,
        (project_id, subject_code),
    ).fetchone()
    if row is None:
        raise RuntimeError(f"Could not load subject {subject_code!r}")
    return int(row[0])


def insert_sample(
    conn: sqlite3.Connection,
    *,
    subject_id: int,
    sample_code: str,
    sample_type: str,
    time_from_treatment_start: int,
) -> int:
    conn.execute(
        """
        INSERT INTO samples (
            subject_id, sample_code, sample_type, time_from_treatment_start
        )
        VALUES (?, ?, ?, ?)
        """,
        (subject_id, sample_code, sample_type, time_from_treatment_start),
    )
    row = conn.execute(
        "SELECT sample_id FROM samples WHERE sample_code = ?",
        (sample_code,),
    ).fetchone()
    if row is None:
        raise RuntimeError(f"Could not load sample {sample_code!r}")
    return int(row[0])


def load_rows(conn: sqlite3.Connection) -> int:
    if not DATA_PATH.exists():
        raise FileNotFoundError(
            f"Expected input CSV at {DATA_PATH}. Copy cell-count.csv into data/ first."
        )

    population_ids = {
        name: population_id
        for population_id, name in conn.execute(
            "SELECT population_id, population_name FROM cell_populations"
        )
    }

    with DATA_PATH.open(newline="", encoding="utf-8") as csv_file:
        reader = csv.DictReader(csv_file)
        missing = REQUIRED_COLUMNS.difference(reader.fieldnames or [])
        if missing:
            missing_list = ", ".join(sorted(missing))
            raise ValueError(f"Input CSV is missing required columns: {missing_list}")

        loaded = 0
        with conn:
            for row in reader:
                response = row["response"].strip() or None
                project_id = get_or_create_project(conn, row["project"].strip())
                subject_id = get_or_create_subject(
                    conn,
                    project_id=project_id,
                    subject_code=row["subject"].strip(),
                    indication=row["condition"].strip(),
                    age=int(row["age"]) if row["age"].strip() else None,
                    gender=row["sex"].strip(),
                    treatment=row["treatment"].strip(),
                    response=response,
                )
                sample_id = insert_sample(
                    conn,
                    subject_id=subject_id,
                    sample_code=row["sample"].strip(),
                    sample_type=row["sample_type"].strip(),
                    time_from_treatment_start=int(row["time_from_treatment_start"]),
                )
                conn.executemany(
                    """
                    INSERT INTO cell_counts (sample_id, population_id, count)
                    VALUES (?, ?, ?)
                    """,
                    [
                        (sample_id, population_ids[population], int(row[population]))
                        for population in POPULATIONS
                    ],
                )
                loaded += 1

    return loaded


def main() -> None:
    if DB_PATH.exists():
        DB_PATH.unlink()

    with connect() as conn:
        initialize_database(conn)
        loaded = load_rows(conn)

    print(f"Created {DB_PATH.name} and loaded {loaded:,} samples.")


if __name__ == "__main__":
    main()
