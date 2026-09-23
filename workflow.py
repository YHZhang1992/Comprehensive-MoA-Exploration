#!/usr/bin/env python3
"""Dependency-free reference workflow for mechanism-of-action exploration."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path


REQUIRED_COLUMNS = {
    "gene",
    "pathway",
    "model_log2fc",
    "treatment_log2fc",
    "human_log2fc",
    "network_degree",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = REQUIRED_COLUMNS - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Input is missing columns: {', '.join(sorted(missing))}")
        rows = list(reader)
    if not rows:
        raise ValueError("Input must contain at least one row")
    if len({row["gene"] for row in rows}) != len(rows):
        raise ValueError("Gene identifiers must be unique in the reference profile")
    return rows


def numeric(row: dict[str, str], field: str) -> float:
    try:
        value = float(row[field])
    except ValueError as error:
        raise ValueError(f"{row.get('gene', '<unknown>')}: {field} must be numeric") from error
    if not math.isfinite(value):
        raise ValueError(f"{row.get('gene', '<unknown>')}: {field} must be finite")
    return value


def score(rows: list[dict[str, str]], config: dict) -> list[dict[str, object]]:
    weights = config.get("weights", {})
    required_weights = {"model_signal", "treatment_reversal", "human_translation", "network_importance"}
    missing = required_weights - set(weights)
    if missing:
        raise ValueError(f"Config is missing weights: {', '.join(sorted(missing))}")
    total_weight = sum(float(weights[name]) for name in required_weights)
    if total_weight <= 0:
        raise ValueError("At least one score weight must be positive")

    minimum_signal = float(config.get("minimum_signal_log2fc", 1.0))
    if minimum_signal <= 0:
        raise ValueError("minimum_signal_log2fc must be positive")
    maximum_degree = max(numeric(row, "network_degree") for row in rows)
    if maximum_degree < 0:
        raise ValueError("network_degree cannot be negative")

    results: list[dict[str, object]] = []
    for row in rows:
        model = numeric(row, "model_log2fc")
        treatment = numeric(row, "treatment_log2fc")
        human = numeric(row, "human_log2fc")
        degree = numeric(row, "network_degree")
        if degree < 0:
            raise ValueError(f"{row['gene']}: network_degree cannot be negative")

        model_signal = min(1.0, abs(model) / minimum_signal)
        treatment_reversal = 0.0
        if model * treatment < 0 and model != 0:
            treatment_reversal = min(1.0, abs(treatment) / abs(model))
        human_translation = 1.0 / (1.0 + abs(model - human))
        network_importance = degree / maximum_degree if maximum_degree else 0.0
        components = {
            "model_signal": model_signal,
            "treatment_reversal": treatment_reversal,
            "human_translation": human_translation,
            "network_importance": network_importance,
        }
        composite = sum(components[name] * float(weights[name]) for name in required_weights) / total_weight
        results.append(
            {
                "gene": row["gene"],
                "pathway": row["pathway"],
                **{name: round(value, 6) for name, value in components.items()},
                "composite_score": round(composite, 6),
            }
        )
    return sorted(results, key=lambda item: (-float(item["composite_score"]), str(item["gene"])))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def pathway_summary(results: list[dict[str, object]]) -> list[dict[str, object]]:
    groups: dict[str, list[float]] = {}
    for row in results:
        groups.setdefault(str(row["pathway"]), []).append(float(row["composite_score"]))
    return [
        {"pathway": pathway, "gene_count": len(values), "mean_composite_score": round(sum(values) / len(values), 6)}
        for pathway, values in sorted(groups.items(), key=lambda item: (-sum(item[1]) / len(item[1]), item[0]))
    ]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    rows = load_rows(args.input)
    config = json.loads(args.config.read_text(encoding="utf-8"))
    results = score(rows, config)
    args.output.mkdir(parents=True, exist_ok=True)
    write_csv(args.output / "results.csv", results)
    write_csv(args.output / "pathway_summary.csv", pathway_summary(results))
    manifest = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": platform.python_version(),
        "input_sha256": sha256(args.input),
        "config_sha256": sha256(args.config),
        "row_count": len(results),
        "limitations": "Reference prioritization only; requires biological and translational validation.",
    }
    (args.output / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
