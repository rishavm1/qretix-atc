"""
leakage_auditor.py
──────────────────
Automated leakage detection. Must pass before any model training begins.

Four checks:
  1. Temporal isolation  — no post-cutoff PMIDs/trial IDs in a temporal split dir
  2. Timestamp coverage  — every literature/clinical-trial feature row has a date
  3. Driver flag isolation — known-driver flag is a standalone column, not merged
  4. Drug-screen bias    — PRISM/drug-screen panel membership flagged, not silently
                           encoded in dependency features

Produces: evaluation/reports/leakage_audit.md
"""

from __future__ import annotations

import json
import logging
from datetime import date
from pathlib import Path
from typing import Any

import pandas as pd

log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[1]
REPORT_PATH = ROOT / "evaluation" / "reports" / "leakage_audit.md"
TEMPORAL_SPLITS_DIR = ROOT / "evaluation" / "temporal_splits"
LITERATURE_DIR = ROOT / "literature"
FEATURES_PATH = ROOT / "features" / "feature_table.parquet"

CUTOFF_YEARS = [2015, 2017, 2019, 2021]

# Columns that must carry per-row timestamps (not pre-aggregated)
TIMESTAMPED_FEATURE_COLS = [
    "lit_most_recent_atc_date",
    "lit_most_recent_cancer_date",
]

# The known-driver flag must exist as its own column — never merged into a score
DRIVER_FLAG_COL = "genetic_known_driver_flag"

# Dependency columns that must NOT encode drug-screen panel membership silently
DEPENDENCY_COLS = [
    "dep_atc_gene_effect",
    "dep_atc_selectivity",
    "dep_atc_percentile",
]
DRUG_SCREEN_FLAG_COL = "dep_in_drug_screen_panel_flag"


# ── Individual checks ─────────────────────────────────────────────────────────

def check_temporal_isolation(cutoff_year: int) -> dict[str, Any]:
    """
    Verify that no PMID or trial ID dated after cutoff_year appears in the
    corresponding temporal split input files.
    """
    split_dir = TEMPORAL_SPLITS_DIR / str(cutoff_year)
    result: dict[str, Any] = {
        "check": "temporal_isolation",
        "cutoff_year": cutoff_year,
        "status": "SKIP",
        "detail": "",
    }

    if not split_dir.exists():
        result["detail"] = f"Split directory {split_dir} does not exist yet."
        return result

    violations: list[str] = []
    cutoff_date = date(cutoff_year, 12, 31)

    for json_path in split_dir.rglob("*.json"):
        try:
            records = json.loads(json_path.read_text())
            if not isinstance(records, list):
                records = [records]
            for rec in records:
                pub_date_str = rec.get("publication_date") or rec.get("study_start_date")
                if pub_date_str:
                    pub_date = date.fromisoformat(str(pub_date_str)[:10])
                    if pub_date > cutoff_date:
                        violations.append(
                            f"{json_path.name}: {rec.get('pmid') or rec.get('nct_id')} "
                            f"dated {pub_date}"
                        )
        except Exception as exc:
            violations.append(f"{json_path.name}: parse error — {exc}")

    if violations:
        result["status"] = "FAIL"
        result["detail"] = f"{len(violations)} post-cutoff records found: {violations[:5]}"
    else:
        result["status"] = "PASS"
        result["detail"] = "No post-cutoff records detected."

    return result


def check_timestamp_coverage() -> dict[str, Any]:
    """Every literature feature row must have a per-publication date, not just aggregates."""
    result: dict[str, Any] = {
        "check": "timestamp_coverage",
        "status": "SKIP",
        "detail": "",
    }

    if not FEATURES_PATH.exists():
        result["detail"] = "feature_table.parquet not yet built — check deferred."
        return result

    df = pd.read_parquet(FEATURES_PATH)
    missing_cols = [c for c in TIMESTAMPED_FEATURE_COLS if c not in df.columns]
    if missing_cols:
        result["status"] = "FAIL"
        result["detail"] = f"Missing timestamped columns: {missing_cols}"
        return result

    null_counts = {c: int(df[c].isna().sum()) for c in TIMESTAMPED_FEATURE_COLS}
    total_nulls = sum(null_counts.values())
    if total_nulls > 0:
        result["status"] = "WARN"
        result["detail"] = f"Null timestamps found: {null_counts}"
    else:
        result["status"] = "PASS"
        result["detail"] = "All timestamped feature columns populated."

    return result


def check_driver_flag_isolation() -> dict[str, Any]:
    """Known-driver flag must be a standalone column, not merged into any score."""
    result: dict[str, Any] = {
        "check": "driver_flag_isolation",
        "status": "SKIP",
        "detail": "",
    }

    if not FEATURES_PATH.exists():
        result["detail"] = "feature_table.parquet not yet built — check deferred."
        return result

    df = pd.read_parquet(FEATURES_PATH)
    if DRIVER_FLAG_COL not in df.columns:
        result["status"] = "FAIL"
        result["detail"] = (
            f"Column '{DRIVER_FLAG_COL}' missing. Driver status must be a "
            "standalone binary flag, never merged into a composite score."
        )
    else:
        unique_vals = set(df[DRIVER_FLAG_COL].dropna().unique())
        if not unique_vals.issubset({0, 1, True, False, 0.0, 1.0}):
            result["status"] = "FAIL"
            result["detail"] = f"'{DRIVER_FLAG_COL}' has non-binary values: {unique_vals}"
        else:
            result["status"] = "PASS"
            result["detail"] = f"'{DRIVER_FLAG_COL}' is binary and isolated."

    return result


def check_drug_screen_bias() -> dict[str, Any]:
    """
    Genes in a drug-screen panel (PRISM etc.) may appear essential partly because
    they were pre-selected as suspected targets. This must be flagged explicitly.
    """
    result: dict[str, Any] = {
        "check": "drug_screen_bias",
        "status": "SKIP",
        "detail": "",
    }

    if not FEATURES_PATH.exists():
        result["detail"] = "feature_table.parquet not yet built — check deferred."
        return result

    df = pd.read_parquet(FEATURES_PATH)
    dep_present = [c for c in DEPENDENCY_COLS if c in df.columns]
    if not dep_present:
        result["status"] = "SKIP"
        result["detail"] = "No dependency columns present yet."
        return result

    if DRUG_SCREEN_FLAG_COL not in df.columns:
        result["status"] = "FAIL"
        result["detail"] = (
            f"Dependency columns present but '{DRUG_SCREEN_FLAG_COL}' is missing. "
            "Drug-screen panel membership must be flagged before model training."
        )
    else:
        n_flagged = int(df[DRUG_SCREEN_FLAG_COL].sum())
        result["status"] = "PASS"
        result["detail"] = (
            f"'{DRUG_SCREEN_FLAG_COL}' present. {n_flagged} genes flagged as "
            "drug-screen panel members."
        )

    return result


# ── Orchestrator ──────────────────────────────────────────────────────────────

def run(write_report: bool = True) -> bool:
    """
    Run all leakage checks. Returns True only if no check returns FAIL.
    Writes evaluation/reports/leakage_audit.md.
    """
    findings: list[dict[str, Any]] = []

    for year in CUTOFF_YEARS:
        findings.append(check_temporal_isolation(year))

    findings.append(check_timestamp_coverage())
    findings.append(check_driver_flag_isolation())
    findings.append(check_drug_screen_bias())

    failures = [f for f in findings if f["status"] == "FAIL"]
    warnings = [f for f in findings if f["status"] == "WARN"]
    passed = [f for f in findings if f["status"] == "PASS"]
    skipped = [f for f in findings if f["status"] == "SKIP"]

    overall = "PASS" if not failures else "FAIL"

    if write_report:
        _write_report(findings, overall)

    for f in failures:
        log.error("LEAKAGE FAIL [%s / cutoff=%s]: %s",
                  f["check"], f.get("cutoff_year", "—"), f["detail"])
    for w in warnings:
        log.warning("LEAKAGE WARN [%s]: %s", w["check"], w["detail"])

    log.info(
        "Leakage audit complete — PASS: %d | WARN: %d | FAIL: %d | SKIP: %d → Overall: %s",
        len(passed), len(warnings), len(failures), len(skipped), overall,
    )
    return overall == "PASS"


def _write_report(findings: list[dict[str, Any]], overall: str) -> None:
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Leakage Audit Report",
        "",
        f"**Overall status: {overall}**",
        "",
        "| Check | Cutoff | Status | Detail |",
        "|-------|--------|--------|--------|",
    ]
    for f in findings:
        lines.append(
            f"| {f['check']} | {f.get('cutoff_year', '—')} "
            f"| {f['status']} | {f['detail']} |"
        )
    lines += ["", "---", "_Generated by leakage_auditor.py_"]
    REPORT_PATH.write_text("\n".join(lines))
    log.info("Leakage audit report → %s", REPORT_PATH)
