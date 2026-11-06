"""
provenance.py
─────────────
Utilities for writing, loading, updating, and validating per-dataset provenance
JSON records stored in data/provenance/.

Every dataset that enters the pipeline must have a provenance record before any
downstream processing. The schema is enforced here — no partial records allowed.
"""

from __future__ import annotations

import hashlib
import json
import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[1]
PROVENANCE_DIR = ROOT / "data" / "provenance"

# ── Required fields (Phase 1 acceptance criterion) ────────────────────────────
REQUIRED_FIELDS = {
    "source",
    "accession",
    "url",
    "date_pulled",
    "sample_count",
    "modality",
    "license",
    "raw_files",
    "checksums",
    "notes",
    "qc_thresholds_applied",
    "is_processed",
}


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return f"sha256:{h.hexdigest()}"


def record_path(accession: str) -> Path:
    return PROVENANCE_DIR / f"{accession}.json"


def write(accession: str, record: dict[str, Any], overwrite: bool = False) -> Path:
    """
    Write a provenance JSON for the given accession.
    Raises ValueError if required fields are missing.
    Raises FileExistsError if record already exists and overwrite=False.
    """
    PROVENANCE_DIR.mkdir(parents=True, exist_ok=True)
    missing = REQUIRED_FIELDS - set(record.keys())
    if missing:
        raise ValueError(f"Provenance record for {accession} missing fields: {missing}")

    path = record_path(accession)
    if path.exists() and not overwrite:
        raise FileExistsError(
            f"Provenance record already exists: {path}. Pass overwrite=True to replace."
        )

    record["accession"] = accession  # ensure accession is embedded
    path.write_text(json.dumps(record, indent=2, default=str))
    log.info("Provenance written → %s", path)
    return path


def load(accession: str) -> dict[str, Any]:
    """Load and return a provenance record. Raises FileNotFoundError if absent."""
    path = record_path(accession)
    if not path.exists():
        raise FileNotFoundError(f"No provenance record for {accession} at {path}")
    return json.loads(path.read_text())


def update(accession: str, updates: dict[str, Any]) -> None:
    """Merge updates into an existing provenance record (e.g., add QC thresholds)."""
    record = load(accession)
    record.update(updates)
    record_path(accession).write_text(json.dumps(record, indent=2, default=str))
    log.info("Provenance updated for %s.", accession)


def add_qc_thresholds(accession: str, thresholds: dict[str, Any]) -> None:
    """Convenience wrapper: set qc_thresholds_applied on an existing record."""
    update(accession, {"qc_thresholds_applied": thresholds})


def checksum_files(file_paths: list[Path]) -> dict[str, str]:
    """Return {filename: sha256_hex} for a list of paths."""
    return {p.name: _sha256(p) for p in file_paths if p.exists()}


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def make_record(
    source: str,
    accession: str,
    url: str,
    sample_count: int,
    modality: str,
    license_: str,
    raw_files: list[str],
    checksums: dict[str, str],
    notes: str,
    is_processed: bool = False,
    qc_thresholds_applied: dict | None = None,
    **extra: Any,
) -> dict[str, Any]:
    """
    Build a provenance dict from named arguments.
    Pass to write() to persist it.
    """
    record: dict[str, Any] = {
        "source": source,
        "accession": accession,
        "url": url,
        "date_pulled": now_iso(),
        "sample_count": sample_count,
        "modality": modality,
        "license": license_,
        "raw_files": raw_files,
        "checksums": checksums,
        "notes": notes,
        "qc_thresholds_applied": qc_thresholds_applied,
        "is_processed": is_processed,
    }
    record.update(extra)
    return record


def validate_all() -> bool:
    """
    Phase 1 acceptance check: every JSON in data/provenance/ has all required fields.
    Returns True only if ALL records pass.
    """
    PROVENANCE_DIR.mkdir(parents=True, exist_ok=True)
    records = list(PROVENANCE_DIR.glob("*.json"))
    if not records:
        log.warning("No provenance records found in %s.", PROVENANCE_DIR)
        return False

    all_pass = True
    for path in records:
        record = json.loads(path.read_text())
        missing = REQUIRED_FIELDS - set(record.keys())
        if missing:
            log.error("FAIL [%s]: missing fields %s", path.name, missing)
            all_pass = False
        else:
            log.info("PASS [%s]", path.name)

    return all_pass


def list_accessions() -> list[str]:
    """Return all accessions that have provenance records."""
    return [p.stem for p in PROVENANCE_DIR.glob("*.json")]
