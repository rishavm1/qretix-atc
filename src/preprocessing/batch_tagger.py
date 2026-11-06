"""
batch_tagger.py
───────────────
Phase 2.4: Verify and enforce platform/batch tags on all processed files.

Every sample in every dataset must carry platform and batch columns
through all downstream steps. This module audits processed files and
patches any missing tags.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

log = logging.getLogger(__name__)

ROOT     = Path(__file__).resolve().parents[2]
PROC_DIR = ROOT / "data" / "processed"

# Expected tags per file
EXPECTED_TAGS: dict[str, dict] = {
    "GSE76039_metadata.parquet": {
        "platform": "Affymetrix_HGU133Plus2_GPL570",
        "batch": "GSE76039",
    },
    "atc_mutations.parquet": {
        "platform": "DepMap_precalled_somatic",
        "batch": "depmap_24Q4",
    },
    "atc_cna.parquet": {
        "platform": "DepMap_precalled_CNA",
        "batch": "depmap_24Q4",
    },
    "depmap_atc_dependency.csv": {
        "platform": None,   # no platform col in this file — noted
        "batch": None,
    },
}


def audit_and_patch() -> dict[str, str]:
    """
    Check each processed file for platform/batch columns.
    Patch if missing. Returns {filename: status}.
    """
    results = {}

    for fname, tags in EXPECTED_TAGS.items():
        fpath = PROC_DIR / fname
        if not fpath.exists():
            results[fname] = "MISSING_FILE"
            log.warning("batch_tagger: %s not found.", fname)
            continue

        if tags["platform"] is None:
            results[fname] = "SKIP_NO_TAG_REQUIRED"
            continue

        try:
            if fname.endswith(".parquet"):
                df = pd.read_parquet(fpath)
            else:
                df = pd.read_csv(fpath, nrows=5)

            patched = False
            for col, val in tags.items():
                if col not in df.columns:
                    log.warning("%s missing column '%s' — patching.", fname, col)
                    # Reload full file, add column, resave
                    if fname.endswith(".parquet"):
                        full = pd.read_parquet(fpath)
                    else:
                        full = pd.read_csv(fpath)
                    full[col] = val
                    if fname.endswith(".parquet"):
                        full.to_parquet(fpath, index=False)
                    else:
                        full.to_csv(fpath, index=False)
                    patched = True

            results[fname] = "PATCHED" if patched else "OK"
            log.info("batch_tagger: %s → %s", fname, results[fname])

        except Exception as exc:
            results[fname] = f"ERROR: {exc}"
            log.error("batch_tagger: %s failed: %s", fname, exc)

    return results


def run() -> None:
    results = audit_and_patch()
    failures = {k: v for k, v in results.items()
                if v.startswith("ERROR") or v == "MISSING_FILE"}
    if failures:
        log.warning("batch_tagger: %d issues: %s", len(failures), failures)
    else:
        log.info("batch_tagger: all files tagged correctly.")
