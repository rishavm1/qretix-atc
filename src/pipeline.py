"""
pipeline.py
───────────
Lightweight DAG orchestrator for the QRETIX-ATC v0.1 pipeline.

Each phase is registered with:
  - run()      : executes the phase
  - validate() : checks acceptance criteria (returns bool)
  - deps       : list of phase numbers that must pass before this phase runs

Usage:
  python -m src.pipeline run    --phase 1
  python -m src.pipeline run    --all
  python -m src.pipeline validate --phase 2
  python -m src.pipeline validate --all
  python -m src.pipeline status
"""

from __future__ import annotations

import argparse
import logging
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[1]


# ── Phase descriptor ──────────────────────────────────────────────────────────

@dataclass
class Phase:
    number: int
    name: str
    deps: list[int] = field(default_factory=list)
    _run_fn: Callable[[], None] | None = field(default=None, repr=False)
    _validate_fn: Callable[[], bool] | None = field(default=None, repr=False)

    def run(self) -> None:
        if self._run_fn is None:
            log.warning("Phase %d (%s): no run() registered yet.", self.number, self.name)
            return
        log.info("── Running Phase %d: %s ──", self.number, self.name)
        self._run_fn()

    def validate(self) -> bool:
        if self._validate_fn is None:
            log.warning("Phase %d (%s): no validate() registered yet.", self.number, self.name)
            return False
        result = self._validate_fn()
        status = "PASS" if result else "FAIL"
        log.info("Phase %d (%s) validation: %s", self.number, self.name, status)
        return result


# ── Registry ──────────────────────────────────────────────────────────────────

_PHASES: dict[int, Phase] = {}


def register(
    number: int,
    name: str,
    deps: list[int] | None = None,
    run_fn: Callable[[], None] | None = None,
    validate_fn: Callable[[], bool] | None = None,
) -> Phase:
    p = Phase(
        number=number,
        name=name,
        deps=deps or [],
        _run_fn=run_fn,
        _validate_fn=validate_fn,
    )
    _PHASES[number] = p
    return p


def get(number: int) -> Phase:
    if number not in _PHASES:
        raise KeyError(f"Phase {number} not registered.")
    return _PHASES[number]


# ── Topological execution ─────────────────────────────────────────────────────

def _topo_order(target: int) -> list[int]:
    """Return phases in dependency-first order up to and including target."""
    visited: list[int] = []
    seen: set[int] = set()

    def visit(n: int) -> None:
        if n in seen:
            return
        seen.add(n)
        for dep in _PHASES.get(n, Phase(n, "")).deps:
            visit(dep)
        visited.append(n)

    visit(target)
    return visited


def run_phase(number: int, skip_validation: bool = False) -> bool:
    """Run a single phase (after validating its deps)."""
    phase = get(number)

    for dep in phase.deps:
        dep_phase = get(dep)
        if not dep_phase.validate():
            log.error(
                "Phase %d (%s) blocked: dependency Phase %d (%s) failed validation.",
                number, phase.name, dep, dep_phase.name,
            )
            return False

    phase.run()

    if not skip_validation:
        return phase.validate()
    return True


def run_all(skip_validation: bool = False) -> bool:
    """Run all registered phases in dependency order."""
    order = sorted(_PHASES.keys())
    for n in order:
        ok = run_phase(n, skip_validation=skip_validation)
        if not ok:
            log.error("Pipeline halted at Phase %d.", n)
            return False
    return True


def validate_phase(number: int) -> bool:
    return get(number).validate()


def validate_all() -> bool:
    return all(get(n).validate() for n in sorted(_PHASES.keys()))


def status() -> None:
    """Print a summary table of all registered phases and their validation state."""
    print(f"\n{'#':<4} {'Phase':<40} {'Deps':<15} {'Valid?'}")
    print("─" * 72)
    for n in sorted(_PHASES.keys()):
        p = _PHASES[n]
        deps_str = str(p.deps) if p.deps else "—"
        try:
            valid = "✓" if p.validate() else "✗"
        except Exception:
            valid = "?"
        print(f"{n:<4} {p.name:<40} {deps_str:<15} {valid}")
    print()


# ── Phase registrations ───────────────────────────────────────────────────────
# Each phase imports its own module lazily to avoid circular imports and to
# allow partial pipeline execution before all modules are implemented.

def _register_all() -> None:

    # Phase 0 — Scaffolding
    def _validate_phase0() -> bool:
        from src import gene_universe, provenance  # noqa: F401
        required_dirs = [
            ROOT / "data" / "raw",
            ROOT / "data" / "processed",
            ROOT / "data" / "provenance",
            ROOT / "features",
            ROOT / "labels",
            ROOT / "knowledge_graph",
            ROOT / "literature",
            ROOT / "models" / "baselines",
            ROOT / "models" / "qretix_ranker",
            ROOT / "evaluation" / "temporal_splits",
            ROOT / "evaluation" / "reports",
            ROOT / "evidence_cards",
            ROOT / "configs",
            ROOT / "docs",
        ]
        missing = [d for d in required_dirs if not d.exists()]
        if missing:
            log.error("Phase 0 FAIL: missing dirs: %s", missing)
            return False
        required_files = [
            ROOT / "requirements.txt",
            ROOT / "src" / "gene_universe.py",
            ROOT / "src" / "provenance.py",
            ROOT / "src" / "leakage_auditor.py",
            ROOT / "src" / "pipeline.py",
            ROOT / "docs" / "QRETIX-ATC-v0.1-spec.md",
            ROOT / "docs" / "decisions_log.md",
        ]
        missing_files = [f for f in required_files if not f.exists()]
        if missing_files:
            log.error("Phase 0 FAIL: missing files: %s", missing_files)
            return False
        log.info("Phase 0 PASS: all directories and core files present.")
        return True

    register(0, "Project Scaffolding", deps=[], validate_fn=_validate_phase0)

    # Phase 1 — Data Acquisition
    def _validate_phase1() -> bool:
        from src import provenance as prov
        prov_dir = ROOT / "data" / "provenance"
        if not any(prov_dir.glob("*.json")):
            log.error("Phase 1 FAIL: no provenance JSON files found.")
            return False
        inventory = ROOT / "docs" / "dataset_inventory.md"
        if not inventory.exists():
            log.error("Phase 1 FAIL: docs/dataset_inventory.md missing.")
            return False
        feature_files = list((ROOT / "features").glob("*.parquet"))
        if feature_files:
            log.error("Phase 1 FAIL: feature files found — Phase 1 must not create features.")
            return False
        return prov.validate_all()

    register(1, "Data Inventory & Acquisition", deps=[0], validate_fn=_validate_phase1)

    # Phase 2 — QC & Preprocessing
    def _validate_phase2() -> bool:
        qc_report = ROOT / "docs" / "qc_report.md"
        if not qc_report.exists():
            log.error("Phase 2 FAIL: docs/qc_report.md missing.")
            return False
        processed = list((ROOT / "data" / "processed").glob("*.parquet")) + \
                    list((ROOT / "data" / "processed").glob("*.h5ad"))
        if not processed:
            log.error("Phase 2 FAIL: no processed files in data/processed/.")
            return False
        log.info("Phase 2 PASS: QC report present, %d processed files.", len(processed))
        return True

    register(2, "QC & Preprocessing", deps=[1], validate_fn=_validate_phase2)

    # Phase 3 — Feature Engineering
    def _validate_phase3() -> bool:
        ft = ROOT / "features" / "feature_table.parquet"
        spec = ROOT / "docs" / "feature_spec.md"
        if not ft.exists():
            log.error("Phase 3 FAIL: features/feature_table.parquet missing.")
            return False
        if not spec.exists():
            log.error("Phase 3 FAIL: docs/feature_spec.md missing.")
            return False
        import pandas as pd
        df = pd.read_parquet(ft)
        if "ensembl_gene_id" not in df.columns:
            log.error("Phase 3 FAIL: feature_table missing ensembl_gene_id primary key.")
            return False
        log.info("Phase 3 PASS: feature_table has %d genes, %d features.", len(df), len(df.columns))
        return True

    register(3, "Feature Engineering", deps=[2], validate_fn=_validate_phase3)

    # Phase 4 — Label Construction
    def _validate_phase4() -> bool:
        labels = ROOT / "labels" / "labels.parquet"
        spec = ROOT / "labels" / "label_spec.md"
        dist = ROOT / "labels" / "class_distribution.md"
        for p in [labels, spec, dist]:
            if not p.exists():
                log.error("Phase 4 FAIL: %s missing.", p)
                return False
        import pandas as pd
        df = pd.read_parquet(labels)
        required_cols = {"ensembl_gene_id", "stage", "earliest_evidence_date",
                         "evidence_sources", "is_true_negative", "is_hard_negative"}
        missing = required_cols - set(df.columns)
        if missing:
            log.error("Phase 4 FAIL: labels.parquet missing columns: %s", missing)
            return False
        pos_frac = (df["stage"] >= 1).mean()
        if pos_frac > 0.02:
            log.error(
                "Phase 4 FAIL: positive class fraction %.3f > 2%% — check for leakage.",
                pos_frac,
            )
            return False
        log.info("Phase 4 PASS: %d labeled genes, positive fraction %.3f%%.",
                 len(df), pos_frac * 100)
        return True

    register(4, "Label Construction", deps=[1, 3], validate_fn=_validate_phase4)

    # Phase 5 — Modular Scoring / Models
    def _validate_phase5() -> bool:
        report = ROOT / "evaluation" / "reports" / "model_comparison.md"
        if not report.exists():
            log.error("Phase 5 FAIL: evaluation/reports/model_comparison.md missing.")
            return False
        configs = list((ROOT / "configs").glob("*.json"))
        if not configs:
            log.error("Phase 5 FAIL: no model config files in configs/.")
            return False
        log.info("Phase 5 PASS: model comparison report present, %d configs.", len(configs))
        return True

    register(5, "Modular Scoring Architecture", deps=[3, 4], validate_fn=_validate_phase5)

    # Phase 6 — Novelty Engine
    def _validate_phase6() -> bool:
        novelty = ROOT / "features" / "novelty_scores.parquet"
        if not novelty.exists():
            log.error("Phase 6 FAIL: features/novelty_scores.parquet missing.")
            return False
        import pandas as pd
        df = pd.read_parquet(novelty)
        required = {
            "novelty_disease", "novelty_target", "novelty_mechanistic",
            "novelty_combination", "novelty_saturation",
        }
        missing = required - set(df.columns)
        if missing:
            log.error("Phase 6 FAIL: novelty_scores missing sub-scores: %s", missing)
            return False
        log.info("Phase 6 PASS: all 5 novelty sub-scores present for %d genes.", len(df))
        return True

    register(6, "Novelty Engine", deps=[1, 3], validate_fn=_validate_phase6)

    # Phase 7 — Temporal Retrospective Validation
    def _validate_phase7() -> bool:
        all_pass = True
        for year in [2015, 2017, 2019, 2021]:
            results = ROOT / "evaluation" / "temporal_splits" / str(year) / "results.md"
            if not results.exists():
                log.error("Phase 7 FAIL: %s missing.", results)
                all_pass = False
        if all_pass:
            log.info("Phase 7 PASS: all 4 temporal split result files present.")
        return all_pass

    register(7, "Temporal Retrospective Validation", deps=[5], validate_fn=_validate_phase7)

    # Phase 8 — Evidence Cards
    def _validate_phase8() -> bool:
        cards = list((ROOT / "evidence_cards").glob("*.md"))
        if not cards:
            log.error("Phase 8 FAIL: no evidence cards generated.")
            return False
        # Spot-check: every card must contain the disclaimer and a Supporting Evidence section
        failures = []
        for card in cards:
            text = card.read_text()
            if "Computational hypothesis" not in text:
                failures.append(f"{card.name}: missing disclaimer")
            if "Supporting Evidence" not in text:
                failures.append(f"{card.name}: missing Supporting Evidence section")
        if failures:
            log.error("Phase 8 FAIL: %d cards failed spot-check: %s", len(failures), failures[:3])
            return False
        log.info("Phase 8 PASS: %d evidence cards, all passed spot-check.", len(cards))
        return True

    register(8, "Evidence Cards", deps=[5, 6], validate_fn=_validate_phase8)

    # Phase 9 — Benchmarking
    def _run_phase9() -> None:
        from src.benchmarking import benchmark_suite
        benchmark_suite.run()
        from src.benchmarking import benchmark_export
        benchmark_export.run()

    def _validate_phase9() -> bool:
        report = ROOT / "evaluation" / "reports" / "benchmark_comparison.md"
        if not report.exists():
            log.error("Phase 9 FAIL: evaluation/reports/benchmark_comparison.md missing.")
            return False
        text = report.read_text()
        required_methods = [
            "DE-only", "DepMap-only", "Open Targets", "Network centrality", "QRETIX"
        ]
        missing = [m for m in required_methods if m not in text]
        if missing:
            log.error("Phase 9 FAIL: benchmark report missing methods: %s", missing)
            return False
        release_dir = ROOT / "evaluation" / "benchmark_release"
        required_release = [
            "benchmark_metadata.json", "ground_truth.parquet",
            "eval_script.py", "q3dc_submission_form.md",
        ]
        missing_release = [f for f in required_release
                           if not (release_dir / f).exists()]
        if missing_release:
            log.error("Phase 9 FAIL: benchmark_release missing: %s", missing_release)
            return False
        log.info("Phase 9 PASS: benchmark report + Q3DC release package present.")
        return True

    register(9, "Benchmarking", deps=[5], run_fn=_run_phase9, validate_fn=_validate_phase9)


# ── CLI ───────────────────────────────────────────────────────────────────────

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m src.pipeline",
        description="QRETIX-ATC v0.1 pipeline orchestrator",
    )
    sub = parser.add_subparsers(dest="command")

    run_p = sub.add_parser("run", help="Execute a phase or all phases")
    run_p.add_argument("--phase", type=int, help="Phase number to run")
    run_p.add_argument("--all", action="store_true", help="Run all phases in order")

    val_p = sub.add_parser("validate", help="Validate a phase or all phases")
    val_p.add_argument("--phase", type=int, help="Phase number to validate")
    val_p.add_argument("--all", action="store_true", help="Validate all phases")

    sub.add_parser("status", help="Print phase status table")

    return parser


def main() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)-8s %(message)s",
        datefmt="%H:%M:%S",
    )
    _register_all()
    parser = _build_parser()
    args = parser.parse_args()

    if args.command == "run":
        if args.all:
            ok = run_all()
        elif args.phase is not None:
            ok = run_phase(args.phase)
        else:
            parser.print_help()
            sys.exit(1)
        sys.exit(0 if ok else 1)

    elif args.command == "validate":
        if args.all:
            ok = validate_all()
        elif args.phase is not None:
            ok = validate_phase(args.phase)
        else:
            parser.print_help()
            sys.exit(1)
        sys.exit(0 if ok else 1)

    elif args.command == "status":
        status()

    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
