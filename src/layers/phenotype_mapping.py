"""phenotype_mapping.py — HP low-level → MP system-level mapping via HP2MP anchor terms.

Implements the approach from Barbitoff et al. (2025):
  HP low-level  →(hp.obo ancestors)→  HP anchor term  →(HP2MP.tsv)→  MP system-level

This script is intentionally decoupled from the export pipeline so that
alternative mapping strategies (SSSOM, semantic similarity, etc.) can be
swapped in without touching the export or adapter layers.

Outputs written to out_dir:
  edge_hp_to_mp_top.tsv   — (hp_id, mp_top_id): every HP term → MP system-level terms
  node_mp_top_names.tsv   — (mp_id, name): names for all MP system-level terms used
  qc_hp_no_anchor.tsv     — HP terms with no anchor in HP2MP (QC)
  qc_hp_no_mapping.tsv    — HP terms with anchor but no MP mapping (should be empty)
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

MORTALITY_AGING_MP = "MP:0010768"


# ── 1. Loaders ────────────────────────────────────────────────────────────────

def load_hp_ontology(data_dir: Path):
    """Load HP ontology from hp.obo. Import is lazy to avoid slow parse at import time."""
    from pronto import Ontology
    hp_path = data_dir / "hp.obo"
    if not hp_path.exists():
        raise FileNotFoundError(hp_path)
    logger.debug("Loading HP ontology from %s", hp_path)
    return Ontology(str(hp_path), encoding="utf-8")


def load_hp2mp(hp2mp_path: Path, exclude_mortality: bool = True) -> dict[str, list[str]]:
    """Load HP2MP.tsv and return dict {hp_anchor_id: [mp_top_id, ...]}.

    HP2MP.tsv uses object_id = HP anchor term, subject_id = MP system-level term.
    Following Barbitoff et al. (2025), mortality/aging (MP:0010768) is excluded
    by default due to strong annotation bias.
    """
    df = pd.read_csv(hp2mp_path, sep="\t", dtype=str)

    required = {"object_id", "subject_id"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"HP2MP missing columns: {missing}. Got: {list(df.columns)}")

    if exclude_mortality:
        before = len(df)
        df = df[df["subject_id"] != MORTALITY_AGING_MP]
        logger.debug("Excluded mortality/aging term: %d rows removed", before - len(df))

    hp2mp: dict[str, list[str]] = {}
    for _, row in df.iterrows():
        hp_anchor = row["object_id"]
        mp_top    = row["subject_id"]
        if not (hp_anchor.startswith("HP:") and mp_top.startswith("MP:")):
            continue
        hp2mp.setdefault(hp_anchor, []).append(mp_top)

    logger.info(
        "HP2MP loaded: %d anchor HP terms → %d unique MP system-level terms",
        len(hp2mp),
        len({mp for mps in hp2mp.values() for mp in mps}),
    )
    return hp2mp


# ── 2. Core lift function ─────────────────────────────────────────────────────

def lift_hp_to_anchors(hp_id: str, hp_ont, anchor_set: set[str]) -> list[str]:
    """Return HP2MP anchor terms that are ancestors (or self) of hp_id.

    The anchor set is the set of HP terms present in HP2MP.tsv as object_id.
    These are not necessarily canonical top-level terms (distance=1 from root)
    — they are MGI-curated pivot points at various ontology levels.
    """
    if hp_id not in hp_ont:
        return []
    if hp_id in anchor_set:
        return [hp_id]
    ancestors = {t.id for t in hp_ont[hp_id].superclasses(with_self=False)}
    return sorted(ancestors & anchor_set)


# ── 3. Main mapping builder ───────────────────────────────────────────────────

def build_hp_to_mp_top(
    hp_ids: list[str],
    hp_ont,
    hp2mp: dict[str, list[str]],
) -> pd.DataFrame:
    """Map a list of HP term IDs to MP system-level terms.

    For each hp_id:
      1. Find HP anchor ancestors via ontology hierarchy.
      2. Map each anchor to MP system-level terms via HP2MP.tsv.

    Returns DataFrame with columns [hp_id, mp_top_id].
    """
    anchor_set = set(hp2mp.keys())
    rows = []

    for hp_id in sorted(set(hp_ids)):
        anchors = lift_hp_to_anchors(hp_id, hp_ont, anchor_set)
        mp_tops = {mp for a in anchors for mp in hp2mp.get(a, [])}
        for mp_top in sorted(mp_tops):
            rows.append({"hp_id": hp_id, "mp_top_id": mp_top})

    df = pd.DataFrame(rows) if rows else pd.DataFrame(columns=["hp_id", "mp_top_id"])
    logger.info(
        "HP→MP top mapping: %d / %d unique HP terms mapped (%d edges total)",
        df["hp_id"].nunique() if not df.empty else 0,
        len(set(hp_ids)),
        len(df),
    )
    return df


# ── 4. End-to-end pipeline ────────────────────────────────────────────────────

def run_mapping(
    genes_to_phenotype_path: str | Path,
    hp2mp_path: str | Path,
    data_dir: str | Path,
    out_dir: str | Path = "./out",
    exclude_mortality: bool = True,
    mp_obo_path: str | Path | None = None,
) -> dict[str, pd.DataFrame]:
    """Run HP→MP system-level mapping and write output TSVs.

    Parameters
    ----------
    genes_to_phenotype_path:
        HPO genes_to_phenotype.txt file.
    hp2mp_path:
        HP2MP.tsv anchor mapping file from MGI.
    data_dir:
        Directory containing hp.obo (and optionally mp.obo for term names).
    out_dir:
        Directory to write output TSVs.
    exclude_mortality:
        Remove MP:0010768 (mortality/aging) following Barbitoff et al. (2025).
    mp_obo_path:
        Optional path to mp.obo for resolving MP term names.
        If None, looks for mp.obo in data_dir.

    Outputs
    -------
    edge_hp_to_mp_top.tsv      hp_id, mp_top_id
    node_mp_top_names.tsv      mp_id, name
    qc_hp_no_anchor.tsv        hp_id  (HP terms with no anchor in HP2MP)
    qc_hp_no_mapping.tsv       hp_id  (HP terms with anchor but no MP mapping — should be empty)
    """
    out_dir  = Path(out_dir)
    data_dir = Path(data_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=== Phenotype mapping pipeline (Barbitoff et al. approach) ===")
    logger.info("genes_to_phenotype : %s", genes_to_phenotype_path)
    logger.info("HP2MP              : %s", hp2mp_path)
    logger.info("data_dir           : %s", data_dir)
    logger.info("out_dir            : %s", out_dir)
    logger.info("exclude_mortality  : %s", exclude_mortality)

    # ── Load ontology and mapping table ──────────────────────────────────────
    hp_ont = load_hp_ontology(data_dir)
    hp2mp  = load_hp2mp(Path(hp2mp_path), exclude_mortality=exclude_mortality)

    # ── Read input gene-phenotype associations ────────────────────────────────
    g2p = pd.read_csv(genes_to_phenotype_path, sep="\t", dtype=str)
    required = {"hpo_id", "gene_symbol"}
    missing = required - set(g2p.columns)
    if missing:
        raise ValueError(f"genes_to_phenotype missing columns: {missing}")

    unique_hp = list(g2p["hpo_id"].dropna().unique())
    logger.info("Unique HP terms in input: %d", len(unique_hp))

    # ── Build HP → MP top mapping ─────────────────────────────────────────────
    hp_to_mp_top = build_hp_to_mp_top(unique_hp, hp_ont, hp2mp)
    hp_to_mp_top.to_csv(out_dir / "edge_hp_to_mp_top.tsv", sep="\t", index=False)
    logger.info("Written: edge_hp_to_mp_top.tsv (%d rows)", len(hp_to_mp_top))

    # ── MP top term names ─────────────────────────────────────────────────────
    mp_obo_path = Path(mp_obo_path) if mp_obo_path else data_dir / "mp.obo"
    if mp_obo_path.exists():
        from pronto import Ontology
        logger.debug("Loading MP ontology for term names from %s", mp_obo_path)
        mp_ont = Ontology(str(mp_obo_path), encoding="utf-8")
        mp_top_ids = sorted(hp_to_mp_top["mp_top_id"].dropna().unique()) if not hp_to_mp_top.empty else []
        names_rows = [
            {"mp_id": mp_id, "name": (mp_ont[mp_id].name or "") if mp_id in mp_ont else ""}
            for mp_id in mp_top_ids
        ]
    else:
        logger.warning("mp.obo not found at %s — names will be empty", mp_obo_path)
        mp_top_ids = sorted(hp_to_mp_top["mp_top_id"].dropna().unique()) if not hp_to_mp_top.empty else []
        names_rows = [{"mp_id": mp_id, "name": ""} for mp_id in mp_top_ids]

    mp_top_names = pd.DataFrame(names_rows)
    mp_top_names.to_csv(out_dir / "node_mp_top_names.tsv", sep="\t", index=False)
    logger.info("Written: node_mp_top_names.tsv (%d terms)", len(mp_top_names))

    # ── QC outputs ────────────────────────────────────────────────────────────
    anchor_set  = set(hp2mp.keys())
    mapped_hp   = set(hp_to_mp_top["hp_id"].dropna()) if not hp_to_mp_top.empty else set()

    # HP terms with no anchor at all in HP2MP (cannot be lifted)
    no_anchor = [
        hp for hp in unique_hp
        if not lift_hp_to_anchors(hp, hp_ont, anchor_set)
    ]
    pd.DataFrame({"hp_id": sorted(no_anchor)}).to_csv(
        out_dir / "qc_hp_no_anchor.tsv", sep="\t", index=False
    )

    # HP terms that have an anchor but still got no MP mapping (should be empty)
    has_anchor_no_mp = [
        hp for hp in unique_hp
        if lift_hp_to_anchors(hp, hp_ont, anchor_set) and hp not in mapped_hp
    ]
    pd.DataFrame({"hp_id": sorted(has_anchor_no_mp)}).to_csv(
        out_dir / "qc_hp_no_mapping.tsv", sep="\t", index=False
    )

    logger.info("QC: HP terms with no anchor     : %d / %d (%.1f%%)",
                len(no_anchor), len(unique_hp), len(no_anchor) / len(unique_hp) * 100)
    logger.info("QC: HP terms anchor→no MP top   : %d", len(has_anchor_no_mp))
    logger.info("HP terms successfully mapped     : %d / %d (%.1f%%)",
                len(mapped_hp), len(unique_hp), len(mapped_hp) / len(unique_hp) * 100)
    logger.info("MP system-level terms used       : %d", len(mp_top_names))

    return {
        "hp_to_mp_top":  hp_to_mp_top,
        "mp_top_names":  mp_top_names,
    }


# ── CLI ───────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import argparse
    import logging
    from src.utils.logging_config import setup_logging

    p = argparse.ArgumentParser(
        description="Map HP low-level terms to MP system-level terms via HP2MP anchor table."
    )
    p.add_argument("--genes-to-phenotype", required=True,
                   help="Path to genes_to_phenotype.txt (HPO)")
    p.add_argument("--hp2mp", required=True,
                   help="Path to HP2MP.tsv (MGI anchor mapping)")
    p.add_argument("--data-dir", required=True,
                   help="Directory containing hp.obo (and optionally mp.obo)")
    p.add_argument("--out", default="./out",
                   help="Output directory (default: ./out)")
    p.add_argument("--keep-mortality", action="store_true",
                   help="Do NOT exclude MP:0010768 (mortality/aging)")
    p.add_argument("--mp-obo", default=None,
                   help="Explicit path to mp.obo (defaults to <data-dir>/mp.obo)")
    p.add_argument("--log-level", default="INFO",
                   choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    args = p.parse_args()

    setup_logging(out_dir=args.out, console_level=getattr(logging, args.log_level))

    run_mapping(
        genes_to_phenotype_path=args.genes_to_phenotype,
        hp2mp_path=args.hp2mp,
        data_dir=args.data_dir,
        out_dir=args.out,
        exclude_mortality=not args.keep_mortality,
        mp_obo_path=args.mp_obo,
    )