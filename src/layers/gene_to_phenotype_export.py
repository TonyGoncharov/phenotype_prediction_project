from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import pandas as pd

# pronto is imported lazily inside load_ontologies() so that importing
# this module does NOT trigger slow OBO parsing at import time.

HP_ROOT = "HP:0000118"
MP_ROOT = "MP:0000001"


# ----------------------------
# 0) Ontology helpers
# ----------------------------

def load_ontologies(data_dir: Path):
    from pronto import Ontology
    hp_path = data_dir / "hp.obo"
    mp_path = data_dir / "mp.obo"
    if not hp_path.exists():
        raise FileNotFoundError(hp_path)
    if not mp_path.exists():
        raise FileNotFoundError(mp_path)
    return Ontology(str(hp_path)), Ontology(str(mp_path))


def top_children(ont, root_id: str) -> set[str]:
    root = ont[root_id]
    tops = {t.id for t in root.subclasses(distance=1)}
    tops.discard(root_id)
    return tops


def ancestors(ont, term_id: str) -> set[str]:
    return {t.id for t in ont[term_id].superclasses()}


def collapse_to_top(ont, term_id: str, top_set: set[str], root_id: str) -> list[str]:
    if term_id not in ont:
        return []
    anc = ancestors(ont, term_id)
    anc.add(term_id)
    return sorted((anc & top_set) - {root_id})


def build_top_name_table(ont, top_set: set[str], id_col: str) -> pd.DataFrame:
    """Return a DataFrame {<id_col>, name} for all top-level terms.

    pronto stores the human-readable label in term.name.
    Terms with no name get an empty string.
    """
    rows = [
        {id_col: tid, "name": (ont[tid].name or "") if tid in ont else ""}
        for tid in sorted(top_set)
    ]
    return pd.DataFrame(rows)


# ----------------------------
# 1) SSSOM mapping loader (HP -> MP)
# ----------------------------

@dataclass(frozen=True)
class Mapping:
    mp_id: str
    predicate_id: Optional[str] = None
    confidence: Optional[float] = None
    source: Optional[str] = None


def load_hp_to_mp_sssom(sssom_path: str | Path) -> dict[str, list[Mapping]]:
    sssom_path = Path(sssom_path)
    df = pd.read_csv(sssom_path, sep="\t", dtype=str, comment="#")

    required = {"subject_id", "object_id"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"SSSOM missing required columns: {missing}. Got: {list(df.columns)}")

    conf_col = "confidence" if "confidence" in df.columns else None
    pred_col = "predicate_id" if "predicate_id" in df.columns else None
    src_col = next(
        (c for c in ("mapping_source", "creator_id", "mapping_provider", "source")
         if c in df.columns),
        None,
    )

    hp2mp: dict[str, list[Mapping]] = defaultdict(list)
    for _, r in df.iterrows():
        subj, obj = r["subject_id"], r["object_id"]
        if subj.startswith("HP:") and obj.startswith("MP:"):
            hp_id, mp_id = subj, obj
        elif subj.startswith("MP:") and obj.startswith("HP:"):
            hp_id, mp_id = obj, subj
        else:
            continue

        hp2mp[hp_id].append(Mapping(
            mp_id=mp_id,
            predicate_id=r[pred_col] if pred_col else None,
            confidence=float(r[conf_col]) if (conf_col and pd.notna(r[conf_col])) else None,
            source=r[src_col] if (src_col and pd.notna(r[src_col])) else None,
        ))

    return dict(hp2mp)


def rank_mappings(mappings: list[Mapping]) -> list[Mapping]:
    def score(m: Mapping) -> tuple[int, float]:
        p = (m.predicate_id or "").lower()
        if "equivalent" in p or "exact" in p:
            s = 0
        elif "close" in p:
            s = 1
        elif "broad" in p:
            s = 2
        elif "narrow" in p:
            s = 3
        else:
            s = 9
        return (s, -(m.confidence or 0.0))
    return sorted(mappings, key=score)


# ----------------------------
# 2) Raw data readers
# ----------------------------

def read_genes_to_phenotype(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    required = {"hpo_id", "gene_symbol"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"genes_to_phenotype file missing columns: {missing}. Got: {list(df.columns)}")
    return df


def read_mgi_phenogenomp(path: str | Path) -> pd.DataFrame:
    cols = ["genotype", "allele", "background", "mp_id", "pmid", "marker_mgi", "genotype_mgi"]
    df = pd.read_csv(path, sep="\t", header=None, names=cols, dtype=str)
    return df


# ----------------------------
# 3) Edge builders
# ----------------------------

def build_hp_to_mp_edges(
    g2p: pd.DataFrame,
    hp2mp: dict[str, list[Mapping]],
    best_only: bool = True,
) -> pd.DataFrame:
    rows = []
    for _, r in g2p.iterrows():
        hp_id = r["hpo_id"]
        gene = r.get("gene_symbol") or r.get("ncbi_gene_id")
        disease_id = r.get("disease_id")

        mappings = hp2mp.get(hp_id, [])
        if not mappings:
            continue

        ranked = rank_mappings(mappings)
        selected = ranked[:1] if best_only else ranked

        for m in selected:
            rows.append({
                "gene_symbol": gene,
                "hp_id": hp_id,
                "mp_id": m.mp_id,
                "predicate_id": m.predicate_id,
                "confidence": m.confidence,
                "source": m.source,
                "disease_id": disease_id,
            })
    return pd.DataFrame(rows)


def build_top_edges(
    terms: pd.Series,
    ont,
    top_set: set[str],
    root_id: str,
    term_col: str,
    top_col: str,
) -> pd.DataFrame:
    edges = []
    for tid in sorted(set(t for t in terms.dropna().unique() if isinstance(t, str))):
        for top in collapse_to_top(ont, tid, top_set, root_id):
            edges.append({term_col: tid, top_col: top})
    return pd.DataFrame(edges)


def build_human_gene_mp_top_edges(
    g2p: pd.DataFrame,
    hp2mp: dict[str, list[Mapping]],
    mp_in_top: pd.DataFrame,
    best_only: bool = True,
) -> pd.DataFrame:
    hp_mp_rows = []
    for _, r in g2p[["gene_symbol", "hpo_id"]].drop_duplicates().iterrows():
        mappings = hp2mp.get(r["hpo_id"], [])
        if not mappings:
            continue
        ranked = rank_mappings(mappings)
        selected = ranked[:1] if best_only else ranked
        for m in selected:
            hp_mp_rows.append({"gene_symbol": r["gene_symbol"], "mp_id": m.mp_id})

    if not hp_mp_rows:
        return pd.DataFrame(columns=["gene_symbol", "mp_top_id"])

    gene_mp = pd.DataFrame(hp_mp_rows).drop_duplicates()
    mp_to_top = mp_in_top[["mp_id", "mp_top_id"]].dropna().drop_duplicates()
    merged = (
        gene_mp
        .merge(mp_to_top, on="mp_id", how="inner")
        [["gene_symbol", "mp_top_id"]]
        .drop_duplicates()
    )
    return merged


def build_mouse_gene_mp_top_edges(
    mouse_gene_mp: pd.DataFrame,
    mp_in_top: pd.DataFrame,
) -> pd.DataFrame:
    mp_to_top = mp_in_top[["mp_id", "mp_top_id"]].dropna().drop_duplicates()
    merged = (
        mouse_gene_mp[["gene_symbol", "mp_id"]]
        .dropna()
        .drop_duplicates()
        .merge(mp_to_top, on="mp_id", how="inner")
        [["gene_symbol", "mp_top_id"]]
        .drop_duplicates()
    )
    return merged


# ----------------------------
# 4) End-to-end pipeline
# ----------------------------

def run_pipeline(
    genes_to_phenotype_path: str | Path,
    mgi_path: str | Path,
    sssom_path: str | Path,
    data_dir: str | Path,
    gene_info_path: str | Path,
    out_dir: str | Path = "./out",
    best_only: bool = True,
    species: str = "both",
) -> dict[str, pd.DataFrame]:
    """Run the full HPO/MGD → system-level phenotype edge pipeline.

    Outputs written to *out_dir*:
      edge_human_gene_has_hp.tsv        - HumanGene → HPO term
      edge_hp_mapped_to_mp.tsv          - HPO term  → MP term (via SSSOM)
      edge_hp_in_top.tsv                - HPO term  → HP system-level term
      edge_mp_in_top_from_hp.tsv        - MP term   → MP system-level term (from HPO)
      edge_mouse_gene_has_mp.tsv        - MouseGene (gene_symbol) → MP term
      edge_mp_in_top_from_mgi.tsv       - MP term   → MP system-level term (from MGD)
      edge_human_gene_has_mp_top.tsv    - HumanGene → MP top term (direct)
      edge_mouse_gene_has_mp_top.tsv    - MouseGene (gene_symbol) → MP top term (direct)
      node_mp_top_names.tsv             - MP top term names (mp_id, name)
      qc_unmapped_hp.tsv                - HP terms with no MP mapping (QC)
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    build_human = species in ("human", "both")
    build_mouse  = species in ("mouse", "both")

    hp_ont, mp_ont = load_ontologies(Path(data_dir))
    hp_top = top_children(hp_ont, HP_ROOT)
    mp_top = top_children(mp_ont, MP_ROOT)

    # Write MP top term names — species-agnostic, written once
    mp_top_names = build_top_name_table(mp_ont, mp_top, id_col="mp_id")
    mp_top_names.to_csv(out_dir / "node_mp_top_names.tsv", sep="\t", index=False)

    # ── Human branch ──────────────────────────────────────────────────
    if build_human:
        g2p   = read_genes_to_phenotype(genes_to_phenotype_path)
        hp2mp = load_hp_to_mp_sssom(sssom_path)

        human_gene_hp = g2p[["gene_symbol", "hpo_id", "disease_id"]].drop_duplicates()
        human_gene_hp.to_csv(out_dir / "edge_human_gene_has_hp.tsv", sep="\t", index=False)

        hp_mp = build_hp_to_mp_edges(g2p, hp2mp, best_only=best_only).drop_duplicates()
        hp_mp.to_csv(out_dir / "edge_hp_mapped_to_mp.tsv", sep="\t", index=False)

        hp_in_top = build_top_edges(
            g2p["hpo_id"], hp_ont, hp_top, HP_ROOT, "hp_id", "hp_top_id"
        )
        hp_in_top.to_csv(out_dir / "edge_hp_in_top.tsv", sep="\t", index=False)

        mp_in_top_from_hp = build_top_edges(
            hp_mp["mp_id"], mp_ont, mp_top, MP_ROOT, "mp_id", "mp_top_id"
        )
        mp_in_top_from_hp.to_csv(out_dir / "edge_mp_in_top_from_hp.tsv", sep="\t", index=False)

        human_gene_mp_top = build_human_gene_mp_top_edges(
            g2p, hp2mp, mp_in_top_from_hp, best_only=best_only
        )
        human_gene_mp_top.to_csv(out_dir / "edge_human_gene_has_mp_top.tsv", sep="\t", index=False)

        mapped_hp = set(hp_mp["hp_id"].dropna()) if not hp_mp.empty else set()
        unmapped = sorted(set(g2p["hpo_id"].dropna()) - mapped_hp)
        pd.DataFrame({"unmapped_hp_id": unmapped}).to_csv(
            out_dir / "qc_unmapped_hp.tsv", sep="\t", index=False
        )
    else:
        human_gene_hp = hp_mp = hp_in_top = mp_in_top_from_hp = pd.DataFrame()
        human_gene_mp_top = pd.DataFrame()
        unmapped = []

    # ── Mouse branch ──────────────────────────────────────────────────
    if build_mouse:
        from src.utils.gene_info import build_mouse_symbol_maps
        _, mgi_to_symbol = build_mouse_symbol_maps(gene_info_path)

        mgi = read_mgi_phenogenomp(mgi_path)

        raw_mouse = (
            mgi[["marker_mgi", "mp_id"]]
            .dropna(subset=["marker_mgi", "mp_id"])
            .drop_duplicates()
            .copy()
        )
        raw_mouse["mgi_single"] = raw_mouse["marker_mgi"].str.split("|")
        raw_mouse = raw_mouse.explode("mgi_single")
        raw_mouse["mgi_single"] = raw_mouse["mgi_single"].str.strip()
        raw_mouse["gene_symbol"] = raw_mouse["mgi_single"].map(mgi_to_symbol)

        n_total  = raw_mouse["mgi_single"].nunique()
        n_mapped = raw_mouse.dropna(subset=["gene_symbol"])["mgi_single"].nunique()
        n_miss   = n_total - n_mapped
        if n_miss:
            print(f"  QC: {n_miss}/{n_total} unique MGI IDs could not be mapped "
                  f"to a gene symbol (alleles, QTLs, transgenes, etc.)")

        mouse_gene_mp = (
            raw_mouse.dropna(subset=["gene_symbol"])
            [["gene_symbol", "mp_id"]]
            .drop_duplicates()
        )
        mouse_gene_mp.to_csv(out_dir / "edge_mouse_gene_has_mp.tsv", sep="\t", index=False)

        mp_in_top_from_mgi = build_top_edges(
            mgi["mp_id"], mp_ont, mp_top, MP_ROOT, "mp_id", "mp_top_id"
        )
        mp_in_top_from_mgi.to_csv(out_dir / "edge_mp_in_top_from_mgi.tsv", sep="\t", index=False)

        mouse_gene_mp_top = build_mouse_gene_mp_top_edges(mouse_gene_mp, mp_in_top_from_mgi)
        mouse_gene_mp_top.to_csv(out_dir / "edge_mouse_gene_has_mp_top.tsv", sep="\t", index=False)
    else:
        mouse_gene_mp = mp_in_top_from_mgi = mouse_gene_mp_top = pd.DataFrame()

    # ── Summary ───────────────────────────────────────────────────────
    print(f"HP system-level terms : {len(hp_top)}")
    print(f"MP system-level terms : {len(mp_top)}")
    if build_human:
        print(f"Human gene→HP edges   : {len(human_gene_hp)}")
        print(f"HP→MP mappings        : {len(hp_mp)}")
        print(f"Unmapped HP terms     : {len(unmapped)}")
        print(f"Human gene→MP top     : {len(human_gene_mp_top)}")
    if build_mouse:
        print(f"Mouse gene→MP edges   : {len(mouse_gene_mp)}")
        print(f"Mouse gene→MP top     : {len(mouse_gene_mp_top)}")

    return {
        "human_gene_hp": human_gene_hp,
        "hp_mp": hp_mp,
        "hp_in_top": hp_in_top,
        "mp_in_top_from_hp": mp_in_top_from_hp,
        "mouse_gene_mp": mouse_gene_mp,
        "mp_in_top_from_mgi": mp_in_top_from_mgi,
        "human_gene_mp_top": human_gene_mp_top,
        "mouse_gene_mp_top": mouse_gene_mp_top,
    }


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("--genes-to-phenotype", required=True)
    p.add_argument("--mgi", required=True)
    p.add_argument("--sssom", required=True)
    p.add_argument("--data-dir", required=True)
    p.add_argument("--out", default="./out")
    p.add_argument("--gene-info", required=True)
    p.add_argument("--all-mappings", action="store_true")
    args = p.parse_args()

    run_pipeline(
        genes_to_phenotype_path=args.genes_to_phenotype,
        mgi_path=args.mgi,
        sssom_path=args.sssom,
        data_dir=args.data_dir,
        gene_info_path=args.gene_info,
        out_dir=args.out,
        best_only=not args.all_mappings,
    )