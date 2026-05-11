"""uberon_export.py - Uberon anatomical ontology pipeline for the gene-phenotype KG.

What this pipeline produces
---------------------------
node_uberon_terms.tsv
    One row per Uberon anatomical term that is referenced by at least one
    GTEx tissue mapping.  Columns: uberon_id, name, definition, synonyms.

edge_tissue_mapped_to_uberon.tsv
    GTEx tissue → Uberon term edges.
    Columns: gtex_tissue_id, uberon_id, gtex_tissue_name.

Data sources
------------
- Uberon OBO file    : http://purl.obolibrary.org/obo/uberon/basic.obo
  (basic.obo includes only Uberon-native terms, no imports — faster to parse)
- GTEx tissue names  : derived from the existing expression layer TSV, OR from
  the curated GTEX_TO_UBERON mapping defined in this file.

Usage
-----
python uberon_export.py \
    --uberon-obo  /data/raw/uberon_basic.obo \
    --expression-tsv  /data/processed/edge_human_gene_expressed_in_tissue.tsv \
    --out  /data/processed

The --expression-tsv argument is optional.  Without it the pipeline uses the
curated GTEX_TO_UBERON mapping directly; with it, only tissues that actually
appear in the expression layer are written to the output.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd

from src.layers.gene_expression_export import _tissue_id

logger = logging.getLogger(__name__)

# ── Curated GTEx-tissue → Uberon mapping ─────────────────────────────────────
#
# Keys  : GTEx tissue IDs as stored in the expression layer (GTEX:<name> format,
#         spaces replaced by underscores to match the graph node IDs).
# Values: canonical Uberon IDs (UBERON:XXXXXXX).
#
# Sources used for curation:
#   - GTEx v8 tissue list: https://gtexportal.org/home/tissueSummaryPage
#   - Uberon term browser : https://www.ebi.ac.uk/ols/ontologies/uberon
#   - BRENDA tissue ontology crossrefs embedded in Uberon OBO
#
GTEX_TO_UBERON: dict[str, str] = {
    # ── Brain regions ────────────────────────────────────────────────── #
    "GTEX:Brain_Amygdala":                          "UBERON:0001876",
    "GTEX:Brain_Anterior_cingulate_cortex_BA24":    "UBERON:0006101",
    "GTEX:Brain_Caudate_basal_ganglia":             "UBERON:0001873",
    "GTEX:Brain_Cerebellar_Hemisphere":             "UBERON:0002037",
    "GTEX:Brain_Cerebellum":                        "UBERON:0002037",
    "GTEX:Brain_Cortex":                            "UBERON:0001851",
    "GTEX:Brain_Frontal_Cortex_BA9":                "UBERON:0001870",
    "GTEX:Brain_Hippocampus":                       "UBERON:0002310",
    "GTEX:Brain_Hypothalamus":                      "UBERON:0001898",
    "GTEX:Brain_Nucleus_accumbens_basal_ganglia":   "UBERON:0001882",
    "GTEX:Brain_Putamen_basal_ganglia":             "UBERON:0001874",
    "GTEX:Brain_Spinal_cord_cervical_c_1":          "UBERON:0002726",
    "GTEX:Brain_Substantia_nigra":                  "UBERON:0002038",
    # ── Heart ────────────────────────────────────────────────────────── #
    "GTEX:Heart_Atrial_Appendage":                  "UBERON:0006618",
    "GTEX:Heart_Left_Ventricle":                    "UBERON:0002084",
    # ── Artery ───────────────────────────────────────────────────────── #
    "GTEX:Artery_Aorta":                            "UBERON:0000947",
    "GTEX:Artery_Coronary":                         "UBERON:0001621",
    "GTEX:Artery_Tibial":                           "UBERON:0007610",
    # ── Digestive system ─────────────────────────────────────────────── #
    "GTEX:Colon_Sigmoid":                           "UBERON:0001159",
    "GTEX:Colon_Transverse":                        "UBERON:0001157",
    "GTEX:Small_Intestine_Terminal_Ileum":          "UBERON:0002116",
    "GTEX:Stomach":                                 "UBERON:0000945",
    "GTEX:Liver":                                   "UBERON:0002107",
    "GTEX:Pancreas":                                "UBERON:0001264",
    "GTEX:Esophagus_Gastroesophageal_Junction":     "UBERON:0007650",
    "GTEX:Esophagus_Mucosa":                        "UBERON:0002469",
    "GTEX:Esophagus_Muscularis":                    "UBERON:0004648",
    # ── Reproductive system ──────────────────────────────────────────── #
    "GTEX:Ovary":                                   "UBERON:0000992",
    "GTEX:Uterus":                                  "UBERON:0000995",
    "GTEX:Vagina":                                  "UBERON:0000996",
    "GTEX:Prostate":                                "UBERON:0002367",
    "GTEX:Testis":                                  "UBERON:0000473",
    # ── Endocrine ────────────────────────────────────────────────────── #
    "GTEX:Adrenal_Gland":                           "UBERON:0002369",
    "GTEX:Thyroid":                                 "UBERON:0002046",
    "GTEX:Pituitary":                               "UBERON:0000007",
    # ── Respiratory ──────────────────────────────────────────────────── #
    "GTEX:Lung":                                    "UBERON:0002048",
    # ── Musculoskeletal ──────────────────────────────────────────────── #
    "GTEX:Muscle_Skeletal":                         "UBERON:0001134",
    "GTEX:Minor_Salivary_Gland":                    "UBERON:0006330",
    # ── Skin ─────────────────────────────────────────────────────────── #
    "GTEX:Skin_Not_Sun_Exposed_Suprapubic":         "UBERON:0000014",
    "GTEX:Skin_Sun_Exposed_Lower_leg":              "UBERON:0000014",
    # ── Nerve ────────────────────────────────────────────────────────── #
    "GTEX:Nerve_Tibial":                            "UBERON:0001323",
    # ── Blood / Immune ───────────────────────────────────────────────── #
    "GTEX:Whole_Blood":                             "UBERON:0000178",
    "GTEX:Spleen":                                  "UBERON:0002106",
    # ── Urinary ──────────────────────────────────────────────────────── #
    "GTEX:Kidney_Cortex":                           "UBERON:0001225",
    "GTEX:Kidney_Medulla":                          "UBERON:0001293",
    "GTEX:Bladder":                                 "UBERON:0001255",
    # ── Adipose / connective tissue ──────────────────────────────────── #
    "GTEX:Adipose_Subcutaneous":                    "UBERON:0002190",
    "GTEX:Adipose_Visceral_Omentum":                "UBERON:0003688",
    # ── Other ────────────────────────────────────────────────────────── #
    "GTEX:Breast_Mammary_Tissue":                   "UBERON:0001911",
    "GTEX:Cells_Cultured_fibroblasts":              "UBERON:0002082",
    "GTEX:Cells_EBV_transformed_lymphocytes":       "UBERON:0000029",
    "GTEX:Cervix_Ectocervix":                       "UBERON:0012249",
    "GTEX:Cervix_Endocervix":                       "UBERON:0000458",
    "GTEX:Fallopian_Tube":                          "UBERON:0003889",
}


# ── Uberon OBO reader ─────────────────────────────────────────────────────────

def _parse_uberon_obo(obo_path: str | Path) -> pd.DataFrame:
    """Parse an Uberon OBO file with the ``obonet`` library.

    Returns a DataFrame with columns:
        uberon_id | name | definition | synonyms
    Only UBERON:* terms are retained (no cross-ontology terms).

    Requires:
        pip install obonet
    """
    try:
        import obonet  # type: ignore
    except ImportError as exc:
        raise ImportError(
            "obonet is required to parse the Uberon OBO file.\n"
            "Install it with:  pip install obonet"
        ) from exc

    logger.info("Parsing Uberon OBO: %s", obo_path)
    graph = obonet.read_obo(str(obo_path))
    logger.info("OBO graph loaded — %d total terms", graph.number_of_nodes())

    rows: list[dict] = []
    for node_id, data in graph.nodes(data=True):
        if not node_id.startswith("UBERON:"):
            continue
        name = data.get("name", "")
        # definition may carry a trailing citation in brackets: strip it
        raw_def: str = data.get("def", "")
        definition = raw_def.split('" [')[0].lstrip('"') if raw_def else ""
        # synonyms: list of strings like '"exact synonym" EXACT []'
        raw_syns: list[str] = data.get("synonym", [])
        synonym_names = [s.split('"')[1] for s in raw_syns if '"' in s]
        rows.append({
            "uberon_id":  node_id,
            "name":       name,
            "definition": definition,
            "synonyms":   "|".join(synonym_names),
        })

    df = pd.DataFrame(rows, columns=["uberon_id", "name", "definition", "synonyms"])
    logger.info("Uberon terms parsed: %d", len(df))
    return df


def _parse_uberon_obo_fallback(uberon_ids: set[str]) -> pd.DataFrame:
    """Return a minimal stub DataFrame for the provided Uberon IDs.

    Used when no OBO file is supplied — only the IDs from GTEX_TO_UBERON are
    kept, with empty name/definition/synonyms fields.  Operators should
    re-run with --uberon-obo to populate the full metadata.
    """
    logger.warning(
        "No Uberon OBO file provided — writing stub nodes with empty metadata. "
        "Re-run with --uberon-obo <path> to populate names and definitions."
    )
    rows = [{"uberon_id": uid, "name": "", "definition": "", "synonyms": ""}
            for uid in sorted(uberon_ids)]
    return pd.DataFrame(rows, columns=["uberon_id", "name", "definition", "synonyms"])


# ── Mapping builders ─────────────────────────────────────────────────────────

def build_tissue_to_uberon_edges(
    mapping: dict[str, str],
    tissue_ids_in_graph: set[str] | None = None,
) -> pd.DataFrame:
    """Build the tissue → Uberon edge table from the curated mapping.

    Args:
        mapping:              GTEX_TO_UBERON or any {gtex_id: uberon_id} dict.
        tissue_ids_in_graph:  If provided, restrict output to tissues that
                              actually appear in the expression layer.

    Returns a DataFrame: gtex_tissue_id | uberon_id | gtex_tissue_name
    """
    rows: list[dict] = []
    for gtex_id, uberon_id in mapping.items():
        if tissue_ids_in_graph is not None and gtex_id not in tissue_ids_in_graph:
            logger.debug("Tissue not in expression layer, skipping: %s", gtex_id)
            continue
        tissue_name = gtex_id.replace("GTEX:", "").replace("_", " ")
        rows.append({
            "gtex_tissue_id":   gtex_id,
            "uberon_id":        uberon_id,
            "gtex_tissue_name": tissue_name,
        })

    df = pd.DataFrame(rows, columns=["gtex_tissue_id", "uberon_id", "gtex_tissue_name"])
    logger.info("Tissue→Uberon edges built: %d", len(df))
    return df


def _read_expression_tissue_ids(expression_tsv: str | Path) -> set[str]:
    """Return the set of tissue node IDs present in the expression edge TSV.

    Expects the TSV to have a ``tissue_id`` column (e.g. GTEX:Brain_Amygdala).
    Falls back to the ``tissue_name`` column if ``tissue_id`` is absent.
    """
    path = Path(expression_tsv)
    if not path.exists():
        logger.warning("Expression TSV not found: %s — using full curated mapping", path)
        return set()

    df = pd.read_csv(path, sep="\t", dtype=str)
    if "tissue_id" in df.columns:
        ids = set(df["tissue_id"].dropna().unique())
    elif "tissue_name" in df.columns:
        ids = {_tissue_id(n) for n in df["tissue_name"].dropna().unique()}
    else:
        logger.warning(
            "Expression TSV has neither 'tissue_id' nor 'tissue_name' column — "
            "using full curated mapping."
        )
        return set()

    logger.debug("Tissue IDs found in expression layer: %d", len(ids))
    return ids


# ── End-to-end pipeline ───────────────────────────────────────────────────────

def run_uberon_pipeline(
    out_dir: str | Path,
    uberon_obo_path: str | Path | None = None,
    expression_tsv: str | Path | None = None,
    mapping: dict[str, str] | None = None,
) -> dict[str, pd.DataFrame]:
    """Run the Uberon anatomical ontology pipeline.

    Args:
        out_dir:          Directory where TSV outputs are written.
        uberon_obo_path:  Path to ``uberon/basic.obo``.  When *None*, a stub
                          table is written with empty metadata fields.
        expression_tsv:   Path to the expression edge TSV (used to filter
                          output to tissues actually present in the graph).
                          When *None*, all curated mappings are written.
        mapping:          Custom {gtex_id: uberon_id} dict.  Defaults to the
                          built-in GTEX_TO_UBERON table.

    Returns:
        {"uberon_terms": <DataFrame>, "tissue_uberon_edges": <DataFrame>}
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    mapping = mapping or GTEX_TO_UBERON

    logger.info("=" * 60)
    logger.info("Starting Uberon pipeline")
    logger.info("Output dir    : %s", out_dir)
    logger.info("OBO file      : %s", uberon_obo_path or "<not provided — stub mode>")
    logger.info("Expression TSV: %s", expression_tsv or "<not provided>")
    logger.info("Mapping size  : %d GTEx tissues", len(mapping))
    logger.info("=" * 60)

    # ── 1. Build tissue → Uberon edge table ──────────────────────────── #
    tissue_ids_in_graph: set[str] | None = None
    if expression_tsv:
        tissue_ids_in_graph = _read_expression_tissue_ids(expression_tsv)
        if not tissue_ids_in_graph:
            logger.warning("Expression TSV yielded 0 tissue IDs — using full mapping.")
            tissue_ids_in_graph = None

    tissue_uberon = build_tissue_to_uberon_edges(mapping, tissue_ids_in_graph)
    tissue_uberon_path = out_dir / "edge_tissue_mapped_to_uberon.tsv"
    tissue_uberon.to_csv(tissue_uberon_path, sep="\t", index=False)
    logger.info("Written: %s  (%d rows)", tissue_uberon_path.name, len(tissue_uberon))

    # ── 2. Build Uberon term node table ──────────────────────────────── #
    referenced_uberon_ids = set(tissue_uberon["uberon_id"].dropna().unique())
    logger.debug("Unique Uberon IDs referenced: %d", len(referenced_uberon_ids))

    if uberon_obo_path and Path(uberon_obo_path).exists():
        all_uberon = _parse_uberon_obo(uberon_obo_path)
        uberon_terms = all_uberon[all_uberon["uberon_id"].isin(referenced_uberon_ids)].copy()
        missing = referenced_uberon_ids - set(uberon_terms["uberon_id"])
        if missing:
            logger.warning(
                "%d referenced Uberon IDs not found in the OBO file: %s",
                len(missing), sorted(missing)[:10],
            )
        logger.info(
            "Uberon term nodes: %d / %d referenced IDs found in OBO",
            len(uberon_terms), len(referenced_uberon_ids),
        )
    else:
        uberon_terms = _parse_uberon_obo_fallback(referenced_uberon_ids)

    uberon_node_path = out_dir / "node_uberon_terms.tsv"
    uberon_terms.to_csv(uberon_node_path, sep="\t", index=False)
    logger.info("Written: %s  (%d rows)", uberon_node_path.name, len(uberon_terms))

    # ── 3. QC report ─────────────────────────────────────────────────── #
    all_mapped = set(tissue_uberon["gtex_tissue_id"])
    all_curated = set(mapping.keys())
    unmapped = sorted(all_curated - all_mapped)
    if unmapped:
        qc_path = out_dir / "qc_uberon_unmapped_tissues.tsv"
        pd.DataFrame({"gtex_tissue_id": unmapped}).to_csv(qc_path, sep="\t", index=False)
        logger.info(
            "QC: %d curated GTEx tissues were not written (filtered out or "
            "missing in expression layer) → %s",
            len(unmapped), qc_path.name,
        )

    logger.info("Uberon pipeline complete.")
    return {"uberon_terms": uberon_terms, "tissue_uberon_edges": tissue_uberon}


# ── CLI entry-point ───────────────────────────────────────────────────────────

if __name__ == "__main__":
    from src.utils.logging_config import setup_logging

    p = argparse.ArgumentParser(
        description="Export Uberon anatomical ontology edges for the gene-phenotype KG."
    )
    p.add_argument(
        "--uberon-obo",
        default=None,
        help=(
            "Path to the Uberon basic OBO file "
            "(uberon/basic.obo, download from http://purl.obolibrary.org/obo/uberon/basic.obo). "
            "When omitted, stub nodes with empty metadata are written."
        ),
    )
    p.add_argument(
        "--expression-tsv",
        default=None,
        help=(
            "Path to edge_human_gene_expressed_in_tissue.tsv (or equivalent). "
            "Used to restrict the output to tissues present in the expression layer. "
            "When omitted, all curated GTEx→Uberon mappings are written."
        ),
    )
    p.add_argument(
        "--out",
        default="./out",
        help="Output directory for TSV files (default: ./out).",
    )
    args = p.parse_args()

    setup_logging(out_dir=args.out)

    run_uberon_pipeline(
        out_dir=args.out,
        uberon_obo_path=args.uberon_obo,
        expression_tsv=args.expression_tsv,
    )
