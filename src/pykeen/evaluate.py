"""src/pykeen/evaluate.py — AUPRC evaluation for the gene↔MP link prediction task.

Usage
-----
    uv run python -m src.pykeen.evaluate --model-dir pykeen_out/rotate_fast/

The model directory must contain:
  trained_model.pkl       (saved by run_training)
  training_triples/       (saved by run_training via result.save_to_directory)
  test_triples/           (saved by run_training via test_tf.to_path_binary)
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import pandas as pd
from sklearn.metrics import average_precision_score

from .config import GENE_PHENOTYPE_RELATION
from .predict import GenePhenotypePredictor

logger = logging.getLogger(__name__)


def evaluate_auprc(
    predictor: GenePhenotypePredictor,
    test_tf,
    relation: str = GENE_PHENOTYPE_RELATION,
    override_positives: "dict[str, set[str]] | None" = None,
) -> dict:
    """Compute per-class AUPRC on the test set.

    For each MP top term that appears in the test set, all genes are ranked by
    the model score (filter_known=False — training positives stay in the candidate
    pool).  y_true=1 for genes that are test-set positives for that MP term.

    Parameters
    ----------
    override_positives
        If provided, use this fixed mapping {mp_term_id → set of gene_ids} as the
        ground-truth test positives instead of extracting them from *test_tf*.
        Pass this in ablation studies so that every condition is evaluated on
        exactly the same gene→phenotype test pairs (derived from the full-graph
        split), making AUPRC scores directly comparable across conditions.

    Returns dict with keys:
      mean_auprc   float
      num_classes  int   (MP terms with ≥1 test positive that could be scored)
      per_class    DataFrame  [mp_term_id, ap, n_pos_test, n_scored]
    """
    # ── 1. Determine positive pairs ───────────────────────────────────────────
    if override_positives is not None:
        # Fixed test set supplied externally (ablation study use-case).
        positives = override_positives
        total_pairs = sum(len(v) for v in positives.values())
        logger.info(
            "Test set (fixed across conditions): %d gene→MP pairs across %d MP terms",
            total_pairs,
            len(positives),
        )
    else:
        # Default: extract positive pairs from the condition-specific test split.
        rel_id = test_tf.relation_to_id.get(relation)
        if rel_id is None:
            raise ValueError(
                f"Relation '{relation}' not found in test TriplesFactory.\n"
                f"Available: {sorted(test_tf.relation_to_id.keys())}"
            )

        mask = test_tf.mapped_triples[:, 1] == rel_id
        pos_triples = test_tf.mapped_triples[mask]

        id_to_entity = {v: k for k, v in test_tf.entity_to_id.items()}

        # mp_term_id → set of gene_ids that are test positives
        positives: dict[str, set[str]] = {}
        for head_id, _, tail_id in pos_triples.tolist():
            gene    = id_to_entity[head_id]
            mp_term = id_to_entity[tail_id]
            positives.setdefault(mp_term, set()).add(gene)

        logger.info(
            "Test set: %d gene→MP pairs across %d MP terms",
            len(pos_triples),
            len(positives),
        )

    # ── 2. Score all genes for each MP term ───────────────────────────────────
    n_genes = len(predictor.gene_entities)
    rows = []
    skipped = 0

    for mp_term_id, pos_genes in sorted(positives.items()):
        if mp_term_id not in predictor.training_tf.entity_to_id:
            skipped += 1
            continue

        # Retrieve scores for ALL genes (no cap, include known associations).
        pred_df = predictor.predict_for_phenotype(
            mp_term_id,
            top_k=n_genes,
            filter_known=False,
            only_genes=True,
        )

        # Genes in the fixed test set that are absent from this condition's
        # entity space get score=-inf so they rank last.  Silently forgiving
        # them would inflate AUPRC for minimal conditions (e.g. pheno_only).
        if override_positives is not None:
            scored_genes = set(pred_df["gene_id"])
            missing_pos  = pos_genes - scored_genes
            if missing_pos:
                missing_df = pd.DataFrame({
                    "gene_id": list(missing_pos),
                    "score":   float("-inf"),
                })
                pred_df = pd.concat([pred_df, missing_df], ignore_index=True)

        y_true  = [1 if g in pos_genes else 0 for g in pred_df["gene_id"]]
        y_score = pred_df["score"].tolist()

        if sum(y_true) == 0:
            # Only reachable in the non-override path (entity mismatch in
            # condition-specific test set — rare).
            skipped += 1
            continue

        ap = float(average_precision_score(y_true, y_score))
        rows.append({
            "mp_term_id":  mp_term_id,
            "ap":          ap,
            "n_pos_test":  len(pos_genes),
            "n_scored":    len(pred_df),
        })

    if skipped:
        logger.warning("Skipped %d MP terms (not in training graph or no scored positives)", skipped)

    per_class = (
        pd.DataFrame(rows, columns=["mp_term_id", "ap", "n_pos_test", "n_scored"])
        .sort_values("ap", ascending=False)
        .reset_index(drop=True)
    )
    mean_auprc = float(per_class["ap"].mean()) if len(per_class) > 0 else float("nan")

    logger.info(
        "AUPRC — mean: %.4f  |  classes evaluated: %d",
        mean_auprc,
        len(per_class),
    )

    return {
        "mean_auprc":  mean_auprc,
        "num_classes": len(per_class),
        "per_class":   per_class,
    }


# ── CLI ───────────────────────────────────────────────────────────────────────

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Compute per-class AUPRC for a trained RotatE model.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--model-dir", required=True,
        help="Directory produced by train.py (must contain test_triples/).",
    )
    p.add_argument(
        "--out-csv", default=None,
        help="Path for per-class AUPRC CSV (default: <model-dir>/auprc_per_phenotype.csv).",
    )
    return p.parse_args()


if __name__ == "__main__":
    from pykeen.triples import TriplesFactory

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-8s  %(name)s  %(message)s",
        datefmt="%H:%M:%S",
    )

    args = _parse_args()
    model_dir = Path(args.model_dir)

    test_triples_path = model_dir / "test_triples"
    if not test_triples_path.exists():
        raise FileNotFoundError(
            f"test_triples/ not found in {model_dir}.\n"
            "Re-run training with the current train.py to save test_triples/."
        )

    predictor = GenePhenotypePredictor.from_directory(model_dir)
    test_tf   = TriplesFactory.from_path_binary(test_triples_path)

    results = evaluate_auprc(predictor, test_tf)

    csv_path = Path(args.out_csv) if args.out_csv else model_dir / "auprc_per_phenotype.csv"
    results["per_class"].to_csv(csv_path, index=False)
    logger.info("Per-class AUPRC saved to: %s", csv_path)

    # Patch summary.json if it exists
    summary_path = model_dir / "summary.json"
    if summary_path.exists():
        summary = json.loads(summary_path.read_text())
        summary["metrics"]["mean_auprc"]   = results["mean_auprc"]
        summary["metrics"]["auprc_classes"] = results["num_classes"]
        summary_path.write_text(json.dumps(summary, indent=2, default=str))
        logger.info("summary.json updated with mean_auprc=%.4f", results["mean_auprc"])

    print(f"\nmean AUPRC : {results['mean_auprc']:.4f}")
    print(f"classes    : {results['num_classes']}")
    print(f"\nTop-10 classes by AP:")
    print(results["per_class"].head(10).to_string(index=False))