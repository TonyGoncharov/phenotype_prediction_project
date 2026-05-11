"""src/pykeen/predict.py — gene ↔ MP top term link prediction using a trained RotatE model.

    predictor = GenePhenotypePredictor.from_directory("pykeen_out/rotate/")
    df = predictor.predict_for_gene("HGNC:TP53", top_k=20)
    df = predictor.predict_for_phenotype("MP:0001265", top_k=30)
    df = predictor.predict_batch_for_phenotypes(["MP:0005387", "MP:0005384"], top_k=50)
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd
import torch
from pykeen.predict import predict_target
from pykeen.triples import TriplesFactory

from .config import GENE_PHENOTYPE_RELATION

logger = logging.getLogger(__name__)


class GenePhenotypePredictor:
    """Wrapper around a trained PyKEEN model for gene ↔ phenotype prediction."""

    # ── Construction ─────────────────────────────────────────────────────────

    def __init__(
        self,
        model: torch.nn.Module,
        training_tf: TriplesFactory,
        relation: str = GENE_PHENOTYPE_RELATION,
    ) -> None:
        if relation not in training_tf.relation_to_id:
            available = sorted(training_tf.relation_to_id.keys())
            raise ValueError(
                f"Relation '{relation}' not found in the TriplesFactory.\n"
                f"Available relations: {available}"
            )

        self.model = model.eval()
        self.training_tf = training_tf
        self.relation = relation

        # Convenience: partition entities by type based on ID prefix
        all_entities = list(training_tf.entity_to_id.keys())
        self.mp_entities:   list[str] = [e for e in all_entities if e.startswith("MP:")]
        self.gene_entities: list[str] = [e for e in all_entities if e.startswith("HGNC:")]

        logger.info(
            "Predictor ready — %d genes | %d MP top-terms | relation: %s",
            len(self.gene_entities),
            len(self.mp_entities),
            relation,
        )

    @classmethod
    def from_directory(
        cls,
        model_dir: str | Path,
        relation: str = GENE_PHENOTYPE_RELATION,
    ) -> "GenePhenotypePredictor":
        """Load a predictor from a directory produced by `run_training()`.

        Loads the trained model (trained_model.pkl) and the training TriplesFactory
        (training_triples/) independently — PipelineResult.from_directory() was
        removed in PyKEEN 1.11.
        """
        import torch
        from pykeen.triples import TriplesFactory

        model_dir = Path(model_dir)
        if not model_dir.exists():
            raise FileNotFoundError(f"Model directory not found: {model_dir}")

        logger.info("Loading model from %s …", model_dir)
        model = torch.load(model_dir / "trained_model.pkl", weights_only=False)
        training_tf = TriplesFactory.from_path_binary(model_dir / "training_triples")

        return cls(model=model, training_tf=training_tf, relation=relation)

    # ── Prediction — tail (gene → phenotype) ─────────────────────────────────

    @torch.no_grad()
    def predict_for_gene(
        self,
        gene_id: str,
        top_k: int = 20,
        filter_known: bool = True,
        only_mp_terms: bool = True,
    ) -> pd.DataFrame:
        """Predict the top-K most likely MP phenotypes for a single gene.

        filter_known=True excludes training-set pairs (novel predictions only).
        only_mp_terms=True restricts candidates to MP: entities.
        Returns DataFrame: gene_id | mp_term_id | score | rank.
        """
        self._check_entity(gene_id)

        pred = predict_target(
            model=self.model,
            head=gene_id,
            relation=self.relation,
            triples_factory=self.training_tf,
        )
        if filter_known:
            pred = pred.filter_triples(self.training_tf)

        pred_df = pred.df
        if only_mp_terms:
            pred_df = pred_df[pred_df["tail_label"].str.startswith("MP:")].copy()

        pred_df = pred_df.sort_values("score", ascending=False).head(top_k)
        pred_df.insert(0, "gene_id", gene_id)
        pred_df = pred_df.rename(columns={"tail_label": "mp_term_id"})
        pred_df["rank"] = range(1, len(pred_df) + 1)

        return pred_df[["gene_id", "mp_term_id", "score", "rank"]].reset_index(drop=True)

    # ── Prediction — head (phenotype → gene) ─────────────────────────────────

    @torch.no_grad()
    def predict_for_phenotype(
        self,
        mp_term_id: str,
        top_k: int = 30,
        filter_known: bool = True,
        only_genes: bool = True,
    ) -> pd.DataFrame:
        """Predict the top-K genes most likely associated with an MP phenotype.

        Inverse query: fixes the tail and scores all heads.
        Returns DataFrame: mp_term_id | gene_id | score | rank.
        """
        self._check_entity(mp_term_id)

        pred = predict_target(
            model=self.model,
            tail=mp_term_id,
            relation=self.relation,
            triples_factory=self.training_tf,
        )
        if filter_known:
            pred = pred.filter_triples(self.training_tf)

        pred_df = pred.df
        if only_genes:
            pred_df = pred_df[pred_df["head_label"].str.startswith("HGNC:")].copy()

        pred_df = pred_df.sort_values("score", ascending=False).head(top_k)
        pred_df.insert(0, "mp_term_id", mp_term_id)
        pred_df = pred_df.rename(columns={"head_label": "gene_id"})
        pred_df["rank"] = range(1, len(pred_df) + 1)

        return pred_df[["mp_term_id", "gene_id", "score", "rank"]].reset_index(drop=True)

    # ── Batch prediction ──────────────────────────────────────────────────────

    def predict_batch_for_genes(
        self,
        gene_ids: list[str],
        top_k: int = 20,
        filter_known: bool = True,
    ) -> pd.DataFrame:
        """Predict phenotypes for a list of genes.

        Genes not found in the graph are skipped with a warning.
        Returns DataFrame: gene_id | mp_term_id | score | rank.
        """
        frames: list[pd.DataFrame] = []
        for gid in gene_ids:
            try:
                df = self.predict_for_gene(gid, top_k=top_k, filter_known=filter_known)
                frames.append(df)
            except ValueError as exc:
                logger.warning("Skipping %s: %s", gid, exc)

        if not frames:
            return pd.DataFrame(columns=["gene_id", "mp_term_id", "score", "rank"])
        return pd.concat(frames, ignore_index=True)

    def predict_batch_for_phenotypes(
        self,
        mp_term_ids: list[str],
        top_k: int = 50,
        filter_known: bool = True,
    ) -> pd.DataFrame:
        """Rank genes for a list of MP top terms.

        Preferred for ablation study and cross-species comparison: produces a ranked
        gene list per phenotype class evaluable with AUPRC.
        MP terms not in the graph are skipped with a warning.
        Returns DataFrame: mp_term_id | gene_id | score | rank.
        """
        frames: list[pd.DataFrame] = []
        for mp_id in mp_term_ids:
            try:
                df = self.predict_for_phenotype(mp_id, top_k=top_k, filter_known=filter_known)
                frames.append(df)
            except ValueError as exc:
                logger.warning("Skipping %s: %s", mp_id, exc)

        if not frames:
            return pd.DataFrame(columns=["mp_term_id", "gene_id", "score", "rank"])
        return pd.concat(frames, ignore_index=True)

    # ── Embedding utilities ───────────────────────────────────────────────────

    def entity_embedding(self, entity_id: str) -> torch.Tensor:
        """Return the learned embedding vector for an entity (detached CPU tensor)."""
        self._check_entity(entity_id)
        idx = self.training_tf.entity_to_id[entity_id]
        emb = self.model.entity_representations[0](
            indices=torch.tensor([idx])
        )
        return emb.detach().cpu().squeeze(0)

    def most_similar_genes(
        self,
        query_gene_id: str,
        top_k: int = 10,
    ) -> pd.DataFrame:
        """Find genes with the most similar embeddings (cosine similarity).

        Useful for exploratory analysis: genes that are geometrically close
        in embedding space have been trained to associate with similar contexts.
        """
        self._check_entity(query_gene_id)

        query_emb = self.entity_embedding(query_gene_id)
        gene_indices = torch.tensor(
            [self.training_tf.entity_to_id[g] for g in self.gene_entities]
        )
        all_embs = self.model.entity_representations[0](
            indices=gene_indices
        ).detach().cpu()

        query_norm = query_emb / query_emb.norm().clamp(min=1e-8)
        all_norms  = all_embs / all_embs.norm(dim=1, keepdim=True).clamp(min=1e-8)
        sims = (all_norms @ query_norm).numpy()

        idx_sorted = sims.argsort()[::-1]
        results = []
        for i in idx_sorted[:top_k + 1]:
            gid = self.gene_entities[i]
            if gid == query_gene_id:
                continue
            results.append({
                "gene_id": gid,
                "cosine_similarity": float(sims[i]),
                "rank": len(results) + 1,
            })
            if len(results) == top_k:
                break

        return pd.DataFrame(results)

    # ── Private helpers ───────────────────────────────────────────────────────

    def _check_entity(self, entity_id: str) -> None:
        if entity_id not in self.training_tf.entity_to_id:
            prefix = entity_id.split(":")[0] + ":" if ":" in entity_id else ""
            similar = [
                e for e in self.training_tf.entity_to_id
                if e.startswith(prefix)
            ][:5]
            raise ValueError(
                f"Entity '{entity_id}' not found in the graph.\n"
                f"Tip — IDs use prefixes like 'HGNC:' or 'MP:'.\n"
                f"Examples with same prefix: {similar}"
            )


# ── CLI entry point ───────────────────────────────────────────────────────────


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Gene ↔ phenotype link prediction with a trained RotatE model.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--model-dir", required=True,
                   help="Directory produced by train.py (contains trained_model.pkl etc.)")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--gene",      help="Gene ID to predict phenotypes for, e.g. HGNC:TP53")
    group.add_argument("--phenotype", help="MP term to predict genes for, e.g. MP:0001265")
    p.add_argument("--top-k", type=int, default=20)
    p.add_argument("--include-known", action="store_true",
                   help="Include already-known associations (filtered out by default).")
    p.add_argument("--out-csv", default=None,
                   help="Optional: save predictions to this CSV file.")
    return p.parse_args()


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-8s  %(message)s",
        datefmt="%H:%M:%S",
    )
    args = _parse_args()

    predictor = GenePhenotypePredictor.from_directory(args.model_dir)

    if args.gene:
        df = predictor.predict_for_gene(
            args.gene,
            top_k=args.top_k,
            filter_known=not args.include_known,
        )
        print(f"\nTop-{args.top_k} predicted phenotypes for {args.gene}:\n")
    else:
        df = predictor.predict_for_phenotype(
            args.phenotype,
            top_k=args.top_k,
            filter_known=not args.include_known,
        )
        print(f"\nTop-{args.top_k} predicted genes for {args.phenotype}:\n")

    print(df.to_string(index=False))

    if args.out_csv:
        df.to_csv(args.out_csv, index=False)
        print(f"\nSaved to {args.out_csv}")