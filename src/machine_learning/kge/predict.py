"""
src/machine_learning/kge/predict.py - gene ↔ MP top term link prediction with a trained KGE model.

Usage
-----
    predictor = GenePhenotypePredictor.from_directory('results/ml/human/kge/all/')
    df = predictor.predict_for_gene('HGNC:TP53', top_k=20)
    df = predictor.predict_for_phenotype('MP:0001265', top_k=30)
    df = predictor.predict_batch_for_phenotypes(['MP:0005387', 'MP:0005384'], top_k=50)
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import pandas as pd
import torch
from pykeen.predict import predict_target
from pykeen.triples import TriplesFactory
from src.machine_learning.kge.config import target_relation_for_species

logger = logging.getLogger(__name__)


class GenePhenotypePredictor:
    """
    Wrapper around a trained PyKEEN model for gene ↔ phenotype prediction.
    """

    # ── Construction ─────────────────────────────────────────────────────────

    def __init__(
        self,
        model: torch.nn.Module,
        training_tf: TriplesFactory,
        relation: str,
        gene_prefix: str = 'HGNC:',
        phenotype_prefix: str = 'MP:',
    ) -> None:
        """
        Initialise the predictor from a trained model and its TriplesFactory.

        Args:
            model: Trained PyKEEN model (set to eval mode on construction).
            training_tf: TriplesFactory used during training (entity/relation mapping).
            relation: Target gene→phenotype relation name; must exist in training_tf.
            gene_prefix: Prefix shared by all gene entity IDs (e.g. 'HGNC:').
            phenotype_prefix: Prefix shared by all MP term IDs (e.g. 'MP:').

        Raises:
            ValueError: If *relation* is not in training_tf.relation_to_id.
        """
        if relation not in training_tf.relation_to_id:
            available = sorted(training_tf.relation_to_id.keys())
            raise ValueError(
                f"Relation '{relation}' not found in the TriplesFactory.\n"
                f'Available relations: {available}'
            )

        self.model = model.eval()
        self.training_tf = training_tf
        self.relation = relation
        self.gene_prefix = gene_prefix
        self.phenotype_prefix = phenotype_prefix

        all_entities = list(training_tf.entity_to_id.keys())
        rel_id = training_tf.relation_to_id[relation]
        target_triples = training_tf.mapped_triples[
            training_tf.mapped_triples[:, 1] == rel_id
        ]
        id_to_entity = {v: k for k, v in training_tf.entity_to_id.items()}

        if len(target_triples) > 0:
            self.gene_entities = sorted(
                {id_to_entity[int(idx)] for idx in target_triples[:, 0].tolist()}
            )
            self.mp_entities = sorted(
                {id_to_entity[int(idx)] for idx in target_triples[:, 2].tolist()}
            )
        else:
            self.mp_entities = [
                e for e in all_entities if e.startswith(phenotype_prefix)
            ]
            self.gene_entities = [
                e for e in all_entities if gene_prefix and e.startswith(gene_prefix)
            ]
        self._gene_entity_set = set(self.gene_entities)
        self._mp_entity_set = set(self.mp_entities)

        logger.info(
            'Predictor ready - %d genes | %d MP top-terms | relation: %s',
            len(self.gene_entities),
            len(self.mp_entities),
            relation,
        )

    @classmethod
    def from_directory(
        cls,
        model_dir: str | Path,
        relation: str | None = None,
        species: str | None = None,
    ) -> GenePhenotypePredictor:
        """
        Load a predictor from a directory produced by run_training().

        Loads trained_model.pkl and training_triples/ independently -
        PipelineResult.from_directory() was removed in PyKEEN 1.11.

        Args:
            model_dir: Directory containing trained_model.pkl and training_triples/.
            relation: Target relation name. When omitted, read from summary.json
                or resolve from species.
            species: Species used when relation is omitted and summary.json does
                not record a target relation.

        Returns:
            Initialised GenePhenotypePredictor.
        """
        model_dir = Path(model_dir)
        if not model_dir.exists():
            raise FileNotFoundError(f'Model directory not found: {model_dir}')

        logger.info('Loading model from %s …', model_dir)
        model = torch.load(model_dir / 'trained_model.pkl', weights_only=False)
        training_tf = TriplesFactory.from_path_binary(model_dir / 'training_triples')
        relation = relation or _relation_from_summary(model_dir, species=species)

        return cls(model=model, training_tf=training_tf, relation=relation)

    # ── Prediction - tail (gene → phenotype) ─────────────────────────────────

    @torch.no_grad()
    def predict_for_gene(
        self,
        gene_id: str,
        top_k: int = 20,
        filter_known: bool = True,
        only_mp_terms: bool = True,
    ) -> pd.DataFrame:
        """
        Predict the top-K most likely MP phenotypes for a single gene.

        Args:
            gene_id: Gene entity ID (e.g. 'HGNC:TP53').
            top_k: Number of predictions to return.
            filter_known: Exclude training-set pairs (novel predictions only).
            only_mp_terms: Restrict candidates to MP: entities.

        Returns:
            DataFrame: gene_id | mp_term_id | score | rank
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
            pred_df = pred_df[pred_df['tail_label'].isin(self._mp_entity_set)].copy()

        pred_df = pred_df.sort_values('score', ascending=False).head(top_k)
        pred_df.insert(0, 'gene_id', gene_id)
        pred_df = pred_df.rename(columns={'tail_label': 'mp_term_id'})
        pred_df['rank'] = range(1, len(pred_df) + 1)
        return pred_df[['gene_id', 'mp_term_id', 'score', 'rank']].reset_index(
            drop=True
        )

    # ── Prediction - head (phenotype → gene) ─────────────────────────────────

    @torch.no_grad()
    def predict_for_phenotype(
        self,
        mp_term_id: str,
        top_k: int = 30,
        filter_known: bool = True,
        only_genes: bool = True,
    ) -> pd.DataFrame:
        """
        Predict the top-K genes most likely associated with an MP phenotype.

        Args:
            mp_term_id: MP term entity ID (e.g. 'MP:0001265').
            top_k: Number of predictions to return.
            filter_known: Exclude training-set pairs.
            only_genes: Restrict candidates to gene entities.

        Returns:
            DataFrame: mp_term_id | gene_id | score | rank
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
            pred_df = pred_df[pred_df['head_label'].isin(self._gene_entity_set)].copy()

        pred_df = pred_df.sort_values('score', ascending=False).head(top_k)
        pred_df.insert(0, 'mp_term_id', mp_term_id)
        pred_df = pred_df.rename(columns={'head_label': 'gene_id'})
        pred_df['rank'] = range(1, len(pred_df) + 1)
        return pred_df[['mp_term_id', 'gene_id', 'score', 'rank']].reset_index(
            drop=True
        )

    # ── Batch prediction ──────────────────────────────────────────────────────

    def predict_batch_for_genes(
        self,
        gene_ids: list[str],
        top_k: int = 20,
        filter_known: bool = True,
    ) -> pd.DataFrame:
        """
        Predict phenotypes for a list of genes.

        Genes not found in the graph are skipped with a warning.

        Args:
            gene_ids: List of gene entity IDs.
            top_k: Predictions per gene.
            filter_known: Exclude training-set pairs.

        Returns:
            DataFrame: gene_id | mp_term_id | score | rank
        """
        frames: list[pd.DataFrame] = []
        for gid in gene_ids:
            try:
                frames.append(
                    self.predict_for_gene(gid, top_k=top_k, filter_known=filter_known)
                )
            except ValueError as exc:
                logger.warning('Skipping %s: %s', gid, exc)

        if not frames:
            return pd.DataFrame(columns=['gene_id', 'mp_term_id', 'score', 'rank'])
        return pd.concat(frames, ignore_index=True)

    def predict_batch_for_phenotypes(
        self,
        mp_term_ids: list[str],
        top_k: int = 50,
        filter_known: bool = True,
    ) -> pd.DataFrame:
        """
        Rank genes for a list of MP top terms.

        MP terms not in the graph are skipped with a warning.

        Args:
            mp_term_ids: List of MP term entity IDs.
            top_k: Genes to return per phenotype.
            filter_known: Exclude training-set pairs.

        Returns:
            DataFrame: mp_term_id | gene_id | score | rank
        """
        frames: list[pd.DataFrame] = []
        for mp_id in mp_term_ids:
            try:
                frames.append(
                    self.predict_for_phenotype(
                        mp_id, top_k=top_k, filter_known=filter_known
                    )
                )
            except ValueError as exc:
                logger.warning('Skipping %s: %s', mp_id, exc)

        if not frames:
            return pd.DataFrame(columns=['mp_term_id', 'gene_id', 'score', 'rank'])
        return pd.concat(frames, ignore_index=True)

    # ── Embedding utilities ───────────────────────────────────────────────────

    def entity_embedding(self, entity_id: str) -> torch.Tensor:
        """
        Return the learned embedding vector for an entity (detached CPU tensor).

        Args:
            entity_id: Entity string ID.

        Returns:
            1-D tensor of shape (embedding_dim,).
        """
        self._check_entity(entity_id)
        idx = self.training_tf.entity_to_id[entity_id]
        emb = self.model.entity_representations[0](indices=torch.tensor([idx]))
        return emb.detach().cpu().squeeze(0)

    def most_similar_genes(
        self,
        query_gene_id: str,
        top_k: int = 10,
    ) -> pd.DataFrame:
        """
        Find genes with the most similar embeddings (cosine similarity).

        Args:
            query_gene_id: Gene entity ID to query.
            top_k: Number of similar genes to return.

        Returns:
            DataFrame: gene_id | cosine_similarity | rank
        """
        self._check_entity(query_gene_id)

        query_emb = self.entity_embedding(query_gene_id)
        gene_indices = torch.tensor(
            [self.training_tf.entity_to_id[g] for g in self.gene_entities]
        )
        all_embs = (
            self.model.entity_representations[0](indices=gene_indices).detach().cpu()
        )

        query_norm = query_emb / query_emb.norm().clamp(min=1e-8)
        all_norms = all_embs / all_embs.norm(dim=1, keepdim=True).clamp(min=1e-8)
        sims = (all_norms @ query_norm).numpy()

        idx_sorted = sims.argsort()[::-1]
        results = []
        for i in idx_sorted[: top_k + 1]:
            gid = self.gene_entities[i]
            if gid == query_gene_id:
                continue
            results.append(
                {
                    'gene_id': gid,
                    'cosine_similarity': float(sims[i]),
                    'rank': len(results) + 1,
                }
            )
            if len(results) == top_k:
                break

        return pd.DataFrame(results)

    # ── Private helpers ───────────────────────────────────────────────────────

    def _check_entity(self, entity_id: str) -> None:
        """
        Raise ValueError if *entity_id* is not present in the training graph.
        """
        if entity_id not in self.training_tf.entity_to_id:
            prefix = entity_id.split(':')[0] + ':' if ':' in entity_id else ''
            similar = [
                e for e in self.training_tf.entity_to_id if e.startswith(prefix)
            ][:5]
            raise ValueError(
                f"Entity '{entity_id}' not found in the graph.\n"
                f"Tip - IDs use prefixes like 'HGNC:' or 'MP:'.\n"
                f'Examples with same prefix: {similar}'
            )


# ── CLI entry point ───────────────────────────────────────────────────────────


def _relation_from_summary(model_dir: Path, species: str | None = None) -> str:
    """
    Read the target relation from summary.json, falling back to species default.
    """
    summary_path = model_dir / 'summary.json'
    if summary_path.exists():
        try:
            summary = json.loads(summary_path.read_text())
        except json.JSONDecodeError:
            logger.warning('Could not parse %s; falling back to species', summary_path)
        else:
            relation = summary.get('target_relation') or summary.get('config', {}).get(
                'target_relation'
            )
            if relation:
                return str(relation)
            species = (
                species
                or summary.get('species')
                or summary.get('config', {}).get('species')
            )
    return target_relation_for_species(species or 'human')


def _parse_args() -> argparse.Namespace:
    """
    Parse and return CLI arguments for the predict entry point.
    """
    p = argparse.ArgumentParser(
        description='Gene ↔ phenotype link prediction with a trained KGE model.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        '--model-dir',
        required=True,
        help='Directory produced by train.py (contains trained_model.pkl).',
    )
    p.add_argument(
        '--species',
        default=None,
        choices=['human', 'mouse'],
        help='Species used to select the relation if summary.json lacks it.',
    )
    p.add_argument(
        '--target-relation',
        default=None,
        help='Explicit target relation. Overrides summary.json and --species.',
    )
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument(
        '--gene', help='Gene ID to predict phenotypes for, e.g. HGNC:TP53'
    )
    group.add_argument(
        '--phenotype', help='MP term to predict genes for, e.g. MP:0001265'
    )
    p.add_argument('--top-k', type=int, default=20)
    p.add_argument(
        '--include-known',
        action='store_true',
        help='Include already-known associations (filtered out by default).',
    )
    p.add_argument('--out-csv', default=None, help='Save predictions to this CSV file.')
    return p.parse_args()


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s  %(levelname)-8s  %(message)s',
        datefmt='%H:%M:%S',
    )
    args = _parse_args()
    predictor = GenePhenotypePredictor.from_directory(
        args.model_dir,
        relation=args.target_relation,
        species=args.species,
    )

    if args.gene:
        df = predictor.predict_for_gene(
            args.gene,
            top_k=args.top_k,
            filter_known=not args.include_known,
        )
        print(f'\nTop-{args.top_k} predicted phenotypes for {args.gene}:\n')
    else:
        df = predictor.predict_for_phenotype(
            args.phenotype,
            top_k=args.top_k,
            filter_known=not args.include_known,
        )
        print(f'\nTop-{args.top_k} predicted genes for {args.phenotype}:\n')

    print(df.to_string(index=False))

    if args.out_csv:
        df.to_csv(args.out_csv, index=False)
        print(f'\nSaved to {args.out_csv}')
