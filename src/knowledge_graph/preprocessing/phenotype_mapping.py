"""
src/knowledge_graph/preprocessing/phenotype_mapping.py - HP→MP top-term mapping pipeline.

Implements the Barbitoff et al. approach: lifts low-level Human Phenotype Ontology
(HP) terms to MP system-level anchor terms via a curated HP2MP bridge table, then
writes the resulting edge TSV consumed by the human phenotype export layer.
"""

from __future__ import annotations

import logging
from collections import deque
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

MORTALITY_AGING_MP = 'MP:0010768'


def load_hp_ontology(data_dir: Path):
    """
    Load the HP ontology from hp.obo in data_dir using pronto.

    Args:
        data_dir: Directory containing hp.obo.

    Returns:
        A pronto.Ontology instance.

    Raises:
        FileNotFoundError: If hp.obo is absent.
    """
    from pronto import Ontology

    hp_path = data_dir / 'hp.obo'
    if not hp_path.exists():
        raise FileNotFoundError(hp_path)
    logger.debug('Loading HP ontology from %s', hp_path)
    return Ontology(str(hp_path), encoding='utf-8')


def load_hp2mp(
    hp2mp_path: Path, exclude_mortality: bool = True
) -> dict[str, list[str]]:
    """
    Load the HP2MP anchor mapping table.

    Args:
        hp2mp_path: Path to HP2MP.tsv with columns object_id (HP) and subject_id (MP).
        exclude_mortality: If True, drop rows mapping to MORTALITY_AGING_MP.

    Returns:
        Dict mapping HP anchor ID → list of MP system-level term IDs.

    Raises:
        ValueError: If required columns are missing.
    """
    df = pd.read_csv(hp2mp_path, sep='\t', dtype=str)

    required = {'object_id', 'subject_id'}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f'HP2MP missing columns: {missing}. Got: {list(df.columns)}')

    if exclude_mortality:
        before = len(df)
        df = df[df['subject_id'] != MORTALITY_AGING_MP]
        logger.debug('Excluded mortality/aging term: %d rows removed', before - len(df))

    valid = df[
        df['object_id'].str.startswith('HP:') & df['subject_id'].str.startswith('MP:')
    ]
    hp2mp: dict[str, list[str]] = (
        valid.groupby('object_id')['subject_id'].apply(list).to_dict()
    )

    logger.info(
        'HP2MP loaded: %d anchor HP terms → %d unique MP system-level terms',
        len(hp2mp),
        len({mp for mps in hp2mp.values() for mp in mps}),
    )
    return hp2mp


def lift_hp_to_anchors(hp_id: str, hp_ont, anchor_set: set[str]) -> dict[str, int]:
    """
    BFS-walk the HP ontology upward from hp_id to find all anchor terms.

    Args:
        hp_id: HP term ID to start from.
        hp_ont: pronto.Ontology instance.
        anchor_set: Set of HP anchor term IDs (keys from hp2mp).

    Returns:
        Dict mapping anchor_hp_id → distance from hp_id, empty if none found.
    """
    if hp_id not in hp_ont:
        return {}
    if hp_id in anchor_set:
        return {hp_id: 0}

    visited: dict[str, int] = {hp_id: 0}
    queue: deque[str] = deque([hp_id])
    found: dict[str, int] = {}

    while queue:
        node = queue.popleft()
        node_dist = visited[node]
        if node not in hp_ont:
            continue
        for parent in hp_ont[node].superclasses(distance=1, with_self=False):
            pid = parent.id
            if pid in visited:
                continue
            visited[pid] = node_dist + 1
            if pid in anchor_set:
                found[pid] = node_dist + 1
            else:
                queue.append(pid)

    return found


def build_hp_to_mp_top(hp_ids, hp_ont, hp2mp):
    """
    Map HP term IDs to MP system-level terms via anchor lifting.

    Args:
        hp_ids: Iterable of HP term IDs to map.
        hp_ont: pronto.Ontology instance.
        hp2mp: Dict mapping HP anchor ID → list of MP top-term IDs.

    Returns:
        DataFrame with columns: hp_id, mp_top_id, anchor_hp_id, hp_to_anchor_distance.
    """
    anchor_set = set(hp2mp.keys())
    rows = []

    for hp_id in sorted(set(hp_ids)):
        anchor_distances = lift_hp_to_anchors(hp_id, hp_ont, anchor_set)

        best: dict[str, tuple[str, int]] = {}
        for anchor, dist in anchor_distances.items():
            for mp_top in hp2mp.get(anchor, []):
                if mp_top not in best or dist < best[mp_top][1]:
                    best[mp_top] = (anchor, dist)

        for mp_top, (anchor, dist) in sorted(best.items()):
            rows.append(
                {
                    'hp_id': hp_id,
                    'mp_top_id': mp_top,
                    'anchor_hp_id': anchor,
                    'hp_to_anchor_distance': dist,
                }
            )

    cols = ['hp_id', 'mp_top_id', 'anchor_hp_id', 'hp_to_anchor_distance']
    df = pd.DataFrame(rows) if rows else pd.DataFrame(columns=cols)
    return df


def run_mapping(
    genes_to_phenotype_path: str | Path,
    hp2mp_path: str | Path,
    data_dir: str | Path,
    out_dir: str | Path = './out',
    exclude_mortality: bool = True,
    mp_obo_path: str | Path | None = None,
) -> dict[str, pd.DataFrame]:
    """
    Run the HP→MP top-term mapping pipeline.

    Args:
        genes_to_phenotype_path: Path to HPO genes_to_phenotype.txt.
        hp2mp_path: Path to HP2MP.tsv anchor mapping.
        data_dir: Directory containing hp.obo (and optionally mp.obo).
        out_dir: Output directory for TSV files.
        exclude_mortality: If True, exclude MP:0010768 (mortality/aging).
        mp_obo_path: Explicit path to mp.obo for term name lookup.

    Returns:
        Dict with keys 'hp_to_mp_top' and 'mp_top_names'.
    """
    out_dir = Path(out_dir)
    data_dir = Path(data_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    logger.info('Phenotype mapping pipeline (Barbitoff et al. approach)')
    logger.info('genes_to_phenotype: %s', genes_to_phenotype_path)
    logger.info('HP2MP: %s', hp2mp_path)
    logger.info('data_dir: %s', data_dir)
    logger.info('out_dir: %s', out_dir)
    logger.info('exclude_mortality: %s', exclude_mortality)

    hp_ont = load_hp_ontology(data_dir)
    hp2mp = load_hp2mp(Path(hp2mp_path), exclude_mortality=exclude_mortality)

    g2p = pd.read_csv(genes_to_phenotype_path, sep='\t', dtype=str)
    required = {'hpo_id', 'gene_symbol'}
    missing = required - set(g2p.columns)
    if missing:
        raise ValueError(f'genes_to_phenotype missing columns: {missing}')

    unique_hp = list(g2p['hpo_id'].dropna().unique())
    logger.info('Unique HP terms in input: %d', len(unique_hp))

    hp_to_mp_top = build_hp_to_mp_top(unique_hp, hp_ont, hp2mp)
    hp_to_mp_top.to_csv(out_dir / 'edge_hp_to_mp_top.tsv', sep='\t', index=False)
    logger.info('Written: edge_hp_to_mp_top.tsv (%d rows)', len(hp_to_mp_top))

    mp_obo_path = Path(mp_obo_path) if mp_obo_path else data_dir / 'mp.obo'
    if mp_obo_path.exists():
        from pronto import Ontology

        logger.debug('Loading MP ontology for term names from %s', mp_obo_path)
        mp_ont = Ontology(str(mp_obo_path), encoding='utf-8')
        mp_top_ids = (
            sorted(hp_to_mp_top['mp_top_id'].dropna().unique())
            if not hp_to_mp_top.empty
            else []
        )
        names_rows = [
            {
                'mp_id': mp_id,
                'name': (mp_ont[mp_id].name or '') if mp_id in mp_ont else '',
            }
            for mp_id in mp_top_ids
        ]
    else:
        logger.warning('mp.obo not found at %s - names will be empty', mp_obo_path)
        mp_top_ids = (
            sorted(hp_to_mp_top['mp_top_id'].dropna().unique())
            if not hp_to_mp_top.empty
            else []
        )
        names_rows = [{'mp_id': mp_id, 'name': ''} for mp_id in mp_top_ids]

    mp_top_names = pd.DataFrame(names_rows)
    mp_top_names.to_csv(out_dir / 'node_mp_top_names.tsv', sep='\t', index=False)
    logger.info('Written: node_mp_top_names.tsv (%d terms)', len(mp_top_names))

    anchor_set = set(hp2mp.keys())
    mapped_hp = set(hp_to_mp_top['hp_id'].dropna()) if not hp_to_mp_top.empty else set()

    no_anchor: list[str] = []
    has_anchor_no_mp: list[str] = []
    for hp in unique_hp:
        if hp in mapped_hp:
            continue
        if lift_hp_to_anchors(hp, hp_ont, anchor_set):
            has_anchor_no_mp.append(hp)
        else:
            no_anchor.append(hp)

    pd.DataFrame({'hp_id': sorted(no_anchor)}).to_csv(
        out_dir / 'qc_hp_no_anchor.tsv', sep='\t', index=False
    )
    pd.DataFrame({'hp_id': sorted(has_anchor_no_mp)}).to_csv(
        out_dir / 'qc_hp_no_mapping.tsv', sep='\t', index=False
    )

    logger.info(
        'QC: HP terms with no anchor: %d / %d (%.1f%%)',
        len(no_anchor),
        len(unique_hp),
        len(no_anchor) / len(unique_hp) * 100,
    )
    logger.info('QC: HP terms anchor→no MP top: %d', len(has_anchor_no_mp))
    logger.info(
        'HP terms successfully mapped: %d / %d (%.1f%%)',
        len(mapped_hp),
        len(unique_hp),
        len(mapped_hp) / len(unique_hp) * 100,
    )
    logger.info('MP system-level terms used: %d', len(mp_top_names))

    return {
        'hp_to_mp_top': hp_to_mp_top,
        'mp_top_names': mp_top_names,
    }


if __name__ == '__main__':
    import argparse
    import logging

    from src.utils.logging_config import setup_logging

    p = argparse.ArgumentParser(
        description='Map HP low-level terms to MP system-level terms via HP2MP anchor table.'
    )
    p.add_argument(
        '--genes-to-phenotype',
        required=True,
        help='Path to genes_to_phenotype.txt (HPO)',
    )
    p.add_argument(
        '--hp2mp', required=True, help='Path to HP2MP.tsv (MGI anchor mapping)'
    )
    p.add_argument(
        '--data-dir',
        required=True,
        help='Directory containing hp.obo (and optionally mp.obo)',
    )
    p.add_argument('--out', default='./out', help='Output directory (default: ./out)')
    p.add_argument(
        '--keep-mortality',
        action='store_true',
        help='Do NOT exclude MP:0010768 (mortality/aging)',
    )
    p.add_argument(
        '--mp-obo',
        default=None,
        help='Explicit path to mp.obo (defaults to <data-dir>/mp.obo)',
    )
    p.add_argument(
        '--log-level', default='INFO', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR']
    )
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
