"""
src/knowledge_graph/build_graph.py - BioCypher graph orchestrator.

Registers all per-species adapter layers, deduplicates nodes and edges across
layers, and writes BioCypher-format CSVs ready for Neo4j import.  Entry point
for programmatic use; CLI users run `run_kg.py` instead.
"""

from __future__ import annotations

import argparse
import itertools
import logging
import time
from pathlib import Path
from typing import Literal

import biocypher
from src.knowledge_graph.adapters.base import BaseAdapter
from src.knowledge_graph.adapters.gene_expression_adapter import (
    HumanExpressionAdapter,
    MouseExpressionAdapter,
)
from src.knowledge_graph.adapters.gene_ontology_adapter import (
    HumanGOAdapter,
    MouseGOAdapter,
)
from src.knowledge_graph.adapters.gene_to_phenotype_adapter import (
    HumanPhenotypeAdapter,
    MousePhenotypeAdapter,
)
from src.knowledge_graph.adapters.ppi_adapter import HumanPPIAdapter, MousePPIAdapter
from src.knowledge_graph.adapters.uberon_adapter import (
    HumanUberonAdapter,
)
from src.utils.logging_config import setup_logging

logger = logging.getLogger(__name__)

Species = Literal['human', 'mouse', 'both']


SPECIES_LAYERS: dict[str, list] = {
    'human': [
        HumanPhenotypeAdapter,
        HumanGOAdapter,
        HumanExpressionAdapter,
        HumanPPIAdapter,
        HumanUberonAdapter,
    ],
    'mouse': [
        MousePhenotypeAdapter,
        MouseGOAdapter,
        MouseExpressionAdapter,
        MousePPIAdapter,
    ],
}


_CONFIG_DIR = Path(__file__).resolve().parents[2] / 'configs' / 'knowledge_graph'
DEFAULT_SCHEMA: dict[str, Path] = {
    'human': _CONFIG_DIR / 'schema_config_human.yaml',
    'mouse': _CONFIG_DIR / 'schema_config_mouse.yaml',
}


def _layer_name(cls) -> str:
    """
    Return the normalised layer name for a given adapter class.

    Args:
        cls: Adapter class with an optional layer_name attribute.

    Returns:
        Stripped, lowercased layer name string.
    """
    return getattr(cls, 'layer_name', cls.__name__).strip().lower()


def build_species(
    species: str,
    data_dir: Path,
    out_dir: Path,
    schema_config_path: Path | None,
    biocypher_config_path: Path | None,
    skip_layers: list[str],
    allow_missing_layers: bool = False,
) -> None:
    """
    Build a BioCypher knowledge graph for a single species.

    Args:
        species: 'human' or 'mouse'.
        data_dir: Directory containing processed per-layer TSVs.
        out_dir: Output directory for BioCypher CSV files.
        schema_config_path: Path to the schema YAML (None → species default).
        biocypher_config_path: Path to biocypher_config.yaml (None → project default).
        skip_layers: Layer names to skip.
        allow_missing_layers: If True, log a warning instead of raising on missing layers.

    Raises:
        FileNotFoundError: If the schema config does not exist.
        RuntimeError: If a required layer fails to initialise and allow_missing_layers is False.
    """
    schema = schema_config_path or DEFAULT_SCHEMA[species]
    if not schema.exists():
        raise FileNotFoundError(f'Schema config not found: {schema}')

    logger.info('=' * 70)
    logger.info('Building %s graph  →  %s', species.upper(), out_dir)
    logger.info('Schema: %s', schema)
    logger.info('Data: %s', data_dir)
    logger.info('=' * 70)

    bc = biocypher.BioCypher(
        schema_config_path=str(schema),
        biocypher_config_path=str(
            biocypher_config_path or _CONFIG_DIR / 'biocypher_config.yaml'
        ),
        output_directory=str(out_dir),
    )

    active_adapters: list[tuple[str, object]] = []
    for LayerCls in SPECIES_LAYERS[species]:
        name = _layer_name(LayerCls)
        if name in skip_layers:
            logger.info('Skipping layer: %s', name)
            continue
        logger.debug('Initialising layer: %s (%s)', name, LayerCls.__name__)
        try:
            active_adapters.append((name, LayerCls(data_dir)))
            logger.info('Layer ready: %s', name)
        except FileNotFoundError as exc:
            if allow_missing_layers:
                logger.warning("Layer '%s' skipped - %s", name, exc)
            else:
                logger.error("Layer '%s' failed to initialise: %s", name, exc)
                raise RuntimeError(
                    f"\nLayer '{name}' failed to initialise:\n  {exc}\n\n"
                    'This usually means an export step did not run or wrote to\n'
                    'a different directory.  Re-run without --skip-export, or\n'
                    'pass --allow-missing-layers to build a partial graph.'
                ) from exc

    if not active_adapters:
        logger.warning(
            "No active adapters for species '%s' - nothing to write.", species
        )
        return

    layer_names = ', '.join(name for name, _ in active_adapters)
    logger.info('Active layers: %s', layer_names)

    def _all_nodes():
        """
        Yield deduplicated nodes from all active adapters.
        """
        all_gen = itertools.chain.from_iterable(
            adapter.get_nodes() for _, adapter in active_adapters
        )
        yield from BaseAdapter._unique_nodes(all_gen)

    logger.info('Writing nodes...')
    t0 = time.perf_counter()
    bc.write_nodes(_all_nodes())
    logger.info('Nodes written in %.1fs', time.perf_counter() - t0)

    def _all_edges():
        """
        Yield deduplicated edges from all active adapters.
        """
        all_gen = itertools.chain.from_iterable(
            adapter.get_edges() for _, adapter in active_adapters
        )
        yield from BaseAdapter._unique_edges(all_gen)

    logger.info('Writing edges...')
    t0 = time.perf_counter()
    bc.write_edges(_all_edges())
    logger.info('Edges written in %.1fs', time.perf_counter() - t0)

    bc.write_import_call()
    logger.info('Done → %s', out_dir)


def build(
    data_dir: str | Path,
    out_dir: str | Path,
    species: Species = 'both',
    schema_config_path: str | Path | None = None,
    biocypher_config_path: str | Path | None = None,
    skip_layers: list[str] | None = None,
    allow_missing_layers: bool = False,
) -> None:
    """
    Build BioCypher knowledge graphs for one or both species.

    Args:
        data_dir: Directory containing processed per-layer TSVs.
        out_dir: Root output directory; species subdirectories are created within.
        species: 'human', 'mouse', or 'both'.
        schema_config_path: Path to schema YAML (None → species default).
        biocypher_config_path: Path to biocypher_config.yaml (None → project default).
        skip_layers: Layer names to skip.
        allow_missing_layers: If True, warn instead of failing on missing layers.
    """
    data_dir = Path(data_dir)
    out_dir = Path(out_dir)
    skip_layers = [s.lower() for s in (skip_layers or [])]
    schema_path = Path(schema_config_path) if schema_config_path else None
    bc_config = Path(biocypher_config_path) if biocypher_config_path else None

    setup_logging(out_dir=out_dir)
    logger.info(
        'Pipeline started - species=%s, data_dir=%s, out_dir=%s',
        species,
        data_dir,
        out_dir,
    )
    if skip_layers:
        logger.info('Skipping layers: %s', skip_layers)

    targets = ['human', 'mouse'] if species == 'both' else [species]

    t_total = time.perf_counter()
    for sp in targets:
        sp_out = out_dir / sp
        build_species(
            species=sp,
            data_dir=data_dir,
            out_dir=sp_out,
            schema_config_path=schema_path,
            biocypher_config_path=bc_config,
            skip_layers=skip_layers,
            allow_missing_layers=allow_missing_layers,
        )

    logger.info('Pipeline complete - total time: %.1fs', time.perf_counter() - t_total)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description='Build BioCypher gene-phenotype graphs')
    p.add_argument('--data-dir', required=True)
    p.add_argument('--out-dir', default='./biocypher_out')
    p.add_argument('--species', choices=['human', 'mouse', 'both'], default='both')
    p.add_argument('--schema-config', default=None)
    p.add_argument('--biocypher-config', default=None)
    p.add_argument('--skip-layers', nargs='*', default=[])
    p.add_argument(
        '--allow-missing-layers',
        action='store_true',
        help="Warn instead of failing when a layer's TSV files are missing.",
    )
    args = p.parse_args()
    build(
        data_dir=args.data_dir,
        out_dir=args.out_dir,
        species=args.species,
        schema_config_path=args.schema_config,
        biocypher_config_path=args.biocypher_config,
        skip_layers=args.skip_layers,
        allow_missing_layers=args.allow_missing_layers,
    )
