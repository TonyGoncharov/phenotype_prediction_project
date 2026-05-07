"""src/pipeline/build_graph.py - assemble species-specific BioCypher graphs.

Adding a new layer
------------------
1. Create src/layers/<n>/adapter.py with classes inheriting BaseAdapter.
2. Set ``layer_name`` as a class attribute in each subclass.
3. Register the classes in SPECIES_LAYERS below — that's the only file to touch.
4. Add schema sections to config/schema_config_human.yaml / _mouse.yaml.
"""

from __future__ import annotations

import argparse
import itertools
import logging
import time
from pathlib import Path
from typing import Literal

import biocypher

from src.adapters.base import BaseAdapter
from src.adapters.gene_to_phenotype_adapter import HumanPhenotypeAdapter, MousePhenotypeAdapter
from src.adapters.gene_ontology_adapter import HumanGOAdapter, MouseGOAdapter
from src.adapters.gene_expression_adapter import HumanExpressionAdapter
from src.adapters.ppi_adapter import HumanPPIAdapter
from src.adapters.uberon_adapter import HumanUberonAdapter, MouseUberonAdapter
from src.utils.logging_config import setup_logging

logger = logging.getLogger(__name__)

Species = Literal["human", "mouse", "both"]

# ── Layer registry ────────────────────────────────────────────────────
#
# To add a layer: import its adapter classes and append them here.
# Order determines the order nodes/edges are written.
#
SPECIES_LAYERS: dict[str, list] = {
    "human": [
        HumanPhenotypeAdapter,
        HumanGOAdapter,
        HumanExpressionAdapter,
        HumanPPIAdapter,
        HumanUberonAdapter,
    ],
    "mouse": [
        MousePhenotypeAdapter,
        MouseGOAdapter,
        MouseUberonAdapter,
        # MouseExpressionAdapter,
    ],
}

# ── Default schema config paths ───────────────────────────────────────
_CONFIG_DIR = Path(__file__).resolve().parent.parent / "config"
DEFAULT_SCHEMA: dict[str, Path] = {
    "human": _CONFIG_DIR / "schema_config_human.yaml",
    "mouse": _CONFIG_DIR / "schema_config_mouse.yaml",
}


def _layer_name(cls) -> str:
    return getattr(cls, "layer_name", cls.__name__).strip().lower()


def build_species(
    species: str,
    data_dir: Path,
    out_dir: Path,
    schema_config_path: Path | None,
    biocypher_config_path: Path | None,
    skip_layers: list[str],
    allow_missing_layers: bool = False,
) -> None:
    """Build the BioCypher graph for a single species.

    Args:
        allow_missing_layers: When False (default), a missing TSV raises
            immediately — almost always a bug or a forgotten export step.
            Pass True only when a partial graph is intentional.
    """
    schema = schema_config_path or DEFAULT_SCHEMA[species]
    if not schema.exists():
        raise FileNotFoundError(f"Schema config not found: {schema}")

    logger.info("=" * 60)
    logger.info("Building %s graph  →  %s", species.upper(), out_dir)
    logger.info("Schema : %s", schema)
    logger.info("Data   : %s", data_dir)
    logger.info("=" * 60)

    bc = biocypher.BioCypher(
        schema_config_path=str(schema),
        biocypher_config_path=str(biocypher_config_path) if biocypher_config_path else None,
        output_directory=str(out_dir),
    )

    # ── Initialise adapters ──────────────────────────────────────────
    active_adapters: list[tuple[str, object]] = []
    for LayerCls in SPECIES_LAYERS[species]:
        name = _layer_name(LayerCls)
        if name in skip_layers:
            logger.info("Skipping layer: %s", name)
            continue
        logger.debug("Initialising layer: %s (%s)", name, LayerCls.__name__)
        try:
            active_adapters.append((name, LayerCls(data_dir)))
            logger.info("Layer ready: %s", name)
        except FileNotFoundError as exc:
            if allow_missing_layers:
                logger.warning("Layer '%s' skipped — %s", name, exc)
            else:
                logger.error("Layer '%s' failed to initialise: %s", name, exc)
                raise RuntimeError(
                    f"\nLayer '{name}' failed to initialise:\n  {exc}\n\n"
                    "This usually means an export step did not run or wrote to\n"
                    "a different directory.  Re-run without --skip-export, or\n"
                    "pass --allow-missing-layers to build a partial graph."
                ) from exc

    if not active_adapters:
        logger.warning("No active adapters for species '%s' — nothing to write.", species)
        return

    layer_names = ", ".join(name for name, _ in active_adapters)
    logger.info("Active layers: %s", layer_names)

    # ── Node writing ─────────────────────────────────────────────────
    def _all_nodes():
        all_gen = itertools.chain.from_iterable(
            adapter.get_nodes() for _, adapter in active_adapters
        )
        yield from BaseAdapter._unique_nodes(all_gen)

    logger.info("Writing nodes...")
    t0 = time.perf_counter()
    bc.write_nodes(_all_nodes())
    logger.info("Nodes written in %.1fs", time.perf_counter() - t0)

    # ── Edge writing ─────────────────────────────────────────────────
    def _all_edges():
        for name, adapter in active_adapters:
            logger.debug("Collecting edges from layer: %s", name)
            yield from adapter.get_edges()

    logger.info("Writing edges...")
    t0 = time.perf_counter()
    bc.write_edges(_all_edges())
    logger.info("Edges written in %.1fs", time.perf_counter() - t0)

    bc.write_import_call()
    logger.info("Done → %s", bc._output_directory)


def build(
    data_dir: str | Path,
    out_dir: str | Path,
    species: Species = "both",
    schema_config_path: str | Path | None = None,
    biocypher_config_path: str | Path | None = None,
    skip_layers: list[str] | None = None,
    allow_missing_layers: bool = False,
) -> None:
    """Build one or both species graphs.

    When species='both', output goes to <out_dir>/human/ and <out_dir>/mouse/
    so the two graphs don't collide.
    """
    data_dir    = Path(data_dir)
    out_dir     = Path(out_dir)
    skip_layers = [s.lower() for s in (skip_layers or [])]
    schema_path = Path(schema_config_path) if schema_config_path else None
    bc_config   = Path(biocypher_config_path) if biocypher_config_path else None

    setup_logging(out_dir=out_dir)
    logger.info("Pipeline started — species=%s, data_dir=%s, out_dir=%s",
                species, data_dir, out_dir)
    if skip_layers:
        logger.info("Skipping layers: %s", skip_layers)

    targets = ["human", "mouse"] if species == "both" else [species]

    t_total = time.perf_counter()
    for sp in targets:
        sp_out = out_dir / sp if species == "both" else out_dir
        build_species(
            species=sp,
            data_dir=data_dir,
            out_dir=sp_out,
            schema_config_path=schema_path,
            biocypher_config_path=bc_config,
            skip_layers=skip_layers,
            allow_missing_layers=allow_missing_layers,
        )

    logger.info("Pipeline complete — total time: %.1fs", time.perf_counter() - t_total)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Build BioCypher gene-phenotype graphs")
    p.add_argument("--data-dir", required=True)
    p.add_argument("--out-dir", default="./biocypher_out")
    p.add_argument("--species", choices=["human", "mouse", "both"], default="both")
    p.add_argument("--schema-config", default=None)
    p.add_argument("--biocypher-config", default=None)
    p.add_argument("--skip-layers", nargs="*", default=[])
    p.add_argument(
        "--allow-missing-layers", action="store_true",
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
