"""src/adapters/base.py - shared base class for all layer adapters."""

from __future__ import annotations

from pathlib import Path
from typing import Generator, Iterator

import pandas as pd

# Concrete type aliases instead of bare `tuple`
NodeTuple = tuple[str, str, dict]        # (id, label, properties)
EdgeTuple = tuple[str, str, str, dict]   # (source_id, target_id, label, properties)


class BaseAdapter:
    """Base class for all layer adapters.

    Provides:
      - _read()         — read a TSV from data_dir into a DataFrame
      - _unique_nodes() — deduplicate node tuples by id

    Each subclass must:
      - Set ``layer_name`` (used by build_graph.py for --skip-layers)
      - Implement get_nodes() and/or get_edges()
      - Accept __init__(self, data_dir) with no other required arguments
        so that build_graph.py can instantiate all adapters uniformly.

    Adding a new layer
    ------------------
    1. Create src/layers/<name>/adapter.py with a class inheriting BaseAdapter.
    2. Set layer_name = "<name>" as a class attribute.
    3. Register the class in src/pipeline/build_graph.py → SPECIES_LAYERS.
    4. Add the schema section to config/schema_config_<species>.yaml.
    """

    #: Stable, lowercase identifier for this layer.  Override in every subclass.
    layer_name: str = ""

    def __init__(self, data_dir: str | Path):
        self.data_dir = Path(data_dir)

    # ── I/O ────────────────────────────────────────────────────────── #

    def _read(self, filename: str) -> pd.DataFrame:
        """Read a tab-separated file from data_dir.  Returns empty DataFrame on empty file."""
        try:
            return pd.read_csv(self.data_dir / filename, sep="\t", dtype=str)
        except pd.errors.EmptyDataError:
            return pd.DataFrame()

    # ── Deduplication ──────────────────────────────────────────────── #

    def _unique_nodes(
        self, tuples: Iterator[NodeTuple]
    ) -> Generator[NodeTuple, None, None]:
        """Yield node tuples, skipping duplicates by id.

        Always pass a lazy iterator (e.g. itertools.chain), not a materialised
        tuple/list — unpacking generators with (*gen1, *gen2) loads everything
        into memory and defeats the purpose of lazy generation.
        """
        seen: set[str] = set()
        for node_id, label, props in tuples:
            if node_id not in seen:
                seen.add(node_id)
                yield node_id, label, props

    # ── Default no-op generators ───────────────────────────────────── #

    def get_nodes(self) -> Generator[NodeTuple, None, None]:
        yield from ()

    def get_edges(self) -> Generator[EdgeTuple, None, None]:
        yield from ()
