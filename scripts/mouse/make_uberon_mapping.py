"""
scripts/mouse/make_uberon_mapping.py - generate gxd_uberon_mapping.tsv.

Extracts unique tissue names from the GXD RNA-seq archive and maps them to
UBERON IDs using basic.obo.  Unmatched rows are written with an empty
uberon_id so they can be curated manually.

Usage:
    uv run python scripts/mouse/make_uberon_mapping.py
    uv run python scripts/mouse/make_uberon_mapping.py \\
        --gxd-tar data/raw/gxdrnaseq.tar.gz \\
        --uberon-obo data/raw/basic.obo \\
        --out data/raw/gxd_uberon_mapping.tsv
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from src.knowledge_graph.layers.uberon_mapping import build_uberon_mapping

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
_RAW_DATA_DIR = _PROJECT_ROOT / 'data' / 'raw'


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        '--gxd-tar',
        default=None,
        help=f'GXD RNA-seq archive (default: {_RAW_DATA_DIR}/gxdrnaseq.tar.gz)',
    )
    p.add_argument(
        '--uberon-obo',
        default=None,
        help=f'Uberon basic.obo (default: {_RAW_DATA_DIR}/basic.obo)',
    )
    p.add_argument(
        '--out',
        default=None,
        help=f'Output TSV path (default: {_RAW_DATA_DIR}/gxd_uberon_mapping.tsv)',
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    tar_path = (
        Path(args.gxd_tar) if args.gxd_tar else _RAW_DATA_DIR / 'gxdrnaseq.tar.gz'
    )
    obo_path = Path(args.uberon_obo) if args.uberon_obo else _RAW_DATA_DIR / 'basic.obo'
    out_path = Path(args.out) if args.out else _RAW_DATA_DIR / 'gxd_uberon_mapping.tsv'
    build_uberon_mapping(tar_path, obo_path, out_path)


if __name__ == '__main__':
    main()
