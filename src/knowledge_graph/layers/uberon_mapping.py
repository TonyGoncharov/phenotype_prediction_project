"""
src/knowledge_graph/layers/uberon_mapping.py - map GXD tissue names to UBERON IDs.

Reads unique Anatomical Structure values from the GXD RNA-seq archive and
attempts to match them against UBERON term names and synonyms from basic.obo.
Writes a two-column TSV (gxd_name, uberon_id); unmatched rows have an empty
uberon_id and can be filled in manually before the full pipeline run.
"""

from __future__ import annotations

import gzip
import io
import logging
import re
import tarfile
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


def parse_uberon_obo(obo_path: Path) -> dict[str, str]:
    """
    Parse basic.obo and return a mapping lowercased name/synonym → UBERON ID.

    Obsolete terms and non-UBERON IDs are excluded.  When the same name maps
    to more than one term, the first occurrence wins (synonyms never overwrite
    primary names due to dict.setdefault).
    """
    name_to_id: dict[str, str] = {}
    current_id: str | None = None
    current_name: str | None = None
    is_obsolete = False
    synonyms: list[str] = []

    def _flush() -> None:
        """
        Commit the current OBO term block to name_to_id if it is a valid UBERON term.
        """
        if (
            current_id
            and current_name
            and not is_obsolete
            and current_id.startswith('UBERON:')
        ):
            name_to_id.setdefault(current_name.lower(), current_id)
            for syn in synonyms:
                name_to_id.setdefault(syn.lower(), current_id)

    with obo_path.open(encoding='utf-8') as fh:
        for raw_line in fh:
            line = raw_line.rstrip('\n')
            if line == '[Term]':
                _flush()
                current_id = current_name = None
                is_obsolete = False
                synonyms = []
            elif line.startswith('id: '):
                current_id = line[4:]
            elif line.startswith('name: '):
                current_name = line[6:]
            elif line == 'is_obsolete: true':
                is_obsolete = True
            elif line.startswith('synonym: '):
                m = re.match(r'synonym: "([^"]+)"', line)
                if m:
                    synonyms.append(m.group(1))
    _flush()

    logger.info('UBERON OBO parsed: %d name/synonym entries', len(name_to_id))
    return name_to_id


def _extract_tissue_names(tar_path: Path) -> list[str]:
    """
    Collect unique Anatomical Structure values from all .rpt.gz files in the archive.

    Args:
        tar_path: Path to gxdrnaseq.tar.gz.

    Returns:
        Sorted list of unique tissue name strings.
    """
    tissues: set[str] = set()
    with tarfile.open(tar_path, 'r:gz') as tar:
        members = [m for m in tar.getmembers() if m.name.endswith('.rpt.gz')]
        logger.info('Scanning %d experiment files for tissue names...', len(members))
        for i, member in enumerate(members, 1):
            fobj = tar.extractfile(member)
            if fobj is None:
                continue
            try:
                with gzip.open(io.BytesIO(fobj.read()), 'rb') as gz:
                    raw = gz.read()
                df = pd.read_csv(io.BytesIO(raw), sep='|', dtype=str, low_memory=False)
                df.columns = df.columns.str.strip()
                if 'Anatomical Structure' in df.columns:
                    tissues.update(
                        df['Anatomical Structure'].dropna().str.strip().unique()
                    )
            except Exception as exc:
                logger.warning('Skipping %s: %s', member.name, exc)
            if i % 100 == 0:
                logger.info(
                    '  %d / %d files, %d tissues so far', i, len(members), len(tissues)
                )

    logger.info('Unique tissue names found: %d', len(tissues))
    return sorted(tissues)


def build_uberon_mapping(
    tar_path: Path,
    obo_path: Path,
    out_path: Path,
) -> Path:
    """
    Generate gxd_uberon_mapping.tsv from the GXD archive and basic.obo.

    Matches each unique GXD tissue name (case-insensitive) against UBERON
    primary names and synonyms.  Unmatched tissues are written with an empty
    uberon_id column so they can be curated manually.

    Args:
        tar_path: Path to gxdrnaseq.tar.gz.
        obo_path: Path to uberon basic.obo.
        out_path: Destination for the generated TSV.

    Returns:
        Path to the written TSV.
    """
    name_to_id = parse_uberon_obo(obo_path)
    tissues = _extract_tissue_names(tar_path)

    rows = []
    n_matched = 0
    for tissue in tissues:
        uberon_id = name_to_id.get(tissue.lower(), '')
        if uberon_id:
            n_matched += 1
        rows.append({'gxd_name': tissue, 'uberon_id': uberon_id})

    out_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=['gxd_name', 'uberon_id']).to_csv(
        out_path, sep='\t', index=False
    )

    logger.info(
        'Mapping written: %s  (%d / %d tissues matched)',
        out_path,
        n_matched,
        len(tissues),
    )
    if n_matched < len(tissues):
        logger.warning(
            '%d tissues unmatched - open %s and fill in uberon_id manually',
            len(tissues) - n_matched,
            out_path,
        )
    return out_path
