from __future__ import annotations

import sys
import urllib.request
from pathlib import Path

ZENODO_RECORD = '20379520'
_PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = _PROJECT_ROOT / 'data' / 'raw'

FILES = [
    'hp.obo',
    'mp.obo',
    'basic.obo',
    'genes_to_phenotype.txt',
    'HP2MP.tsv',
    'MGI_PhenoGenoMP.rpt',
    'gene2go.gz',
    'gene_info.gz',
    'GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz',
    'gxdrnaseq.tar.gz',
    'gxd_uberon_mapping.tsv',
    'BIOGRID-ALL-5.0.256.tab3.txt',
]


def _progress(filename: str):
    last = [-1]

    def hook(block: int, block_size: int, total: int) -> None:
        if total <= 0:
            return
        pct = min(100, block * block_size * 100 // total)
        if pct != last[0]:
            last[0] = pct
            bar = '█' * (pct // 5) + '░' * (20 - pct // 5)
            print(f'\r  [{bar}] {pct:3d}%  {filename}', end='', flush=True)
        if pct == 100:
            print()

    return hook


def main() -> None:
    if ZENODO_RECORD == 'TODO':
        print('Error: set ZENODO_RECORD to your Zenodo record ID before running.')
        sys.exit(1)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    base_url = f'https://zenodo.org/record/{ZENODO_RECORD}/files'

    failed: list[str] = []
    for filename in FILES:
        dest = OUT_DIR / filename
        if dest.exists():
            size_mb = dest.stat().st_size / 1e6
            print(f'  ✓  {filename} ({size_mb:.1f} MB)')
            continue

        url = f'{base_url}/{filename}?download=1'
        print(f'  ↓  {filename}')
        tmp = dest.with_suffix(dest.suffix + '.part')
        try:
            urllib.request.urlretrieve(url, tmp, reporthook=_progress(filename))
            tmp.rename(dest)
        except Exception as exc:
            tmp.unlink(missing_ok=True)
            print(f'\n  ✗  {filename}: {exc}')
            failed.append(filename)

    print()
    if failed:
        print(f'Failed ({len(failed)}): {", ".join(failed)}')
        sys.exit(1)
    else:
        print(f'All {len(FILES)} files ready in {OUT_DIR.resolve()}')
        print('Run the pipeline: uv run python run_kg.py')


if __name__ == '__main__':
    main()
