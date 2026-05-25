"""
src/machine_learning/common/results.py - structured result persistence for both ML tracks.

Both KGE and R-GCN write outputs under results/ml/{species}/{track}/{condition}/.
This module standardises summary JSON, per-class CSV, and comparison CSV saving.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import pandas as pd
from src.machine_learning.common.metrics import (
    STANDARD_DELTA_METRICS,
    compute_delta_metrics,
)

logger = logging.getLogger(__name__)


def save_json(data: dict, path: Path | str) -> Path:
    """
    Write *data* as indented JSON to *path*, creating parent directories as needed.

    Args:
        data: Serialisable dict.
        path: Destination file path.

    Returns:
        Resolved Path of the written file.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    # default=str handles numpy scalar types (float32, int64) that json cannot
    # serialise natively.  Any other unexpected type will also be coerced to its
    # str() representation - callers should pass only dicts of Python primitives,
    # numpy scalars, or lists thereof.
    path.write_text(json.dumps(data, indent=2, default=str))
    logger.info('Saved JSON → %s', path)
    return path


def patch_json(path: Path | str, updates: dict) -> None:
    """
    Merge *updates* into an existing JSON file (in-place).

    If the file does not exist, it is created with *updates* as its content.

    Args:
        path: Path to the JSON file.
        updates: Key-value pairs to merge into the existing content.
    """
    path = Path(path)
    if not path.exists():
        save_json(updates, path)
        return
    data = json.loads(path.read_text())
    data.update(updates)
    path.write_text(json.dumps(data, indent=2, default=str))


def save_csv(df: pd.DataFrame, path: Path | str) -> Path:
    """
    Save *df* to CSV, creating parent directories as needed.

    Args:
        df: DataFrame to save.
        path: Destination file path.

    Returns:
        Resolved Path of the written file.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    logger.info('Saved CSV → %s', path)
    return path


def comparison_dataframe(
    rows: list[dict] | pd.DataFrame,
    baseline_condition: str = 'all',
    metric_keys: list[str] | None = None,
    group_keys: tuple[str, ...] = ('species',),
) -> pd.DataFrame:
    """
    Build an ablation comparison DataFrame with delta metrics vs baseline.

    Delta columns are named delta_<metric>_vs_<baseline_condition> and are
    computed only for requested metric columns that are present in the results.
    If species is present, each species gets its own baseline row.
    """
    df = rows.copy() if isinstance(rows, pd.DataFrame) else pd.DataFrame(rows)
    if df.empty or 'condition' not in df.columns:
        return df

    metric_keys = metric_keys or STANDARD_DELTA_METRICS
    present_metrics = [key for key in metric_keys if key in df.columns]
    if not present_metrics:
        return df

    delta_cols = [f'delta_{key}_vs_{baseline_condition}' for key in present_metrics]
    for col in delta_cols:
        if col not in df.columns:
            df[col] = None

    group_cols = [key for key in group_keys if key in df.columns]
    groups = (
        df.groupby(group_cols, dropna=False, sort=False) if group_cols else [((), df)]
    )

    found_baseline = False
    for _, group in groups:
        baseline_rows = group[group['condition'] == baseline_condition]
        if baseline_rows.empty:
            logger.warning(
                "Baseline condition '%s' not found for comparison group %s - "
                'delta metrics skipped for this group.',
                baseline_condition,
                group_cols or '<all>',
            )
            continue

        found_baseline = True
        baseline = baseline_rows.iloc[0].to_dict()
        for idx, row in group.iterrows():
            deltas = compute_delta_metrics(
                baseline,
                row.to_dict(),
                metric_keys=present_metrics,
                suffix=f'_vs_{baseline_condition}',
            )
            for key, value in deltas.items():
                df.at[idx, key] = value

    if not found_baseline:
        return df.drop(columns=delta_cols)

    return df


def save_comparison_csv(
    rows: list[dict] | pd.DataFrame,
    out_dir: Path | str,
    filename: str = 'comparison_summary.csv',
) -> Path:
    """
    Save ablation comparison rows to a CSV and return the path.

    Args:
        rows: List of per-condition result dicts.
        out_dir: Output directory.
        filename: Output file name.

    Returns:
        Path to the written CSV.
    """
    return save_csv(comparison_dataframe(rows), Path(out_dir) / filename)
