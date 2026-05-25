"""
src/utils/logging_config.py - shared logging setup for the pipeline.

Configures root-level DEBUG file logging and INFO (or configurable)
console logging with a consistent format.
"""

from __future__ import annotations

import logging
from pathlib import Path


def setup_logging(
    out_dir: str | Path = 'biocypher-log',
    log_filename: str = 'pipeline.log',
    console_level: int = logging.INFO,
) -> None:
    """
    Configure root logger with a file handler (DEBUG) and a console handler.

    Args:
        out_dir: Directory for the log file; created if absent.
        log_filename: Log file name within out_dir.
        console_level: Logging level for the console handler.
    """
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)s [%(name)s] %(message)s',
        handlers=[
            logging.FileHandler(out_path / log_filename),
            logging.StreamHandler(),
        ],
        force=True,
    )
    for handler in logging.getLogger().handlers:
        if isinstance(handler, logging.StreamHandler) and not isinstance(
            handler, logging.FileHandler
        ):
            handler.setLevel(console_level)
    logging.getLogger('chardet').setLevel(logging.WARNING)
