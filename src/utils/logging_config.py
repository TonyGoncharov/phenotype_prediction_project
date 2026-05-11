"""src/utils/logging_config.py — centralised logging setup for the KG pipeline.

Call setup_logging(out_dir=...) once before any export step.
Console: INFO.  File (<out_dir>/pipeline.log): DEBUG.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

# Module-level sentinel so setup_logging() is idempotent.
_LOGGING_CONFIGURED = False

PIPELINE_LOGGER_NAME = "kg_pipeline"


def setup_logging(
    out_dir: str | Path = ".",
    log_filename: str = "pipeline.log",
    console_level: int = logging.INFO,
    file_level: int = logging.DEBUG,
) -> logging.Logger:
    """Configure root logger with console and file handlers.

    Idempotent — subsequent calls are no-ops and return the already-configured logger.
    out_dir is created if absent; log file written to <out_dir>/pipeline.log.
    """
    global _LOGGING_CONFIGURED
    if _LOGGING_CONFIGURED:
        return logging.getLogger(PIPELINE_LOGGER_NAME)

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = out_dir / log_filename

    # ── Formatters ──────────────────────────────────────────────────── #

    _DETAILED_FMT = "%(asctime)s | %(levelname)-8s | %(name)s | %(message)s"
    _CONSOLE_FMT  = "%(asctime)s | %(levelname)-8s | %(message)s"
    _DATE_FMT     = "%Y-%m-%d %H:%M:%S"

    file_formatter    = logging.Formatter(_DETAILED_FMT, datefmt=_DATE_FMT)
    console_formatter = logging.Formatter(_CONSOLE_FMT,  datefmt=_DATE_FMT)

    # ── Handlers ────────────────────────────────────────────────────── #

    file_handler = logging.FileHandler(log_path, mode="a", encoding="utf-8")
    file_handler.setLevel(file_level)
    file_handler.setFormatter(file_formatter)

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(console_level)
    console_handler.setFormatter(console_formatter)

    # ── Root logger ─────────────────────────────────────────────────── #
    # Set to DEBUG so that handlers can filter independently.

    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    root.addHandler(file_handler)
    root.addHandler(console_handler)

    # Silence overly verbose third-party loggers at WARNING level.
    for noisy in ("neo4j", "urllib3", "httpx"):
        logging.getLogger(noisy).setLevel(logging.WARNING)

    # Route BioCypher WARNING+ into pipeline.log regardless of import order
    # (BioCypher attaches its own FileHandler on first import before our root
    # handler exists, so we add file_handler directly rather than relying on propagate).
    bc_logger = logging.getLogger("biocypher")
    bc_logger.setLevel(logging.WARNING)
    bc_logger.addHandler(file_handler)

    _LOGGING_CONFIGURED = True

    logger = logging.getLogger(PIPELINE_LOGGER_NAME)
    logger.info("Logging initialised — file: %s", log_path)
    return logger