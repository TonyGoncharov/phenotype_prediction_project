"""src/utils/logging_config.py — centralised logging setup for the KG pipeline.

Usage (call once, as early as possible — in build_graph.py or a CLI entry point):

    from src.utils.logging_config import setup_logging
    setup_logging(out_dir="/path/to/output")

Every other module then simply does:

    import logging
    logger = logging.getLogger(__name__)

and uses logger.info / logger.debug / logger.warning / logger.error.

Log levels
----------
  Console  : INFO  — progress, counts, warnings
  File     : DEBUG — everything above + internal details useful for debugging

Log file location
-----------------
  <out_dir>/pipeline.log

  When building both species (species="both"), build_graph.py passes the
  shared parent out_dir so both species share one log file for the run.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

# Module-level sentinel so setup_logging() is idempotent.
_LOGGING_CONFIGURED = False

# The logger that build_graph.py and the pipeline itself use.
PIPELINE_LOGGER_NAME = "kg_pipeline"


def setup_logging(
    out_dir: str | Path = ".",
    log_filename: str = "pipeline.log",
    console_level: int = logging.INFO,
    file_level: int = logging.DEBUG,
) -> logging.Logger:
    """Configure root logger with a console handler and a rotating file handler.

    Safe to call multiple times — subsequent calls are no-ops and return the
    already-configured pipeline logger.

    Args:
        out_dir:       Directory where the log file is written.
                       Created if it does not exist.
        log_filename:  Name of the log file (default: pipeline.log).
        console_level: Minimum level for console output (default: INFO).
        file_level:    Minimum level for file output (default: DEBUG).

    Returns:
        The pipeline logger (``logging.getLogger("kg_pipeline")``).
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

    # ── BioCypher → our file ─────────────────────────────────────────── #
    # BioCypher creates its own FileHandler on first import, writing to
    # biocypher-log/biocypher-<timestamp>.log.  We additionally route its
    # WARNING+ messages into our pipeline.log so both sources are visible
    # in one place.  propagate=True (the default) means root already sees
    # these records, but the root handler may not have been attached yet
    # when BioCypher initialised.  Adding file_handler directly guarantees
    # delivery regardless of import order.
    bc_logger = logging.getLogger("biocypher")
    bc_logger.setLevel(logging.WARNING)
    bc_logger.addHandler(file_handler)

    _LOGGING_CONFIGURED = True

    logger = logging.getLogger(PIPELINE_LOGGER_NAME)
    logger.info("Logging initialised — file: %s", log_path)
    return logger