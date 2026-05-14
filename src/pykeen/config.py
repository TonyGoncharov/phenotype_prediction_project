import logging
import os

GENE_PHENOTYPE_RELATION = "HumanGeneHasMpTopTerm"

CSV_SEP = "\t"

logger = logging.getLogger(__name__)

# ── Default config (laptop / single-GPU workstation) ─────────────────────────

_DEFAULT_CONFIG: dict = {
    "model":                "RotatE",
    "embedding_dim":        128,
    "batch_size":           2048,
    "lr":                   1e-3,
    "num_negs_per_pos":     64,
    "train_frac":           0.80,
    "val_frac":             0.10,
    "random_seed":          42,
    "eval_batch_size":      128,
    "checkpoint_frequency": 25,
    "target_steps":         80_000,
    "es_patience":          10,
    "es_frequency":         10,
    "es_delta":             2e-3,
}

# ── Named overrides selectable via PYKEEN_CONFIG ──────────────────────────────

_SERVER_CONFIG: dict = {
    "model":                "RotatE",
    "embedding_dim":        128,
    "num_epochs":           200,
    "batch_size":           8192,
    "lr":                   3e-3,
    "num_negs_per_pos":     32,
    "train_frac":           0.80,
    "val_frac":             0.10,
    "random_seed":          42,
    "eval_batch_size":      4096,
    "checkpoint_frequency": 50,
}

_NAMED_CONFIGS: dict[str, dict] = {
    "default": _DEFAULT_CONFIG,
    "server":  _SERVER_CONFIG,
}


def _load_config() -> dict:
    """Return the active config.

    Set PYKEEN_CONFIG=server to use high-throughput server settings.
    Omit or set to 'default' for the standard laptop/workstation config.
    """
    env = os.environ.get("PYKEEN_CONFIG", "").strip().lower()
    if not env:
        return _DEFAULT_CONFIG
    if env not in _NAMED_CONFIGS:
        raise ValueError(
            f"PYKEEN_CONFIG='{env}' is not a valid config name. "
            f"Choose from: {sorted(_NAMED_CONFIGS)}"
        )
    logger.debug("Using config '%s' (PYKEEN_CONFIG)", env)
    return _NAMED_CONFIGS[env]


DEFAULT_CONFIG: dict = _load_config()
