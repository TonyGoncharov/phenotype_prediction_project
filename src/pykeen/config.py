import os
import socket

GENE_PHENOTYPE_RELATION = "HumanGeneHasMpTopTerm"

CSV_SEP = "\t"

# ── Per-machine configs ───────────────────────────────────────────────────────

_PROT1_CONFIG: dict = {
    "model":            "RotatE",
    "embedding_dim":    128,
    "num_epochs":       200,
    "batch_size":       8192,
    "lr":               3e-3,
    "num_negs_per_pos": 32,
    "train_frac":       0.80,
    "val_frac":         0.10,
    "random_seed":      42,
    "eval_batch_size":  4096,
    "checkpoint_frequency": 50,
}

_MACBOOK_CONFIG: dict = {
    "model":            "RotatE",
    "embedding_dim":    128,
    "num_epochs":       200,
    "batch_size":       2048,
    "lr":               1e-3,
    "num_negs_per_pos": 32,
    "train_frac":       0.80,
    "val_frac":         0.10,
    "random_seed":      42,
    "eval_batch_size":  512,
    "checkpoint_frequency": 50,
}

# ── Named configs for env-var selection ──────────────────────────────────────

_NAMED_CONFIGS: dict[str, dict] = {
    "prot1":   _PROT1_CONFIG,
    "macbook": _MACBOOK_CONFIG,
}

# ── Autodetect ────────────────────────────────────────────────────────────────

def _detect_config() -> dict:
    """Return the config for the current machine.

    Override via PYKEEN_CONFIG env var: set to 'prot1' or 'macbook' to bypass
    hostname detection (useful in Docker, CI, or remote environments).
    """
    env = os.environ.get("PYKEEN_CONFIG", "").strip().lower()
    if env:
        if env not in _NAMED_CONFIGS:
            raise ValueError(
                f"PYKEEN_CONFIG='{env}' is not a valid config name. "
                f"Choose from: {sorted(_NAMED_CONFIGS)}"
            )
        return _NAMED_CONFIGS[env]

    hostname = socket.gethostname()
    if "prot1" in hostname:
        return _PROT1_CONFIG
    if "TonyMacBook" in hostname:
        return _MACBOOK_CONFIG
    return _MACBOOK_CONFIG

DEFAULT_CONFIG: dict = _detect_config()
