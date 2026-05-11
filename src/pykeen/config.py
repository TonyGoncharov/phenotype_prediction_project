GENE_PHENOTYPE_RELATION = "HumanGeneHasMpTopTerm"

# ── BioCypher CSV format ──────────────────────────────────────────────────────
# Configured in biocypher_config.yaml:  neo4j.delimiter = '\t'
CSV_SEP = "\t"

# ── Model and training defaults ───────────────────────────────────────────────
DEFAULT_CONFIG: dict = {
    "model": "RotatE",
    # RotatE requires embedding_dim to be EVEN (complex-space rotation).
    "embedding_dim": 256,

    "num_epochs": 300,
    "batch_size": 512,
    "lr": 1e-3,

    # sLCWA — stochastic Local Closed World Assumption
    "num_negs_per_pos": 64,

    # must sum to 1.0; test_frac is derived: 1 - train_frac - val_frac
    "train_frac": 0.80,
    "val_frac":   0.10,

    "random_seed": 42,

    # separate from training batch size
    "eval_batch_size": 256,

    # 0 = disabled
    "checkpoint_frequency": 50,
}
