GENE_PHENOTYPE_RELATION = "HumanGeneHasMpTopTerm"

# ── BioCypher CSV format ──────────────────────────────────────────────────────
# Configured in biocypher_config.yaml:  neo4j.delimiter = '\t'
CSV_SEP = "\t"

# ── Model and training defaults ───────────────────────────────────────────────
DEFAULT_CONFIG: dict = {
    # Model
    "model": "RotatE",
    # RotatE requires embedding_dim to be EVEN (complex-space rotation).
    # 256 gives a good trade-off between capacity and memory.
    "embedding_dim": 256,

    # Training
    "num_epochs": 300,
    "batch_size": 512,
    "lr": 1e-3,

    # Negative sampling (sLCWA — stochastic Local Closed World Assumption)
    "num_negs_per_pos": 64,

    # Train / val / test split ratios (must sum to 1.0)
    "train_frac": 0.80,
    "val_frac":   0.10,
    # test_frac is derived: 1 - train_frac - val_frac

    # Reproducibility
    "random_seed": 42,

    # Evaluation batch size (separate from training batch size)
    "eval_batch_size": 256,

    # Checkpoint: save every N epochs (0 = disabled)
    "checkpoint_frequency": 50,
}
