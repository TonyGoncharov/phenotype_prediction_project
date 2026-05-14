"""
scripts/train_multilabel.py

Multi-label classifier on R-GCN embeddings.

R-GCN trains WITHOUT phenotype edges (expression + GO + PPI + encodes + uberon).
Phenotype edges (HumanGeneHasMpTopTerm) are used ONLY as labels for the
downstream LogisticRegression — no data leakage.

Usage:
    uv run python scripts/train_multilabel.py --data-dir biocypher_out/human
    uv run python scripts/train_multilabel.py --data-dir biocypher_out/human --skip-rgcn
    uv run python scripts/train_multilabel.py --data-dir biocypher_out/human --fast
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score
from sklearn.multiclass import OneVsRestClassifier
from sklearn.preprocessing import MultiLabelBinarizer
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from src.pykeen.loader import load_triples  # noqa: E402
from src.pykeen.rgcn import RGCNEncoder     # noqa: E402

logger = logging.getLogger(__name__)

# ── Relation names — exact strings from BioCypher *-header.csv :TYPE columns ─
PHENOTYPE_REL = "HumanGeneHasMpTopTerm"        # labels only, NOT in R-GCN graph
ENCODES_REL   = "HumanGeneEncodesProtein"
EXPR_REL      = "HumanGeneExpressedInTissue"
GO_REL        = "HumanGeneHasGoTerm"
PPI_REL       = "HumanProteinInteractsWith"
UBERON_REL    = "TissueMappedToUberon"          # tissue→uberon ontology context

CONTEXT_RELS = {ENCODES_REL, EXPR_REL, GO_REL, PPI_REL, UBERON_REL}

_FAST_CONFIG = {
    "rgcn_epochs": 20,
    "dim":         64,
}


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Train R-GCN encoder + multi-label phenotype classifier.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--data-dir", required=True,
        help="BioCypher output directory (e.g. biocypher_out/human)",
    )
    p.add_argument(
        "--out-dir", default="results/multilabel",
        help="Directory to write results and the saved encoder",
    )
    p.add_argument(
        "--skip-rgcn", action="store_true",
        help="Load existing rgcn_encoder.pt instead of retraining",
    )
    p.add_argument("--rgcn-epochs", type=int, default=100)
    p.add_argument("--dim",         type=int, default=128)
    p.add_argument("--seed",        type=int, default=42)
    p.add_argument(
        "--fast", action="store_true",
        help="Smoke-test preset: dim=64, epochs=20. Completes in ~5 min on CPU.",
    )
    return p.parse_args()


# ── Graph construction ────────────────────────────────────────────────────────

def build_graph(data_dir: str) -> tuple:
    """
    Load BioCypher triples via loader.load_triples() and build the R-GCN graph.

    The phenotype relation is separated out as labels — it is NEVER added to
    the graph, preventing label leakage into the encoder.

    Returns
    -------
    context_df   : DataFrame of context edges fed to R-GCN
    ent2id       : entity string → integer index
    rel2id       : relation string → integer index
    edge_index   : (2, 2E) bidirectional edge tensor
    edge_type    : (2E,) relation-type tensor
    num_rel_inv  : total number of relation types (forward + inverse)
    pheno_df     : DataFrame of HumanGeneHasMpTopTerm edges (labels)
    all_genes    : set of gene IDs that have at least one phenotype annotation
    """
    logger.info("Loading triples from BioCypher output...")
    raw    = load_triples(data_dir)
    all_df = pd.DataFrame(raw, columns=["head", "relation", "tail"])

    # ── Phenotype edges → classification labels only ──────────────────────────
    pheno_df  = all_df[all_df["relation"] == PHENOTYPE_REL].copy()
    all_genes = set(pheno_df["head"].unique())

    logger.info("  Genes with phenotype labels : %d", len(all_genes))
    logger.info("  Gene→Phenotype triples      : %d", len(pheno_df))
    logger.info("  Unique MP top-terms         : %d", pheno_df["tail"].nunique())

    # ── Context graph (no phenotype) ──────────────────────────────────────────
    # 1. Encodes: gene→protein  (restrict to phenotype-annotated genes)
    enc      = all_df[all_df["relation"] == ENCODES_REL]
    enc      = enc[enc["head"].isin(all_genes)]
    proteins = set(enc["tail"].unique())

    # 2. PPI: protein→protein  (restrict to proteins encoded by the above genes)
    ppi = all_df[all_df["relation"] == PPI_REL]
    ppi = ppi[ppi["head"].isin(proteins) | ppi["tail"].isin(proteins)]

    # 3. Expression: gene→tissue
    expr    = all_df[all_df["relation"] == EXPR_REL]
    expr    = expr[expr["head"].isin(all_genes)]
    tissues = set(expr["tail"].unique())

    # 4. GO annotations: gene→go_term
    go = all_df[all_df["relation"] == GO_REL]
    go = go[go["head"].isin(all_genes)]

    # 5. Tissue→Uberon: adds ontology context for tissues
    uberon = all_df[all_df["relation"] == UBERON_REL]
    uberon = uberon[uberon["head"].isin(tissues)]

    context_df = pd.concat(
        [enc, ppi, expr, go, uberon], ignore_index=True
    ).drop_duplicates()

    entities  = sorted(set(context_df["head"]) | set(context_df["tail"]))
    relations = sorted(context_df["relation"].unique())
    ent2id = {e: i for i, e in enumerate(entities)}
    rel2id = {r: i for i, r in enumerate(relations)}

    heads  = torch.tensor([ent2id[h] for h in context_df["head"]],     dtype=torch.long)
    tails  = torch.tensor([ent2id[t] for t in context_df["tail"]],     dtype=torch.long)
    rtypes = torch.tensor([rel2id[r] for r in context_df["relation"]], dtype=torch.long)

    # Bidirectional: forward + inverse (standard R-GCN practice)
    edge_index  = torch.stack([
        torch.cat([heads, tails]),
        torch.cat([tails, heads]),
    ])
    edge_type   = torch.cat([rtypes, rtypes + len(relations)])
    num_rel_inv = len(relations) * 2

    logger.info("R-GCN context graph (NO phenotype edges):")
    logger.info("  Entities         : %d", len(entities))
    logger.info("  Relation types   : %d × 2 (+ inverse) = %d", len(relations), num_rel_inv)
    logger.info("  Context triples  : %d", len(context_df))
    logger.info("    Encodes        : %d", len(enc))
    logger.info("    PPI            : %d", len(ppi))
    logger.info("    Expression     : %d", len(expr))
    logger.info("    GO             : %d", len(go))
    logger.info("    Uberon         : %d", len(uberon))

    return (
        context_df, ent2id, rel2id,
        edge_index, edge_type, num_rel_inv,
        pheno_df, all_genes,
    )


# ── R-GCN training ────────────────────────────────────────────────────────────

def train_rgcn(
    edge_index:    torch.Tensor,
    edge_type:     torch.Tensor,
    num_entities:  int,
    num_relations: int,
    dim:           int,
    epochs:        int,
    device:        torch.device,
) -> RGCNEncoder:
    """
    Train the R-GCN encoder via link-reconstruction on context edges.

    Uses binary cross-entropy with negative sampling (5 negatives per positive,
    batch size 4096). This is the same self-supervised objective as in the
    original R-GCN link prediction paper.
    """
    model     = RGCNEncoder(num_entities, num_relations, dim).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)

    ei    = edge_index.to(device)
    et    = edge_type.to(device)
    pos_h = ei[0]
    pos_t = ei[1]

    logger.info("Training R-GCN for %d epochs, dim=%d ...", epochs, dim)
    for epoch in range(1, epochs + 1):
        model.train()
        optimizer.zero_grad()

        z = model(ei, et)

        # Random batch of positive edges
        idx = torch.randperm(len(pos_h), device=device)[:4096]
        ph, pt = pos_h[idx], pos_t[idx]

        pos_scores = model.decode(z, ph, pt)

        # Corrupt tail (5× negatives per positive)
        neg_t      = torch.randint(0, num_entities, (len(ph) * 5,), device=device)
        neg_h      = ph.repeat_interleave(5)
        neg_scores = model.decode(z, neg_h, neg_t)

        loss = (
            F.binary_cross_entropy_with_logits(
                pos_scores, torch.ones_like(pos_scores)
            )
            + F.binary_cross_entropy_with_logits(
                neg_scores, torch.zeros_like(neg_scores)
            )
        )
        loss.backward()
        optimizer.step()

        if epoch % 20 == 0:
            logger.info("  Epoch %3d/%d  loss=%.4f", epoch, epochs, loss.item())

    return model


def get_embeddings(
    model:      RGCNEncoder,
    edge_index: torch.Tensor,
    edge_type:  torch.Tensor,
    device:     torch.device,
) -> np.ndarray:
    model.eval()
    with torch.no_grad():
        z = model(edge_index.to(device), edge_type.to(device))
    return z.cpu().numpy()


# ── Classifier ────────────────────────────────────────────────────────────────

def train_classifier(
    embeddings: np.ndarray,
    pheno_df:   pd.DataFrame,
    ent2id:     dict[str, int],
    all_genes:  set[str],
    seed:       int = 42,
) -> tuple[float, pd.DataFrame, np.ndarray]:
    """
    Train a OneVsRest LogisticRegression on gene embeddings.

    Returns mean AUPRC, per-phenotype metrics DataFrame, and class labels.
    Genes absent from the R-GCN graph (ent2id) are silently skipped —
    this can happen if a gene has no context edges after filtering.
    """
    gene_to_pheno: dict[str, list[str]] = {}
    for _, row in pheno_df.iterrows():
        gene_to_pheno.setdefault(row["head"], []).append(row["tail"])

    genes  = [g for g in all_genes if g in ent2id]
    n_skip = len(all_genes) - len(genes)
    if n_skip:
        logger.warning("%d genes skipped (no context edges in R-GCN graph)", n_skip)

    X       = np.array([embeddings[ent2id[g]] for g in genes])
    y_lists = [gene_to_pheno.get(g, [])       for g in genes]

    mlb = MultiLabelBinarizer()
    Y   = mlb.fit_transform(y_lists)

    rng     = np.random.default_rng(seed)
    idx     = np.arange(len(genes))
    rng.shuffle(idx)
    n_train = int(len(idx) * 0.8)

    X_train, X_test = X[idx[:n_train]], X[idx[n_train:]]
    Y_train, Y_test = Y[idx[:n_train]], Y[idx[n_train:]]

    logger.info("  Train / test genes : %d / %d", len(X_train), len(X_test))
    logger.info("  Phenotype classes  : %d", len(mlb.classes_))

    clf = OneVsRestClassifier(
        LogisticRegression(max_iter=1000, random_state=seed, C=1.0),
        n_jobs=-1,
    )
    clf.fit(X_train, Y_train)
    Y_prob = clf.predict_proba(X_test)

    rows = []
    for i, mp_term in enumerate(mlb.classes_):
        y_true  = Y_test[:, i]
        y_score = Y_prob[:, i]
        if y_true.sum() == 0:
            continue          # skip classes with no positives in test split
        rows.append({
            "mp_term_id":  mp_term,
            "ap":          round(float(average_precision_score(y_true, y_score)), 4),
            "n_pos_test":  int(y_true.sum()),
            "n_test":      int(len(y_true)),
            "baseline":    round(float(y_true.mean()), 4),
        })

    per_class  = (
        pd.DataFrame(rows)
        .sort_values("ap", ascending=False)
        .reset_index(drop=True)
    )
    mean_auprc = float(per_class["ap"].mean())
    return mean_auprc, per_class, mlb.classes_


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    args = parse_args()
    if args.fast:
        args.dim         = _FAST_CONFIG["dim"]
        args.rgcn_epochs = _FAST_CONFIG["rgcn_epochs"]
        logger.info("--fast preset: %s", _FAST_CONFIG)

    device  = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Device      : %s", device)
    logger.info("Data dir    : %s", args.data_dir)
    logger.info("Output dir  : %s", out_dir)
    logger.info("R-GCN trains WITHOUT phenotype edges → leakage-free evaluation")

    (context_df, ent2id, rel2id, edge_index, edge_type,
     num_rel, pheno_df, all_genes) = build_graph(args.data_dir)

    num_entities = len(ent2id)
    rgcn_path    = out_dir / "rgcn_encoder.pt"

    if args.skip_rgcn and rgcn_path.exists():
        logger.info("Loading R-GCN encoder from %s ...", rgcn_path)
        rgcn = RGCNEncoder(num_entities, num_rel, args.dim).to(device)
        rgcn.load_state_dict(
            torch.load(rgcn_path, map_location=device, weights_only=True)
        )
    else:
        logger.info("Training R-GCN (%d epochs, dim=%d) ...", args.rgcn_epochs, args.dim)
        rgcn = train_rgcn(
            edge_index, edge_type,
            num_entities, num_rel,
            args.dim, args.rgcn_epochs, device,
        )
        torch.save(rgcn.state_dict(), rgcn_path)
        logger.info("Encoder saved: %s", rgcn_path)

    logger.info("Extracting embeddings ...")
    embeddings = get_embeddings(rgcn, edge_index, edge_type, device)
    logger.info("  Shape: %s", embeddings.shape)

    logger.info("Training multi-label classifier ...")
    mean_auprc, per_class, classes = train_classifier(
        embeddings, pheno_df, ent2id, all_genes, seed=args.seed,
    )

    per_class.to_csv(out_dir / "auprc_per_phenotype.csv", index=False)

    mean_baseline = float(per_class["baseline"].mean())
    improvement   = mean_auprc / mean_baseline if mean_baseline > 0 else 0.0

    logger.info("=" * 52)
    logger.info("  Mean AUPRC            : %.4f", mean_auprc)
    logger.info("  Random baseline       : %.4f", mean_baseline)
    logger.info("  Improvement vs random : %.2f×", improvement)
    logger.info("  Top-10 phenotypes by AUPRC:\n%s", per_class.head(10).to_string(index=False))
    logger.info("=" * 52)

    n_genes_in_graph = sum(1 for g in all_genes if g in ent2id)
    summary = {
        "method":          "R-GCN (leakage-free) + LogisticRegression",
        "rgcn_dim":        args.dim,
        "rgcn_epochs":     args.rgcn_epochs,
        "leakage_free":    True,
        "mean_auprc":      round(mean_auprc, 4),
        "mean_baseline":   round(mean_baseline, 4),
        "improvement_x":   round(improvement, 3),
        "n_genes":         n_genes_in_graph,
        "n_genes_skipped": len(all_genes) - n_genes_in_graph,
        "n_phenotypes":    len(classes),
        "n_entities":      num_entities,
        "n_context_edges": len(context_df),
    }
    pd.DataFrame([summary]).to_csv(out_dir / "summary.csv", index=False)
    logger.info("Saved to %s/", out_dir)


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-8s  %(name)s  %(message)s",
        datefmt="%H:%M:%S",
    )
    main()
