"""
scripts/ablation_rgcn.py

Ablation study for the R-GCN multi-label phenotype classifier.

Trains one R-GCN encoder per condition, removing one context layer at a time.
Phenotype edges are NEVER in the R-GCN graph — always labels only.

Conditions
----------
  all          : encodes + expression + GO + PPI + uberon  (full context)
  no_go        : encodes + expression + PPI + uberon
  no_ppi       : encodes + expression + GO + uberon
  no_expr      : encodes + GO + PPI
  encodes_only : gene→protein only  (minimal baseline)

Test-set integrity
------------------
A fixed 80/20 gene split is derived ONCE before the ablation loop (seed=42).
Every condition evaluates on exactly the same held-out genes. Genes absent from
a condition's graph (no context edges after filtering) receive zero embeddings —
they are NOT dropped, ensuring an honest penalty for missing information.

MultiLabelBinarizer is fitted once on the full gene set so the class list is
identical across all conditions.

Parallel execution
------------------
Add --parallel to run all conditions simultaneously.  Each worker process is
launched with multiprocessing start_method="spawn" so that OMP_NUM_THREADS
and related env vars are set BEFORE PyTorch/OpenMP are imported — the only
correct approach (torch.set_num_threads() cannot retract the OpenMP pool after
first import).  Thread budget is divided automatically:
    threads_per_worker = max(1, cpu_count // n_workers)
Override with --threads-per-worker.

Usage
-----
    # Sequential
    uv run python scripts/ablation_rgcn.py --data-dir biocypher_out/human

    # Parallel (all conditions at once)
    uv run python scripts/ablation_rgcn.py --data-dir biocypher_out/human --parallel

    # Parallel, explicit worker count and thread budget
    uv run python scripts/ablation_rgcn.py --data-dir biocypher_out/human \\
        --parallel --max-workers 3 --threads-per-worker 42

    # Quick smoke test
    uv run python scripts/ablation_rgcn.py --data-dir biocypher_out/human --fast

    # Single condition
    uv run python scripts/ablation_rgcn.py --data-dir biocypher_out/human \\
        --conditions all no_ppi
"""

from __future__ import annotations

import argparse
import json
import multiprocessing
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

# ── Thread environment variables — must be set before torch import ─────────────
_THREAD_ENV_VARS = [
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "PYTORCH_NUM_THREADS",
]


def _set_thread_env(n_threads: int) -> None:
    """ProcessPoolExecutor initializer: set all thread caps before torch import."""
    val = str(n_threads)
    for var in _THREAD_ENV_VARS:
        os.environ[var] = val


# ── Deferred heavy imports (only after thread env is set in worker) ────────────

def _import_torch():
    import torch
    import torch.nn.functional as F
    from torch_geometric.nn import RGCNConv
    return torch, F, RGCNConv


# ── Relation names (exact BioCypher :TYPE strings) ────────────────────────────

PHENOTYPE_REL = "HumanGeneHasMpTopTerm"    # labels only — NEVER in R-GCN graph
ENCODES_REL   = "HumanGeneEncodesProtein"
EXPR_REL      = "HumanGeneExpressedInTissue"
GO_REL        = "HumanGeneHasGoTerm"
PPI_REL       = "HumanProteinInteractsWith"
UBERON_REL    = "TissueMappedToUberon"

_ALL_CONTEXT = {ENCODES_REL, EXPR_REL, GO_REL, PPI_REL, UBERON_REL}

# ── Ablation conditions ───────────────────────────────────────────────────────

_ENCODES = {ENCODES_REL}
_GO      = {GO_REL}
_EXPR    = {EXPR_REL, UBERON_REL}
_PPI     = {PPI_REL}

CONDITIONS: dict[str, set[str] | None] = {
    "all":          None,                     # None → _ALL_CONTEXT at runtime
    "no_go":        _ENCODES | _EXPR | _PPI,
    "no_ppi":       _ENCODES | _GO   | _EXPR,
    "no_expr":      _ENCODES | _GO   | _PPI,
    "encodes_only": _ENCODES,
}

_FAST_OVERRIDE = {"rgcn_epochs": 20, "dim": 64}


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="R-GCN ablation study — layer-removal experiment.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--data-dir", required=True,
                   help="BioCypher output directory (biocypher_out/human)")
    p.add_argument("--out-dir", default="results/ablation_rgcn",
                   help="Root output dir; each condition gets a subdirectory")
    p.add_argument("--fast", action="store_true",
                   help="Smoke test: 20 epochs, dim=64")
    p.add_argument("--conditions", nargs="+", default=None,
                   choices=list(CONDITIONS.keys()),
                   help="Subset of conditions to run (default: all)")
    p.add_argument("--rgcn-epochs", type=int, default=100)
    p.add_argument("--dim",         type=int, default=128)
    p.add_argument("--seed",        type=int, default=42)

    # ── Parallel options ──────────────────────────────────────────────────────
    p.add_argument("--parallel", action="store_true",
                   help="Run conditions in parallel (spawn-based subprocesses)")
    p.add_argument("--max-workers", type=int, default=None,
                   help="Max parallel workers. Default: number of selected conditions.")
    p.add_argument("--threads-per-worker", type=int, default=None,
                   help="OMP_NUM_THREADS per worker. Default: cpu_count // n_workers.")
    return p.parse_args()


# ── Data helpers (called in main process to prepare shared inputs) ─────────────

def _load_all_df(data_dir: str) -> pd.DataFrame:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from src.pykeen.loader import load_triples
    raw = load_triples(data_dir)
    return pd.DataFrame(raw, columns=["head", "relation", "tail"])


def _get_pheno_data(
    all_df: pd.DataFrame,
) -> tuple[pd.DataFrame, set[str]]:
    pheno_df  = all_df[all_df["relation"] == PHENOTYPE_REL].copy()
    all_genes = set(pheno_df["head"].unique())
    return pheno_df, all_genes


def _make_fixed_gene_split(
    all_genes: set[str],
    train_frac: float = 0.8,
    seed: int = 42,
) -> tuple[list[str], list[str]]:
    genes   = sorted(all_genes)
    rng     = np.random.default_rng(seed)
    idx     = np.arange(len(genes))
    rng.shuffle(idx)
    n_train = int(len(idx) * train_frac)
    return [genes[i] for i in idx[:n_train]], [genes[i] for i in idx[n_train:]]



# ── Single-condition worker (runs in subprocess) ───────────────────────────────

def _run_condition(
    cond_name:    str,
    include_rels: set[str] | None,
    all_df_dict:  list[list],          # list of [head, relation, tail] rows
    gene_to_pheno: dict[str, list[str]],
    train_genes:  list[str],
    test_genes:   list[str],
    mlb_classes:  list[str],           # class list (not the fitted object — not picklable)
    all_genes:    set[str],
    out_dir:      str,
    rgcn_epochs:  int,
    dim:          int,
    seed:         int,
    log_file:     str,
    n_jobs:       int = -1,
) -> dict:
    """
    Train R-GCN for one condition and return a result dict.

    This function runs inside a worker subprocess (spawn start method).
    Thread env vars are already set by _set_thread_env() initializer before
    this function is called.  We import torch here — after env vars are set.
    """
    import logging

    # ── Per-condition log file ────────────────────────────────────────────────
    # Use an explicit FileHandler on the named logger instead of basicConfig so
    # that sequential calls in the same process each get their own log file.
    # (basicConfig is a no-op if the root logger already has handlers.)
    Path(log_file).parent.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.INFO)
    fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    log = logging.getLogger(cond_name)
    log.setLevel(logging.INFO)
    log.addHandler(fh)
    log.propagate = False

    # Redirect print() to the log file as well
    import builtins
    _orig_print = builtins.print
    def _log_print(*args, **kwargs):
        msg = " ".join(str(a) for a in args)
        log.info(msg)
        _orig_print(f"[{cond_name}]", *args, **kwargs)
    builtins.print = _log_print

    t0 = time.time()
    log.info("Condition '%s' started (pid=%d)", cond_name, os.getpid())
    log.info("OMP_NUM_THREADS=%s", os.environ.get("OMP_NUM_THREADS", "unset"))

    try:
        torch, F, RGCNConv = _import_torch()
        from sklearn.linear_model import LogisticRegression
        from sklearn.metrics import average_precision_score
        from sklearn.multiclass import OneVsRestClassifier
        from sklearn.preprocessing import MultiLabelBinarizer

        all_df = pd.DataFrame(all_df_dict, columns=["head", "relation", "tail"])

        active_rels = include_rels if include_rels is not None else _ALL_CONTEXT
        ctx = all_df[all_df["relation"].isin(active_rels)]

        # ── Build context graph ───────────────────────────────────────────────
        enc      = ctx[ctx["relation"] == ENCODES_REL]
        enc      = enc[enc["head"].isin(all_genes)]
        proteins = set(enc["tail"].unique())

        ppi     = ctx[ctx["relation"] == PPI_REL]
        ppi     = ppi[ppi["head"].isin(proteins) | ppi["tail"].isin(proteins)]

        expr    = ctx[ctx["relation"] == EXPR_REL]
        expr    = expr[expr["head"].isin(all_genes)]
        tissues = set(expr["tail"].unique())

        go     = ctx[ctx["relation"] == GO_REL]
        go     = go[go["head"].isin(all_genes)]

        uberon = ctx[ctx["relation"] == UBERON_REL]
        uberon = uberon[uberon["head"].isin(tissues)]

        context_df = pd.concat(
            [enc, ppi, expr, go, uberon], ignore_index=True
        ).drop_duplicates()

        if context_df.empty:
            raise RuntimeError(f"No context edges for condition '{cond_name}'")

        entities  = sorted(set(context_df["head"]) | set(context_df["tail"]))
        relations = sorted(context_df["relation"].unique())
        ent2id    = {e: i for i, e in enumerate(entities)}
        rel2id    = {r: i for i, r in enumerate(relations)}

        heads  = torch.tensor([ent2id[h] for h in context_df["head"]],     dtype=torch.long)
        tails  = torch.tensor([ent2id[t] for t in context_df["tail"]],     dtype=torch.long)
        rtypes = torch.tensor([rel2id[r] for r in context_df["relation"]], dtype=torch.long)

        edge_index  = torch.stack([torch.cat([heads, tails]), torch.cat([tails, heads])])
        edge_type   = torch.cat([rtypes, rtypes + len(relations)])
        num_rel_inv = len(relations) * 2
        num_entities = len(entities)

        n_missing = sum(1 for g in all_genes if g not in ent2id)
        log.info(
            "Graph: %d entities | %d rel types | %d triples | %d genes missing",
            num_entities, len(relations), len(context_df), n_missing,
        )
        print(f"  Graph: {num_entities:,} entities | "
              f"{len(relations)} rel types | "
              f"{len(context_df):,} triples | "
              f"{n_missing} genes missing")

        # ── R-GCN encoder ─────────────────────────────────────────────────────

        # num_entities, num_rel_inv, dim captured from enclosing scope
        class RGCNEncoder(torch.nn.Module):
            def __init__(self):
                super().__init__()
                self.entity_emb = torch.nn.Embedding(num_entities, dim)
                self.conv1 = RGCNConv(dim, dim, num_rel_inv)
                self.conv2 = RGCNConv(dim, dim, num_rel_inv)
                torch.nn.init.xavier_uniform_(self.entity_emb.weight)

            def forward(self, ei, et):
                x = self.entity_emb.weight
                x = F.relu(self.conv1(x, ei, et))
                x = self.conv2(x, ei, et)
                return x

            def decode(self, z, h, t):
                return (z[h] * z[t]).sum(dim=-1)

        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        model  = RGCNEncoder().to(device)
        optim  = torch.optim.Adam(model.parameters(), lr=1e-3)

        ei = edge_index.to(device)
        et = edge_type.to(device)
        ph_all = ei[0]
        pt_all = ei[1]

        print(f"  Training R-GCN: {rgcn_epochs} epochs, dim={dim} ...")
        for epoch in range(1, rgcn_epochs + 1):
            model.train()
            optim.zero_grad()

            z   = model(ei, et)
            idx = torch.randperm(len(ph_all), device=device)[:4096]
            ph, pt = ph_all[idx], pt_all[idx]

            pos_s = model.decode(z, ph, pt)
            neg_t = torch.randint(0, num_entities, (len(ph) * 5,), device=device)
            neg_h = ph.repeat_interleave(5)
            neg_s = model.decode(z, neg_h, neg_t)

            loss = (
                F.binary_cross_entropy_with_logits(pos_s, torch.ones_like(pos_s))
                + F.binary_cross_entropy_with_logits(neg_s, torch.zeros_like(neg_s))
            )
            loss.backward()
            optim.step()

            if epoch % 20 == 0:
                log.info("Epoch %d/%d  loss=%.4f", epoch, rgcn_epochs, loss.item())
                print(f"    Epoch {epoch:3d}/{rgcn_epochs}  loss={loss.item():.4f}")

        # ── Extract embeddings ────────────────────────────────────────────────
        model.eval()
        with torch.no_grad():
            embeddings = model(ei, et).cpu().numpy()   # (num_entities, dim)

        # ── Save encoder ──────────────────────────────────────────────────────
        cond_path = Path(out_dir) / cond_name
        cond_path.mkdir(parents=True, exist_ok=True)
        torch.save(model.state_dict(), cond_path / "rgcn_encoder.pt")

        # ── Classify (fixed split, fixed class list) ──────────────────────────
        mlb = MultiLabelBinarizer(classes=mlb_classes)
        mlb.fit([])   # fit with explicit classes — no need for data

        def _embed(genes: list[str]) -> np.ndarray:
            rows = []
            for g in genes:
                rows.append(embeddings[ent2id[g]] if g in ent2id
                            else np.zeros(dim, dtype=np.float32))
            return np.array(rows, dtype=np.float32)

        X_train = _embed(train_genes)
        X_test  = _embed(test_genes)
        Y_train = mlb.transform([gene_to_pheno.get(g, []) for g in train_genes])
        Y_test  = mlb.transform([gene_to_pheno.get(g, []) for g in test_genes])

        clf = OneVsRestClassifier(
            LogisticRegression(max_iter=1000, random_state=seed, C=1.0),
            n_jobs=n_jobs,
        )
        clf.fit(X_train, Y_train)
        Y_prob = clf.predict_proba(X_test)

        rows_out = []
        for i, mp_term in enumerate(mlb.classes_):
            y_true  = Y_test[:, i]
            y_score = Y_prob[:, i]
            if y_true.sum() == 0:
                continue
            rows_out.append({
                "mp_term_id":  mp_term,
                "ap":          round(float(average_precision_score(y_true, y_score)), 4),
                "n_pos_test":  int(y_true.sum()),
                "n_test":      int(len(y_true)),
                "baseline":    round(float(y_true.mean()), 4),
            })

        per_class  = (
            pd.DataFrame(rows_out)
            .sort_values("ap", ascending=False)
            .reset_index(drop=True)
        )
        mean_auprc    = float(per_class["ap"].mean())
        mean_baseline = float(per_class["baseline"].mean())
        improvement   = mean_auprc / mean_baseline if mean_baseline > 0 else 0.0
        elapsed       = time.time() - t0

        per_class.to_csv(cond_path / "auprc_per_phenotype.csv", index=False)

        result = {
            "condition":        cond_name,
            "relations":        "all" if include_rels is None else "|".join(sorted(include_rels)),
            "mean_auprc":       round(mean_auprc, 4),
            "mean_baseline":    round(mean_baseline, 4),
            "improvement_x":    round(improvement, 3),
            "n_entities":       num_entities,
            "n_context_edges":  len(context_df),
            "n_genes_in_graph": sum(1 for g in all_genes if g in ent2id),
            "n_genes_missing":  n_missing,
            "rgcn_epochs":      rgcn_epochs,
            "rgcn_dim":         dim,
            "elapsed_s":        round(elapsed, 1),
        }
        (cond_path / "summary.json").write_text(json.dumps(result, indent=2))

        log.info(
            "Condition '%s' done — AUPRC=%.4f  improvement=%.2fx  (%.0fs)",
            cond_name, mean_auprc, improvement, elapsed,
        )
        print(f"  Done: AUPRC={mean_auprc:.4f}  "
              f"baseline={mean_baseline:.4f}  "
              f"improvement={improvement:.2f}×  ({elapsed:.0f}s)")

        log.removeHandler(fh)
        fh.close()
        builtins.print = _orig_print
        return result

    except Exception as exc:
        elapsed = time.time() - t0
        log.exception("Condition '%s' FAILED", cond_name)
        log.removeHandler(fh)
        fh.close()
        builtins.print = _orig_print
        return {
            "condition": cond_name,
            "error":     str(exc),
            "elapsed_s": round(elapsed, 1),
        }


# ── Shared data preparation ────────────────────────────────────────────────────

def _prepare_shared_data(
    data_dir: str,
    seed: int,
) -> tuple[list, dict, list[str], list[str], set[str], list[str]]:
    """
    Load triples and derive everything that is shared across conditions.

    Returns data in pickle-safe forms (lists/dicts, not DataFrames/sklearn objects)
    so they can be passed to worker processes via ProcessPoolExecutor.

    Returns
    -------
    all_df_rows    : list of [head, rel, tail] rows (pickle-safe DataFrame substitute)
    gene_to_pheno  : gene → list of MP terms
    train_genes    : fixed train split (list)
    test_genes     : fixed test split (list)
    all_genes      : all phenotype-annotated genes (set)
    mlb_classes    : sorted MP class list
    """
    print("Loading BioCypher triples ...")
    all_df    = _load_all_df(data_dir)
    pheno_df, all_genes = _get_pheno_data(all_df)

    print(f"  Phenotype-annotated genes : {len(all_genes):,}")
    print(f"  Gene→Phenotype triples    : {len(pheno_df):,}")
    print(f"  Unique MP top-terms       : {pheno_df['tail'].nunique():,}")

    train_genes, test_genes = _make_fixed_gene_split(all_genes, seed=seed)
    print(f"\nFixed gene split (seed={seed}): "
          f"{len(train_genes):,} train / {len(test_genes):,} test")

    gene_to_pheno: dict[str, list[str]] = {}
    for _, row in pheno_df.iterrows():
        gene_to_pheno.setdefault(row["head"], []).append(row["tail"])

    # Derive class list (sorted MP terms that appear in any gene's labels)
    all_mp = sorted({mp for mps in gene_to_pheno.values() for mp in mps})
    print(f"Phenotype classes (shared)  : {len(all_mp)}\n")

    # Convert DataFrame to plain list for pickling
    all_df_rows = all_df.values.tolist()

    return all_df_rows, gene_to_pheno, train_genes, test_genes, all_genes, all_mp


# ── Parallel runner ────────────────────────────────────────────────────────────

def _run_parallel(
    conditions:          dict[str, set[str] | None],
    shared:              tuple,
    out_dir:             str,
    rgcn_epochs:         int,
    dim:                 int,
    seed:                int,
    max_workers:         int,
    threads_per_worker:  int,
) -> list[dict]:
    """
    Launch all conditions in parallel using ProcessPoolExecutor (spawn).

    spawn is mandatory: it starts a fresh Python interpreter per worker so that
    OMP_NUM_THREADS set by _set_thread_env() takes effect before torch is
    imported.  With 'fork', torch would already be imported in the parent and
    the OpenMP thread pool would already be initialised — setting the env var
    afterwards has no effect.
    """
    all_df_rows, gene_to_pheno, train_genes, test_genes, all_genes, mlb_classes = shared

    futures = {}
    results = []

    print(f"Parallel mode: {max_workers} workers × {threads_per_worker} threads each")
    print(f"Total CPU budget: {max_workers * threads_per_worker} / {os.cpu_count()} cores\n")

    with ProcessPoolExecutor(
        max_workers=max_workers,
        mp_context=multiprocessing.get_context("spawn"),
        initializer=_set_thread_env,
        initargs=(threads_per_worker,),
    ) as pool:
        for cond_name, include_rels in conditions.items():
            log_file = str(Path(out_dir) / cond_name / "worker.log")
            future = pool.submit(
                _run_condition,
                cond_name     = cond_name,
                include_rels  = include_rels,
                all_df_dict   = all_df_rows,
                gene_to_pheno = gene_to_pheno,
                train_genes   = train_genes,
                test_genes    = test_genes,
                mlb_classes   = mlb_classes,
                all_genes     = all_genes,
                out_dir       = out_dir,
                rgcn_epochs   = rgcn_epochs,
                dim           = dim,
                seed          = seed,
                log_file      = log_file,
                n_jobs        = 1,   # each worker already occupies its CPU budget
            )
            futures[future] = cond_name
            print(f"  Submitted: {cond_name}  (log → {log_file})")

        print()
        for future in as_completed(futures):
            cond_name = futures[future]
            try:
                result = future.result()
                results.append(result)
                if "error" in result:
                    print(f"[FAIL] {cond_name}: {result['error']}")
                else:
                    print(f"[ OK ] {cond_name}: "
                          f"AUPRC={result['mean_auprc']:.4f}  "
                          f"({result['elapsed_s']:.0f}s)")
            except Exception as exc:
                print(f"[FAIL] {cond_name}: {exc}")
                results.append({"condition": cond_name, "error": str(exc)})

    return results


# ── Sequential runner ─────────────────────────────────────────────────────────

def _run_sequential(
    conditions:  dict[str, set[str] | None],
    shared:      tuple,
    out_dir:     str,
    rgcn_epochs: int,
    dim:         int,
    seed:        int,
) -> list[dict]:
    all_df_rows, gene_to_pheno, train_genes, test_genes, all_genes, mlb_classes = shared
    results = []

    for cond_name, include_rels in conditions.items():
        rels_desc = "all" if include_rels is None else "|".join(sorted(include_rels))
        print("=" * 60)
        print(f"  Condition : {cond_name}")
        print(f"  Relations : {rels_desc}")
        print("=" * 60)

        log_file = str(Path(out_dir) / cond_name / "worker.log")
        result = _run_condition(
            cond_name     = cond_name,
            include_rels  = include_rels,
            all_df_dict   = all_df_rows,
            gene_to_pheno = gene_to_pheno,
            train_genes   = train_genes,
            test_genes    = test_genes,
            mlb_classes   = mlb_classes,
            all_genes     = all_genes,
            out_dir       = out_dir,
            rgcn_epochs   = rgcn_epochs,
            dim           = dim,
            seed          = seed,
            log_file      = log_file,
            n_jobs        = -1,   # sequential: use all cores for the classifier
        )
        results.append(result)

    return results


# ── Entry point ───────────────────────────────────────────────────────────────

def main() -> None:
    args   = parse_args()
    out_dir = args.out_dir
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    epochs = _FAST_OVERRIDE["rgcn_epochs"] if args.fast else args.rgcn_epochs
    dim    = _FAST_OVERRIDE["dim"]         if args.fast else args.dim
    if args.fast:
        print("Fast mode: epochs=20, dim=64\n")

    selected = (
        {k: CONDITIONS[k] for k in args.conditions}
        if args.conditions else CONDITIONS
    )

    # ── Load data once in main process ───────────────────────────────────────
    shared = _prepare_shared_data(args.data_dir, seed=args.seed)

    # ── Run conditions ────────────────────────────────────────────────────────
    if args.parallel:
        n_workers = args.max_workers or len(selected)
        n_threads = args.threads_per_worker or max(1, (os.cpu_count() or 1) // n_workers)
        results = _run_parallel(
            conditions         = selected,
            shared             = shared,
            out_dir            = out_dir,
            rgcn_epochs        = epochs,
            dim                = dim,
            seed               = args.seed,
            max_workers        = n_workers,
            threads_per_worker = n_threads,
        )
    else:
        results = _run_sequential(
            conditions  = selected,
            shared      = shared,
            out_dir     = out_dir,
            rgcn_epochs = epochs,
            dim         = dim,
            seed        = args.seed,
        )

    # ── Summary table ─────────────────────────────────────────────────────────
    df = pd.DataFrame(results)
    csv_path = Path(out_dir) / "comparison_summary.csv"
    df.to_csv(csv_path, index=False)

    print("\n" + "=" * 66)
    print("  R-GCN ABLATION SUMMARY")
    print("=" * 66)
    cols = ["condition", "mean_auprc", "mean_baseline", "improvement_x",
            "n_entities", "n_context_edges", "elapsed_s"]
    show = [c for c in cols if c in df.columns]
    print(df[show].to_string(index=False, float_format="%.4f"))
    print("=" * 66)
    print(f"\nResults → {out_dir}/")
    print(f"Per-condition logs → {out_dir}/<condition>/worker.log")


if __name__ == "__main__":
    main()