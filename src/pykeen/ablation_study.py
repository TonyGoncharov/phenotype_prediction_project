from src.pykeen.train import run_training

layers = {
    "all":         None,
    "no_go":       {"HumanGeneHasMpTopTerm", "HumanGeneExpressedInTissue", "HumanGeneEncodesProtein", "HumanProteinInteractsWith"},
    "no_ppi":      {"HumanGeneHasMpTopTerm", "HumanGeneHasGoTerm", "HumanGeneExpressedInTissue"},
    "no_expr":     {"HumanGeneHasMpTopTerm", "HumanGeneHasGoTerm", "HumanGeneEncodesProtein", "HumanProteinInteractsWith"},
    "pheno_only":  {"HumanGeneHasMpTopTerm"},
}

if __name__ == "__main__":
    for name, rels in layers.items():
        run_training(
            data_dir="biocypher_out/human/",
            out_dir=f"pykeen_out/ablation/{name}",
            config={"include_relations": rels},
        )