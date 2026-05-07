#!/bin/bash

DIRS=("biocypher_out" "biocypher-log" "out" "graph_tsv")

echo "🧹 Cleaning up output directories..."

for dir in "${DIRS[@]}"; do
    if [ -d "$dir" ]; then
        rm -rf "$dir"
        echo "  ✅ Removed: $dir"
    else
        echo "  ⚠️  Skipped (not found): $dir"
    fi
done

echo "Done."