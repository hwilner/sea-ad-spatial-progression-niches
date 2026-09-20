#!/usr/bin/env python
"""Run the spatial pipeline on real staged SEA-AD MTG MERFISH cells.

Per section: k-NN neighbor graph on real (x, y) coordinates, subclass
neighborhood-enrichment z-scores (label-permutation null), and niche
detection on neighbor-composition vectors (k-means). Aggregated results
are written to ``reports/`` as small CSV/JSON files (safe to commit).

Usage
-----
    python scripts/analyze_seaad.py --cache-dir data/seaad
    python scripts/analyze_seaad.py --cache-dir data/seaad --max-sections 12
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from seaad_niches.io import load_merfish_cells  # noqa: E402
from seaad_niches.niches import assign_niches  # noqa: E402
from seaad_niches.spatial import (  # noqa: E402
    build_neighbor_graph,
    neighbor_composition,
    neighborhood_enrichment,
)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache-dir", default="data/seaad")
    ap.add_argument("--reports-dir", default="reports")
    ap.add_argument("--max-sections", type=int, default=None,
                    help="limit the number of sections analyzed (default: all)")
    ap.add_argument("--k-neighbors", type=int, default=6)
    ap.add_argument("--n-permutations", type=int, default=100)
    ap.add_argument("--n-niches", type=int, default=8)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    cells = load_merfish_cells(args.cache_dir)
    sections = list(pd.unique(cells["section"]))
    if args.max_sections:
        sections = sections[: args.max_sections]
    print(f"Loaded {len(cells)} cells, {cells['donor'].nunique()} donors, "
          f"{cells['section'].nunique()} sections; analyzing {len(sections)} sections")

    z_mats, niche_rows = [], []
    for i, sec in enumerate(sections):
        sub = cells[cells["section"] == sec]
        coords = sub[["x", "y"]].to_numpy(dtype=float)
        adj = build_neighbor_graph(coords, k=args.k_neighbors)
        z = neighborhood_enrichment(
            adj, sub["cell_type"],
            n_permutations=args.n_permutations, seed=args.seed + i,
        )
        z_mats.append(z.to_numpy())
        comp = neighbor_composition(adj, sub["cell_type"])
        niches = assign_niches(comp, n_clusters=args.n_niches,
                               method="kmeans", seed=args.seed)
        sizes = np.bincount(niches, minlength=args.n_niches)
        niche_rows.append(
            {
                "section": sec,
                "donor": sub["donor"].iloc[0],
                "n_cells": len(sub),
                "n_niches": int(len(np.unique(niches))),
                "largest_niche_fraction": float(sizes.max() / sizes.sum()),
                "smallest_niche_fraction": float(sizes.min() / sizes.sum()),
            }
        )
        print(f"  [{i + 1}/{len(sections)}] {sec}: {len(sub)} cells done")

    cats = list(z.index)
    mean_z = pd.DataFrame(np.mean(z_mats, axis=0), index=cats, columns=cats)
    sem_z = pd.DataFrame(
        np.std(z_mats, axis=0) / np.sqrt(len(z_mats)), index=cats, columns=cats
    )

    reports = Path(args.reports_dir)
    reports.mkdir(parents=True, exist_ok=True)
    mean_z.to_csv(reports / "seaad_merfish_enrichment_mean_z.csv")
    sem_z.to_csv(reports / "seaad_merfish_enrichment_sem_z.csv")
    niche_df = pd.DataFrame(niche_rows)
    niche_df.to_csv(reports / "seaad_merfish_niche_summary.csv", index=False)

    # Top enriched / depleted subclass pairs (upper triangle, mean z).
    iu = np.triu_indices(len(cats), k=1)
    pairs = pd.DataFrame(
        {
            "type_a": [cats[i] for i in iu[0]],
            "type_b": [cats[j] for j in iu[1]],
            "mean_z": mean_z.to_numpy()[iu],
            "sem_z": sem_z.to_numpy()[iu],
        }
    ).sort_values("mean_z", ascending=False)
    pairs.to_csv(reports / "seaad_merfish_top_pairs.csv", index=False)

    run_summary = {
        "dataset": "SEA-AD MTG MERFISH (open access, AWS Open Data Registry)",
        "n_cells_loaded": int(len(cells)),
        "n_donors": int(cells["donor"].nunique()),
        "n_sections_analyzed": int(len(sections)),
        "n_cell_types": int(len(cats)),
        "k_neighbors": args.k_neighbors,
        "n_permutations": args.n_permutations,
        "n_niches": args.n_niches,
        "seed": args.seed,
        "top_enriched_pair": pairs.iloc[0].to_dict(),
        "top_depleted_pair": pairs.iloc[-1].to_dict(),
    }
    (reports / "seaad_merfish_run_summary.json").write_text(
        json.dumps(run_summary, indent=2)
    )
    print(f"Wrote reports to {reports}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
