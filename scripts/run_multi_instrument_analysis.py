#!/usr/bin/env python3
"""Run the I/L sibling spectral-similarity analysis across multiple instruments.

This is the runnable end-to-end form of the new multi-instrument extension. It
digests a FASTA file, builds reproducible I/L sibling pairs, fetches predicted
MS2 spectra from Koina for each requested instrument configuration, scores
original-vs-switched similarity with every metric, writes the tidy results to
parquet, prints a summary table, and saves an overlay figure.

The only network dependency is Koina (``koina.wilhelmlab.org``); no bulk raw-data
download is required because the predicted branch does not need observed spectra.

Examples
--------
    # default 2000 peptides, three instruments
    python scripts/run_multi_instrument_analysis.py

    # custom selection and sample size
    python scripts/run_multi_instrument_analysis.py \
        --n 5000 \
        --instruments unispec_lumos,unispec_qe,alphapept_timstof,prosit_cid \
        --out temp_data/multi_instrument_similarity.parquet \
        --fig temp_data/figure7_multi_instrument.png

Run ``--list`` to see the available instrument presets.
"""

import argparse
import os
import sys

# Make the repo root importable when run from anywhere.
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from make_predictions.instruments import INSTRUMENT_PRESETS, get_config  # noqa: E402
from make_predictions.il_pipeline import (  # noqa: E402
    load_peptides,
    build_il_pairs,
    run_comparison,
    summarize,
)

DEFAULT_INSTRUMENTS = ["unispec_lumos", "astral_dda", "alphapept_timstof"]
CORE_METRICS = ["spectral_angle", "mse", "weighted_dot_product"]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--fasta",
        default=os.path.join(REPO_ROOT, "fasta", "UP000005640_9606.fasta"),
        help="FASTA file to digest (default: bundled human proteome).",
    )
    parser.add_argument("--n", type=int, default=2000, help="Number of peptides to sample (default 2000).")
    parser.add_argument("--seed", type=int, default=42, help="Random seed (default 42).")
    parser.add_argument(
        "--instruments",
        default=",".join(DEFAULT_INSTRUMENTS),
        help="Comma-separated instrument preset keys (see --list).",
    )
    parser.add_argument(
        "--out",
        default=os.path.join(REPO_ROOT, "temp_data", "multi_instrument_similarity.parquet"),
        help="Output parquet path for tidy per-pair scores.",
    )
    parser.add_argument(
        "--fig",
        default=os.path.join(REPO_ROOT, "temp_data", "figure7_multi_instrument.png"),
        help="Output figure path.",
    )
    parser.add_argument("--list", action="store_true", help="List available instrument presets and exit.")
    return parser.parse_args()


def plot_comparison(scores, fig_path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    metrics = [m for m in CORE_METRICS if m in scores.columns]
    fig, axes = plt.subplots(1, len(metrics), figsize=(5 * len(metrics), 4))
    if len(metrics) == 1:
        axes = [axes]
    for ax, metric in zip(axes, metrics):
        for instrument, sub in scores.groupby("instrument"):
            values = sub[metric].dropna()
            if len(values) > 1:
                sns.kdeplot(values, ax=ax, label=instrument, fill=True, alpha=0.25)
        ax.set_title(metric)
        ax.set_xlabel(metric)
        ax.legend(fontsize=8)
    fig.suptitle("I/L sibling spectral similarity by instrument (predicted spectra)")
    fig.tight_layout()
    fig.savefig(fig_path, dpi=150)
    print(f"Wrote figure: {fig_path}")


def main():
    args = parse_args()

    if args.list:
        print("Available instrument presets:")
        for key, cfg in sorted(INSTRUMENT_PRESETS.items()):
            print(
                f"  {key:20s} model={cfg.intensity_model:30s} "
                f"instrument={cfg.instrument_type} frag={cfg.fragmentation_type} ce={cfg.collision_energy}"
            )
        return

    configs = [get_config(name.strip()) for name in args.instruments.split(",") if name.strip()]
    print(f"Instruments: {[c.name for c in configs]}")

    print(f"Loading and sampling peptides from {args.fasta} (n={args.n}, seed={args.seed}) ...")
    peptides = load_peptides(args.fasta, sample_size=args.n, seed=args.seed)
    originals, switched = build_il_pairs(peptides, seed=args.seed)
    print(f"Built {len(originals)} I/L sibling pairs.")

    scores = run_comparison(originals, switched, configs)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    scores.to_parquet(args.out, index=False)
    print(f"Wrote scores: {args.out}  (rows={len(scores)})")

    print("\nMean similarity per metric per instrument:")
    print(summarize(scores, CORE_METRICS).round(4).to_string())

    plot_comparison(scores, args.fig)


if __name__ == "__main__":
    main()
