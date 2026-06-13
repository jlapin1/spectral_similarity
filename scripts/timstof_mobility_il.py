#!/usr/bin/env python3
"""timsTOF observed ion-mobility separation of I/L siblings (model-free).

Uses the observed ion mobility (1/K0) in a ProteoBench diaPASEF DIA-NN report. No
prediction and no mobility/CCS unit conversion are needed, because I/L siblings are
isobaric (identical precursor m/z and charge), so their observed 1/K0 values are
directly comparable. This is Level A (physics) evidence for the ion-mobility modality:
does an I/L swap change observed mobility beyond measurement noise?

Three distributions of |delta(1/K0)| are compared:
  1. replicate null   - same peptide measured in different runs (measurement noise),
  2. I/L siblings     - pairs where both the sequence and an I/L-swapped sibling are
                        independently identified (the I/L effect),
  3. random reference - random m/z-matched non-sibling peptide pairs (the scale of a
                        genuine mobility difference between unrelated peptides).

Usage:
    python scripts/timstof_mobility_il.py \
        --parquet temp_data/pb_diann_inspect/input_file.parquet \
        --out-fig temp_data/figure_timstof_mobility_il.png \
        --out-csv temp_data/timstof_mobility_il.csv
"""

import argparse
import os
import re
import random

import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--parquet", required=True, help="diaPASEF DIA-NN report parquet with IM column.")
    p.add_argument("--out-fig", required=True)
    p.add_argument("--out-csv", required=True)
    p.add_argument("--charge", type=int, default=2)
    p.add_argument("--qvalue", type=float, default=0.01)
    p.add_argument("--max-len", type=int, default=30)
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


def all_swaps(seq):
    out = []
    for m in re.finditer("[IL]", seq):
        i = m.start()
        out.append(seq[:i] + ("L" if seq[i] == "I" else "I") + seq[i + 1:])
    return out


def main():
    args = parse_args()
    rng = random.Random(args.seed)

    df = pd.read_parquet(
        args.parquet,
        columns=["Stripped.Sequence", "Precursor.Charge", "IM", "Precursor.Mz", "Q.Value", "Run"],
    )
    df = df[(df["Q.Value"] < args.qvalue) & (df["Precursor.Charge"] == args.charge)]
    df = df[df["Stripped.Sequence"].str.contains("[IL]", regex=True)]
    df = df[df["Stripped.Sequence"].str.len() <= args.max_len]

    # replicate null: per-peptide spread of observed IM across runs
    rep_diffs = []
    per_seq_im = {}
    per_seq_mz = {}
    for seq, sub in df.groupby("Stripped.Sequence"):
        ims = sub["IM"].to_numpy()
        per_seq_im[seq] = float(np.median(ims))
        per_seq_mz[seq] = float(np.median(sub["Precursor.Mz"]))
        if len(ims) >= 2:
            # split-half difference of medians: noise on the same statistic used for siblings
            idx = list(range(len(ims)))
            rng.shuffle(idx)
            half = len(idx) // 2
            a = np.median(ims[idx[:half]]) if half else ims[idx[0]]
            b = np.median(ims[idx[half:]])
            rep_diffs.append(abs(a - b))
    rep_diffs = np.array(rep_diffs)

    seqset = set(per_seq_im)

    # I/L siblings: both members observed
    pairs = set()
    for s in seqset:
        for sw in all_swaps(s):
            if sw in seqset:
                pairs.add(tuple(sorted((s, sw))))
    sib_diffs = np.array([abs(per_seq_im[a] - per_seq_im[b]) for a, b in pairs])

    # random m/z-matched non-sibling reference (same isobaric constraint as siblings)
    seqs = sorted(seqset)
    mzs = np.array([per_seq_mz[s] for s in seqs])
    order = np.argsort(mzs)
    seqs_sorted = [seqs[i] for i in order]
    mzs_sorted = mzs[order]
    rand_diffs = []
    n_rand = max(2000, 20 * len(sib_diffs))
    for _ in range(n_rand):
        i = rng.randrange(len(seqs_sorted))
        # nearest-in-mz window, pick a different peptide
        lo = np.searchsorted(mzs_sorted, mzs_sorted[i] - 1.0)
        hi = np.searchsorted(mzs_sorted, mzs_sorted[i] + 1.0)
        if hi - lo < 2:
            continue
        j = rng.randrange(lo, hi)
        if j == i:
            continue
        rand_diffs.append(abs(per_seq_im[seqs_sorted[i]] - per_seq_im[seqs_sorted[j]]))
    rand_diffs = np.array(rand_diffs)

    # significance: are sibling diffs larger than replicate noise?
    try:
        from scipy.stats import mannwhitneyu
        u_stat, p_val = mannwhitneyu(sib_diffs, rep_diffs, alternative="greater")
    except Exception:
        p_val = float("nan")

    summary = pd.DataFrame([
        {"group": "replicate_noise", "n": len(rep_diffs), "median_dIM": np.median(rep_diffs)},
        {"group": "il_siblings", "n": len(sib_diffs), "median_dIM": np.median(sib_diffs)},
        {"group": "random_mz_matched", "n": len(rand_diffs), "median_dIM": np.median(rand_diffs)},
    ])
    os.makedirs(os.path.dirname(args.out_csv), exist_ok=True)
    summary.to_csv(args.out_csv, index=False)
    print(summary.to_string(index=False))
    print(f"\nMann-Whitney (siblings > replicate noise) p = {p_val:.3g}")
    print(f"I/L sibling median |dIM| = {np.median(sib_diffs):.5f} vs replicate noise "
          f"{np.median(rep_diffs):.5f} vs random {np.median(rand_diffs):.5f} (1/K0)")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 5))
    bins = np.linspace(0, np.percentile(rand_diffs, 95), 60)
    ax.hist(rand_diffs, bins=bins, density=True, alpha=0.4, label=f"random m/z-matched (n={len(rand_diffs)})", color="grey")
    ax.hist(rep_diffs, bins=bins, density=True, alpha=0.5, label=f"replicate noise (n={len(rep_diffs)})", color="steelblue")
    ax.hist(sib_diffs, bins=bins, density=True, alpha=0.6, label=f"I/L siblings (n={len(sib_diffs)})", color="firebrick")
    ax.axvline(np.median(sib_diffs), color="firebrick", ls="--", lw=1)
    ax.set_xlabel("|delta (1/K0)| between paired precursors")
    ax.set_ylabel("density")
    ax.set_title("timsTOF observed ion mobility: I/L siblings vs noise vs random (diaPASEF)")
    ax.legend()
    fig.tight_layout()
    os.makedirs(os.path.dirname(args.out_fig), exist_ok=True)
    fig.savefig(args.out_fig, dpi=150)
    print(f"Wrote figure: {args.out_fig}")
    print(f"Wrote table:  {args.out_csv}")


if __name__ == "__main__":
    main()
