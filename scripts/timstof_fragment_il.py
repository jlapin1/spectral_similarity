#!/usr/bin/env python3
"""timsTOF observed fragment-intensity I/L discriminative power (Sage outputs).

The timsTOF DDA-PASEF data (Van Puyvelde PXD028735) has no published search results, and
the Bruker .d cannot be converted to mzML in this environment (the Bruker SDK aborts), so
identifications are produced with Sage, which reads .d natively via timsrust. Sage is run
with --annotate-matches, which writes the matched b/y fragment ions and their OBSERVED
intensities per PSM. That is exactly the observed spectrum needed here, so no Python .d
reader is required.

This is the timsTOF analogue of scripts/observed_discriminative_power.py (Astral): for each
identified, unmodified, charge-matched I/L PSM, score observed-vs-predicted(original)
against observed-vs-predicted(I/L swap) and report the fraction correct per metric. The
scoring helpers are reused from observed_discriminative_power.

Usage:
    python scripts/timstof_fragment_il.py \
        --sage-dir temp_data/timstof_sage_out \
        --instrument alphapept_timstof \
        --out-fig temp_data/figure_timstof_fragment_il.png \
        --out-csv temp_data/timstof_fragment_il.csv
"""

import argparse
import os
import re
import sys
import random

import numpy as np
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
for _p in (REPO_ROOT, SCRIPT_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from make_predictions.instruments import get_config, predict_intensities  # noqa: E402
from observed_discriminative_power import (  # noqa: E402
    HIGHER_IS_SIMILAR, pred_map, score_pair, _metric,
)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sage-dir", required=True, help="Directory with results.sage.tsv and matched_fragments.sage.tsv.")
    p.add_argument("--instrument", default="alphapept_timstof")
    p.add_argument("--out-fig", required=True)
    p.add_argument("--out-csv", required=True)
    p.add_argument("--qvalue", type=float, default=0.01)
    p.add_argument("--charge", type=int, default=2)
    p.add_argument("--max-len", type=int, default=30)
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


def bare_sequence(peptide):
    return re.sub(r"\[[^\]]*\]", "", str(peptide)).replace(".", "").strip()


def load_results(path, qvalue, charge, max_len):
    df = pd.read_csv(path, sep="\t")
    if "label" in df.columns:
        df = df[df["label"] == 1]
    if "rank" in df.columns:
        df = df[df["rank"] == 1]
    if "spectrum_q" in df.columns:
        df = df[df["spectrum_q"] < qvalue]
    df = df[~df["peptide"].astype(str).str.contains(r"\[")]  # unmodified
    df["seq"] = df["peptide"].map(bare_sequence)
    df = df[df["charge"] == charge]
    df = df[df["seq"].str.len() <= max_len]
    df = df[df["seq"].str.contains("[IL]", regex=True)]
    return df[["psm_id", "seq"]].drop_duplicates("psm_id").reset_index(drop=True)


def load_matched_fragments(path, psm_ids):
    """{psm_id: {(ion, charge): observed_intensity}} restricted to wanted psm_ids."""
    wanted = set(psm_ids)
    out = {}
    frag = pd.read_csv(path, sep="\t")
    frag = frag[frag["psm_id"].isin(wanted)]
    for psm_id, g in frag.groupby("psm_id"):
        d = {}
        for ftype, ordinal, fcharge, inten in zip(
            g["fragment_type"], g["fragment_ordinals"], g["fragment_charge"], g["fragment_intensity"]
        ):
            if ftype not in ("b", "y"):
                continue
            inten = float(inten)
            if inten <= 0:
                continue
            d[(f"{ftype}{int(ordinal)}", int(fcharge))] = inten
        if len(d) >= 2:
            out[psm_id] = d
    return out


def main():
    args = parse_args()
    rng = random.Random(args.seed)
    config = get_config(args.instrument)

    results_tsv = os.path.join(args.sage_dir, "results.sage.tsv")
    frag_tsv = os.path.join(args.sage_dir, "matched_fragments.sage.tsv")

    print(f"Loading Sage PSMs from {results_tsv} ...")
    res = load_results(results_tsv, args.qvalue, args.charge, args.max_len)
    print(f"{len(res)} unmodified charge-{args.charge} I/L PSMs")

    print(f"Loading observed matched fragments from {frag_tsv} ...")
    obs_by_psm = load_matched_fragments(frag_tsv, res["psm_id"].tolist())
    res = res[res["psm_id"].isin(obs_by_psm)].reset_index(drop=True)
    print(f"{len(res)} PSMs with >=2 matched b/y fragments")

    seqs = sorted(res["seq"].unique())

    def swap_il(seq):
        pos = [i for i, c in enumerate(seq) if c in "IL"]
        if not pos:
            return None
        i = pos[0] if len(pos) == 1 else rng.choice(pos[1:])
        return seq[:i] + ("L" if seq[i] == "I" else "I") + seq[i + 1:]

    swapped = {s: swap_il(s) for s in seqs}
    swapped = {s: v for s, v in swapped.items() if v and v != s}
    seqs = [s for s in seqs if s in swapped]
    sw_list = sorted(set(swapped.values()))
    print(f"{len(seqs)} unique sequences with an I/L swap")

    print(f"Predicting original spectra via {config.intensity_model} ...")
    pred_orig = pred_map(predict_intensities(np.array(seqs), config))
    print("Predicting I/L-swapped spectra ...")
    pred_sw = pred_map(predict_intensities(np.array(sw_list), config))

    metric_keys = list(HIGHER_IS_SIMILAR)
    correct = {k: 0 for k in metric_keys}
    total = {k: 0 for k in metric_keys}
    n_scored = 0

    psm_to_seq = dict(zip(res["psm_id"], res["seq"]))
    for psm_id, seq in psm_to_seq.items():
        sw = swapped.get(seq)
        if sw is None or seq not in pred_orig or sw not in pred_sw:
            continue
        obs = obs_by_psm.get(psm_id)
        if not obs:
            continue
        vo = score_pair(obs, pred_orig[seq])
        vs = score_pair(obs, pred_sw[sw])
        if vo is None or vs is None:
            continue
        n_scored += 1
        for k in metric_keys:
            try:
                so = _metric(k, vo)
                ss = _metric(k, vs)
            except Exception:
                continue
            if np.isnan(so) or np.isnan(ss):
                continue
            total[k] += 1
            correct[k] += int(so > ss) if HIGHER_IS_SIMILAR[k] else int(so < ss)

    rows = [{"metric": k, "ratio_correct": (correct[k] / total[k] if total[k] else float("nan")), "n": total[k]}
            for k in metric_keys]
    res_df = pd.DataFrame(rows).sort_values("ratio_correct", ascending=False)
    os.makedirs(os.path.dirname(args.out_csv), exist_ok=True)
    res_df.to_csv(args.out_csv, index=False)
    print(f"\nScored {n_scored} PSMs. Ratio correct per metric ({config.name}):")
    print(res_df.to_string(index=False))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(res_df["metric"], res_df["ratio_correct"], color="darkorange")
    ax.axhline(0.5, color="k", ls="--", lw=1, label="chance")
    ax.set_ylim(0, 1)
    ax.set_ylabel("fraction correct (original > swap)")
    ax.set_title(f"I/L fragment-intensity discrimination, timsTOF observed ({config.name}, n={n_scored})")
    ax.tick_params(axis="x", rotation=60)
    ax.legend()
    fig.tight_layout()
    os.makedirs(os.path.dirname(args.out_fig), exist_ok=True)
    fig.savefig(args.out_fig, dpi=150)
    print(f"Wrote figure: {args.out_fig}")
    print(f"Wrote table:  {args.out_csv}")


if __name__ == "__main__":
    main()
