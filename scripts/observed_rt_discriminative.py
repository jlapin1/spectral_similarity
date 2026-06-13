#!/usr/bin/env python3
"""Observed retention-time I/L discriminative power (RT panel), any instrument.

Calibrates DeepLC on a subset of the dataset's observed retention times, then for an
evaluation subset predicts the RT of the identified (original) sequence and of its
I/L-swapped sibling. A PSM is "correct" if the observed RT is closer to the predicted
original than to the predicted swap. Reports the fraction correct and the magnitude of
the predicted I/L RT shift.

Runs under the py312_deeplc conda environment:
    /home/robbin/anaconda3/envs/py312_deeplc/bin/python scripts/observed_rt_discriminative.py ...

Inputs:
  --format maxquant : Sequence / Modifications / (Calibrated) Retention time / Charge
  --format diann    : Stripped.Sequence / Precursor.Charge / RT / Q.Value (parquet)
"""

import argparse
import os
import random

import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True)
    p.add_argument("--format", choices=["maxquant", "diann"], required=True)
    p.add_argument("--instrument", required=True, help="Label for outputs.")
    p.add_argument("--out-fig", required=True)
    p.add_argument("--combined-csv", required=True, help="Shared CSV; one row appended per instrument.")
    p.add_argument("--charge", type=int, default=2)
    p.add_argument("--qvalue", type=float, default=0.01)
    p.add_argument("--max-len", type=int, default=30)
    p.add_argument("--n-cal", type=int, default=1500)
    p.add_argument("--n-eval", type=int, default=4000)
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


def load(path, fmt, max_len, charge, qval):
    if fmt == "maxquant":
        df = pd.read_csv(path, sep="\t", low_memory=False)
        df = df[df["Sequence"].str.len() <= max_len]
        if "Modifications" in df.columns:
            df = df[df["Modifications"].astype(str).str.lower() == "unmodified"]
        if charge and "Charge" in df.columns:
            df = df[df["Charge"] == charge]
        rtcol = "Calibrated retention time" if "Calibrated retention time" in df.columns else "Retention time"
        df = df.rename(columns={"Sequence": "seq", rtcol: "rt"})
    else:
        df = pd.read_parquet(path, columns=["Stripped.Sequence", "Precursor.Charge", "RT", "Q.Value"])
        df = df[df["Q.Value"] < qval]
        if charge:
            df = df[df["Precursor.Charge"] == charge]
        df = df[df["Stripped.Sequence"].str.len() <= max_len]
        df = df.rename(columns={"Stripped.Sequence": "seq", "RT": "rt"})
    df = df.dropna(subset=["seq", "rt"])
    df = df[df["seq"].str.contains("[IL]", regex=True)]
    return df.groupby("seq")["rt"].median().reset_index()


def swap_il(seq, rng):
    pos = [i for i, c in enumerate(seq) if c in "IL"]
    if not pos:
        return None
    i = pos[0] if len(pos) == 1 else rng.choice(pos[1:])
    return seq[:i] + ("L" if seq[i] == "I" else "I") + seq[i + 1:]


def deeplc_preds(seqs, charge, dlc):
    from psm_utils import PSM, PSMList
    psms = PSMList(psm_list=[
        PSM(peptidoform=f"{s}/{charge}", spectrum_id=str(i)) for i, s in enumerate(seqs)
    ])
    return np.asarray(dlc.make_preds(psm_list=psms))


def main():
    args = parse_args()
    rng = random.Random(args.seed)

    g = load(args.input, args.format, args.max_len, args.charge, args.qvalue)
    g["sw"] = g["seq"].map(lambda s: swap_il(s, rng))
    g = g[g["sw"].notna() & (g["sw"] != g["seq"])].reset_index(drop=True)
    print(f"[{args.instrument}] {len(g)} unique I/L peptides with observed RT")

    g = g.sample(frac=1.0, random_state=args.seed).reset_index(drop=True)
    n_cal = min(args.n_cal, len(g) // 3)
    cal = g.iloc[:n_cal]
    ev = g.iloc[n_cal:n_cal + args.n_eval].reset_index(drop=True)
    print(f"  calibration={len(cal)}  eval={len(ev)}")

    from psm_utils import PSM, PSMList
    from deeplc import DeepLC
    dlc = DeepLC()
    cal_psms = PSMList(psm_list=[
        PSM(peptidoform=f"{s}/{args.charge}", retention_time=float(rt), spectrum_id=str(i))
        for i, (s, rt) in enumerate(zip(cal["seq"], cal["rt"]))
    ])
    dlc.calibrate_preds(psm_list=cal_psms)

    pred_orig = deeplc_preds(ev["seq"].tolist(), args.charge, dlc)
    pred_sw = deeplc_preds(ev["sw"].tolist(), args.charge, dlc)
    obs = ev["rt"].to_numpy()

    err_orig = np.abs(obs - pred_orig)
    err_sw = np.abs(obs - pred_sw)
    valid = np.isfinite(err_orig) & np.isfinite(err_sw)
    correct = (err_orig[valid] < err_sw[valid])
    ratio = float(correct.mean())
    n = int(valid.sum())
    dpred = (pred_orig - pred_sw)[valid]
    median_abs_dpred = float(np.median(np.abs(dpred)))
    median_obs_err = float(np.median(err_orig[valid]))

    # Wilson 95% CI for the proportion
    from math import sqrt
    z = 1.96
    phat = ratio
    denom = 1 + z * z / n
    centre = (phat + z * z / (2 * n)) / denom
    half = z * sqrt(phat * (1 - phat) / n + z * z / (4 * n * n)) / denom
    ci_lo, ci_hi = centre - half, centre + half

    print(f"  ratio_correct={ratio:.4f}  (95% CI {ci_lo:.4f}-{ci_hi:.4f}, n={n})")
    print(f"  median |predicted RT shift orig vs swap| = {median_abs_dpred:.4f} min")
    print(f"  median observed-vs-predicted RT error    = {median_obs_err:.4f} min")

    row = pd.DataFrame([{
        "instrument": args.instrument, "ratio_correct": ratio, "ci_lo": ci_lo, "ci_hi": ci_hi,
        "n": n, "median_abs_dRT_pred": median_abs_dpred, "median_obs_pred_error": median_obs_err,
    }])
    if os.path.exists(args.combined_csv):
        prev = pd.read_csv(args.combined_csv)
        prev = prev[prev["instrument"] != args.instrument]
        row = pd.concat([prev, row], ignore_index=True)
    os.makedirs(os.path.dirname(args.combined_csv) or ".", exist_ok=True)
    row.to_csv(args.combined_csv, index=False)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4))
    a1.bar(["original closer", "swap closer"], [correct.mean(), 1 - correct.mean()], color=["seagreen", "lightgrey"])
    a1.axhline(0.5, color="k", ls="--", lw=1)
    a1.set_ylim(0, 1)
    a1.set_ylabel("fraction of PSMs")
    a1.set_title(f"RT recovery ({args.instrument})\nratio correct={ratio:.3f} [{ci_lo:.3f}-{ci_hi:.3f}], n={n}")
    a2.hist(dpred, bins=60, color="steelblue")
    a2.axvline(0, color="k", lw=1)
    a2.set_xlabel("predicted RT(original) - RT(I/L swap)  [min]")
    a2.set_ylabel("count")
    a2.set_title(f"predicted I/L RT shift (median |.|={median_abs_dpred:.3f} min)")
    fig.tight_layout()
    os.makedirs(os.path.dirname(args.out_fig) or ".", exist_ok=True)
    fig.savefig(args.out_fig, dpi=150)
    print(f"  wrote {args.out_fig}")


if __name__ == "__main__":
    main()
