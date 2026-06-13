#!/usr/bin/env python3
"""Observed-spectrum I/L discriminative power for any instrument (figure3 ported).

For each identified, unmodified, charge-matched PSM with an I/L residue:
  1. annotate the observed MS2 spectrum against the theoretical b/y ions,
  2. predict the MS2 spectrum of the identified (original) sequence and of its
     I/L-swapped sibling via Koina (instrument preset),
  3. on the ions common to all three, score similarity(observed, predicted_original)
     vs similarity(observed, predicted_swap) for every metric,
  4. mark the PSM "correct" for a metric if the observed spectrum is more similar to the
     original than to the swap, and report the fraction correct per metric.

This is the model-coupled "recovery" experiment (Level B in MANUSCRIPT_IMPROVEMENTS.md):
the FASTA/database sequence is the ground truth; the question is whether the evidence
recovers it over the I/L sibling. Inputs are instrument agnostic, so it runs for the
Orbitrap, Astral, or timsTOF DDA data given the matching identifications and MGFs.

Usage
-----
    python scripts/observed_discriminative_power.py \
        --input temp_data/pb_astral_mq/input_file.txt \
        --mgf-dir temp_data/astral_mgf \
        --instrument astral_dda \
        --out-fig temp_data/figure_astral_ratio_correct.png \
        --out-csv temp_data/astral_ratio_correct.csv
"""

import argparse
import os
import re
import sys

import numpy as np
import pandas as pd
from pyteomics import mgf, mass

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

import random  # noqa: E402
import metrics.metrics as M  # noqa: E402
from make_predictions.instruments import get_config, predict_intensities  # noqa: E402

# Direction of each metric: True if a larger value means more similar.
HIGHER_IS_SIMILAR = {
    "pearson_correlation": True, "spearman_correlation": True, "kendall_tau": True,
    "spectral_angle": True, "weighted_dot_product": True, "fit": True,
    "ruzicka_similarity_1": True, "ruzicka_similarity_2": True,
    "x_corr": True, "mutual_information": True,
    "cosine_similarity": False,  # scipy cosine is a distance (1 - cos)
    "mse": False, "canberra_distance": False, "wasserstein": False, "bray_curtis": False,
}
ANN_RE = re.compile(r"([a-zA-Z]+\d+)\+(\d+)")
SCAN_RE = re.compile(r"[Ss]can[=:](\d+)")
DOTTED_RE = re.compile(r"\.(\d+)\.\d+\.\d+")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, help="MaxQuant-style identification table (input_file.txt).")
    p.add_argument("--mgf-dir", required=True, help="Directory of MGF files (basenames must match the 'Raw file' column).")
    p.add_argument("--instrument", default="astral_dda", help="Instrument preset for prediction.")
    p.add_argument("--out-fig", required=True)
    p.add_argument("--out-csv", required=True)
    p.add_argument("--max-len", type=int, default=30)
    p.add_argument("--charge", type=int, default=2)
    p.add_argument("--tol-ppm", type=float, default=20.0)
    p.add_argument("--max-charge", type=int, default=2, help="Max fragment charge to annotate.")
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


def load_psms(path, max_len, charge):
    df = pd.read_csv(path, sep="\t", low_memory=False)
    df = df[df["Sequence"].str.len() <= max_len]
    # unmodified only (the intensity models here are run without PTMs)
    if "Modifications" in df.columns:
        df = df[df["Modifications"].astype(str).str.lower() == "unmodified"]
    df = df[df["Charge"] == charge]
    df = df.dropna(subset=["MS/MS scan number"])
    df = df[df["Sequence"].str.contains("[IL]", regex=True)]
    df = df.drop_duplicates(subset=["Sequence", "Raw file", "MS/MS scan number"])
    return df.reset_index(drop=True)


def switch_il(seq, rng):
    pos = [m.start() for m in re.finditer(r"[IL]", seq)]
    if not pos:
        return None
    i = pos[0] if len(pos) == 1 else rng.choice(pos[1:])
    return seq[:i] + ("L" if seq[i] == "I" else "I") + seq[i + 1:]


def parse_mgf_dir(mgf_dir):
    """Return {raw_basename: {scan_number: (mz_array, intensity_array)}}."""
    out = {}
    for fn in os.listdir(mgf_dir):
        if not fn.lower().endswith(".mgf"):
            continue
        base = fn[:-4]
        scans = {}
        for spec in mgf.read(os.path.join(mgf_dir, fn)):
            params = spec.get("params", {})
            scan = None
            if "scans" in params:
                try:
                    scan = int(params["scans"])
                except ValueError:
                    scan = None
            if scan is None:
                title = str(params.get("title", ""))
                m = SCAN_RE.search(title) or DOTTED_RE.search(title)
                if m:
                    scan = int(m.group(1))
            if scan is not None:
                scans[scan] = (np.asarray(spec["m/z array"]), np.asarray(spec["intensity array"]))
        out[base] = scans
        print(f"  parsed {fn}: {len(scans)} MS2 scans")
    return out


def theoretical_by(seq, max_charge):
    frags = {}
    n = len(seq)
    for i in range(1, n):
        for c in range(1, max_charge + 1):
            frags[(f"b{i}", c)] = mass.fast_mass(seq[:i], ion_type="b", charge=c)
            frags[(f"y{n - i}", c)] = mass.fast_mass(seq[i:], ion_type="y", charge=c)
    return frags


def annotate(mz_arr, int_arr, frags, tol_ppm):
    """Match observed peaks to theoretical b/y ions; return {(ion,charge): intensity}."""
    out = {}
    if len(mz_arr) == 0:
        return out
    order = np.argsort(mz_arr)
    mz_sorted = mz_arr[order]
    int_sorted = int_arr[order]
    for key, theo in frags.items():
        lo = theo * (1 - tol_ppm * 1e-6)
        hi = theo * (1 + tol_ppm * 1e-6)
        l = np.searchsorted(mz_sorted, lo)
        r = np.searchsorted(mz_sorted, hi)
        if r > l:
            # pick the most intense peak in the window
            j = l + int(np.argmax(int_sorted[l:r]))
            out[key] = float(int_sorted[j])
    return out


def pred_map(df_pred):
    """{peptide: {(ion,charge): (intensity, mz)}} from a Koina prediction frame."""
    out = {}
    for seq, g in df_pred.groupby("peptide_sequences"):
        d = {}
        for ann, inten, mz in zip(g["annotation"], g["intensities"], g["mz"]):
            m = ANN_RE.match(str(ann))
            if not m:
                continue
            inten = float(inten)
            if inten <= 0:
                continue
            d[(m.group(1), int(m.group(2)))] = (inten, float(mz))
        out[seq] = d
    return out


def score_pair(obs, pred):
    """Aligned (obs_int, pred_int, mz) vectors over common ions, plus annotations."""
    keys = sorted(set(obs) & set(pred))
    if len(keys) < 2:
        return None
    oi = np.array([obs[k] for k in keys], dtype=float)
    pi = np.array([pred[k][0] for k in keys], dtype=float)
    mz = np.array([pred[k][1] for k in keys], dtype=float)
    ann = np.array([f"{k[0]}+{k[1]}" for k in keys])
    return oi, pi, mz, ann


def _metric(key, vec):
    """similarity(observed, predicted) for one metric; vec = (obs_int, pred_int, mz, ann)."""
    oi, pi, mz, ann = vec
    return float(getattr(M, key)(
        intensity1=oi, intensity2=pi, mz1=mz, mz2=mz,
        annotation1=ann, annotation2=ann, diagnostic_mz=np.array([]), mz=mz))


def main():
    args = parse_args()
    rng = random.Random(args.seed)
    config = get_config(args.instrument)

    print(f"Loading PSMs from {args.input} ...")
    df = load_psms(args.input, args.max_len, args.charge)
    print(f"{len(df)} unmodified charge-{args.charge} I/L PSMs")

    seqs = sorted(df["Sequence"].unique())
    swapped = {s: switch_il(s, rng) for s in seqs}
    swapped = {s: v for s, v in swapped.items() if v and v != s}
    seqs = [s for s in seqs if s in swapped]
    sw_list = sorted(set(swapped.values()))
    print(f"{len(seqs)} unique sequences with an I/L swap")

    print(f"Predicting original spectra via {config.intensity_model} ...")
    pred_orig = pred_map(predict_intensities(np.array(seqs), config))
    print(f"Predicting I/L-swapped spectra ...")
    pred_sw = pred_map(predict_intensities(np.array(sw_list), config))

    print(f"Parsing MGFs in {args.mgf_dir} ...")
    mgf_dict = parse_mgf_dir(args.mgf_dir)

    metric_keys = [k for k in HIGHER_IS_SIMILAR]
    correct = {k: 0 for k in metric_keys}
    total = {k: 0 for k in metric_keys}
    n_scored = 0

    psm_iter = df[["Sequence", "Raw file", "MS/MS scan number"]].itertuples(index=False, name=None)
    for seq, raw, scan in psm_iter:
        scan = int(scan)
        sw = swapped.get(seq)
        if sw is None or seq not in pred_orig or sw not in pred_sw:
            continue
        spec = mgf_dict.get(raw, {}).get(scan)
        if spec is None:
            continue
        obs = annotate(spec[0], spec[1], theoretical_by(seq, args.max_charge), args.tol_ppm)
        if len(obs) < 2:
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

    rows = []
    for k in metric_keys:
        ratio = correct[k] / total[k] if total[k] else float("nan")
        rows.append({"metric": k, "ratio_correct": ratio, "n": total[k]})
    res = pd.DataFrame(rows).sort_values("ratio_correct", ascending=False)
    os.makedirs(os.path.dirname(args.out_csv), exist_ok=True)
    res.to_csv(args.out_csv, index=False)
    print(f"\nScored {n_scored} PSMs. Ratio correct per metric ({config.name}):")
    print(res.to_string(index=False))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(res["metric"], res["ratio_correct"], color="steelblue")
    ax.axhline(0.5, color="k", ls="--", lw=1, label="chance")
    ax.set_ylim(0, 1)
    ax.set_ylabel("fraction correct (original > swap)")
    ax.set_title(f"I/L discriminative power, observed vs predicted ({config.name}, n={n_scored})")
    ax.tick_params(axis="x", rotation=60)
    ax.legend()
    fig.tight_layout()
    os.makedirs(os.path.dirname(args.out_fig), exist_ok=True)
    fig.savefig(args.out_fig, dpi=150)
    print(f"Wrote figure: {args.out_fig}")
    print(f"Wrote table:  {args.out_csv}")


if __name__ == "__main__":
    main()
