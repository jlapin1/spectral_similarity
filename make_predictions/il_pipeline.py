"""Reproducible I/L sibling similarity pipeline across instrument configs.

This module turns the notebook glue (sample peptides, build I/L pairs, predict,
align by annotation, score) into reusable functions, and fixes the two
reproducibility hazards documented in ``CLAUDE.md``:

1. ``seq_utils.fasta_to_peptides.create_tryptic_peptides`` returns a
   ``list(set(...))`` whose order is not stable across runs. We sort before
   sampling so the sampled set is reproducible for a given seed.
2. ``seq_utils.peptide.switch_random_il`` uses the global ``random`` module.
   ``build_il_pairs`` seeds it once, and we carry an explicit
   original->switched mapping (a ``pair_id``) into the prediction frames so no
   downstream code ever has to replay the FASTA sampling to recover the link
   (the fragility in ``figure2.ipynb``).
"""

import random

import numpy as np
import pandas as pd

from seq_utils.peptide import (
    switch_random_il,
    remove_non_il,
    remove_ux_containing,
)
from seq_utils.fasta_to_peptides import create_tryptic_peptides
from metrics.get_metrics import metrics_comparison
from make_predictions.instruments import InstrumentConfig, predict_intensities


def load_peptides(fasta_file, min_length=7, max_length=30, sample_size=None, seed=42):
    """Digest a FASTA file and return a reproducible, I/L-containing peptide list.

    The result is sorted (for stable ordering) and, if ``sample_size`` is set,
    sampled with ``numpy.random.default_rng(seed)``.
    """
    peptides = create_tryptic_peptides(
        fasta_file, min_length=min_length, max_length=max_length
    )
    peptides = remove_ux_containing(remove_non_il(peptides))
    peptides = sorted(set(peptides))

    if sample_size is not None and sample_size < len(peptides):
        rng = np.random.default_rng(seed)
        idx = rng.choice(len(peptides), size=sample_size, replace=False)
        peptides = [peptides[i] for i in sorted(idx)]
    return peptides


def build_il_pairs(peptides, seed=42):
    """Build deterministic (original, switched) I/L pairs.

    Seeds the global ``random`` module (used internally by ``switch_random_il``)
    so the swaps are reproducible. Peptides whose swap is a no-op (no I/L outside
    brackets) are dropped. Returns two parallel lists ``(originals, switched)``.
    """
    random.seed(seed)
    originals, switched = [], []
    seen = set()
    for peptide in sorted(set(peptides)):
        if peptide in seen:
            continue
        swap = switch_random_il(peptide)
        if swap != peptide:
            originals.append(peptide)
            switched.append(swap)
            seen.add(peptide)
    return originals, switched


_CORE_METRICS = ["spectral_angle", "mse", "weighted_dot_product"]


def similarity_for_instrument(originals, switched, config: InstrumentConfig, charges=None):
    """Predict and score original-vs-switched spectra for one instrument config.

    Returns a tidy DataFrame: one row per peptide pair, columns are the metrics
    from ``metrics.get_metrics.metric_keys`` plus ``instrument`` and
    ``pair_key``. Metric columns are coerced to numeric (broken metrics such as
    the known ``hyper_score`` bug become ``NaN`` rather than raising).
    """
    pairs = pd.DataFrame(
        {"orig_seq": originals, "sw_seq": switched}
    ).drop_duplicates("orig_seq").reset_index(drop=True)
    pairs["ID"] = np.arange(len(pairs))

    orig_pred = predict_intensities(pairs["orig_seq"].to_numpy(), config, charges=charges)
    sw_pred = predict_intensities(pairs["sw_seq"].to_numpy(), config, charges=charges)

    # Attach the shared pair ID so metrics_comparison can align the two spectra.
    orig_pred = orig_pred.merge(
        pairs[["ID", "orig_seq"]], left_on="peptide_sequences", right_on="orig_seq", how="inner"
    )
    sw_pred = sw_pred.merge(
        pairs[["ID", "sw_seq"]], left_on="peptide_sequences", right_on="sw_seq", how="inner"
    )

    scores = metrics_comparison(orig_pred, sw_pred)
    scores = scores.apply(pd.to_numeric, errors="coerce")
    scores["instrument"] = config.name
    return scores.reset_index(names="pair_key")


def run_comparison(originals, switched, configs, charges=None):
    """Run ``similarity_for_instrument`` for each config and concatenate results."""
    frames = []
    for config in configs:
        print(
            f"[{config.name}] predicting {len(originals)} pairs "
            f"via {config.intensity_model} ..."
        )
        frames.append(similarity_for_instrument(originals, switched, config, charges=charges))
    return pd.concat(frames, ignore_index=True)


def summarize(scores, metrics=None):
    """Return mean per metric per instrument (a quick comparison table)."""
    metrics = metrics or _CORE_METRICS
    cols = [m for m in metrics if m in scores.columns]
    return scores.groupby("instrument")[cols].mean()
