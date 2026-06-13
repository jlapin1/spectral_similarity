"""Instrument configuration for multi-instrument I/L sibling analysis.

This module centralizes the instrument and acquisition parameters that were
previously hardcoded inside ``make_predictions/intensity_predictions.py`` and the
notebooks (``instrument_types="LUMOS"``, ``fragmentation_types="HCD"``,
``collision_energies=28``, ``precursor_charges=2``). Adding a new instrument is
now a single entry in ``INSTRUMENT_PRESETS`` rather than an edit scattered across
notebook cells.

Each ``InstrumentConfig`` names a Koina intensity model and the acquisition
parameters that model expects. Models differ in which input columns they accept:

* Prosit HCD/CID models take ``collision_energies`` only (no instrument or
  fragmentation column). Set ``instrument_type`` and ``fragmentation_type`` to
  ``None`` for those.
* UniSpec and AlphaPeptDeep are instrument aware and take ``collision_energies``,
  ``instrument_types`` and ``fragmentation_types``. Vary ``instrument_type`` to
  compare instruments with the same model.

The accepted enum values for ``instrument_type`` and ``fragmentation_type`` are
defined by each Koina model card (https://koina.wilhelmlab.org). If Koina rejects
a value, adjust the preset to a value the model accepts.
"""

from dataclasses import dataclass
from typing import Optional
import time

import numpy as np
import pandas as pd
from koinapy import Koina

KOINA_SERVER = "koina.wilhelmlab.org:443"


@dataclass(frozen=True)
class InstrumentConfig:
    """A named instrument/acquisition setting for spectrum prediction.

    Parameters
    ----------
    name:
        Short identifier used as a label in result tables and figures.
    intensity_model:
        Koina MS2 intensity model (e.g. ``"UniSpec"``,
        ``"Prosit_2020_intensity_HCD"``, ``"AlphaPeptDeep_ms2_generic"``).
    fragmentation_type:
        Fragmentation enum for instrument-aware models (e.g. ``"HCD"``,
        ``"CID"``). ``None`` for models that do not take this column.
    instrument_type:
        Instrument enum for instrument-aware models (e.g. ``"LUMOS"``, ``"QE"``,
        ``"TIMSTOF"``). ``None`` for models that do not take this column.
    collision_energy:
        Normalized collision energy. ``None`` for models that ignore it.
    rt_model, ccs_model:
        Koina models for the retention time and CCS branches (used by callers
        that extend the comparison to RT/CCS).
    default_charge:
        Precursor charge applied when no per-peptide charge is provided.
    """

    name: str
    intensity_model: str
    fragmentation_type: Optional[str] = "HCD"
    instrument_type: Optional[str] = "LUMOS"
    collision_energy: Optional[float] = 28.0
    rt_model: str = "Deeplc_hela_hf"
    ccs_model: str = "AlphaPept_ccs_generic"
    default_charge: int = 2


# Preset registry. Add a new instrument by adding an entry here.
# The UniSpec and AlphaPeptDeep families are instrument aware, so varying
# instrument_type with one model gives a controlled instrument comparison.
INSTRUMENT_PRESETS = {
    # UniSpec across Orbitrap instruments (same model, different instrument).
    "unispec_lumos": InstrumentConfig(
        "unispec_lumos", "UniSpec", "HCD", "LUMOS", 28.0
    ),
    "unispec_qe": InstrumentConfig(
        "unispec_qe", "UniSpec", "HCD", "QE", 28.0
    ),
    "unispec_elite": InstrumentConfig(
        "unispec_elite", "UniSpec", "HCD", "ELITE", 28.0
    ),
    # AlphaPeptDeep across an Orbitrap and a TOF instrument.
    "alphapept_lumos": InstrumentConfig(
        "alphapept_lumos", "AlphaPeptDeep_ms2_generic", "HCD", "LUMOS", 28.0
    ),
    # timsTOF Pro DDA-PASEF benchmark instrument (Van Puyvelde Generation Alpha,
    # PXD028735). Activation on Bruker timsTOF is collision-induced; AlphaPeptDeep is
    # instrument aware. Verify the accepted instrument/fragmentation enum on the model
    # card if Koina rejects a value. This instrument also feeds the observed-CCS arm
    # via the IM2Deep / AlphaPept CCS models, because timsTOF yields observed ion
    # mobility (the current CCS analysis is predicted only).
    "alphapept_timstof": InstrumentConfig(
        "alphapept_timstof", "AlphaPeptDeep_ms2_generic", "HCD", "TIMSTOF", 30.0
    ),
    # Orbitrap Astral DDA benchmark instrument (Generation Beta, PXD070049/PXD071205):
    # HCD at NCE 30, 15 min gradient. Astral is not in the UniSpec/AlphaPeptDeep
    # instrument enums, so Prosit (which needs only the collision energy) is the robust
    # choice. The NCE may need calibration to Prosit's internal scale.
    "astral_dda": InstrumentConfig(
        "astral_dda", "Prosit_2020_intensity_HCD",
        fragmentation_type=None, instrument_type=None, collision_energy=30.0,
    ),
    # Prosit fragmentation comparison (no instrument/fragmentation column).
    "prosit_hcd": InstrumentConfig(
        "prosit_hcd", "Prosit_2020_intensity_HCD",
        fragmentation_type=None, instrument_type=None, collision_energy=28.0,
    ),
    "prosit_cid": InstrumentConfig(
        "prosit_cid", "Prosit_2020_intensity_CID",
        fragmentation_type=None, instrument_type=None, collision_energy=35.0,
    ),
}


def get_config(name: str) -> InstrumentConfig:
    """Return a preset by key, raising a clear error listing valid keys."""
    try:
        return INSTRUMENT_PRESETS[name]
    except KeyError as exc:
        valid = ", ".join(sorted(INSTRUMENT_PRESETS))
        raise KeyError(f"Unknown instrument preset '{name}'. Valid: {valid}") from exc


def _build_input_frame(peptides, charges, config: InstrumentConfig) -> pd.DataFrame:
    """Assemble the Koina input frame, including only the columns the model needs."""
    peptides = np.asarray(peptides)
    n = len(peptides)
    if charges is None:
        charges = np.full(n, config.default_charge)
    charges = np.asarray(charges)

    inputs = pd.DataFrame(
        {"peptide_sequences": peptides, "precursor_charges": charges}
    )
    if config.collision_energy is not None:
        inputs["collision_energies"] = np.full(n, config.collision_energy)
    if config.instrument_type is not None:
        inputs["instrument_types"] = np.full(n, config.instrument_type)
    if config.fragmentation_type is not None:
        inputs["fragmentation_types"] = np.full(n, config.fragmentation_type)
    return inputs.drop_duplicates().reset_index(drop=True)


def _predict_once(model, inputs, max_retries, delay):
    """Single Koina predict with retry on transient errors; raises on final failure."""
    for attempt in range(max_retries):
        try:
            return model.predict(inputs, debug=True)
        except Exception:  # noqa: BLE001 - Koina raises generic errors
            if attempt < max_retries - 1:
                time.sleep(delay)
            else:
                raise


def predict_intensities(
    peptides,
    config: InstrumentConfig,
    charges=None,
    max_retries: int = 3,
    delay: float = 1.0,
    chunk_size: int = 256,
) -> pd.DataFrame:
    """Predict MS2 fragment intensities for ``peptides`` under ``config``.

    Predictions are made in chunks. ``koinapy`` makes each request atomic, so a single
    peptide Koina rejects (an unsupported length or composition) would otherwise abort
    the whole run. If a chunk fails it is retried peptide by peptide and the offending
    peptides are skipped and counted. Returns the concatenated Koina frame with
    ``peptide_sequences``, ``precursor_charges``, ``intensities``, ``mz`` and
    ``annotation`` (decoded). Raises ``RuntimeError`` only if every peptide fails.
    """
    inputs = _build_input_frame(peptides, charges, config)
    model = Koina(config.intensity_model, KOINA_SERVER)

    frames = []
    dropped = 0
    for start in range(0, len(inputs), chunk_size):
        chunk = inputs.iloc[start:start + chunk_size]
        try:
            frames.append(_predict_once(model, chunk, max_retries, delay))
        except Exception:  # noqa: BLE001 - isolate the failing peptide(s), keep the rest
            for i in range(len(chunk)):
                try:
                    frames.append(_predict_once(model, chunk.iloc[[i]], 1, delay))
                except Exception:  # noqa: BLE001
                    dropped += 1

    if dropped:
        print(f"  [{config.name}] skipped {dropped} peptide(s) Koina could not predict")
    if not frames:
        raise RuntimeError(
            f"Koina intensity prediction failed for '{config.name}' "
            f"(model {config.intensity_model}): no peptides could be predicted."
        )

    predictions = pd.concat(frames, ignore_index=True)
    # Koina returns the annotation column as bytes for several models.
    if "annotation" in predictions.columns and len(predictions):
        first = predictions["annotation"].iloc[0]
        if isinstance(first, (bytes, bytearray)):
            predictions["annotation"] = predictions["annotation"].map(
                lambda x: x.decode("utf-8")
            )
    predictions["instrument"] = config.name
    return predictions
