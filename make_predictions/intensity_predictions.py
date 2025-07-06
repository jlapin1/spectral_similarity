import numpy as np
import pandas as pd
from koinapy import Koina


def obtain_predictions_pairs(
    peptides,
    charges=[2],
    collision_energie=28,
    instrument_types="LUMOS",
    fragmentation_types="HCD",
    switched=False,
    model="UniSpec",
):
    """
    Function to obtain intensity predictions for a set of peptides.
    """
    num_peptides = peptides.shape[0]

    inputs = pd.DataFrame()
    inputs["peptide_sequences"] = np.array(peptides)
    inputs["precursor_charges"] = np.array(num_peptides * charges)
    inputs["collision_energies"] = np.array(num_peptides * [collision_energie])
    inputs["instrument_types"] = np.array(num_peptides * [instrument_types])
    inputs["fragmentation_types"] = np.array(num_peptides * [fragmentation_types])

    model = Koina(model, "koina.wilhelmlab.org:443")
    predictions = model.predict(inputs, debug = True)

    predictions["annotation"] = predictions["annotation"].map(
        lambda x: x.decode("utf-8")
    )

    predictions["non_switched"] = switched

    return predictions
