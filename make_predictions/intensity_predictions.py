import numpy as np
import pandas as pd
from koinapy import Koina


def obtain_predictions_pairs(
    peptides,
    charges=[2],
    collision_energie=28,
    instrument_types="LUMOS",
    switched=False,
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

    model = Koina("UniSpec", "koina.wilhelmlab.org:443")
    predictions = model.predict(inputs)

    predictions["annotation"] = predictions["annotation"].map(
        lambda x: x.decode("utf-8")
    )

    predictions["ID"] = (
        predictions.groupby(
            ["peptide_sequences", "precursor_charges"], sort=False
        ).ngroup()
        + 1
    )

    predictions["non_switched"] = switched

    return predictions
