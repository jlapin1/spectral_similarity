import numpy as np
import pandas as pd
import metrics.metrics as M

metric_keys = [
    m for m in dir(M) if ((m[:2] != "__") & (m != "binarize") & (m != "normalize"))
]
metric_keys = [
    "mse",
    "sequest_score",
    "andromeda_score",
    "pearson_correlation",
    "spearman_correlation",
    "dot_product",
    "mara_similarity",
    "modified_dot_product",
    "massbank_score",
    "gnps_score",
    "stein_scott_score",
    "wasserstein",
    "kendall_tau",
    "mutual_information",
    "bray_curtis",
    "canberra_distance",
    "mara_weighted_similarity",
    "diagnostic_weighted_similarity",
]


def add_gaussian_noise(arr, mean=0, std_dev=0.01):
    noise = np.random.normal(mean, std_dev, arr.shape)
    return arr + noise


def swap_two(arr):
    idx1, idx2 = np.random.choice(len(arr), 2, replace=False)
    arr[idx1], arr[idx2] = arr[idx2], arr[idx1]
    return arr


def metrics_comparison(
    peptides_predictions,
    peptides_switch_predictions,
    num_randomization_rounds=1,
    noise_mean=0,
    noise_std_dev=0.001,
    num_randomizations=1,
    randomize_gaussian=False,
    randomize_switched=False,
):
    peptide_dict = {}

    # Iterate through unique peptide IDs
    for i in set(peptides_predictions["ID"]):
        # Filter rows for the given ID
        selected_peptide = peptides_predictions[peptides_predictions["ID"] == i]
        selected_peptide_switched = peptides_switch_predictions[
            peptides_switch_predictions["ID"] == i
        ]

        # Identify the first annotation in selected_peptide that is also in selected_peptide_switched
        first_common_annotation = next(
            (
                ann
                for ann in selected_peptide["annotation"]
                if ann in set(selected_peptide_switched["annotation"])
            ),
            None,
        )

        # If a common annotation exists, filter on it; otherwise, return empty structures
        if first_common_annotation is not None:
            selected_peptide = selected_peptide[
                selected_peptide["annotation"] == first_common_annotation
            ].copy()
            selected_peptide_switched = selected_peptide_switched[
                selected_peptide_switched["annotation"] == first_common_annotation
            ].copy()
            original_intensities = selected_peptide["intensities"].to_numpy()
        else:
            selected_peptide = selected_peptide.iloc[0:0].copy()
            selected_peptide_switched = selected_peptide_switched.iloc[0:0].copy()
            original_intensities = np.array([])

        """
        selected_peptide = peptides_predictions[peptides_predictions["ID"] == i]
        selected_peptide_switched = peptides_switch_predictions[
            peptides_switch_predictions["ID"] == i
        ]

        common_annotations = set(selected_peptide["annotation"]).intersection(
            selected_peptide_switched["annotation"]
        )

        selected_peptide = selected_peptide[
            selected_peptide["annotation"].isin(common_annotations)
        ].copy()

        selected_peptide_switched = selected_peptide_switched[
            selected_peptide_switched["annotation"].isin(common_annotations)
        ].copy()

        original_intensities = selected_peptide["intensities"].to_numpy()"
        """

        for j in range(num_randomization_rounds):
            score_dict = {}

            # Maybe this should be renamed as these are not always noisy
            noisy_intensities = selected_peptide["intensities"].to_numpy()
            if randomize_gaussian:
                # Add Gaussian noise instead of swapping
                noisy_intensities = add_gaussian_noise(
                    original_intensities, mean=noise_mean, std_dev=noise_std_dev
                )
            if randomize_switched:
                noisy_intensities = selected_peptide["intensities"].to_numpy()
                for _ in range(num_randomizations):
                    try:
                        noisy_intensities = swap_two(noisy_intensities)
                    except:
                        noisy_intensities = []

            noisy_intensities = np.clip(noisy_intensities, 0, None)

            for key in metric_keys:
                inp = {
                    "intensity1": noisy_intensities,
                    "intensity2": selected_peptide_switched["intensities"].to_numpy(),
                    "mz1": selected_peptide["mz"].to_numpy(),
                    "mz2": selected_peptide_switched["mz"].to_numpy(),
                    "diagnostic_mz": np.array([]),
                    "mz": selected_peptide["mz"].to_numpy(),
                }

                try:
                    score = getattr(M, key)(**inp)
                except:
                    score = np.nan

                score_dict[key] = score

            peptide_dict[
                selected_peptide["peptide_sequences"].iloc[0]
                + "|"
                + selected_peptide_switched["peptide_sequences"].iloc[0]
                + "|"
                + str(j)
            ] = score_dict

    score_df = pd.DataFrame(peptide_dict).T
    return score_df
