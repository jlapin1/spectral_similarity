import numpy as np
import pandas as pd
import metrics.metrics as M

metric_keys = [
    m for m in dir(M) if ((m[:2] != "__") & (m != "binarize") & (m != "normalize"))
]
metric_keys = [
    "x_corr",
    "pearson_correlation",
    "spearman_correlation",
    "kendall_tau",
    "cosine_similarity",
    "spectral_angle",
    "weighted_dot_product",
    "fit",
    "ruzicka_similarity_1",
    "ruzicka_similarity_2",
    "mse",
    "canberra_distance",
    "wasserstein",
    "bray_curtis",
    "mutual_information",
    "hyper_score",
]


def add_gaussian_noise(arr, mean=0, std_dev=0.01):
    noise = np.random.normal(mean, std_dev, arr.shape)
    return arr + noise


def swap_two(arr):
    idx1, idx2 = np.random.choice(len(arr), 2, replace=True)
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
        selected_peptide = peptides_predictions[peptides_predictions["ID"] == i].copy()
        selected_peptide_switched = peptides_switch_predictions[
            peptides_switch_predictions["ID"] == i
        ].copy()

        # Find all common annotations between the two DataFrames
        common_annotations = set(selected_peptide["annotation"]) & set(
            selected_peptide_switched["annotation"]
        )

        if len(common_annotations) == 0:
            # No common annotations, skip this peptide
            continue

        # Sort common annotations to ensure consistent ordering
        common_annotations_sorted = sorted(common_annotations)

        # Filter both DataFrames to only include common annotations
        selected_peptide = selected_peptide[
            selected_peptide["annotation"].isin(common_annotations_sorted)
        ].copy()
        selected_peptide_switched = selected_peptide_switched[
            selected_peptide_switched["annotation"].isin(common_annotations_sorted)
        ].copy()

        # Set annotation as index and reindex to align by common annotations
        selected_peptide = selected_peptide.set_index("annotation").loc[
            common_annotations_sorted
        ]
        selected_peptide_switched = selected_peptide_switched.set_index(
            "annotation"
        ).loc[common_annotations_sorted]

        # Now both DataFrames are aligned by annotation order
        original_intensities = selected_peptide["intensities"].to_numpy()
        switched_intensities = selected_peptide_switched["intensities"].to_numpy()
        mz_values = selected_peptide["mz"].to_numpy()
        mz_values_switched = selected_peptide_switched["mz"].to_numpy()

        annotation = selected_peptide.index.to_numpy()
        annotation_switched = selected_peptide_switched.index.to_numpy()

        # Get peptide sequence info (from the first row since they're all the same ID)
        peptide_seq = selected_peptide["peptide_sequences"].iloc[0]
        peptide_seq_switched = selected_peptide_switched["peptide_sequences"].iloc[0]

        for j in range(num_randomization_rounds):
            score_dict = {}

            # Maybe this should be renamed as these are not always noisy
            noisy_intensities = original_intensities.copy()
            if randomize_gaussian:
                # Add Gaussian noise instead of swapping
                noisy_intensities = add_gaussian_noise(
                    original_intensities, mean=noise_mean, std_dev=noise_std_dev
                )
            if randomize_switched:
                noisy_intensities = original_intensities.copy()
                for _ in range(num_randomizations):
                    noisy_intensities = swap_two(noisy_intensities)
            noisy_intensities = np.clip(noisy_intensities, 0, None)

            for key in metric_keys:
                inp = {
                    "intensity1": noisy_intensities,
                    "intensity2": switched_intensities,
                    "mz1": mz_values,
                    "mz2": mz_values_switched,
                    "annotation1": annotation,
                    "annotation2": annotation_switched,
                    "diagnostic_mz": np.array([]),
                    "mz": mz_values,
                }

                try:
                    score = getattr(M, key)(**inp)
                except Exception as e:
                    score = np.nan

                score_dict[key] = score

            try:
                peptide_dict[
                    peptide_seq + "|" + peptide_seq_switched + "|" + str(j)
                ] = score_dict
            except Exception as e:
                continue

    score_df = pd.DataFrame(peptide_dict).T
    return score_df
