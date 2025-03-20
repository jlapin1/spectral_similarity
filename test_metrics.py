import numpy as np
from metrics import (
    spectral_angle,
    pearson_correlation,
    spearman_correlation,
    mse,
    sequest_score,
    andromeda_score,
    dot_product,
    mara_similarity,
    modified_dot_product,
    massbank_score,
    gnps_score,
    stein_scott_score,
    wasserstein,
    kendall_tau,
    mutual_information,
    bray_curtis,
    canberra_distance,
    mara_weighted_similarity,
    diagnostic_weighted_similarity
)

if __name__ == "__main__":
    # Example test spectra (aligned)
    mz = np.array([100, 150, 200, 250, 300, 350, 400, 450, 500, 550])
    intensity1 = np.array([0.2, 0.0, 0.4, 0.6, 0.1, 0.3, 0.5, 0.7, 0.2, 0.4])
    intensity2 = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.0, 0.3, 0.2])

    # Normalize intensities if required
    intensity1_norm = intensity1 / np.linalg.norm(intensity1)
    intensity2_norm = intensity2 / np.linalg.norm(intensity2)

    print("\nSpectral Similarity Metrics Results:")
    print("="*50)

    print(f"1. Spectral Angle (Cosine Similarity): {spectral_angle(intensity1_norm, intensity2_norm):.4f}")
    print(f"2. Pearson Correlation: {pearson_correlation(intensity1, intensity2):.4f}")
    print(f"3. Spearman Correlation: {spearman_correlation(intensity1, intensity2):.4f}")
    print(f"4. Mean Squared Error: {mse(intensity1, intensity2):.4f}")
    print(f"5. Sequest Score: {sequest_score(intensity1, intensity2):.4f}")
    print(f"6. Andromeda Score: {andromeda_score(intensity1, intensity2):.4f}")
    print(f"7. Dot Product: {dot_product(intensity1_norm, intensity2_norm):.4f}")
    print(f"8. Mara Cluster Similarity: {mara_similarity(intensity1, intensity2):.4f}")
    print(f"9. Modified Dot Product: {modified_dot_product(mz, intensity1_norm, mz, intensity2_norm):.4f}")
    print(f"10. MASSBANK Score: {massbank_score(intensity1, intensity2):.4f}")
    print(f"11. GNPS Score: {gnps_score(intensity1, intensity2):.4f}")
    print(f"12. Stein-Scott Similarity Score: {stein_scott_score(intensity1, intensity2):.4f}")
    print(f"13. Wasserstein Distance: {wasserstein(mz, intensity1, mz, intensity2):.4f}")
    print(f"14. Kendall's Tau Correlation: {kendall_tau(intensity1, intensity2):.4f}")
    print(f"15. Mutual Information: {mutual_information(intensity1, intensity2):.4f}")
    print(f"16. Bray-Curtis Dissimilarity: {bray_curtis(intensity1, intensity2):.4f}")
    print(f"17. Canberra Distance: {canberra_distance(intensity1, intensity2):.4f}")
    print(f"18. Mara Weighted Similarity (mz-weighted): {mara_weighted_similarity(mz, intensity1, mz, intensity2):.4f}")