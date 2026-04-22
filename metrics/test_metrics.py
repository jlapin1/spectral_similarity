import numpy as np
from metrics.metrics import (
    spectral_angle,
    pearson_correlation,
    spearman_correlation,
    mse,
    x_corr,
    cosine_similarity,
    weighted_dot_product,
    fit,
    ruzicka_similarity_1,
    ruzicka_similarity_2,
    wasserstein,
    kendall_tau,
    mutual_information,
    bray_curtis,
    canberra_distance,
    hyper_score,
)

if __name__ == "__main__":
    # Example test spectra (aligned)
    mz = np.array([100, 150, 200, 250, 300, 350, 400, 450, 500, 550])
    intensity1 = np.array([0.2, 0.0, 0.4, 0.6, 0.1, 0.3, 0.5, 0.7, 0.2, 0.4])
    intensity2 = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.0, 0.3, 0.2])

    annotation1 = np.array(['b2+1', 'b3+1', 'b4+1', 'b4+2', 'b7+1', 'y1+1', 'y10+1', 'y10+2', 'y2+1', 'y3+1'])
    annotation2 = np.array(['b2+1', 'b3+1', 'b4+1', 'b4+2', 'b7+1', 'y1+1', 'y10+1', 'y10+2', 'y2+1', 'y3+1'])

    # Normalize intensities if required
    intensity1_norm = intensity1 / np.linalg.norm(intensity1)
    intensity2_norm = intensity2 / np.linalg.norm(intensity2)

    print("\nSpectral Similarity Metrics Results:")

    print(f"1. Spectral Angle: {spectral_angle(intensity1_norm, intensity2_norm):.4f}")
    print(f"2. Pearson Correlation: {pearson_correlation(intensity1, intensity2):.4f}")
    print(
        f"3. Spearman Correlation: {spearman_correlation(intensity1, intensity2):.4f}"
    )
    print(f"4. Mean Squared Error: {mse(intensity1, intensity2):.4f}")
    print(f"5. Sequest Score: {x_corr(intensity1, intensity2):.4f}")
    print(f"6. Cosine Similarity: {cosine_similarity(intensity1, intensity2):.4f}")
    print(
        f"7. Weighted Dot Product: {weighted_dot_product(intensity1_norm, intensity2_norm):.4f}"
    )
    print(f"8. Fit: {fit(intensity1, intensity2):.4f}")
    print(
        f"9. Ruzicka Similarity (L1 norm): {ruzicka_similarity_1(mz, intensity1_norm, mz, intensity2_norm):.4f}"
    )

    print(
        f"10. Ruzicka Similarity (L2 norm): {ruzicka_similarity_2(mz, intensity1_norm, mz, intensity2_norm):.4f}"
    )
    print(f"11. Hyper_score: {hyper_score(annotation1, intensity1, annotation2, intensity2):.4f}")

    print(
        f"12. Wasserstein Distance: {wasserstein(mz, intensity1, mz, intensity2):.4f}"
    )
    print(f"13. Kendall's Tau Correlation: {kendall_tau(intensity1, intensity2):.4f}")
    print(f"14. Mutual Information: {mutual_information(intensity1, intensity2):.4f}")
    print(f"15. Bray-Curtis Dissimilarity: {bray_curtis(intensity1, intensity2):.4f}")
    print(f"16. Canberra Distance: {canberra_distance(intensity1, intensity2):.4f}")
