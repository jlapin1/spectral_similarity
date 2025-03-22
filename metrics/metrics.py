import numpy as np
from scipy.spatial.distance import cosine, braycurtis, canberra
from scipy.stats import pearsonr, spearmanr, kendalltau, wasserstein_distance
from sklearn.metrics import mean_squared_error, mutual_info_score

# Utility functions
def normalize(intensity):
    norm = np.linalg.norm(intensity, **kwargs)
    return intensity / norm if norm > 0 else intensity

def binarize(intensity, threshold=0.01, **kwargs):
    return (intensity > threshold).astype(int)

# Spectral angle (Cosine similarity)
def spectral_angle(intensity1, intensity2, **kwargs):
    return 1 - cosine(intensity1, intensity2)

# Pearson correlation
def pearson_correlation(intensity1, intensity2, **kwargs):
    corr, _ = pearsonr(intensity1, intensity2)
    return corr

# Spearman correlation
def spearman_correlation(intensity1, intensity2, **kwargs):
    corr, _ = spearmanr(intensity1, intensity2)
    return corr

# Mean Squared Error (MSE)
def mse(intensity1, intensity2, **kwargs):
    return mean_squared_error(intensity1, intensity2)

# Sequest-scoring (sanity; returns identical input)
def sequest_score(intensity1, intensity2, **kwargs):
    return np.dot(intensity1, intensity2)

# Andromeda-scoring (sanity; returns identical input)
def andromeda_score(intensity1, intensity2, **kwargs):
    return np.dot(intensity1, intensity2)

# Dot product (same as spectral angle)
def dot_product(intensity1, intensity2, **kwargs):
    return spectral_angle(intensity1, intensity2)

# Mara Cluster similarity (simple form)
def mara_similarity(intensity1, intensity2, **kwargs):
    sim = np.sum(np.minimum(intensity1, intensity2))
    return sim / np.sum(np.maximum(intensity1, intensity2))

# Modified dot product (weighted cosine similarity)
def modified_dot_product(mz1, intensity1, mz2, intensity2, mz_weight=1, **kwargs):
    weights = (mz1 ** mz_weight) * (mz2 ** mz_weight)
    return np.sum(weights * intensity1 * intensity2) / (
        np.sqrt(np.sum(weights * intensity1**2)) * np.sqrt(np.sum(weights * intensity2**2))
    )

# MASSBANK score (simplified)
def massbank_score(intensity1, intensity2, **kwargs):
    score = np.sum(intensity1 * intensity2) / (np.sum(intensity1**2) + np.sum(intensity2**2) - np.sum(intensity1 * intensity2))
    return score

# GNPS score (simplified)
def gnps_score(intensity1, intensity2, **kwargs):
    dot = np.dot(intensity1, intensity2)
    return dot / (np.sum(intensity1**2) + np.sum(intensity2**2) - dot)

# Stein-Scott similarity score
def stein_scott_score(intensity1, intensity2, **kwargs):
    shared_peaks = intensity1 > 0
    intensity1 = intensity1[shared_peaks]
    intensity2 = intensity2[shared_peaks]
    prod = intensity1 * intensity2
    score = np.sum(prod) / (np.sqrt(np.sum(intensity1**2)) * np.sqrt(np.sum(intensity2**2)))
    return score

# Wasserstein distance
def wasserstein(mz1, intensity1, mz2, intensity2, **kwargs):
    return wasserstein_distance(mz1, mz2, intensity1, intensity2)

# Kendall's Tau Rank correlation
def kendall_tau(intensity1, intensity2, **kwargs):
    tau, _ = kendalltau(intensity1, intensity2)
    return tau

# Mutual Information
def mutual_information(intensity1, intensity2, bins=20, **kwargs):
    c_intensity1 = np.digitize(intensity1, bins=np.histogram_bin_edges(intensity1, bins))
    c_intensity2 = np.digitize(intensity2, bins=np.histogram_bin_edges(intensity2, bins))
    return mutual_info_score(c_intensity1, c_intensity2)

# Bray-Curtis dissimilarity
def bray_curtis(intensity1, intensity2, **kwargs):
    return braycurtis(intensity1, intensity2)

# Canberra distance
def canberra_distance(intensity1, intensity2, **kwargs):
    return canberra(intensity1, intensity2)

# Weighted with m/z (Mara cluster style)
def mara_weighted_similarity(mz1, intensity1, mz2, intensity2, mz_scale=1, **kwargs):
    mz_diff = np.abs(mz1 - mz2)
    weight = np.exp(-mz_diff * mz_scale)
    sim = np.sum(weight * np.minimum(intensity1, intensity2))
    return sim / np.sum(weight * np.maximum(intensity1, intensity2))

# Weighted with diagnostic ions
def diagnostic_weighted_similarity(mz, intensity1, intensity2, diagnostic_mz, diagnostic_weight=2, **kwargs):
    weights = np.ones_like(mz)
    for d_mz in diagnostic_mz:
        weights[np.abs(mz - d_mz) < 0.1] *= diagnostic_weight
    sim = np.sum(weights * intensity1 * intensity2)
    return sim / (np.sqrt(np.sum(weights * intensity1**2)) * np.sqrt(np.sum(weights * intensity2**2)))

