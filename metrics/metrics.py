import numpy as np
from scipy.spatial.distance import cosine, braycurtis, canberra
from scipy.stats import pearsonr, spearmanr, kendalltau, wasserstein_distance
from sklearn.metrics import mean_squared_error
from math import factorial


# Utility functions
def normalize(intensity, **kwargs):
    norm = np.linalg.norm(intensity, **kwargs)
    return intensity / norm if norm > 0 else intensity


def binarize(intensity, threshold=0.01, **kwargs):
    return (intensity > threshold).astype(int)


def x_corr(intensity1, intensity2, **kwargs):
    tau_max = 75
    n = len(intensity1)
    y_padded = np.pad(
        intensity2, (tau_max, tau_max), mode="constant", constant_values=0
    )

    R = np.zeros(2 * tau_max + 1)
    taus = np.arange(-tau_max, tau_max + 1)

    for k, tau in enumerate(taus):
        y_shifted = y_padded[tau_max + tau : tau_max + tau + n]
        R[k] = np.sum(intensity1 * y_shifted)

    return R[tau_max] - np.mean(R)


def pearson_correlation(intensity1, intensity2, **kwargs):
    corr, _ = pearsonr(intensity1, intensity2)
    return corr


def spearman_correlation(intensity1, intensity2, **kwargs):
    corr, _ = spearmanr(intensity1, intensity2)
    return corr


def kendall_tau(intensity1, intensity2, **kwargs):
    tau, _ = kendalltau(intensity1, intensity2)
    return tau


def cosine_similarity(intensity1, intensity2, **kwargs):
    return cosine(intensity1, intensity2)


def spectral_angle(intensity1, intensity2, **kwargs):
    cos_theta = cosine(intensity1, intensity2)
    angle = np.arccos(cos_theta)
    return 1 - (2 * angle / np.pi)


def weighted_dot_product(mz1, intensity1, mz2, intensity2, k=1, m=1, **kwargs):
    w1 = (mz1**k) * (intensity1**m)
    w2 = (mz2**k) * (intensity2**m)
    return (np.sum(w1 * w2)) ** 2 / (np.sum(w1**2) * np.sum(w2**2))


def fit(intensity1, intensity2, **kwargs):
    mask = intensity2 != 0
    return (np.sum(intensity1 * intensity2)) ** 2 / (
        np.sum(intensity1[mask] ** 2) * np.sum(intensity2**2)
    )


def ruzicka_similarity_1(intensity1, intensity2, **kwargs):
    # Test L1 norm + Ruzicka
    intensity1 / np.linalg.norm(intensity1, 1)
    intensity2 / np.linalg.norm(intensity2, 1)
    return np.sum(np.minimum(intensity1, intensity2)) / np.sum(
        np.maximum(intensity1, intensity2)
    )

def ruzicka_similarity_2(intensity1, intensity2, **kwargs):
    # Test L2 norm + Ruzicka
    intensity1 / np.linalg.norm(intensity1, 2)
    intensity2 / np.linalg.norm(intensity2, 2)
    return np.sum(np.minimum(intensity1, intensity2)) / np.sum(
        np.maximum(intensity1, intensity2)
    )


def mse(intensity1, intensity2, **kwargs):
    return mean_squared_error(intensity1, intensity2)


def canberra_distance(intensity1, intensity2, **kwargs):
    return canberra(intensity1, intensity2)


def wasserstein(mz1, intensity1, mz2, intensity2, **kwargs):
    return wasserstein_distance(mz1, mz2, intensity1, intensity2)


def bray_curtis(intensity1, intensity2, **kwargs):
    return braycurtis(intensity1, intensity2)


def mutual_information(intensity1, intensity2, **kwargs):
    eps = 1e-15
    px = intensity1 / np.linalg.norm(intensity1, 1)
    py = intensity2 / np.linalg.norm(intensity2, 1)

    P = np.diag(px * py)  # the joint
    P /= P.sum()

    # marginal
    Pi = P.sum(axis=1)
    Pj = P.sum(axis=0)

    return np.sum(P * np.log((P + eps) / (Pi[:, None] * Pj[None, :] + eps)))


def hyper_score(annotation1, intensity1, annotation2, intensity2, **kwargs):
    b_ions = {ion.split("+")[0][1:] for ion in annotation1 if ion.startswith("b")}
    y_ions = {ion.split("+")[0][1:] for ion in annotation1 if ion.startswith("y")}

    b_ions_switched = {ion.split("+")[0][1:] for ion in annotation2 if ion.startswith("b")}
    y_ions_switched = {ion.split("+")[0][1:] for ion in annotation2 if ion.startswith("y")}

    nb = b_ions.intersection(b_ions_switched)
    ny = y_ions.intersection(y_ions_switched)

    dot = np.dot(intensity1, intensity2)

    return dot * factorial(nb) * factorial(ny)
