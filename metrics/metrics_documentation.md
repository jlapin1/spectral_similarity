# Overview of different similarity metrics

Here we describe the different metrics provided in this project (and in this folder).

## Data preparation

To calculate similarity metrics between peptide-spectrum matches (PSMs), one begins by extracting the relevant spectral representations from each PSM. Each PSM consists of an associated peptide sequence and its corresponding MS/MS spectrum, typically represented as a list of fragment ion m/z values and their intensities. For metric-based comparisons, these spectra must first be converted into aligned vectors, usually by binning the m/z values into a common reference grid and associating the corresponding intensities. Once aligned, the result is a pair of intensity vectors that represent the signal distribution across m/z space for two peptide-spectrum matches. These vectors can then be compared using a range of similarity or distance metrics, which quantify how alike the spectra are in terms of ion presence, intensity patterns, or overall distribution. 

We represent the intensity vectors by $\mathbf{x} = (x_1, x_2, \dots, x_n)$ and $\mathbf{y} = (y_1, y_2, \dots, y_n)$

## Metrics

### 1. Correlation-based Metrics

#### 1.1 Pearson Correlation

The Pearson correlation coefficient is defined as:

$$
r(\mathbf{x}, \mathbf{y}) = \frac{\sum_{i=1}^n (x_i - \bar{x})(y_i - \bar{y})}{\sqrt{\sum_{i=1}^n (x_i - \bar{x})^2} \sqrt{\sum_{i=1}^n (y_i - \bar{y})^2}}
$$

where $\bar{x}$ and $\bar{y}$ are the means of vectors $\mathbf{x}$ and $\mathbf{y}$.

**Description:**  
Measures the linear correlation between two intensity vectors. A value of 1 means perfect positive linear correlation, 0 means no linear correlation, and -1 means perfect negative linear correlation.

**Reference:**  
Pearson, K. (1895). *Note on regression and inheritance in the case of two parents*. Proceedings of the Royal Society of London, 58, 240–242.

**Properties:**  
- Output range: [-1, 1]  
- Sensitive to outliers  
- Invariant to scaling and shifting of intensity vectors  
- Emphasizes linear relationships

**Parameters:**  
- `intensity1`, `intensity2`: aligned intensity vectors (typically normalized or mean-centered)

#### 1.2 Spearman Correlation

The Spearman correlation coefficient is defined as the Pearson correlation between the rank values of the two vectors:

$$
\rho(\mathbf{x}, \mathbf{y}) = \text{Pearson}(\text{rank}(\mathbf{x}), \text{rank}(\mathbf{y}))
$$

where $\text{rank}(\cdot)$ assigns ranks to the elements of the vector.

**Description:**  
Measures the monotonic relationship between two intensity vectors by comparing the ranks of their values. It is less sensitive to outliers and non-linear relationships than Pearson correlation.

**Reference:**  
Spearman, C. (1904). *The proof and measurement of association between two things*. The American Journal of Psychology, 15(1), 72–101.

**Properties:**  
- Output range: [-1, 1]  
- Sensitive to tied ranks  
- Invariant to monotonic transformations  
- Emphasizes rank order rather than absolute values

**Parameters:**  
- `intensity1`, `intensity2`: aligned intensity vectors (typically normalized or mean-centered)

#### 1.3 Kendall's Tau Rank Correlation

Kendall's Tau is defined as:

$$
\tau = \frac{(\text{number of concordant pairs}) - (\text{number of discordant pairs})}{\frac{1}{2} n(n-1)}
$$

where a pair of observations is concordant if the ranks for both elements agree, and discordant if they disagree.

**Description:**  
Measures the ordinal association between two intensity vectors by comparing the number of concordant and discordant pairs. It is robust to outliers and non-linear relationships, and is especially useful for assessing monotonic relationships.

**Reference:**  
Kendall, M. G. (1938). *A New Measure of Rank Correlation*. Biometrika, 30(1/2), 81–93.

**Properties:**  
- Output range: [-1, 1]  
- Robust to outliers  
- Invariant to monotonic transformations  
- Emphasizes rank order rather than absolute values

**Parameters:**  
- `intensity1`, `intensity2`: aligned intensity vectors (typically normalized or mean-centered)

---

### 2. Dot Product and Cosine-based Metrics

#### 2.1 Spectral Angle (also Cosine Similarity or dot product)
The spectral angle is defined as:

$$
\text{SpectralAngle}(\mathbf{x}, \mathbf{y}) = 1 - \cos(\theta) = 1 - \frac{\mathbf{x} \cdot \mathbf{y}}{\|\mathbf{x}\| \|\mathbf{y}\|},
$$

where $\cdot$ denotes the dot product and $\|\cdot\|$ the Euclidean norm.

**Description:**  
This metric measures the cosine of the angle between two vectors in high-dimensional space. In mass spectrometry, it quantifies how similar two spectra are in terms of their shape (relative distribution of intensity), regardless of absolute magnitude. A smaller spectral angle (i.e., a cosine closer to 1) indicates greater similarity.

**Reference:**  
Stein, S. E., & Scott, D. R. (1994). *Optimization and Testing of Mass Spectral Library Search Algorithms for Compound Identification*. J. Am. Soc. Mass Spectrom., 5(9), 859–866.

**Properties:**  
- Output range: [0, 1] for non-negative intensities.  
- Invariant to scaling of intensity vectors  
- Emphasizes **relative** ion intensities, not their absolute levels  
- Commonly used in spectral library searches and proteome annotation tools  
- Sensitive to vector sparsity and noise if normalization is poor

**Parameters:**  
- `intensity1`, `intensity2`: aligned intensity vectors (typically normalized); both inputs must be nonzero to avoid division by zero.  
- Optional preprocessing may include normalization (e.g., total ion current, L2 norm) or binning of spectra into common m/z grids.  

#### 2.2 Sequest/Andromeda Score

The Sequest score is defined as the dot product between two vectors:

$$
\text{SequestScore}(\mathbf{x}, \mathbf{y}) = \sum_{i=1}^n x_i y_i
$$

**Description:**  
Measures the overall similarity between two spectra by summing the products of their corresponding intensities. Used in the Sequest search algorithm for peptide identification.

**Reference:**  
Eng, J. K., McCormack, A. L., & Yates, J. R. (1994). *An approach to correlate tandem mass spectral data of peptides with amino acid sequences in a protein database*. J. Am. Soc. Mass Spectrom., 5(11), 976–989.

**Properties:**  
- Output range: $[0, \infty)$ for non-negative intensities  
- Sensitive to absolute intensity scale  
- Not invariant to scaling  
- Simple and fast to compute

**Parameters:**  
- `intensity1`, `intensity2`: aligned intensity vectors (typically normalized or mean-centered)

#### 2.3 Modified Dot Product (Weighted Cosine Similarity)

The modified dot product introduces m/z-dependent weights to the cosine similarity:

$$
\text{ModifiedDotProduct}(\mathbf{mz}_1, \mathbf{x}, \mathbf{mz}_2, \mathbf{y}, \mathbf{w}) = 
\frac{\sum_{i=1}^n (mz_{1,i}^{w} \cdot mz_{2,i}^{w}) \cdot x_i y_i}
{\sqrt{\sum_{i=1}^n (mz_{1,i}^{w} \cdot x_i^2)} \cdot \sqrt{\sum_{i=1}^n (mz_{2,i}^{w} \cdot y_i^2)}}
$$

where $w$ is the m/z weighting exponent (default: $w=1$).

**Description:**  
This metric weights the contribution of each peak to the similarity score by its m/z value, emphasizing higher-mass fragments if desired. It is a generalization of the cosine similarity.

**Reference:**  
Horai, H., et al. (2010). *MassBank: a public repository for sharing mass spectral data for life sciences*. J. Mass Spectrom., 45(7), 703–714.

**Properties:**  
- Output range: [0, 1] for non-negative intensities  
- Allows tuning of m/z influence via the `mz_weight` parameter  
- Invariant to scaling of intensity vectors  
- Emphasizes peaks with higher m/z if `mz_weight > 0`

**Parameters:**  
- `mz1`, `mz2`: m/z vectors for each spectrum  
- `intensity1`, `intensity2`: aligned intensity vectors  
- `mz_weight`: exponent for m/z weighting (default: 1)

#### 2.4 Diagnostic Weighted Similarity

This metric is defined as:

$$
\text{DiagnosticWeighted}(\mathbf{mz}, \mathbf{x}, \mathbf{y}, D, w) = 
\frac{\sum_{i=1}^n \omega_i x_i y_i}
{\sqrt{\sum_{i=1}^n \omega_i x_i^2} \cdot \sqrt{\sum_{i=1}^n \omega_i y_i^2}}
$$

where $\omega_i$ is a weight that is increased by `diagnostic_weight` for each m/z value within 0.1 of a diagnostic ion in $D$.

**Description:**  
A cosine-like similarity that increases the contribution of peaks near specified diagnostic m/z values. Useful for emphasizing ions of particular biological or analytical interest.

**Reference:**  
Adapted from approaches in targeted mass spectrometry and diagnostic ion scoring.

**Properties:**  
- Output range: [0, 1] for non-negative intensities  
- Emphasizes peaks near diagnostic m/z values  
- Invariant to scaling  
- Sensitive to the choice of diagnostic ions and weighting

**Parameters:**  
- `mz`: m/z vector for the spectra  
- `intensity1`, `intensity2`: aligned intensity vectors  
- `diagnostic_mz`: list of diagnostic m/z values  
- `diagnostic_weight`: multiplicative weight for diagnostic ions (default: 2)

---

### 3. Overlap/Shared Peak Metrics

#### 3.1 Mara Cluster Similarity

The Mara cluster similarity is defined as:

$$
\text{MaraSimilarity}(\mathbf{x}, \mathbf{y}) = \frac{\sum_{i=1}^n \min(x_i, y_i)}{\sum_{i=1}^n \max(x_i, y_i)}
$$

**Description:**  
Measures the overlap between two spectra by comparing the sum of minimum intensities to the sum of maximum intensities at each position. Values closer to 1 indicate greater similarity.

**Reference:**  
Horai, H., et al. (2010). *MassBank: a public repository for sharing mass spectral data for life sciences*. J. Mass Spectrom., 45(7), 703–714.

**Properties:**  
- Output range: [0, 1]  
- Invariant to scaling  
- Emphasizes shared peaks  
- Sensitive to missing or extra peaks

**Parameters:**  
- `intensity1`, `intensity2`: aligned intensity vectors (typically normalized or mean-centered)

#### 3.2 Stein-Scott Similarity Score

The Stein-Scott similarity score is defined as:

$$
\text{SteinScott}(\mathbf{x}, \mathbf{y}) = \frac{\sum_{i \in S} x_i y_i}{\sqrt{\sum_{i \in S} x_i^2} \cdot \sqrt{\sum_{i \in S} y_i^2}}
$$

where $S$ is the set of indices where $x_i > 0$ (shared peaks).

**Description:**  
Measures the similarity between two spectra by considering only the peaks present in the first spectrum. It is a normalized dot product restricted to shared peaks, as used in the Stein-Scott algorithm for spectral library searching.

**Reference:**  
Stein, S. E., & Scott, D. R. (1994). *Optimization and Testing of Mass Spectral Library Search Algorithms for Compound Identification*. J. Am. Soc. Mass Spectrom., 5(9), 859–866.

**Properties:**  
- Output range: [0, 1] for non-negative intensities  
- Focuses on peaks present in the first spectrum  
- Invariant to scaling  
- Emphasizes shared, high-intensity peaks

**Parameters:**  
- `intensity1`, `intensity2`: aligned intensity vectors (typically normalized or mean-centered)


#### 3.3 Mara Weighted Similarity (m/z-weighted Mara cluster style)

This metric is defined as:

$$
\text{MaraWeighted}(\mathbf{mz}_1, \mathbf{x}, \mathbf{mz}_2, \mathbf{y}, \mathbf{mz\_scale}) = 
\frac{\sum_{i=1}^n w_i \cdot \min(x_i, y_i)}{\sum_{i=1}^n w_i \cdot \max(x_i, y_i)}
$$

where $w_i = \exp(-|\text{mz}_{1,i} - \text{mz}_{2,i}| \cdot \text{mz\_scale})$.

**Description:**  
A variant of the Mara cluster similarity that applies an exponential weight based on the m/z difference between corresponding peaks. Peaks with closer m/z values contribute more to the similarity score.

**Reference:**  
Horai, H., et al. (2010). *MassBank: a public repository for sharing mass spectral data for life sciences*. J. Mass Spectrom., 45(7), 703–714.

**Properties:**  
- Output range: [0, 1]  
- Emphasizes peaks with similar m/z  
- Invariant to scaling  
- Sensitive to the `mz_scale` parameter

**Parameters:**  
- `mz1`, `mz2`: m/z vectors for each spectrum  
- `intensity1`, `intensity2`: aligned intensity vectors  
- `mz_scale`: scaling factor for m/z weighting (default: 1)



---

### 4. Distance-based Metrics

#### 4.1 Mean Squared Error (MSE)

The mean squared error between two vectors is defined as:

$$
\text{MSE}(\mathbf{x}, \mathbf{y}) = \frac{1}{n} \sum_{i=1}^n (x_i - y_i)^2
$$

**Description:**  
Measures the average squared difference between corresponding elements of two intensity vectors. Lower values indicate greater similarity.

**Reference:**  
Willmott, C. J., & Matsuura, K. (2005). *Advantages of the mean absolute error (MAE) over the root mean square error (RMSE) in assessing average model performance*. Climate Research, 30(1), 79–82.

**Properties:**  
- Output range: $[0, \infty)$  
- Sensitive to outliers  
- Penalizes large differences more heavily  
- Not invariant to scaling or shifting

**Parameters:**  
- `intensity1`, `intensity2`: aligned intensity vectors (typically normalized or mean-centered)

#### 4.2 Bray-Curtis Dissimilarity

The Bray-Curtis dissimilarity is defined as:

$$
\text{BrayCurtis}(\mathbf{x}, \mathbf{y}) = \frac{\sum_{i=1}^n |x_i - y_i|}{\sum_{i=1}^n (x_i + y_i)}
$$

**Description:**  
Measures the dissimilarity between two intensity vectors by comparing the sum of absolute differences to the sum of all intensities. Values closer to 0 indicate greater similarity, while values closer to 1 indicate greater dissimilarity.

**Reference:**  
Bray, J. R., & Curtis, J. T. (1957). *An Ordination of the Upland Forest Communities of Southern Wisconsin*. Ecological Monographs, 27(4), 325–349.

**Properties:**  
- Output range: [0, 1]  
- Sensitive to differences in both presence/absence and abundance  
- Not invariant to scaling  
- Emphasizes relative differences

**Parameters:**  
- `intensity1`, `intensity2`: aligned intensity vectors (typically normalized or mean-centered)

#### 4.3 Canberra Distance

The Canberra distance is defined as:

$$
\text{Canberra}(\mathbf{x}, \mathbf{y}) = \sum_{i=1}^n \frac{|x_i - y_i|}{|x_i| + |y_i|}
$$

**Description:**  
Measures the distance between two intensity vectors by summing the ratio of absolute differences to the sum of absolute values for each element. It is more sensitive to small changes near zero than other metrics.

**Reference:**  
Lance, G. N., & Williams, W. T. (1966). *A General Theory of Classificatory Sorting Strategies: 1. Hierarchical Systems*. The Computer Journal, 9(4), 373–380.

**Properties:**  
- Output range: $[0, \infty)$  
- Sensitive to small values and zeros  
- Not invariant to scaling  
- Emphasizes proportional differences

**Parameters:**  
- `intensity1`, `intensity2`: aligned intensity vectors (typically normalized or mean-centered)

#### 4.4 Wasserstein Distance

The Wasserstein (or Earth Mover's) distance between two distributions is defined as:

$$
W(\mu, \nu) = \inf_{\gamma \in \Gamma(\mu, \nu)} \int |x - y| \, d\gamma(x, y)
$$

where $\mu$ and $\nu$ are the distributions (here, spectra), and $\Gamma(\mu, \nu)$ is the set of all joint distributions with marginals $\mu$ and $\nu$.

**Description:**  
Measures the minimum "cost" required to transform one spectrum into another, treating the spectra as distributions over m/z. Sensitive to both the location and magnitude of peaks.

**Reference:**  
Rubner, Y., Tomasi, C., & Guibas, L. J. (2000). *The Earth Mover's Distance as a Metric for Image Retrieval*. International Journal of Computer Vision, 40(2), 99–121.

**Properties:**  
- Output range: $[0, \infty)$  
- Sensitive to shifts in m/z  
- Takes both intensity and m/z into account  
- Useful for comparing spectra with different peak locations

**Parameters:**  
- `mz1`, `intensity1`: m/z and intensity vectors for the first spectrum  
- `mz2`, `intensity2`: m/z and intensity vectors for the second spectrum

---

### 5. Other Metrics

#### 5.1 MASSBANK Score (simplified, also GNPS score)

The MASSBANK score is defined as:

$$
\text{MASSBANK}(\mathbf{x}, \mathbf{y}) = \frac{\sum_{i=1}^n x_i y_i}{\sum_{i=1}^n x_i^2 + \sum_{i=1}^n y_i^2 - \sum_{i=1}^n x_i y_i}
$$

**Description:**  
Measures the similarity between two spectra by normalizing the dot product with respect to the sum of squared intensities and their overlap. This metric is used in the MassBank database for spectral matching.

**Reference:**  
Horai, H., et al. (2010). *MassBank: a public repository for sharing mass spectral data for life sciences*. J. Mass Spectrom., 45(7), 703–714.

**Properties:**  
- Output range: [0, 1] for non-negative intensities  
- Invariant to scaling  
- Emphasizes shared peaks  
- Sensitive to missing or extra peaks

**Parameters:**  
- `intensity1`, `intensity2`: aligned intensity vectors (typically normalized or mean-centered)


#### 5.2 Mutual Information

Mutual information is defined as:

$$
I(X; Y) = \sum_{x \in X} \sum_{y \in Y} p(x, y) \log \left( \frac{p(x, y)}{p(x)p(y)} \right)
$$

where $p(x, y)$ is the joint probability distribution of the binned intensities, and $p(x)$, $p(y)$ are the marginal distributions.

**Description:**  
Measures the amount of information shared between two intensity vectors. It captures both linear and non-linear dependencies by quantifying how much knowing the intensity in one spectrum reduces uncertainty about the other.

**Reference:**  
Cover, T. M., & Thomas, J. A. (2006). *Elements of Information Theory* (2nd ed.). Wiley-Interscience.

**Properties:**  
- Output range: $[0, \infty)$  
- Sensitive to binning and normalization  
- Captures both linear and non-linear relationships  
- Not invariant to scaling or shifting

**Parameters:**  
- `intensity1`, `intensity2`: aligned intensity vectors  
- `bins`: number of bins for discretizing intensities (default: 20)


