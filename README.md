# spectral_similarity

## Install

```shell
conda env create -f environment.yml
conda activate lorentz_center_spectrum_similarity
```


## Development

* Python dependencies go into `requirements.txt`
* Libraries and tools go into `environment.yml`


## Ambigous peptides
URL: https://ruhr-uni-bochum.sciebo.de/s/wpGoJpJJwEQNrfq

### Columns
| `sequence`                         | Seqeuence |
| `ambigous_sequence`                | Ambiguous sequence |
| `sequence_raw_files`               | Comma separated list of raw files + scan number (separarted by colon) where the seqeunce was identified |
| `ambiguous_sequence_raw_files`     | Same as `sequence_raw_files` but for ambigous_sequence  |
| `sequence_mz`                      | m/z of the first raw file in `sequence_raw_files` |
| `sequence_intensity`               | intensities of the first raw file in `sequence_raw_files` |
| `ambiguous_sequence_mz`            | m/z of the first raw file in `ambiguous_sequence_raw_files` |
| `ambiguous_sequence_intensity`     | intwensities of the first raw file in `ambiguous_sequence_raw_files` |
