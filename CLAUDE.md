# CLAUDE.md — spectral_similarity

Guidance for Claude Code (and developers) working in this repository. This file
focuses on the architecture, the data contracts between stages, and **how to extend
the project to new datasets and mass-spectrometry instruments**, which is the
current development goal.

> A separate `CODEBASE_NOTES.md` exists with an earlier module-level overview. This
> file supersedes it and corrects two stale points (see *Known bugs and drift*).

---

## 1. What the project does

The project studies the **spectral similarity of I/L sibling peptides**: peptide
pairs that differ only by Isoleucine (I) ↔ Leucine (L) substitutions. I and L are
isobaric, so they are nearly indistinguishable by MS. The central question is how
distinguishable such siblings are in practice, using:

1. **Predicted** spectra (Koina-hosted models: Prosit / AlphaPeptDeep / UniSpec for
   MS2 intensities, DeepLC for retention time, AlphaPeptDeep for CCS), and
2. **Observed** spectra (MaxQuant identifications plus raw peaks from mzML or MGF).

The output figures quantify per-metric "ratio correct" (does the true sequence score
higher than its I/L-swapped decoy) for fragment intensity, retention time, and CCS.

## 2. Repository layout

```
spectral_similarity/
├── seq_utils/          # peptide sequence handling (digest, I/L switching, ProForma)
├── metrics/            # spectral-similarity metric functions + comparison driver
├── make_predictions/   # Koina API client for RT / CCS / fragment-intensity prediction
├── find_siblings/      # FASTA digestion + I/L sibling-group discovery (CLI)
├── ambiguity_search/   # observed-data ingestion: MaxQuant msms.txt + mzML → sibling pairs
├── fasta/              # human proteome FASTA (UniProt UP000005640)
├── notebooks/          # figure pipelines + exploratory analysis (most pipeline glue lives here)
├── temp_data/          # intermediate + input data (CSV, parquet, MGF, MaxQuant TSV)
└── manuscript/         # draft manuscript and figures
```

## 3. Environment

- Conda env `lorentz_center_spectrum_similarity` (`environment.yml`, Python 3.12).
- `pip` deps in `requirements.txt`: `koinapy`, `pyteomics`, `biopython`, `scipy`,
  `scikit-learn`, `pandas`, `numpy`, `pyarrow`, `matplotlib`, `seaborn`, `tqdm`.
- Notebooks additionally import `deeplc` and `psm_utils` (and one uses `rustyms`);
  these are **not** in `requirements.txt` and must be installed separately.
- Setup: `conda env create -f environment.yml && conda activate lorentz_center_spectrum_similarity`.

## 4. End-to-end data flow

```
                         ┌─────────────── PREDICTED branch ───────────────┐
 fasta/UP000005640_9606  → create_tryptic_peptides → remove_non_il/ux
   (or MaxQuant peptides)   → sample (seed 42, n=20000) → switch_random_il
                            → Koina: intensity (Prosit/AlphaPeptDeep/UniSpec)
                            →        RT (DeepLC), CCS (AlphaPeptDeep)
                            → CSV caches in temp_data/

                         ┌─────────────── OBSERVED branch ────────────────┐
 (a) ambiguity_search/maxquant.py: msms.txt + mzML → matched I/L pairs → PXD004732.parquet
 (b) notebooks figure3-6:          input_file.txt (MaxQuant) + *.mgf → annotate b/y ions

 metrics.metrics_comparison(predicted_orig, predicted_switched) → 16 metrics per pair
 notebooks: observed-vs-predicted overlay, ratio-correct bar charts, flanking heatmaps
```

There are **two observed-data ingestion paths** and they are not unified:
`ambiguity_search/maxquant.py` reads `msms.txt` + **mzML**, while the figure notebooks
read `input_file.txt` (a MaxQuant evidence/msms export) + **MGF**.

## 5. Module reference (concise)

### `seq_utils/peptide.py`
- `remove_non_il`, `remove_ux_containing`: peptide-list filters.
- `switch_first_il(pep)`: switch the first I/L.
- `switch_random_il(pep)`: switch one random I/L outside `[...]` brackets. **Uses the
  global `random` module**, not a seeded NumPy RNG — this is the root reproducibility
  hazard (see §9).
- `has_il_outside_brackets(pep)`: ProForma-aware I/L presence check.
- `get_proforma_bracketed(...)`: converts modified sequences to ProForma. Default
  `modification_dict` is **mass-delta keyed** (`+57.0215`→Carbamidomethyl, etc.), but
  the notebooks override it with **name-keyed** dicts (`(ox)`→`UNIMOD:35`). These
  conventions differ and must be matched per dataset.

### `seq_utils/fasta_to_peptides.py`
- `tryptic_digest`, `create_tryptic_peptides(fasta, min_length=7, max_length=30)`.
  **Trypsin only; length 7–30.**

### `metrics/metrics.py`
16 metric functions, all accepting `intensity1`, `intensity2`, and optionally `mz1/mz2`,
`annotation1/annotation2`. Each takes `**kwargs` so `metrics_comparison` can call them
uniformly. Metrics: `x_corr`, `pearson_correlation`, `spearman_correlation`,
`kendall_tau`, `cosine_similarity`, `spectral_angle`, `weighted_dot_product`, `fit`,
`ruzicka_similarity_1` (L1), `ruzicka_similarity_2` (L2), `mse`, `canberra_distance`,
`wasserstein`, `bray_curtis`, `mutual_information`, `hyper_score`. See §9 for caveats
(`cosine_similarity` returns a *distance*; `hyper_score` is broken).

### `metrics/get_metrics.py`
- `metric_keys`: hardcoded ordered list of the 16 metrics above (a dynamic `dir(M)`
  version is computed then immediately overwritten by the literal list).
- `metrics_comparison(peptides_predictions, peptides_switch_predictions, ...)`: the core
  driver. Input contract (both DataFrames): columns **`ID`, `annotation`, `intensities`,
  `mz`, `peptide_sequences`**. It groups by `ID`, aligns the two spectra on the
  intersection of `annotation` labels (sorted), and computes every metric. Output: a
  DataFrame indexed by `"seq|switched_seq|round"`, one column per metric. Optional
  Gaussian-noise / intensity-swap randomization for null distributions.
  **Note:** the prediction CSVs in `temp_data/` do **not** contain an `ID` column;
  callers must add it (figure1/2 build `ID = sequence + charge`). `theoretical_dist.ipynb`
  fails with `KeyError: 'ID'` because it skips this step.

### `make_predictions/intensity_predictions.py`
Koina client (server hardcoded `koina.wilhelmlab.org:443`). This is the **most
instrument-sensitive module**.
- `obtain_predictions_pairs(peptides, charges=[2], collision_energies=28,
  instrument_types="LUMOS", fragmentation_types="HCD", switched=False, model="UniSpec")`
  — fragment-intensity prediction. **All instrument parameters are defaulted here.**
- `obtain_ccs_predictions(..., model="AlphaPept_ccs_generic")`, charges default 2.
- `obtain_rt_predictions(..., model="Deeplc_hela_hf")`; alternatives `Prosit_2019_irt`,
  `AlphaPept_rt_generic`.
- `safe_obtain_*` wrappers add retry logic (3 attempts, 1 s delay). `safe_obtain_predictions`
  forces `charges = [2]*N` and defaults `model="Prosit_2020_intensity_HCD"`.

### `ambiguity_search/maxquant.py` — `MaxQuantAmbiguitySearch`
Observed-data ingestion (CLI: `python -m ambiguity_search.maxquant <out> <mzml_folder> <mq_folder>...`).
- Reads `msms.txt` (hardcoded filename) with `usecols=["Sequence","Raw file","Score","Scan number"]`.
- Builds a peptide → `{rawfile:scan}` index, finds I/L-swapped counterparts (`FIRST_N=1`,
  only the first I/L position; `AMBIGUOUS_AA={I,L}`).
- Loads peaks from **mzML only** (`raw_file + ".mzML"`, pyteomics, Thermo spectrum-id
  format `controllerType=0 controllerNumber=1 scan={n}`).
- Output columns: `sequence`, `ambigous_sequence` (sic), `sequence_raw_files`,
  `ambiguous_sequence_raw_files`, `sequence_mz`, `sequence_intensity`,
  `ambiguous_sequence_mz`, `ambiguous_sequence_intensity`. `ambiguity_search/__init__.py`
  is empty.

### `find_siblings/`
- `digest_find_siblings.py` (CLI): trypsin digest, group peptides by length and by
  I/L-normalized key (I and L → J), report groups with >1 member. **Length filter 6–60
  here**, inconsistent with `fasta_to_peptides.py` (7–30).
- `siblings_in_uniprot_proteomes.py`: downloads reference proteomes from the ExPASy FTP
  and runs the sibling search per proteome. Currently sliced to `eukaryota[:5]` for testing.
- `plot_pred_rt.py`: plotting helper.

## 6. Canonical data contracts (`temp_data/`)

Conform any new dataset to these schemas. **Long format** = one row per fragment ion.

**Fragment-intensity prediction CSV** (`peptides_predictions*.csv`, `*_puyvelde.csv`), long:
`peptide_sequences, precursor_charges, collision_energies, instrument_types,
fragmentation_types, intensities, mz, annotation, non_switched`
(e.g. `LILKPRPHVQ, 2, 28, LUMOS, HCD, 0.186, 147.076, y1+1, True`). `annotation` parsed by
regex `([a-zA-Z]+\d+)\+(\d+)` → `(ion, charge)`. `non_switched` flags original (True) vs
I/L-swapped (False).

**RT prediction CSV** (`*_rt.csv`): `peptide_sequences, irt, non_switched`.
**CCS prediction CSV** (`*_ccs.csv`, `ccs_predictions*.csv`): `peptide_sequences, precursor_charges, ccs, non_switched`.

**`PXD004732.parquet` / `.tsv`** (matched observed pairs, 8 cols):
`sequence, ambigous_sequence, sequence_raw_files, ambiguous_sequence_raw_files,
sequence_mz, sequence_intensity, ambiguous_sequence_mz, ambiguous_sequence_intensity`.
The `*_mz`/`*_intensity` are list columns (parquet) or text (tsv). The misspelling
`ambigous_sequence` is load-bearing in this table.

**`similarity_results_recomputed.parquet`** (20 cols): `sequence1, sequence2,
sequence_raw_files, ambiguous_sequence_raw_files, n_matched_peaks` + 15 metric columns.

**`similarity_results_only_exact_match.parquet`** (34 cols): ProForma/UNIMOD `sequence1/2`
with full per-pair metadata (`raw_file*`, `scan*`, `ce*`, `charge*`, `method*`) and an
**older, different** metric vocabulary (`dot_product`, `mara_similarity`, `massbank_score`,
`gnps_score`, `stein_scott_score`, `sequest_score`, `andromeda_score`, …) — generated by a
prior version of `metrics.py` (see drift note in §9).

**`input_file.txt`** (MaxQuant evidence/msms TSV). Columns consumed: `Sequence`,
`Modified sequence`, `Charge`, `MS/MS scan number`, `Raw file`, `Calibrated retention time`.

**MGF** (`LFQ_Orbitrap_DDA_Condition_{A,B}_Sample_Alpha_{01,02,03}.mgf`): per spectrum
`TITLE` (`...scan=N`), `SCANS`, `RTINSECONDS`, `PEPMASS`, `CHARGE`. **Instrument metadata
(Orbitrap, DDA, condition, replicate) is encoded only in the filename, not the header.**
Lookup key is `mgf_dict[Raw file][MS/MS scan number]`, so the MaxQuant `Raw file` value must
equal the MGF basename and `MS/MS scan number` must equal the parsed scan.

## 7. Instrument / dataset coupling inventory

Every place a specific instrument, dataset, model, charge, CE, organism, or path is baked
in. **This is the checklist to change when adding a dataset or instrument.**

| Location | Hardcoded assumption | Why it matters for a new instrument/dataset |
|---|---|---|
| `make_predictions/intensity_predictions.py:200-208` | `collision_energies=28, instrument_types="LUMOS", fragmentation_types="HCD", model="UniSpec"`, `charges=[2]` | Intensity model + CE + fragmentation must match the instrument. timsTOF/CID/EThcD need a different model and CE. |
| `…:9, 27` (`safe_obtain_predictions`) | model `Prosit_2020_intensity_HCD`, charges forced to 2 | HCD-specific; ignores real precursor charge. |
| `…:46, 88` | CCS model `AlphaPept_ccs_generic`, RT model `Deeplc_hela_hf` | RT model is HeLa-specific; CCS only meaningful for ion-mobility instruments. |
| `…:146,188,225` | Koina server `koina.wilhelmlab.org:443` | Should be configurable / env-var for self-hosted Koina. |
| `ambiguity_search/maxquant.py:58,66` | filename `msms.txt`, `usecols` MaxQuant column names | Other search engines (FragPipe/Sage/Spectronaut) use different filenames/columns. |
| `…:144,154-156` | mzML only; Thermo spectrum-id `controllerType=0 …scan=N` | Bruker/Sciex raw or MGF-only datasets won't match this id scheme. |
| notebooks figure3-6 cell-1 | `df["Charge"] == 2`, `Sequence` length `< 30` | Drops all non-2+ precursors; higher-charge species lost. |
| notebooks figure3-6 cell-3 | `model="AlphaPeptDeep_ms2_generic", collision_energies=27` (figure3-5); `Prosit_2024_intensity_PTMs_gl` (fragment_intensity); `model="UniSpec", instrument='LUMOS'` (notebook1) | Three different intensity models / CEs across notebooks; no single source of truth. |
| notebooks `get_theoretical_fragments` / `annotate_spectrum` | b/y ions only, `max_charge=2`, `tol_ppm=20` | ETD/ECD need c/z ions; ion-trap/TOF may need different tolerance; high charges excluded. |
| figure4 CCS cells | `ccs_charges = [2]*N` | CCS compared predicted-vs-predicted only (no observed mobility in this Orbitrap dataset). |
| figure1/2 sampling | `np.random.default_rng(42)`, `sample_size=20000`; figure2 *replays* this to rebuild the original↔switched map | Fragile: depends on byte-identical FASTA order + RNG + `switch_random_il`; see §9. |
| `seq_utils/fasta_to_peptides.py:8,14,27` | trypsin, length 7–30, FASTA `UP000005640_9606` | Human + trypsin assumed; other enzymes/organisms need parameters. |
| `find_siblings/digest_find_siblings.py:9-10` | length 6–60 | Inconsistent with the 7–30 used elsewhere. |
| All notebooks | relative paths `../temp_data/...`, dataset names `PXD004732`, `puyvelde`, `LFQ_Orbitrap_DDA_*` | No config layer; paths and dataset identity are literals scattered across cells. |

## 8. How to add a new dataset / instrument

### Minimum manual path (works today, no refactor)
1. **Observed data.** Produce search results + raw peaks in one of the supported forms:
   - MaxQuant `msms.txt` + **mzML** → run `ambiguity_search/maxquant.py` to emit a
     `PXD-style` matched-pairs parquet (the `figure2` input format), **or**
   - MaxQuant `input_file.txt` (evidence export) + **MGF** files → the `figure3-6` input
     format. Ensure `Raw file` equals the MGF basename and scan numbers align.
2. **Predicted data.** Pick the intensity model + acquisition parameters that match the
   new instrument and pass them explicitly (do not rely on the LUMOS/HCD/CE=28 defaults):
   - Orbitrap HCD → `Prosit_2020_intensity_HCD` / `AlphaPeptDeep_ms2_generic`, CE from data.
   - timsTOF → an AlphaPeptDeep timsTOF model **and** use the observed CCS for a real
     CCS comparison (figure4 currently only compares predicted CCS).
   - CID / EThcD → a matching fragmentation model **and** extend the annotator to the
     correct ion series (c/z for ETD); b/y-only annotation will silently undercount.
   - Adjust `precursor_charges` to the real charges and relax the `Charge == 2` /
     `max_charge=2` filters if the dataset has higher-charge precursors.
   - Set fragment-match `tol_ppm` to the instrument resolution (20 ppm assumes Orbitrap).
3. **Conform schemas** to §6 (column names, `non_switched` flag, `ID` column for
   `metrics_comparison`).
4. **Paths/identity.** Add the new files under `temp_data/` and update the path/dataset
   literals in the relevant notebook cells.

### Scientific caveats when changing instrument
- The intensity model **must** match fragmentation and instrument; a mismatched model
  produces plausible-but-wrong intensities and invalidates the similarity comparison.
- CCS discrimination is only meaningful with an ion-mobility instrument that yields
  **observed** CCS; otherwise figure4's CCS panel measures only prediction divergence.
- RT comparison uses DeepLC trained/calibrated on each dataset's own observed RT, so it
  transfers across instruments as long as observed RT is available; the Koina
  `Deeplc_hela_hf` model in `make_predictions` is HeLa-specific and not used by the
  notebooks (they train DeepLC locally).

## 9. Known bugs, drift, and risks

- **`hyper_score` is broken** (`metrics/metrics.py:133`): `factorial(nb)` is called with
  `nb`/`ny` as **sets**, not counts → `TypeError`; `metrics_comparison` swallows it as
  `NaN`. Fix to `factorial(len(nb)) * factorial(len(ny))`.
- **Koina error handler hangs / crashes** (`intensity_predictions.py:226-231`): on a
  prediction error it `print`s and calls `input()` (blocks in non-interactive runs) and
  never re-raises, so `predictions` is undefined → `NameError`. Also `obtain_predictions_pairs`
  with default `charges=[2]` and >1 peptide produces a length-mismatched `precursor_charges`
  column; the working path only succeeds because `safe_obtain_predictions` passes
  `[2]*N`.
- **`cosine_similarity` returns a distance**, not a similarity (`scipy…cosine` = 1 − cos).
  Higher = more different. `spectral_angle` is the correct similarity form.
- **figure2 seed reconstruction is fragile** (`figure2.ipynb`): it rebuilds the original↔
  switched mapping by replaying the FASTA parse, `np.random.default_rng(42)`, the 20000-row
  `rng.choice`, and `switch_random_il`. But `switch_random_il` uses the **global `random`
  module** (`seq_utils/peptide.py:52`), which figure2 does not seed, so swapped sequences
  need not match the saved CSVs. The mapping is also lossy (`dict(zip(...))` collapses
  collisions), and most of the 20000 rows drop during alignment. **Recommendation: persist
  the original-sequence/`ID` join key in the prediction CSVs so no replay is needed.**
- **Three metric vocabularies have drifted**: current `metrics.py` / `get_metrics.metric_keys`
  (16 metrics ending in `hyper_score`); the figure4/5 notebooks (15-metric list, no
  `hyper_score`); and `notebook1`/`fragment_intensity_il_swap` plus
  `similarity_results_only_exact_match.parquet` (an older richer set: `dot_product`,
  `mara_similarity`, `massbank_score`, `gnps_score`, `stein_scott_score`, …). Unify on one
  canonical module before adding datasets, or comparisons across datasets will use
  different metric sets.
- **Stale note correction:** `CODEBASE_NOTES.md` claims `ruzicka_similarity_1/2`
  normalization is a no-op. In the **current working tree** (`git status` shows
  `metrics/metrics.py` modified) lines 75-76 and 83-84 **do** assign back, so L1 vs L2
  Ruzicka now differ correctly. That note is out of date.
- **Inconsistent length filters** (7–30 vs 6–60) and **inconsistent figure output dirs**
  (`temp_data/` vs `figs/`).
- Notebooks duplicate large blocks (MaxQuant loading + ProForma parsing, MGF parsing,
  theoretical-fragment annotation, observed↔predicted merge, scoring loop, ratio-correct,
  flanking/distance features) across 5+ files. These should be extracted before adding
  datasets multiplies the copies.

## 10. Recommended refactor for multi-dataset support

To make new instruments a config change rather than a code edit:

1. **Dataset config object** (dataclass or one YAML per dataset) carrying: `input_file`,
   `mgf_files` / `mzml_folder`, prediction-cache paths, `instrument`, `fragmentation`,
   `collision_energy`, `charge_filter`, `intensity_model`, `rt_model`, `ccs_model`,
   `tol_ppm`, `column_map`, and the ProForma `modification_dict`.
2. **`io` module**: `load_search_results(path, column_map)` (generalize beyond MaxQuant
   column names) and a single `parse_spectra(...)` supporting both MGF and mzML with a
   configurable scan-id strategy (replace the two duplicate `parse_mgf` variants).
3. **Extend `seq_utils`**: `build_il_pairs(peptides) -> (switched, orig_to_swap)` and a
   single `get_swap_info`/flanking/distance helper; switch `switch_random_il` to an
   injectable seeded RNG.
4. **Extend `make_predictions`**: one configurable predictor taking `model`,
   `instrument_types`, `fragmentation_types`, `collision_energies`, `charges`; remove the
   inline `safe_obtain_predictions` redefinitions in the notebooks; add the local-DeepLC RT
   path the notebooks actually use; fix the error handler.
5. **`annotation` module**: `annotate_spectrum(mz, intensity, peptidoform, charge,
   tol_ppm, ion_types)` supporting charges > 2 and ion series beyond b/y.
6. **`scoring` module**: one `metric_keys`, `score_pairs(...)`, and `ratio_correct(...)`
   with a single sign/orientation convention.

After this, each notebook reduces to: load config → call pipeline functions → plot, and a
new instrument is a new config plus (if needed) a new model/ion-series entry.

## 11. Common commands

```shell
# environment
conda env create -f environment.yml
conda activate lorentz_center_spectrum_similarity

# find I/L sibling peptides in a FASTA
python find_siblings/digest_find_siblings.py <fasta> <out>

# build observed I/L pairs from MaxQuant + mzML
python -m ambiguity_search.maxquant <out.parquet> <mzml_folder> <maxquant_folder>...

# figures: run the notebooks in notebooks/ (figure1 = predicted; figure2 = observed vs
# predicted; figure3-6 = discriminative power on the van Puyvelde / Orbitrap dataset)
```

## 12. Multi-instrument extension (added)

A first, additive extension to multiple instruments is implemented in the predicted
branch (no observed raw data or bulk download required). It does not modify the existing
notebooks or modules, so the original pipeline is unaffected.

**New files**
- `make_predictions/instruments.py` — `InstrumentConfig` dataclass and `INSTRUMENT_PRESETS`
  registry. This is where the previously hardcoded `LUMOS`/`HCD`/`CE`/charge values now
  live. Adding an instrument is one entry here. `predict_intensities(peptides, config)`
  is a corrected Koina wrapper (sends only the columns each model needs; never calls
  `input()`; always re-raises on failure).
- `make_predictions/il_pipeline.py` — reproducible pipeline: `load_peptides` (sorts before
  sampling), `build_il_pairs` (seeds the global RNG used by `switch_random_il` and carries
  an explicit `ID` so no FASTA replay is needed), `similarity_for_instrument`,
  `run_comparison`, `summarize`.
- `scripts/run_multi_instrument_analysis.py` — runnable end-to-end CLI: digest, pair,
  predict per instrument, score, write parquet + summary + overlay figure.
- `scripts/download_dataset.py` — `ppx`-based downloader for an observed dataset from a
  ProteomeXchange/PRIDE accession (for extending the observed branch later).
- `notebooks/figure7_multi_instrument.ipynb` — interactive version of the analysis.

**Run it** (requires the conda env active and network access to Koina):
```shell
conda activate lorentz_center_spectrum_similarity
python scripts/run_multi_instrument_analysis.py --list           # show instrument presets
python scripts/run_multi_instrument_analysis.py --n 2000 \
    --instruments unispec_lumos,astral_dda,alphapept_timstof
# outputs: temp_data/multi_instrument_similarity.parquet, temp_data/figure7_multi_instrument.png
```

**Add a new instrument**: add an `InstrumentConfig` to `INSTRUMENT_PRESETS` with the Koina
`intensity_model` and the `instrument_type`/`fragmentation_type`/`collision_energy` that
model accepts (check the model card at https://koina.wilhelmlab.org), then pass its key
via `--instruments`.

**Extend the observed branch to a new instrument** (heavier): `scripts/download_dataset.py
<accession> --pattern .mgf --pattern msms.txt` to fetch processed peak lists + search
results, conform them to the schemas in §6, then reuse the annotation/scoring path from
`figure3`-`figure6`. Note the feasibility caveats in §8 (raw->search tooling is not bundled).

**Status note:** the modules, script, and notebook are written and syntax-checked, and the
deterministic pair-building was validated offline. The live Koina run and any dataset
download were not executed in the authoring session because command execution
(network / conda env) was blocked by the session permission mode.

## 13. Benchmark datasets for new instruments (Astral, timsTOF)

To extend the observed-data figures to the Orbitrap Astral and Bruker timsTOF, use the
DDA arms of the multi-species LFQ benchmarks. They share the human/yeast/E. coli design
and the Condition A/B replicate structure of the existing Orbitrap data.

| Instrument | Acquisition | Accession | Public on PXD | Fragmentation | Use |
|---|---|---|---|---|---|
| Orbitrap QE HF-X (current) | DDA | PXD028735 (Van Puyvelde "Gen Alpha") | yes | HCD | existing `LFQ_Orbitrap_DDA_*` |
| Bruker timsTOF Pro | DDA-PASEF + ion mobility | PXD028735 | yes | CID/PASEF | DDA figures + observed CCS |
| timsTOF SCP | diaPASEF (DIA) | PXD062685 (ProteoBench) | not yet | CID | DIA only |
| Orbitrap Astral | DDA (and DIA 2 Th) | PXD070049 / PXD071205 ("Gen Beta") | yes | HCD (DDA NCE 30) | DDA figures |

Identification results: pull ProteoBench curated outputs from the per-module results
repositories under https://github.com/Proteobench (`Results_quant_ion_DDA_Astral`,
`Results_quant_ion_DIA_diaPASEF`, `Results_quant_ion_DDA`). The figure3-6
spectrum-annotation path needs per-PSM scan info (`Raw file`, `MS/MS scan number`),
present in MaxQuant evidence/msms and FragPipe `psm.tsv` but NOT in DIA-NN precursor
reports; conform whichever table you use to the §6 schemas.

Which instrument feeds which figure:
- Predicted-only (figure1, figure6, figure7) and RT: both instruments now, no download.
  Astral uses the `astral_dda` preset (Prosit HCD, NCE 30); timsTOF uses
  `alphapept_timstof` (AlphaPeptDeep, instrument TIMSTOF).
- Observed-spectrum discriminative power (figure2-5): Astral DDA and timsTOF Pro
  DDA-PASEF (b/y under HCD/CID). The DIA arms (Astral 2 Th, diaPASEF) do not fit the
  per-spectrum b/y figures without pseudo-spectrum extraction.
- Observed CCS: timsTOF only; it converts the predicted-only CCS panel into a real
  observed-vs-predicted comparison.

Ingestion: `scripts/download_benchmark_data.py <astral_dda|timstof_dda>` downloads the
matching DDA files (ppx), optionally converts vendor raw to MGF (ThermoRawFileParser for
`.raw`; msconvert/tdf2mzml for timsTOF `.d`, preserving ion mobility for CCS), and can
clone the ProteoBench identification repo. Converters are external and not bundled; the
script detects them and prints exact commands if missing. Run with `--list` first.

## 14. New-instrument figures: environments, scripts, and results

Two conda environments are used; do not mix them:
- `base` (anaconda3): Koina prediction work. `koinapy` and `ppx` are installed here.
  Installing koinapy upgraded protobuf, which breaks DeepLC/TensorFlow in base.
- `py312_deeplc`: DeepLC retention-time work (DeepLC 3.1.8, psm_utils; matplotlib and
  pyarrow added). Run the RT script with
  `/home/robbin/anaconda3/envs/py312_deeplc/bin/python`.

New scripts (parameterized, additive):
- `scripts/run_multi_instrument_analysis.py` (base): predicted multi-instrument similarity (figure 7).
- `scripts/observed_discriminative_power.py` (base): observed fragment-intensity ratio-correct (figure 3 port).
- `scripts/timstof_mobility_il.py` (base): observed ion-mobility I/L separation (model-free).
- `scripts/observed_rt_discriminative.py` (py312_deeplc): observed RT ratio-correct (DeepLC).
- `scripts/download_benchmark_data.py` (base): ppx download plus ThermoRawFileParser / tdf2mzml conversion.

Data obtained: Astral DDA raw (PXD070049, Condition A/B REP1) converted to MGF via
ThermoRawFileParser in `temp_data/astral_mgf/`; ProteoBench MaxQuant Astral identifications
in `temp_data/pb_astral_mq/`; ProteoBench diaPASEF DIA-NN report (timsTOF) in
`temp_data/pb_diann_inspect/`. ThermoRawFileParser v2 (Linux) is in `temp_data/trfp/`.

Results so far (multi-modal I/L discrimination, fraction correct = observed closer to the
original than to the I/L swap):

| modality | Orbitrap | Astral | timsTOF |
|---|---|---|---|
| fragment intensity (best metric) | original analysis | ~0.61 (n=19,180) | DIA export lacks fragments |
| retention time (DeepLC, n=4000) | 0.603 [0.587-0.618] | 0.576 [0.560-0.591] | 0.607 [0.592-0.622] |
| ion mobility (observed) | n/a (no IM) | n/a (no IM) | no separation (siblings = noise, p=0.22) |

Interpretation: fragment intensity and retention time each carry a weak but statistically
significant I/L signal (about 0.58-0.61, 95% CIs exclude 0.5); observed ion mobility does
not resolve I/L at this resolution. Figures: `temp_data/figure7_multi_instrument.png`,
`figure_astral_ratio_correct.png`, `figure_rt_{orbitrap,astral,timstof}.png`,
`figure_rt_combined.png`, `figure_timstof_mobility_il.png`.

timsTOF fragment-intensity figure status (blocked in this environment): the timsTOF Pro
DDA-PASEF raw `.d` files were downloaded (Van Puyvelde PXD028735,
`LFQ_timsTOFPro_PASEF_Condition_A/B_Sample_Alpha_01`, in `temp_data/genalpha_PXD028735/`
and unzipped in `temp_data/timstof_d/`). Two reading routes were tried and both failed:
(1) `tdf2mzml` (Bruker proprietary SDK) aborts with `std::bad_cast` plus
`boost::interprocess::intermodule_singleton initialization failed`, a low-level SDK/sandbox
incompatibility not resolved by locale (C.UTF-8 fixes the bad_cast but the boost error
remains), /dev/shm, or temp-dir fixes. (2) Sage reads the `.d` natively via timsrust and
runs, but yields no confident identifications (median 4 matched fragments per spectrum, no
b/y ladders, all spectrum_q > 0.12), consistent with uncalibrated m/z from the open-source
reader for this 2021 data. To complete this figure, convert the `.d` to mzML on a system
where ProteoWizard msconvert or the Bruker SDK works (a short step), then run Sage on the
mzML and `scripts/timstof_fragment_il.py`, or run `scripts/observed_discriminative_power.py`
on the mzML with the resulting identifications. The Sage config, HYE FASTA, and analysis
scripts are all in place. timsTOF still contributes three modalities (predicted similarity,
observed RT, observed ion mobility).

