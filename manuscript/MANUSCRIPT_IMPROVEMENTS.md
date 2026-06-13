# Manuscript improvement analysis: I/L sibling spectral similarity

A deep assessment of how to strengthen this work for publication, organized from the
most fundamental issues down to specific analyses and engineering.

Revision note (2026-06-12): revised after author feedback and a completed literature
search. A substantial part of the original critique was anchored to
`manuscript_draft.md`, which the authors have confirmed is outdated and superseded;
those manuscript-dependent points have been withdrawn and are flagged where they
occurred (sections 6, 8, and 9). The ground-truth analysis in section 2 was corrected:
the FASTA database, not any search-engine assignment, defines the ground truth. The
literature search is complete and its verified results are in section 10.

---

## 1. Three evidence levels, used as a deliberate narrative

The study spans three questions, and the planned narrative addresses them in order:
first show that predictions differ for I vs L, then show the effect in experimental
data alone, then show recovery (ranking the correct sibling). This staged structure is
sound. The only requirement is to keep the claim at each level scoped to what that
level supports:

- Level C, prediction only: predicted(original) vs predicted(swap). This shows the
  models assign different intensities to I vs L. Scope the claim to "predictions
  differ", not "siblings are distinguishable in reality".
- Level A, model-free: observed spectrum vs observed spectrum, for pairs where both
  siblings are independently observed. This is the model-independent evidence about the
  physics.
- Level B, recovery: observed spectrum vs predicted(original) vs predicted(swap),
  scored against the FASTA ground truth (section 2). This is not circular, because the
  ground truth is genomic, not model-derived. Its performance is jointly determined by
  the real physical difference and the model's fidelity in predicting it.

The gap between Level A and Level B is a useful quantity: it is the extraction
efficiency of the predictors for I/L, that is, how much of the real, observable
difference the models actually capture. It is worth reporting explicitly.

---

## 2. Ground truth is the FASTA database, not a search-engine assignment

The ground truth for I vs L is set by the sequence database, which is built from
genomic and other non-MS evidence and fixes each residue as either I or L. The
mass-spectrometry assignment does not define the truth; it inherits it from the
database. There is therefore no arbitrary coin-flip in the label and no circularity in
the recovery experiment. (An earlier version of this section argued the opposite; that
argument was wrong and has been withdrawn.)

The research question is then clean and well posed: if the database constraint is
removed, so that I and L are treated as a priori equally likely, can the spectral, RT,
and CCS evidence recover the genomic residue? ProteomeTools (Zolg et al., 2017,
DOI 10.1038/nmeth.4153) remains useful as an optional strengthening, because synthetic
peptides have certain identity and no proteoform ambiguity, but it is not required to
fix a ground-truth problem.

One residual subtlety supports the prevalence and use-case framing rather than
undercutting it: the database fixes each protein's residues, but when the I/L sibling
is itself a tryptic peptide of a different database protein, which protein a given
spectrum originated from is genuinely ambiguous at the proteoform level. Those cases,
where two database proteins form an I/L sibling pair, are the practically hard ones and
are the right target to quantify in the prevalence analysis.

---

## 3. The scientific crux (already covered in the manuscript introduction)

The authors confirm this is established in the introduction. Recorded here for
completeness: for I and L the residue mass is identical, so all backbone b and y ions
and the immonium ion at m/z 86 share identical m/z regardless of isomer. Under HCD with
b/y annotation the siblings differ only in intensities, never in peak positions, which
makes HCD I/L discrimination a pure intensity-pattern inference problem. The physical
origin is beta-branching of isoleucine modulating neighboring backbone cleavage
(mobile-proton and steric effects; the mechanistic citations such as Wysocki and Paizs
and Suhai still need verification, see section 10).

ETD, EThcD, and UVPD generate side-chain fragments (w and d ions for the radical
methods) whose masses differ between I and L, giving deterministic, mass-level
discrimination. This is the established reliable route (section 10). Note one
categorization fix: Wu et al. 2023 is intensity-based statistical prior art, not a
deterministic w-ion method, and should be cited accordingly (an earlier draft of this
document grouped it with the w-ion papers, which was incorrect).

---

## 4. Recommended paper framing (choose the strongest story)

Three coherent stories are possible.

- Framing A, predictor benchmark (methods, defensible, modest data needs): "Do
  fragment-intensity, RT, and ion-mobility predictors encode a real isoleucine/leucine
  difference, and how faithfully?" Validate on ProteomeTools synthetic ground truth.
  Quantify per-model I/L sensitivity and correctness.
- Framing B, information budget (comprehensive, ambitious): decompose I/L
  discriminability across evidence levels (A/B/C in section 1), activation methods
  (HCD vs EThcD w-ions), and modalities (MS2 intensity, RT, CCS). Needs EThcD and a
  timsTOF/ion-mobility dataset.
- Framing C, practical tool (applications): proteome-wide prevalence of confusable
  I/L siblings plus a calibrated scoring tool that flags or resolves ambiguous
  assignments.

Recommended synthesis: motivate with prevalence (C), make the predictor and
similarity-metric comparison the core, present multimodal fusion as the practical
payoff, and scope to HCD while positioning against the deterministic w-ion route
(B-lite, ideally with a small EThcD contrast).

---

## 5. Concrete analyses and figures to add or change

- Prevalence figure (motivation): using the proteome-wide digestion in
  `find_siblings/`, report how many identified-or-identifiable tryptic peptides have a
  same-mass I/L sibling that is also tryptic and in-database, per organism. Be precise
  that database co-existence is an upper bound on practical ambiguity; refine toward
  "confusable given typical detectability". This is also where the proteoform-ambiguity
  case from section 2 is quantified.
- Predictor I/L-sensitivity figure: the distribution of predicted similarity between
  a sequence and its swap, per model, and whether the predicted difference is
  calibrated against the observed difference.
- Decomposition figure: Level A (physics) vs Level B (achieved) discriminability;
  the gap is the model extraction efficiency.
- Mechanistic figure: where on the peptide and with which neighboring residues the I/L
  intensity fingerprint appears, interpreted through beta-branching. This is the most
  original science.
- Deterministic contrast (if EThcD data is reachable): w-ion-based discrimination as
  the upper bound, to contextualize the HCD numbers.

---

## 6. Statistics (optional, to back the "relatively good separation" claim)

The earlier version of this section criticized the draft Methods for promising
bootstrap, confidence intervals, and significance tests that were not performed. That
critique was anchored to the outdated `manuscript_draft.md` and is withdrawn.

The only forward-looking point, and it is optional, supports the positive claim in
section 11. To state that the separation is "relatively good", give the adjective a
number: report the recovery rate with a confidence interval and a test against the
chance rate. Each peptide pair is a Bernoulli trial, so a binomial proportion interval
(Wilson or Clopper-Pearson) and an exact binomial test against chance are the natural,
lightweight choices. Nothing heavier is required for the claim to hold.

---

## 7. Metric breadth (suggestion withdrawn)

The authors prefer a comprehensive comparison, and a broad survey of similarity metrics
is a legitimate contribution. The earlier suggestion to reduce to a single principled
discriminant is withdrawn.

One narrow caveat remains, and only if it applies: if the manuscript reports the
accuracy of a single selected best metric or best metric combination (for example a
best-performing trio), estimate that headline number on data not used for the
selection, otherwise it will be optimistic. The broad descriptive comparison across all
metrics does not have this problem.

---

## 8. CCS and multimodal combination (mostly withdrawn)

The earlier version critiqued specific abstract claims (a combined accuracy above 75
percent and an ion-mobility "validation"). Those came from the outdated draft and are
withdrawn.

One manuscript-independent fact is worth keeping, purely as information: in the current
code and data, CCS is predicted only, on Orbitrap data with no observed ion mobility,
so any CCS comparison reflects prediction divergence rather than an observed
measurement. This matters only if and when a timsTOF or other ion-mobility dataset is
added, at which point CCS could become a genuine observed signal.

---

## 9. Engineering credibility (a methods paper will be checked against its repo)

The earlier "abstract numbers do not match the code" item was based on the outdated
draft and is withdrawn. The remaining items are facts in the current repository and
stand:

- Metric definition bugs: `cosine_similarity` returns a distance (1 minus cosine), and
  `spectral_angle` then takes arccos of that distance rather than of the cosine, so
  both are mis-oriented. `hyper_score` calls factorial on a set and silently fails.
  `fit` masks the denominator but not the numerator. Fix or remove each, add unit tests
  (the existing `test_metrics.py` does not catch these), and regenerate figures.
- Metric-vocabulary drift across the module, the figure notebooks, and an orphaned
  notebook should be unified into one canonical scoring module.
- Reproducibility: figure2 rebuilds the original-to-switched mapping by replaying a
  seed, but `switch_random_il` uses the unseeded global random module, so the
  regeneration is not deterministic and the join can silently mispair. Persist the
  original-sequence join key instead, and seed the swap RNG. (A reproducible pipeline
  with a fixed join key is implemented in `make_predictions/il_pipeline.py`.)
- The swap function excludes position 0 from random selection, an undocumented bias in
  which I/L position is tested. Document and justify or remove it, and report the
  position distribution tested.
- FAIR deposit: provide a one-command reproducible pipeline, deposit data and predicted
  spectra (mzSpecLib where appropriate), and report exact model versions and Koina
  settings.

---

## 10. Literature positioning and novelty (search completed 2026-06-12)

The literature search is complete (source: PubMed plus web). The verified results below
supersede the earlier "to verify" list.

NOVELTY VERDICT: the core idea is not new, and a "first systematic evaluation" framing
is not defensible as written. It must be reframed.

1. pDeep already demonstrated the central premise. Zhou et al. (pDeep, Anal Chem 2017,
   DOI 10.1021/acs.analchem.7b02566) state explicitly that pDeep "showed the potential
   to distinguish extremely similar peptides (peptides that contain isobaric amino
   acids, for example, GG = N, AG = Q, or even I = L)". That is the manuscript's exact
   idea in the predicted domain. Cite pDeep and state precisely what is added beyond it.
2. Wu et al. (Anal Chem 2023, DOI 10.1021/acs.analchem.3c00495; already in main.bib)
   provide a statistical framework, with a calculated probability, for distinguishing
   spectra that differ only in peak intensity, validated on Leu/Ile across CID, HCD,
   ETD, and RDD. This is direct prior art for the intensity-based discrimination and the
   statistical angle.

Defensible novelty: a systematic, large-scale comparison of many spectral-similarity
metrics for I/L, combined with orthogonal RT (and, with appropriate data, CCS)
predictors, anchored to observed data, with calibrated practical performance. Reframe
the contribution around "systematic, multi-metric, multi-modal, observed-data
quantification", and note that the precedents above are predicted-domain or framework
results whereas the systematic observed-data demonstration is the new part.

Method-database note: PubMed indexes almost none of this field. The PubMed I/L queries
returned only amino-acid metabolomics noise; every relevant method paper came from
analytical-chemistry, CS, thesis, or preprint venues. Any further search must prioritize
ACS journals (Anal Chem, JASMS, JPR), Google Scholar, bioRxiv, ChemRxiv, and arXiv.

Still-open verification items (not confirmed by this search; do not treat as settled):
a specific Brodbelt UVPD I/L reference; the mobile-proton mechanism citations (Wysocki;
Paizs and Suhai; Tabb); and the quantitative RT and CCS resolvability of I/L
(SSRCalc/Krokhin and ion-mobility evidence).

Verified must-cite references (attribution: PubMed; DOIs below):

| Role | Citation | DOI |
|---|---|---|
| Direct precedent (predicted intensities distinguish I/L) | Zhou et al., pDeep, Anal Chem 89(23):12690, 2017 | 10.1021/acs.analchem.7b02566 |
| Predictor lineage (modified peptides, ProteomeTools-trained) | Zeng et al., pDeep2, Anal Chem 91(15):9724, 2019 | 10.1021/acs.analchem.9b01262 |
| Prior art, statistical intensity-based isomer (incl. Leu/Ile) discrimination | Wu et al., Anal Chem 95(17):6996, 2023 | 10.1021/acs.analchem.3c00495 |
| Deterministic route (ETD/HCD MS3 w-ions, Orbitrap Fusion, Leu/Ile) | Lebedev et al., Anal Chem 86(14):7017, 2014 | 10.1021/ac501200h |
| Deterministic route (w-ion I/L in antibody de novo) | Schulte and Snijder, J Proteome Res 23(8):3552, 2024 | 10.1021/acs.jproteome.4c00188 |
| Synthetic ground truth | Zolg et al., ProteomeTools, Nat Methods 14(3):259, 2017 | 10.1038/nmeth.4153 |
| Intensity predictor | Gessulat et al., Prosit, Nat Methods, 2019 | 10.1038/s41592-019-0426-7 |
| Intensity predictor | Lapin et al., UniSpec, Anal Chem, 2024 | 10.1021/acs.analchem.3c02321 |
| Property predictor framework | Zeng et al., AlphaPeptDeep, Nat Commun, 2022 | 10.1038/s41467-022-34904-3 |
| RT predictor | Bouwmeester et al., DeepLC, Nat Methods, 2021 | 10.1038/s41592-021-01301-5 |
| Ion mobility / CCS predictor | Devreese et al., IM2Deep, Anal Chem, 2025 | 10.1021/acs.analchem.5c01142 |
| Fragment intensity predictor (group tool) | Declercq et al., MS2PIP web server, NAR, 2023 | 10.1093/nar/gkad335 |

Additional competing or adjacent work (verify; not all peer-reviewed):
- Shen, B. "Discrimination of Leucine and Isoleucine in De Novo Peptide Sequencing
  Using Deep Neural Networks", MSc thesis, University of Western Ontario (deep neural
  networks on raw spectra for L/I; references EThcD). A conceptual competitor in the
  deep-learning-for-I/L space, de novo focused and not peer-reviewed.
- ILDIFF, "Algorithm for Isoleucine and Leucine Discrimination ... Based on
  Electron-Activated Dissociation", Proc. ICBBS 2025 (EAD w-ions on a ZenoTOF).
- "Accurate discrimination of leucine and isoleucine ... combining continuous digestion
  with multiple MS3 spectra integration", Talanta, 2022.
- These confirm the deterministic w-ion / EAD route is the established, reliable method;
  position the HCD intensity-inference approach explicitly against it.

These entries have been added to `manuscript/main.bib`.

---

## 11. Scope and claims

The positive framing is warranted. Against the common practice of treating I and L as
indistinguishable by HCD, demonstrating a clear and relatively strong separation is a
noteworthy result. The earlier "embrace the negative result" recommendation is
withdrawn.

Two supports and one condition:
- Precedent: pDeep (Zhou et al., 2017, DOI 10.1021/acs.analchem.7b02566) already showed
  predicted spectra can distinguish I=L, and Wu et al. (2023,
  DOI 10.1021/acs.analchem.3c00495) showed intensity-only statistical discrimination
  for Leu/Ile. The positive claim therefore has precedent and is not a lone outlier.
- The new part: those precedents are predicted-domain or framework results, whereas the
  practitioner belief concerns observed HCD in routine workflows. The systematic
  observed-data demonstration is the genuinely new contribution, so frame it as
  extending pDeep and Wu to observed data at scale.
- The condition: quantify "relatively good" with the above-chance number and interval
  from section 6, so the adjective is backed rather than asserted.

On the title: scoping it to the HCD intensity setting (unless an EThcD/UVPD or
multi-instrument arm is added) keeps the claim aligned with the evidence while still
making the positive point.

---

## 12. Prioritized roadmap

Durable, manuscript-independent items:
1. Engineering (section 9): fix the metric bugs, unify the metric vocabulary, rebuild
   figures from one seeded pipeline, and keep numbers traceable to committed code.
2. Quantify the separation (section 6, optional but recommended): recovery rate with a
   confidence interval and an above-chance test, to back the section 11 claim.
3. Literature and novelty (section 10): cite pDeep and Wu et al. as precedent, position
   against the deterministic w-ion route, and frame the contribution as the systematic,
   observed-data, multi-metric, multi-modal quantification.

Optional strengthening:
4. ProteomeTools as a certain-identity check, and the Level A vs Level B extraction-
   efficiency decomposition (sections 1, 2).
5. Mechanistic interpretation of the positional and flanking fingerprint via
   beta-branching (section 5), and, if data allow, an EThcD/w-ion deterministic contrast
   and a timsTOF observed-CCS arm.

The core idea is publishable and the positive result is defensible. The remaining work
is to make the quantitative claims traceable and backed, and to position the
contribution precisely against pDeep, Wu et al., and the deterministic w-ion route.
