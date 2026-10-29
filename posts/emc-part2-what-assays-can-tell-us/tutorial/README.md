# EMC Part 2: reproducible teaching tutorial

The article's teaching example follows the original research division of labour:
`simCoDaR` in R for hypothetical mineral compositions and theoretical assays;
Python for constrained fitting, distributional scores and mineral-order search.
This 60-specimen demonstration is **not** the six-scenario research study.

## 1. Install dependencies

R >= 4.0:

```r
install.packages("remotes")
remotes::install_github(
  "ClausonGeomet/simCoDaR",
  ref = "d9a0b27676a0cb3a4acb300a3f6b2ebc845d205d",
  upgrade = "never"
)
```

Python >= 3.10:

```sh
python -m pip install -r requirements.txt
```

## 2. Run the example

Run commands from the directory containing the downloaded tutorial scripts:

```sh
Rscript generate_data.R .
python explore_emc.py .
```

This creates `emc_data.csv` (60 synthetic specimens) and `emc_C.csv`
(ideal mineral chemistry), and prints forward-moment checks, mineral
estimates, allocation-order contrasts, score components and the chemical
non-identifiability witness.

The Python script imports the unmodified frozen research source files:
`baselines.py` and `emc_prototypes.py`, included in this download folder. These
are archived as they were in
[ClausonGeomet/writing at research commit d2f69fe](https://github.com/ClausonGeomet/writing/tree/d2f69fe3b4375a4a3dd03ea1ac2cb40149ee8b79/research).

The main research study used six synthetic scenarios × three replicate cohorts
and a fixed 30/10/20 train/validation/test design in each cohort. This tutorial
uses only **one newly generated 60-row cohort**, so its printed errors need
not match the published aggregate results. The exact research provenance,
source settings, output manifests and historical corrections are preserved in
the linked [Part 2A](https://github.com/ClausonGeomet/writing/tree/d2f69fe3b4375a4a3dd03ea1ac2cb40149ee8b79/research/emc-part2a)
and [Part 2B](https://github.com/ClausonGeomet/writing/tree/d2f69fe3b4375a4a3dd03ea1ac2cb40149ee8b79/research/emc-part2b)
branches.

## Limits of what is reproduced

- `simCoDaR` uses hypothetical compositional targets and ideal mineral formulae;
  the XRF/XRD-like outputs are not independent measurements.
- The computed population moments use generated mineral truth only to verify
  the forward statistical identity; a field EMC model cannot assume that truth.
- The Gaussian and finite experimental designs in the study used separate
  assumptions and cohorts, not reproduced by this compact example.
- The fixed four-phase GA searches **orderings**, not geological assemblage
  probabilities. All reported scores are exploratory chemical scores.
- Results can vary slightly with R, NumPy or SciPy versions. For the exact
  historical environment and data, consult the linked freeze manifests.
