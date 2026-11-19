# Little Functions — public article files

This directory contains the public supporting files for [the Little Functions article](https://clausongeomet.com/writing/posts/little-functions/). It is a deliberately limited subset of the original **private** review package.

The five-mineral synthetic feed and response are in [copper_v2.csv](materials/examples/copper_v2/copper_v2.csv); values are positive closed fractions, and recovery is percent of *contained copper*. The three main article figures have summary input CSVs in `assets/`. The older 240-row illustrative study is separate.

To explore R's `easyCODA::STEPR`, install easyCODA and run the article's R example. The log-ratio random-forest selector is in `materials/src/R/lr_ml.R` (requires `ranger`) and `materials/src/python/lr_ml.py` (requires `numpy`, `pandas`, `scikit-learn`). The restricted `materials/src/python/stepr.py` is **not** a complete parity-correct Python replacement for R STEPR. Five of eight historical parity tests passed; three strict third-pair selections did not. See the [scope note](materials/docs/sc023/HISTORICAL_FAILURE_SCOPE.md).

**Reproducibility boundary:** The 14 original saved models, admitted runtime environment, and complete frozen scientific claim register remain private. The linked archived scripts `examples/read_saved.R`, `examples/read_saved.py` and `run-examples.sh` require those materials for exact replay; they are provided for inspection, not as a standalone fresh reproduction. Article values are from previously accepted synthetic evidence, not new modelling or a plant study.

**Licences:** Website article prose is CC BY 4.0. The supplied R/Python source and the Python adaptation of upstream STEPR are **GPL-3.0**, with [the licence](materials/LICENSE-GPL-3.0.txt) included and upstream method attribution retained. These licences must not be conflated. The project has not made its entire private code or model archive public.

The article was editorially approved on 2026-11-19. LF-02 remains blocked; LF-04 remains design-only and unexecuted.
