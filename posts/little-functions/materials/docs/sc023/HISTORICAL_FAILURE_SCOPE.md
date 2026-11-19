# Python STEPR: known limits

`stepr_py` is a limited Gaussian method-1 translation of the stepwise log-ratio selection from R `easyCODA::STEPR`. It is **not** a full port.

Eight strict historical R/Python parity checks used a separate 240-row synthetic fixture. Five passed and three failed. All three failures involve the selected **third ratio**, although the resulting OLS predictions differed by no more than about 1.07 × 10^-13. Near-identical predictions do not make the selected ordered pairs identical.

In the main 1,024-row synthetic copper study, Python STEPR comparator models **reuse R-selected pairs**. They do not provide evidence of independently matching Python STEPR selection on those rows. The R/Python `select_lr_ml` forest selector is a *different method*.

Other historical full-suite issues remain unresolved and are outside this bounded blog result. No causal plant response, general algorithm ranking or full software acceptance is claimed.
