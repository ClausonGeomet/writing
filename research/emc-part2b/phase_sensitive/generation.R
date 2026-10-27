#!/usr/bin/env Rscript
# Sole stochastic path. Invoked only after the Python execution gate.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("usage: generation.R CONFIG OUTPUT ISOLATED_LIBRARY")
config <- jsonlite::fromJSON(args[[1L]], simplifyVector = FALSE)
out <- args[[2L]]
lib <- normalizePath(args[[3L]], mustWork = TRUE)
.libPaths(c(lib, .libPaths()))
if (normalizePath(find.package("simCoDaR"), mustWork = TRUE) !=
    normalizePath(file.path(lib, "simCoDaR"), mustWork = TRUE)) stop("wrong simCoDaR installation")
library(simCoDaR, lib.loc = lib)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
options(digits = 17)
phases <- c("Quartz", "Anatase", "Rutile")
C <- rbind(c(1, 0, 0), c(0, 1, 1))
dimnames(C) <- list(c("SiO2", "TiO2"), phases)
lib_local <- mineral_assay_library()
lib_local$Anatase <- list(formula = "TiO2", atoms = c(Ti = 1, O = 2))
lib_local$Rutile <- list(formula = "TiO2", atoms = c(Ti = 1, O = 2))
eye <- diag(3)
colnames(eye) <- phases
pure <- minerals_to_assay(eye, library = lib_local,
                          basis = "fraction", iron_reporting = "total_Fe2O3")
if (max(abs(t(pure$oxides_reported[, c("SiO2", "TiO2"), drop = FALSE]) / 100 - C)) > 1e-12)
  stop("pure-mineral chemistry mismatch")
Rb <- matrix(c(1e-4, 3e-5, 3e-5, 1e-4), 2)
Rq <- diag(4e-4, 3)
Corr <- matrix(c(1,.3,0,0, .3,1,.5,0, 0,.5,1,0, 0,0,0,1), 4, byrow = TRUE)
R4 <- diag(c(.01,.01,.01,.02)) %*% Corr %*% diag(c(.01,.01,.01,.02))
noise <- function(n, cov, seed) {
  set.seed(seed)
  matrix(rnorm(n * nrow(cov)), nrow = n) %*% chol(cov)
}
truth <- function(center, n, seed) {
  endpoints <- matrix(as.numeric(unlist(center)), nrow = 1,
                      dimnames = list("domain1", phases))
  spec <- new_assay_scenario(endpoints, weights = 1, logratio_sd = .35,
                             structural_parts = NULL, composition_bounds = NULL)
  d <- simulate_mineral_assay_scenario(spec, n = n, seed = seed,
          library = lib_local, iron_reporting = "total_Fe2O3", max_draw_attempts = 1000L)
  raw <- d$truth$mineral_fraction[, phases, drop = FALSE]
  if (!is.data.frame(raw) || !all(vapply(raw, is.numeric, logical(1))))
    stop("truth phase columns must be numeric data frame")
  x <- as.matrix(raw)
  storage.mode(x) <- "double"
  if (any(!is.finite(x)) || max(abs(rowSums(x) - 1)) > 1e-12 || any(x < 0))
    stop("invalid generated truth")
  x
}
root <- file.path(out, "data")
dir.create(root, recursive = TRUE, showWarnings = FALSE)
kinds <- c("calibration", "direct", "pool", "evaluation")
for (kind in kinds) for (suffix in c("truth", "observed")) {
  f <- file.path(root, paste0(kind, "_", suffix, ".csv"))
  if (file.exists(f)) stop("refusing existing data")
}
append <- function(kind, suffix, meta, values) {
  f <- file.path(root, paste0(kind, "_", suffix, ".csv"))
  tab <- data.frame(meta, values, check.names = FALSE)
  write.table(tab, f, sep = ",", row.names = FALSE, col.names = !file.exists(f),
              append = file.exists(f), quote = TRUE, na = "")
}
master <- 731900000L
for (cell in config$cells) {
  s <- as.integer(cell$cell_id)
  for (r in seq_len(as.integer(config$replicates))) {
    base <- master + 100000L*s + 100L*r
    cal <- truth(cell$calibration_center, 20L, base+1L)
    direct <- truth(config$dgp$target_center, 20L, base+2L)
    pool <- truth(config$dgp$target_center, 10L, base+3L)
    eval <- truth(config$dgp$target_center, 30L, base+4L)
    obs <- list(
      calibration = cal + noise(20L, Rq, base+5L),
      direct = direct %*% t(rbind(C,C[2,],c(0,1,0))) + noise(20L,R4,base+6L),
      pool = pool %*% t(C) + noise(10L,Rb,base+7L),
      evaluation = eval %*% t(C) + noise(30L,Rb,base+8L))
    rep_bulk <- pool[1:5,,drop=FALSE] %*% t(C) + noise(5L,Rb,base+9L)
    qxrd <- pool[1:5,,drop=FALSE] + noise(5L,Rq,base+10L)
    colnames(obs$calibration) <- c("QXRD_Q","QXRD_A","QXRD_R")
    colnames(obs$direct) <- c("SiO2","TiO2","redundant_Ti","H_A")
    colnames(obs$pool) <- c("SiO2","TiO2")
    colnames(obs$evaluation) <- c("SiO2","TiO2")
    extra <- matrix(NA_real_, 10, 5,
                    dimnames=list(NULL,c("rep_SiO2","rep_TiO2","QXRD_Q","QXRD_A","QXRD_R")))
    extra[1:5,] <- cbind(rep_bulk, qxrd)
    obs$pool <- cbind(obs$pool, extra)
    for (kind in kinds) {
      x <- get(if (kind == "evaluation") "eval" else kind)
      ids <- sprintf("ps1_s%d_r%d_%s_%d", s, r,
                     switch(kind, calibration="cal",direct="direct",pool="pool",evaluation="eval"),
                     seq_len(nrow(x)))
      meta <- data.frame(cell_id=s,cohort_id=r,sample_id=ids)
      colnames(x) <- c("Q","A","R")
      append(kind,"truth",meta,x)
      append(kind,"observed",meta,obs[[kind]])
    }
  }
}
cat(sprintf("generation_complete cohorts=%d calls=%d rng=%s\n",
            length(config$cells)*as.integer(config$replicates),
            4L*length(config$cells)*as.integer(config$replicates),
            paste(RNGkind(),collapse="/")))
