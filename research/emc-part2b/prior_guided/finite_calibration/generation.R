#!/usr/bin/env Rscript
# Sole declared synthetic draw path. The package is loaded from the exact pinned export.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("usage: generation.R CONFIG OUTPUT ISOLATED_R_LIBRARY")
config <- jsonlite::fromJSON(args[[1L]], simplifyVector = FALSE)
out <- args[[2L]]
lib <- normalizePath(args[[3L]], mustWork = TRUE)
.libPaths(c(lib, .libPaths()))
if (normalizePath(find.package("simCoDaR"), mustWork = TRUE) !=
    normalizePath(file.path(lib, "simCoDaR"), mustWork = TRUE)) stop("wrong package library")
library(simCoDaR, lib.loc = lib)
if (as.character(packageVersion("simCoDaR")) != "0.1.0.9000") stop("unexpected simCoDaR version")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
options(digits = 17)

phases <- c("Quartz", "Anatase", "Rutile")
channels <- c("SiO2", "TiO2")
C <- rbind(c(1, 0, 0), c(0, 1, 1))
dimnames(C) <- list(channels, phases)
lib_local <- mineral_assay_library()
lib_local$Anatase <- list(formula = "TiO2", atoms = c(Ti = 1, O = 2))
lib_local$Rutile <- list(formula = "TiO2", atoms = c(Ti = 1, O = 2))
eye <- diag(3)
colnames(eye) <- phases
pure <- minerals_to_assay(eye, library = lib_local,
                          basis = "fraction", iron_reporting = "total_Fe2O3")
converted <- t(pure$oxides_reported[, channels, drop = FALSE]) / 100
chemistry_error <- max(abs(converted - C))
if (!is.finite(chemistry_error) || chemistry_error > 1e-12) stop("chemistry conversion mismatch")

dir.create(file.path(out, "data"), recursive = TRUE, showWarnings = FALSE)
f_cal_obs <- file.path(out, "data", "calibration_observed.csv")
f_cal_truth <- file.path(out, "data", "calibration_truth.csv")
f_eval_obs <- file.path(out, "data", "evaluation_observed.csv")
f_eval_truth <- file.path(out, "data", "evaluation_truth.csv")
if (any(file.exists(c(f_cal_obs, f_cal_truth, f_eval_obs, f_eval_truth))))
  stop("refusing to overwrite generated data")
head_phase <- c("cell_id", "cohort_id", "sample_id", phases)
head_assay <- c("cell_id", "cohort_id", "sample_id", channels)
write.table(t(head_phase), f_cal_obs, sep = ",", row.names = FALSE,
            col.names = FALSE, quote = TRUE)
write.table(t(head_phase), f_cal_truth, sep = ",", row.names = FALSE,
            col.names = FALSE, quote = TRUE)
write.table(t(head_assay), f_eval_obs, sep = ",", row.names = FALSE,
            col.names = FALSE, quote = TRUE)
write.table(t(head_phase), f_eval_truth, sep = ",", row.names = FALSE,
            col.names = FALSE, quote = TRUE)

append_rows <- function(path, meta, values) {
  tab <- data.frame(meta, values, check.names = FALSE)
  write.table(tab, path, sep = ",", row.names = FALSE, col.names = FALSE,
              append = TRUE, quote = TRUE)
}
draw_truth <- function(center, n, seed, fixed_q) {
  endmembers <- matrix(as.numeric(unlist(center)), nrow = 1,
                       dimnames = list("domain1", phases))
  spec <- new_assay_scenario(endmembers, weights = 1,
                             logratio_sd = 0.35,
                             structural_parts = NULL,
                             composition_bounds = NULL)
  d <- simulate_mineral_assay_scenario(spec, n = n, seed = seed,
          library = lib_local, iron_reporting = "total_Fe2O3",
          max_draw_attempts = 1000L)
  phase_data <- d$truth$mineral_fraction[, phases, drop = FALSE]
  if (!all(vapply(phase_data, is.numeric, logical(1))))
    stop("non-numeric synthetic truth phases")
  x <- as.matrix(phase_data)
  storage.mode(x) <- "double"
  if (!is.null(fixed_q)) {
    ar <- x[, "Anatase"] + x[, "Rutile"]
    if (any(ar <= 0)) stop("invalid fixed-Quartz denominator")
    x[, "Anatase"] <- (1 - fixed_q) * x[, "Anatase"] / ar
    x[, "Rutile"] <- (1 - fixed_q) * x[, "Rutile"] / ar
    x[, "Quartz"] <- fixed_q
  }
  if (any(!is.finite(x)) || max(abs(rowSums(x) - 1)) > 1e-12)
    stop("invalid synthetic truth")
  x
}
draw_noise <- function(n, p, sd, seed) {
  set.seed(seed)
  matrix(rnorm(n * p, sd = sd), nrow = n, ncol = p)
}

master <- as.integer(config$rng$master_seed)
for (cell in config$cells) {
  s <- as.integer(cell$cell_id)
  n_cal <- as.integer(cell$calibration_n)
  q_sd <- as.numeric(cell$qxrd_sd)
  fixed_q <- if (is.null(cell$fixed_quartz)) NULL else as.numeric(cell$fixed_quartz)
  for (r in seq_len(as.integer(config$replicates))) {
    seeds <- master + 100000L * s + 100L * r + 1:4
    cal <- draw_truth(cell$calibration_center, n_cal, seeds[[1L]], fixed_q)
    eval <- draw_truth(config$evaluation_center, as.integer(config$evaluation_n),
                       seeds[[2L]], fixed_q)
    y <- cal + draw_noise(n_cal, 3L, q_sd, seeds[[3L]])
    a <- eval %*% t(C) + draw_noise(nrow(eval), 2L,
                                   as.numeric(config$assay_sd), seeds[[4L]])
    colnames(y) <- phases
    colnames(a) <- channels
    cal_meta <- data.frame(cell_id = s, cohort_id = r,
       sample_id = sprintf("s%d_r%d_cal_%d", s, r, seq_len(n_cal)))
    eval_meta <- data.frame(cell_id = s, cohort_id = r,
       sample_id = sprintf("s%d_r%d_eval_%d", s, r, seq_len(nrow(eval))))
    append_rows(f_cal_obs, cal_meta, y)
    append_rows(f_cal_truth, cal_meta, cal)
    append_rows(f_eval_obs, eval_meta, a)
    append_rows(f_eval_truth, eval_meta, eval)
  }
}
cat(sprintf("generation_complete cells=%d cohorts=%d chemistry_error=%.17g rng=%s\n",
            length(config$cells), length(config$cells) * as.integer(config$replicates),
            chemistry_error, paste(RNGkind(), collapse = "/")))
