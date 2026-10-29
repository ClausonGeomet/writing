#!/usr/bin/env Rscript
# Self-contained teaching example, not the six-scenario research execution.
# simCoDaR: d9a0b27676a0cb3a4acb300a3f6b2ebc845d205d (0.1.0.9000).
suppressPackageStartupMessages(library(simCoDaR))
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[[1]] else "."
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

phases <- c("Quartz", "Hematite", "Kaolinite", "Forsterite")
analytes <- c("SiO2", "Fe2O3_total", "Al2O3", "MgO")
centres <- rbind(
  siliceous = c(Quartz = .55, Hematite = .15, Kaolinite = .15, Forsterite = .15),
  mixed = c(Quartz = .25, Hematite = .25, Kaolinite = .25, Forsterite = .25)
)
scenario <- new_assay_scenario(centres, weights = c(1, 1), logratio_sd = .25)
draw <- simulate_mineral_assay_scenario(
  scenario, n = 60L, seed = 20261026L, iron_reporting = "total_Fe2O3"
)
M <- as.matrix(draw$truth$mineral_fraction[, phases, drop = FALSE])
A <- as.matrix(draw$observed$xrf_oxide_wt_percent[, analytes, drop = FALSE]) / 100

# Pure-mineral responses give the same assumed chemistry as the generator.
pure <- diag(length(phases))
dimnames(pure) <- list(phases, phases)
reported <- minerals_to_assay(pure, library = mineral_assay_library(),
                              basis = "fraction", iron_reporting = "total_Fe2O3")
C <- t(as.matrix(reported$oxides_reported[, analytes, drop = FALSE])) / 100
stopifnot(max(abs(rowSums(M) - 1)) < 1e-10)
stopifnot(max(abs(M %*% t(C) - A)) < 1e-10)

rows <- data.frame(
  sample_id = draw$sample_id,
  domain = as.character(draw$domains[draw$sample_id, "domain"]),
  check.names = FALSE
)
for (j in seq_along(phases)) rows[[paste0("m_", phases[j])]] <- M[, j]
for (i in seq_along(analytes)) rows[[paste0("a_", analytes[i])]] <- A[, i]
write.csv(rows, file.path(out_dir, "emc_data.csv"), row.names = FALSE, quote = FALSE)
write.csv(data.frame(analyte = analytes, C, check.names = FALSE),
          file.path(out_dir, "emc_C.csv"), row.names = FALSE, quote = FALSE)
cat(sprintf("Generated %d hypothetical specimens; max forward error = %.2e\n",
            nrow(M), max(abs(M %*% t(C) - A))))
