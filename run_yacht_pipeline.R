#!/usr/bin/env Rscript

rm(list = ls())

pipeline_scripts <- c(
  "R/00_setup.R",
  "R/01_loglog_prior.R",
  "R/02_bas_bma.R",
  "R/03_gamma_model.R",
  "R/04_posterior_predictive_validation.R"
)

for (script in pipeline_scripts) {
  message("\n--- Running ", script, " ---")
  source(script, echo = FALSE, chdir = FALSE)
}

message("\nPipeline completed successfully.")
