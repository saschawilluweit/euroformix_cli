#!/usr/bin/env Rscript

# nicely load required packages
for (package_name in c("docopt", "euroformix")) {
  if (!require(package_name, character.only = TRUE)) {
    print(paste0("Packege '", package_name, "' is needed. Please install using: R CMD INSTALL ", package_name))
    exit(1)
  }
  suppressPackageStartupMessages(library(package_name, quietly = TRUE, character.only = TRUE))
}

cmd <- "Rscript"
args <- c("./euroformix_cli.R", "--quiet", "--kit_database_file", "tests/one_marker_kit_database.txt", "--kit", "TESTKIT1", "--fst", "0.00", "--analysis_threshold", "50", "--queried_file", "tests/one_marker_queried.txt", "--db_file", "tests/one_marker_allele_frequencies.txt", "--mle_steps", "10", "--mle_iterations", "1000")

additional_args <- list(
  none = c(),
  stutter = c("--stutter"),
  stutter_deg = c("--stutter", "--degradation"),
  stutter_dropin = c("--stutter", "--dropin", "--dropin_probability", "0.01"),
  stutter_deg_dropin = c("--stutter", "--degradation", "--dropin", "--dropin_probability", "0.01")
)

knowns <- list(
  none = c("--unknowns", "1"),
  k10 = c("--known_file", "tests/one_marker_known_K10.txt", "--unknowns", "0"),
  k20 = c("--known_file", "tests/one_marker_known_K20.txt", "--unknowns", "0")
)

for (sample in c("S001", "S002", "S003")) {
  results_file <- file(paste0("results/one_marker_", sample, ".txt"), "w")

  write(paste("mode", "known", paste0(read.table("tests/one_marker_queried.txt", sep = ",", header = TRUE)[,1], collapse = ";"), sep = ";"), file = results_file)
  for (known in names(knowns)) {
    for (additional_arg in names(additional_args)) {
      lrs <- system2(cmd, c(args, additional_args[[additional_arg]], knowns[[known]], "--sample_file", paste0("tests/one_marker_sample_", sample, ".txt")), stdout = TRUE)
      write(paste(additional_arg, known, paste0(unlist(lrs), collapse = ";"), sep = ";"), file = results_file, append = TRUE)
    }
  }
}