#!/usr/bin/env Rscript
# =============================================================================
# CellMorphR test suite entry point.
#
#   Rscript tests/run.R              # run everything
#   Rscript tests/run.R stats        # run only test files matching "stats"
#
# The type-I error simulations are the slow part (a few minutes). Set
# CELLMORPHR_FAST=1 to run them with fewer simulations for a quick check --
# but the committed numbers are the full-length ones.
# =============================================================================

if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("The test suite needs the 'testthat' package: install.packages(\"testthat\")")
}

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])), ".."))
if (!file.exists(file.path(root, "CellMorphR_v3.R"))) root <- normalizePath(".")

args <- commandArgs(trailingOnly = TRUE)
filter <- if (length(args) > 0) args[1] else NULL

cat("CellMorphR test suite\n")
cat("  repo:   ", root, "\n")
cat("  R:      ", R.version.string, "\n")
cat("  lme4:   ", requireNamespace("lme4", quietly = TRUE),
    " lmerTest: ", requireNamespace("lmerTest", quietly = TRUE),
    " emmeans: ", requireNamespace("emmeans", quietly = TRUE), "\n")
cat("  fast:   ", !identical(Sys.getenv("CELLMORPHR_FAST"), ""), "\n\n")

res <- testthat::test_dir(file.path(root, "tests", "testthat"),
                          filter = filter, stop_on_failure = TRUE,
                          reporter = "summary")
invisible(res)
