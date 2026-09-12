# =============================================================================
# Test harness for CellMorphR
#
# The statistical core of CellMorphR_v3.R (everything above the
# "END OF STATISTICAL CORE" marker) is plain R with no Shiny reactivity, so it
# can be sourced standalone and driven directly. This helper does that.
# =============================================================================

cellmorphr_root <- function() {
  # Works whether tests are run from the repo root, from tests/, or by testthat.
  candidates <- c(".", "..", "../..", dirname(sys.frame(1)$ofile %||% "."))
  for (p in candidates) {
    if (file.exists(file.path(p, "CellMorphR_v3.R"))) return(normalizePath(p))
  }
  stop("Could not locate CellMorphR_v3.R -- run the tests from the repository root.")
}

if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (!is.null(a)) a else b

MARKER <- "^# ---- END OF STATISTICAL CORE"

#' Source the statistical core of the app into the calling environment.
load_cellmorphr_core <- function(envir = globalenv()) {
  app_file <- file.path(cellmorphr_root(), "CellMorphR_v3.R")
  lines <- readLines(app_file, warn = FALSE)
  end <- grep(MARKER, lines)
  if (length(end) != 1L) {
    stop("Expected exactly one 'END OF STATISTICAL CORE' marker in CellMorphR_v3.R, found ",
         length(end), ". The test harness and the app have drifted apart.")
  }
  suppressWarnings(suppressMessages(
    eval(parse(text = paste(lines[seq_len(end - 1L)], collapse = "\n")), envir = envir)
  ))
  invisible(end)
}

# -----------------------------------------------------------------------------
# Synthetic data generators used by several test files.
# -----------------------------------------------------------------------------

#' Balanced time-course dataset.
#'
#' @param cond_effect true additive Condition effect (0 => null data)
#' @param rep_sd      SD of the per-dish (replicate) offset, constant over time
#' @param reptime_sd  SD of the per-dish-per-timepoint batch effect (focus drift,
#'                    evaporation, per-session segmentation quality)
#' @param cell_sd     SD of the within-dish cell noise
make_cellmorph_data <- function(cond_effect = 0, rep_sd = 1.0, reptime_sd = 0,
                                cell_sd = 1.0, ncell = 100,
                                conds = c("Ctrl", "Trt"),
                                reps = c("R1", "R2", "R3"),
                                times = c(0, 30, 60, 90, 120),
                                interaction_effect = 0) {
  g <- expand.grid(Cell = seq_len(ncell), Time_min = times, Replicate = reps,
                   Condition = conds, stringsAsFactors = FALSE)
  key <- paste(g$Condition, g$Replicate)
  rtk <- paste(key, g$Time_min)
  offs <- stats::setNames(stats::rnorm(length(unique(key)), 0, rep_sd), unique(key))
  rt   <- stats::setNames(stats::rnorm(length(unique(rtk)), 0, reptime_sd), unique(rtk))
  g$Area <- 10 + offs[key] + rt[rtk] +
    ifelse(g$Condition == conds[2], cond_effect, 0) +
    ifelse(g$Condition == conds[2], interaction_effect * (g$Time_min / max(times, 1)), 0) +
    stats::rnorm(nrow(g), 0, cell_sd)
  g$Cell_ID <- seq_len(nrow(g))
  rownames(g) <- NULL
  g
}

#' Assert a measured type-I error rate is compatible with the nominal alpha.
#'
#' Uses a binomial tolerance so the suite is not flaky: with n_sims replicates
#' of a Bernoulli(alpha) the observed rate has SD sqrt(a(1-a)/n). We allow
#' `sd_tol` standard deviations, plus an absolute floor.
type1_tolerance <- function(n_sims, alpha = 0.05, sd_tol = 3.5, floor_abs = 0.02) {
  max(floor_abs, sd_tol * sqrt(alpha * (1 - alpha) / n_sims))
}

# Load the app's statistical core once, for every test file in this directory.
load_cellmorphr_core(globalenv())
