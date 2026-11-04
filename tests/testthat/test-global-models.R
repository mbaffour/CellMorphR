# =============================================================================
# The two global models. Both were pseudoreplicated before v3.2.0:
#   * the two-way ANOVA counted the same dish once per timepoint  (CM-1)
#   * the LMM tested Time and Condition:Time on cell-level df     (CM-2)
# These tests re-run the auditor's null simulations against the corrected code.
# =============================================================================

fast <- nzchar(Sys.getenv("CELLMORPHR_FAST"))

test_that("fit_rm_anova tests Condition against dish-to-dish variation", {
  set.seed(11)
  d <- make_cellmorph_data(ncell = 50, times = c(0, 30, 60))
  fit <- fit_rm_anova(compute_replicate_summaries(d, "Area"))
  expect_true(fit$ok)
  tb <- fit$table
  expect_setequal(tb$Term, c("Condition", "Time_fac", "Condition:Time_fac"))
  # 6 dishes, 2 conditions -> Condition has 4 denominator df, NOT 12 (the number
  # of replicate-summary rows minus the model terms, which is what aov() used).
  expect_equal(tb$DenDF[tb$Term == "Condition"], 4)
  expect_equal(tb$Error_Stratum[tb$Term == "Condition"], "Rep_unique")
  expect_equal(tb$Error_Stratum[tb$Term == "Condition:Time_fac"], "Within")
})

test_that("fit_rm_anova refuses to invent a global model it cannot fit", {
  set.seed(12)
  d <- make_cellmorph_data(ncell = 20, reps = "R1", times = c(0, 30))
  fit <- fit_rm_anova(compute_replicate_summaries(d, "Area"))
  expect_false(fit$ok)
  expect_match(fit$msg, "at least two replicates")
})

test_that("fit_rm_anova collapses to a one-way ANOVA with a single timepoint", {
  set.seed(13)
  d <- make_cellmorph_data(ncell = 30, times = 0)
  fit <- fit_rm_anova(compute_replicate_summaries(d, "Area"))
  expect_true(fit$ok)
  expect_false(fit$repeated)
  expect_equal(fit$table$Term, "Condition")
  expect_equal(fit$table$DenDF, 4)
})

test_that("run_anova_for_feature keeps its contract for the Multi-Feature tab", {
  set.seed(14)
  d <- make_cellmorph_data(ncell = 30, times = c(0, 30, 60))
  r <- run_anova_for_feature(d, "Area")
  expect_equal(names(r), c("Feature", "Condition_p", "Time_p", "Interaction_p",
                           "Condition_partial_eta2", "Time_partial_eta2",
                           "Interaction_partial_eta2"))
  expect_true(all(is.finite(c(r$Condition_p, r$Time_p, r$Interaction_p))))
})

test_that("repeated-measures ANOVA holds its type-I error on null data (CM-1)", {
  # Before the fix this simulation gave 41.0% for the Condition term.
  n_sims <- if (fast) 120 else 400
  set.seed(101)
  pC <- pT <- pI <- numeric(n_sims)
  for (i in seq_len(n_sims)) {
    d <- make_cellmorph_data(cond_effect = 0, rep_sd = 1, cell_sd = 1, ncell = 100)
    a <- run_anova_for_feature(d, "Area")
    pC[i] <- a$Condition_p; pT[i] <- a$Time_p; pI[i] <- a$Interaction_p
  }
  rate <- function(p) mean(p < 0.05, na.rm = TRUE)
  cat(sprintf("\n  [type-I] RM-ANOVA Condition %.1f%%  Time %.1f%%  Condition:Time %.1f%%\n",
              100 * rate(pC), 100 * rate(pT), 100 * rate(pI)))
  tol <- type1_tolerance(n_sims)
  expect_lt(abs(rate(pC) - 0.05), tol)
  expect_lt(abs(rate(pT) - 0.05), tol)
  expect_lt(abs(rate(pI) - 0.05), tol)
})

test_that("repeated-measures ANOVA still holds when dishes drift over time", {
  # The realistic time-lapse case: each dish has its own per-session offset.
  n_sims <- if (fast) 120 else 300
  set.seed(102)
  pC <- pT <- pI <- numeric(n_sims)
  for (i in seq_len(n_sims)) {
    d <- make_cellmorph_data(cond_effect = 0, rep_sd = 1, reptime_sd = 0.7,
                             cell_sd = 1, ncell = 100)
    a <- run_anova_for_feature(d, "Area")
    pC[i] <- a$Condition_p; pT[i] <- a$Time_p; pI[i] <- a$Interaction_p
  }
  rate <- function(p) mean(p < 0.05, na.rm = TRUE)
  cat(sprintf("  [type-I] RM-ANOVA + dish drift: Condition %.1f%%  Time %.1f%%  Condition:Time %.1f%%\n",
              100 * rate(pC), 100 * rate(pT), 100 * rate(pI)))
  tol <- type1_tolerance(n_sims)
  expect_lt(abs(rate(pC) - 0.05), tol)
  expect_lt(abs(rate(pT) - 0.05), tol)
  expect_lt(abs(rate(pI) - 0.05), tol)
})

test_that("the repeated-measures ANOVA can still find a real effect", {
  set.seed(103)
  n_sims <- if (fast) 40 else 100
  hit <- vapply(seq_len(n_sims), function(i) {
    d <- make_cellmorph_data(cond_effect = 0, interaction_effect = 4,
                             rep_sd = 0.5, cell_sd = 1, ncell = 100)
    run_anova_for_feature(d, "Area")$Interaction_p < 0.05
  }, logical(1))
  expect_gt(mean(hit), 0.8)
})

# -----------------------------------------------------------------------------

test_that("fit_lmm_safe adds the dish x timepoint random term", {
  skip_if_not(has_lme4, "lme4 not installed")
  set.seed(21)
  d <- make_cellmorph_data(ncell = 40, times = c(0, 30, 60))
  r <- fit_lmm_safe(d, "Area")
  expect_true(r$fitted)
  expect_true(r$has_reptime_term)
  expect_match(r$random_structure, "Rep_unique:Time_fac", fixed = TRUE)
})

test_that("a singular LMM fit is reported, not swapped for a phantom fallback (CM-7)", {
  skip_if_not(has_lme4, "lme4 not installed")
  set.seed(22)
  d <- make_cellmorph_data(ncell = 40, rep_sd = 0, reptime_sd = 0, times = c(0, 30, 60))
  r <- fit_lmm_safe(d, "Area")
  expect_true(r$fitted)
  expect_false(is.null(r$model))
  expect_true(r$singular)
  expect_false(r$fallback_used)
  expect_true(any(grepl("Singular fit", r$warnings)))
  # The old code promised an ANOVA that was never computed.
  expect_false(any(grepl("Falling back", r$warnings)))
})

test_that("fit_lmm_safe declines rather than promising a model it will not fit", {
  skip_if_not(has_lme4, "lme4 not installed")
  set.seed(23)
  d <- make_cellmorph_data(ncell = 20, reps = "R1", times = c(0, 30))  # 2 dishes
  r <- fit_lmm_safe(d, "Area")
  expect_null(r$model)
  expect_false(r$fitted)
  expect_match(r$method, "No global model")
  expect_match(paste(r$warnings, collapse = " "), "Pairwise results above are still valid")
})

test_that("the LMM holds its type-I error on null data (CM-2)", {
  skip_if_not(has_lme4 && has_lmerTest, "lme4 + lmerTest needed for LMM p-values")
  # Before the fix: Time 99.2% and Condition:Time 100.0% under dish drift.
  n_sims <- if (fast) 60 else 150
  for (scen in list(list(nm = "dish intercept only", rt = 0),
                    list(nm = "+ dish x time drift", rt = 0.7))) {
    set.seed(201)
    M <- matrix(NA_real_, n_sims, 3)
    for (i in seq_len(n_sims)) {
      d <- make_cellmorph_data(cond_effect = 0, rep_sd = 1, reptime_sd = scen$rt,
                               cell_sd = 1, ncell = 100)
      r <- fit_lmm_safe(d, "Area")
      if (is.null(r$model)) next
      a <- tryCatch(anova(r$model), error = function(e) NULL)
      if (is.null(a) || is.null(a[["Pr(>F)"]])) next
      M[i, ] <- a[["Pr(>F)"]][1:3]
    }
    ok <- stats::complete.cases(M)
    rate <- function(j) mean(M[ok, j] < 0.05)
    cat(sprintf("  [type-I] LMM (%s): Condition %.1f%%  Time %.1f%%  Condition:Time %.1f%%  (n=%d)\n",
                scen$nm, 100 * rate(1), 100 * rate(2), 100 * rate(3), sum(ok)))
    tol <- type1_tolerance(sum(ok))
    expect_lt(abs(rate(1) - 0.05), tol)
    expect_lt(abs(rate(2) - 0.05), tol)
    expect_lt(abs(rate(3) - 0.05), tol)
  }
})
