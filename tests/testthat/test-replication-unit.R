# =============================================================================
# The unit of replication is the biological replicate, not the cell.
# These tests pin down the behaviour the tool is named for.
# =============================================================================

test_that("compute_replicate_summaries returns one row per condition x time x replicate", {
  set.seed(1)
  d <- make_cellmorph_data(ncell = 20, times = c(0, 30, 60))
  rs <- compute_replicate_summaries(d, "Area")
  expect_equal(nrow(rs), 2 * 3 * 3)
  expect_true(all(rs$n_cells == 20))
})

test_that("compute_flexible_replicate_summary collapses cells to replicates", {
  set.seed(2)
  d <- make_cellmorph_data(ncell = 50, times = c(0, 30, 60))
  rf <- compute_flexible_replicate_summary(d, "Area", c("Condition"))
  expect_equal(nrow(rf), 6)              # 2 conditions x 3 replicates
  expect_true(all(rf$n_cells == 50 * 3)) # 3 timepoints pooled
  expect_equal(rf$Stat_Value, rf$median_val)
})

test_that("compute_condition_summaries computes SEM over replicates, not cells", {
  set.seed(3)
  d <- make_cellmorph_data(ncell = 200, times = c(0, 30))
  cs <- compute_condition_summaries(compute_replicate_summaries(d, "Area"))
  expect_true(all(cs$n_reps == 3))
})

test_that("min_wilcox_p reports the minimum attainable two-sided p", {
  expect_equal(round(min_wilcox_p(2, 2), 4), 0.3333)
  expect_equal(round(min_wilcox_p(3, 3), 4), 0.1000)
  expect_equal(round(min_wilcox_p(4, 4), 4), 0.0286)
})

test_that("run_pairwise_tests reports biological replicate n alongside cell n", {
  set.seed(4)
  d <- make_cellmorph_data(ncell = 100, times = c(0, 30, 60))
  rf <- compute_flexible_replicate_summary(d, "Area", c("Condition"))
  ct <- build_comparison_table(rf, "vs_reference", "Condition=Ctrl")
  pw <- run_pairwise_tests(rf, ct, correction = "holm")
  expect_equal(nrow(pw), 1)
  expect_equal(pw$N_BioRep_Group_1, 3)
  expect_equal(pw$N_BioRep_Group_2, 3)
  expect_equal(pw$N_Cells_Group_1, 900)  # 3 replicates x 3 timepoints x 100 cells
  expect_true(pw$df < 5)  # Welch df from n=3 vs n=3, not from 300 cells
})

test_that("run_pairwise_tests holds its nominal type-I error on null data", {
  # The reference point for the whole tool. A naive cell-level test on the same
  # data rejects ~90% of the time; the replicate-level test must be near 5%.
  n_sims <- if (nzchar(Sys.getenv("CELLMORPHR_FAST"))) 100 else 400
  set.seed(7)
  p_app <- numeric(n_sims); p_naive <- numeric(n_sims)
  for (i in seq_len(n_sims)) {
    d <- make_cellmorph_data(cond_effect = 0, rep_sd = 1, cell_sd = 1, ncell = 100)
    p_naive[i] <- stats::t.test(Area ~ Condition, data = d)$p.value
    rf <- compute_flexible_replicate_summary(d, "Area", c("Condition"))
    ct <- build_comparison_table(rf, "vs_reference", "Condition=Ctrl")
    p_app[i] <- run_pairwise_tests(rf, ct, correction = "holm")$p_adj[1]
  }
  rate_app <- mean(p_app < 0.05, na.rm = TRUE)
  rate_naive <- mean(p_naive < 0.05, na.rm = TRUE)
  cat(sprintf("\n  [type-I] pairwise (replicate-level): %.1f%%   naive (cell-level): %.1f%%\n",
              100 * rate_app, 100 * rate_naive))
  # Welch on n=3 vs n=3 is slightly conservative (the auditor measured 2.8%), so
  # the requirement is "not anti-conservative", not "exactly 5%".
  expect_lt(rate_app, 0.05 + type1_tolerance(n_sims))
  expect_gt(rate_naive, 0.5)  # documents what the tool is protecting against
})
