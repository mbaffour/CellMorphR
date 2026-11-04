# =============================================================================
# Comparison modes must report what was built, not what was requested (CM-9).
# =============================================================================

test_that("build_comparison_table reports the mode it actually built (CM-9)", {
  set.seed(36)
  d <- make_cellmorph_data(ncell = 20, times = c(0, 30, 60))

  with_time <- compute_flexible_replicate_summary(d, "Area", c("Condition", "Time_min"))
  ct <- build_comparison_table(with_time, "within_time")
  expect_equal(realised_comparison_mode(ct), "within_time")
  expect_null(attr(ct, "mode_message"))
  expect_equal(nrow(ct), 3)                       # one Ctrl-vs-Trt pair per time

  # Time_min unchecked: the stratified mode is not available. It must say so,
  # not silently report the mode the user asked for.
  no_time <- compute_flexible_replicate_summary(d, "Area", c("Condition"))
  ct2 <- build_comparison_table(no_time, "within_time")
  expect_equal(attr(ct2, "requested_mode"), "within_time")
  expect_equal(realised_comparison_mode(ct2), "vs_reference")
  expect_match(attr(ct2, "mode_message"), "needs Time_min")
})

test_that("build_comparison_table can be asked to error instead of degrading", {
  set.seed(37)
  d <- make_cellmorph_data(ncell = 20, times = c(0, 30))
  no_time <- compute_flexible_replicate_summary(d, "Area", c("Condition"))
  expect_error(build_comparison_table(no_time, "within_time", strict = TRUE),
               "needs Time_min")
})
