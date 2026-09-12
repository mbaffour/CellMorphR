# =============================================================================
# Assumption diagnostics: Levene's test and the variance ratio (CM-12, CM-18).
# They are advisory -- nothing in the app switches test on the strength of them.
# =============================================================================

test_that("levene_test matches car::leveneTest where car is available", {
  skip_if_not_installed("car")
  set.seed(31)
  y <- c(rnorm(10, 0, 1), rnorm(10, 0, 3), rnorm(10, 0, 1))
  g <- factor(rep(c("a", "b", "c"), each = 10))
  ours <- levene_test(y, g)
  theirs <- car::leveneTest(y, g, center = median)
  expect_equal(ours$statistic, theirs[1, "F value"], tolerance = 1e-8)
  expect_equal(ours$p.value, theirs[1, "Pr(>F)"], tolerance = 1e-8)
  expect_equal(ours$df1, theirs[1, "Df"])
})

test_that("levene_test detects real heteroscedasticity and ignores none", {
  set.seed(32)
  hetero <- levene_test(c(rnorm(40, 0, 1), rnorm(40, 0, 6)), rep(c("a", "b"), each = 40))
  expect_lt(hetero$p.value, 0.01)
  homo <- levene_test(c(rnorm(40, 0, 1), rnorm(40, 0, 1)), rep(c("a", "b"), each = 40))
  expect_gt(homo$p.value, 0.05)
})

test_that("levene_test declines rather than erroring on degenerate input", {
  expect_null(levene_test(c(1, 2, 3), rep("a", 3)))          # one group
  expect_null(levene_test(c(1, 2, 3), c("a", "b", "b")))      # a group of 1
  expect_null(levene_test(rep(5, 10), rep(c("a", "b"), 5)))   # no spread
})

test_that("variance_ratio reports the largest/smallest group SD", {
  set.seed(33)
  v <- variance_ratio(c(rnorm(30, 0, 1), rnorm(30, 0, 5)), rep(c("a", "b"), each = 30))
  expect_gt(v$ratio, 3)
  expect_equal(v$max_group, "b")
})

test_that("assumption_diagnostics states that it is advisory", {
  set.seed(34)
  d <- make_cellmorph_data(ncell = 30, times = c(0, 30))
  rf <- compute_flexible_replicate_summary(d, "Area", c("Condition"))
  diag <- assumption_diagnostics(rf)
  expect_false(is.null(diag$levene))
  expect_match(diag$advisory, "ADVISORY")
  txt <- format_assumption_diagnostics(diag)
  expect_match(txt, "Shapiro-Wilk")
  expect_match(txt, "Brown & Forsythe", fixed = TRUE)
  expect_match(txt, "does not switch test")
})
