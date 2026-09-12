# =============================================================================
# Standardised effect size (CM-3, CM-4).
# =============================================================================

test_that("hedges_g uses the (n-1)-weighted pooled SD, valid for unequal n", {
  set.seed(5)
  ctrl <- rnorm(3, 10, 1); trt <- rnorm(6, 11, 3)
  r <- hedges_g(ctrl, trt)
  n1 <- 3; n2 <- 6
  expected_sd <- sqrt(((n1 - 1) * var(ctrl) + (n2 - 1) * var(trt)) / (n1 + n2 - 2))
  expect_equal(r$pooled_sd, expected_sd)
  # The old formula -- the unweighted average of the variances -- differs here.
  old_sd <- sqrt((var(ctrl) + var(trt)) / 2)
  expect_false(isTRUE(all.equal(old_sd, expected_sd)))
})

test_that("hedges_g applies the small-sample correction", {
  set.seed(5)
  ctrl <- rnorm(3, 10, 1); trt <- rnorm(3, 11, 1)
  r <- hedges_g(ctrl, trt)
  expect_equal(r$J, 1 - 3 / (4 * 4 - 1))          # J = 0.8 at n = 3 vs 3
  expect_equal(r$g, r$J * r$d)
  expect_lt(abs(r$g), abs(r$d))                    # always shrinks towards zero
  expect_equal(round(r$J, 3), 0.800)
})

test_that("the confidence interval brackets the estimate and is wide at n = 3", {
  set.seed(5)
  ctrl <- rnorm(3, 10, 1); trt <- rnorm(3, 11, 1)
  r <- hedges_g(ctrl, trt)
  expect_lt(r$ci_low, r$g)
  expect_gt(r$ci_high, r$g)
  # With three dishes per group the interval must not let a user claim
  # "large effect" from the point estimate alone.
  expect_gt(r$ci_high - r$ci_low, 2)
  expect_match(r$note, "very wide")
})

test_that("the confidence interval agrees with a known reference case", {
  # n = 30 per group, d = 1.0 exactly: the noncentral-t 95% interval for d is
  # approximately (0.45, 1.54) (Cumming 2012, Table 11.1 worked example).
  set.seed(1)
  ctrl <- scale(rnorm(30))[, 1]              # mean 0, sd 1 exactly
  trt  <- scale(rnorm(30))[, 1] + 1          # mean 1, sd 1 exactly
  r <- hedges_g(ctrl, trt)
  expect_equal(round(r$d, 3), 1.000)
  expect_lt(abs(r$ci_low - 0.45), 0.05)
  expect_lt(abs(r$ci_high - 1.54), 0.05)
})

test_that("hedges_g returns NA with a message instead of crashing (CM-4)", {
  # An absent condition at a timepoint -- the ordinary "drug added after the
  # first frame" design -- used to reach `if (NA > 0)` and stop the reactive.
  empty <- hedges_g(c(1, 2, 3), numeric(0))
  expect_true(is.na(empty$g))
  expect_match(empty$note, "at least 2 biological replicates")

  single <- hedges_g(c(1, 2, 3), 5)
  expect_true(is.na(single$g))
  expect_match(single$note, "at least 2 biological replicates")

  identical_reps <- hedges_g(c(2, 2, 2), c(2, 2, 2))
  expect_true(is.na(identical_reps$g))
  expect_match(identical_reps$note, "pooled standard deviation is zero")
})

test_that("hedges_g is signed treatment-minus-control", {
  expect_gt(hedges_g(c(1, 2, 3), c(5, 6, 7))$g, 0)
  expect_lt(hedges_g(c(5, 6, 7), c(1, 2, 3))$g, 0)
})
