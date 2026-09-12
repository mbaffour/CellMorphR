# =============================================================================
# One transform pipeline, always reported (CM-13, CM-14).
# =============================================================================

test_that("apply_measure_transform reports the log shift it applied (CM-13)", {
  none <- apply_measure_transform(c(1, 2, 3), "none")
  expect_equal(none$values, c(1, 2, 3))
  expect_equal(none$shift, 0)

  pos <- apply_measure_transform(c(1, 10, 100), "log")
  expect_equal(pos$values, log(c(1, 10, 100)))
  expect_equal(pos$shift, 0)
  expect_equal(pos$label, "log(x)")

  # A data-dependent shift is the part users have to be told about: two datasets
  # get different shifts, so the effect sizes are not comparable between them.
  neg <- apply_measure_transform(c(-5, 0, 10), "log")
  expect_equal(neg$shift, 6)
  expect_equal(neg$values, log(c(1, 6, 16)))
  expect_match(neg$label, "log(x + 6)", fixed = TRUE)
  expect_true(all(is.finite(neg$values)))   # the v3.1.1 bug: log of a negative
})

test_that("the sqrt transform floors at zero", {
  s <- apply_measure_transform(c(-4, 0, 9), "sqrt")
  expect_equal(s$values, c(0, 0, 3))
  expect_equal(s$label, "sqrt(max(x, 0))")
})

test_that("apply_transform_to_columns records every column it touched", {
  df <- data.frame(Condition = "A", Area = c(1, 10, 100), Perimeter = c(-1, 0, 3),
                   Replicate = "R1", stringsAsFactors = FALSE)
  out <- apply_transform_to_columns(df, c("Area", "Perimeter"), "log")
  info <- attr(out, "transform_info")
  expect_equal(sort(names(info)), c("Area", "Perimeter"))
  expect_equal(info$Area$shift, 0)
  expect_equal(info$Perimeter$shift, 2)
  expect_equal(describe_transform(out, "Perimeter"), info$Perimeter$label)
  expect_equal(describe_transform(out, "Area"), "log(x)")
})

test_that("fit_lmm_safe reports the transform it used, shift included (CM-13)", {
  skip_if_not(has_lme4, "lme4 not installed")
  set.seed(35)
  d <- make_cellmorph_data(ncell = 30, times = c(0, 30, 60))
  d$Area <- d$Area - 12          # force non-positive values
  r <- fit_lmm_safe(d, "Area", transform = "log")
  expect_equal(r$transform_used, "log")
  expect_gt(r$transform_shift, 0)
  expect_match(r$transform_label, "log(x + ", fixed = TRUE)
})
