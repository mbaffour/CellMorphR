# =============================================================================
# Input validation. Real ImageJ / CellProfiler exports are messy: "n/a" in a
# measurement column, "30 min" in the time column, a single replicate.
# =============================================================================

base_df <- function() {
  set.seed(99)
  data.frame(Condition = rep(c("A", "B"), each = 6),
             Time_min  = rep(c(0, 30, 60), 4),
             Replicate = rep(c("R1", "R2"), 6),
             Area      = rnorm(12, 10),
             stringsAsFactors = FALSE)
}

test_that("a clean data frame validates", {
  v <- validate_data(base_df())
  expect_true(v$ok)
  expect_equal(v$measures, "Area")
})

test_that("a measurement column that parsed as text reports, never crashes (CM-5)", {
  b <- base_df(); b$Area <- as.character(b$Area); b$Area[3] <- "n/a"
  v <- validate_data(b)
  expect_false(v$ok)
  expect_match(v$msg, "No numeric measurement columns found")
  expect_match(v$msg, "look numeric but were read as text")
  expect_match(v$msg, "Area")
})

test_that("a text column alongside a good one is flagged but does not block loading", {
  b <- base_df(); b$Perimeter <- as.character(round(rnorm(12, 40), 2))
  v <- validate_data(b)
  expect_true(v$ok)
  expect_equal(v$measures, "Area")
  expect_match(v$parse_warning, "Perimeter")
})

test_that("a Time_min column with unit suffixes is caught, not left as text (CM-17)", {
  b <- base_df(); b$Time_min <- c("0 min", "30 min", "60 min")[(seq_len(12) - 1) %% 3 + 1]
  v <- validate_data(b)
  expect_false(v$ok)
  expect_true(v$needs_mapping)
  expect_match(v$msg, "not numeric")
  expect_match(v$msg, "30 min", fixed = TRUE)
})

test_that("a numeric-but-differently-named time column is still coerced", {
  b <- base_df(); names(b)[2] <- "Timepoint"
  v <- validate_data(b)
  expect_true(v$ok)
  expect_true(is.numeric(v$df$Time_min))
})

test_that("data with no time column at all loads with a synthetic timepoint", {
  b <- base_df()[, c("Condition", "Replicate", "Area")]
  v <- validate_data(b)
  expect_true(v$ok)
  expect_true(v$synthetic_time)
  expect_equal(unique(v$df$Time_min), 0)
})

test_that("an unrecognisable condition column asks for mapping", {
  b <- base_df(); names(b)[1] <- "genotype"
  v <- validate_data(b)
  expect_false(v$ok)
  expect_true(v$needs_mapping)
  expect_match(v$msg, "Condition/Group")
})

test_that("constant and ID-like columns are excluded from measurements", {
  b <- base_df(); b$Object_ID <- seq_len(12); b$Constant <- 7
  v <- validate_data(b)
  expect_true(v$ok)
  expect_equal(v$measures, "Area")
})

test_that("the shipped example data still validates", {
  csv <- file.path(cellmorphr_root(), "example_data.csv")
  skip_if_not(file.exists(csv))
  v <- validate_data(utils::read.csv(csv))
  expect_true(v$ok)
  expect_gt(length(v$measures), 0)
  expect_true(is.numeric(v$df$Time_min))
})
