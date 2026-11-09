# =============================================================================
# Sample colour map (3.6.0): hex normalisation, reading a map file, and the
# transfer rule (match by name, then remaining rows in order).
# =============================================================================

load_cellmorphr_core()

test_that("normalise_hex accepts hex, short hex and R colour names", {
  expect_equal(normalise_hex(c("#0072b2", "0072B2", "#abc", "red", "grey30", "", NA, "notacolour")),
               c("#0072B2", "#0072B2", "#AABBCC", "#FF0000", "#4D4D4D", NA, NA, NA))
})

test_that("parse_colour_map finds the sample and colour columns under several names", {
  a <- parse_colour_map(data.frame(Sample = c("WT", "Mut"), Hex = c("#0072B2", "#D55E00")))
  expect_equal(a$Sample, c("WT", "Mut")); expect_equal(a$Hex, c("#0072B2", "#D55E00"))
  b <- parse_colour_map(data.frame(Condition = c("A", "B", ""), Colour = c("blue", "bad colour", "#FFFFFF")))
  expect_equal(b$Sample, "A"); expect_equal(b$Hex, "#0000FF")      # bad colour and blank sample dropped
  c <- parse_colour_map(data.frame(x = c("P", "Q"), y = c("#111111", "#222222")))   # positional fallback
  expect_equal(c$Hex, c("#111111", "#222222"))
  expect_equal(nrow(parse_colour_map(NULL)), 0)
  expect_equal(nrow(parse_colour_map(data.frame(only = "one"))), 0)
})

test_that("assign_colour_map matches by name first, then fills the rest in order", {
  map <- data.frame(Sample = c("Uninfected", "Phage_MOI_1", "Extra1", "Extra2"),
                    Hex = c("#111111", "#222222", "#333333", "#444444"), stringsAsFactors = FALSE)
  cur <- c(WT = "#AAAAAA", Uninfected = "#BBBBBB", Phage_MOI_1 = "#CCCCCC", Novel = "#DDDDDD")
  r <- assign_colour_map(c("WT", "Uninfected", "Phage_MOI_1", "Novel"), map, current = cur)
  expect_equal(unname(r$colours), c("#333333", "#111111", "#222222", "#444444"))
  expect_equal(r$matched_by_name, 2); expect_equal(r$matched_by_order, 2)
  # a map from a differently named experiment transfers purely by order
  r2 <- assign_colour_map(c("A", "B"), map[3:4, ], current = c(A = "#AAAAAA", B = "#BBBBBB"))
  expect_equal(unname(r2$colours), c("#333333", "#444444")); expect_equal(r2$matched_by_name, 0)
  # fewer rows than samples: unmatched samples keep their current colour
  r3 <- assign_colour_map(c("A", "B", "C"), map[1, ], current = c(A = "#AAAAAA", B = "#BBBBBB", C = "#CCCCCC"))
  expect_equal(unname(r3$colours), c("#111111", "#BBBBBB", "#CCCCCC"))
  # no map: current colours, NA where none
  r4 <- assign_colour_map(c("A", "Z"), NULL, current = c(A = "#AAAAAA"))
  expect_equal(unname(r4), c("#AAAAAA", NA))
})
