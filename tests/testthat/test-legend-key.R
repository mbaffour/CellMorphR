# =============================================================================
# Legend key (3.5.0): colour / shape / line-type words, the text starter, and
# the legend-grob extraction, on plain ggplots. Pure functions above the core
# marker; the per-tab wiring is exercised headlessly through shiny::testServer.
# =============================================================================

load_cellmorphr_core()
library(ggplot2)

test_that("colour, shape and line-type words are human and cover the defaults", {
  expect_equal(colour_word(c("#0072B2", "#D55E00", "#009E73", "#E69F00")), c("blue", "vermillion", "bluish green", "orange"))
  expect_equal(colour_word("black"), "black")
  expect_equal(colour_word("#FFFFFF"), "white")
  expect_equal(colour_word("grey30"), "dark grey")
  expect_equal(shape_word(c(21, 22, 24, 16, 17)), c("filled circles", "filled squares", "filled triangles", "circles", "triangles"))
  expect_equal(shape_word(99), "shape 99")
  expect_equal(linetype_word(c("dashed", 1)), c("dashed lines", "solid lines"))
})

test_that("legend_text_starter writes one clause per legend from the plot's own scales", {
  d <- data.frame(x = rep(1:3, 2), y = 1:6, Condition = rep(c("Ctrl", "Trt"), each = 3), Replicate = factor(rep(1:3, 2)))
  p <- ggplot(d, aes(x, y, colour = Condition, shape = Replicate)) + geom_point() +
    scale_colour_manual(values = c(Ctrl = "#0072B2", Trt = "#D55E00")) +
    scale_shape_manual(values = c(21, 22, 24))
  txt <- legend_text_starter(p)
  expect_match(txt, "Condition: Ctrl \\(blue\\); Trt \\(vermillion\\)\\.", fixed = FALSE)
  expect_match(txt, "Replicate: 1 \\(filled circles\\); 2 \\(filled squares\\); 3 \\(filled triangles\\)\\.")
  expect_match(legend_text_starter(p, cvd = "deuteranopia"), "simulated for deuteranopia")
  # colour and fill of the same variable collapse to one clause
  p2 <- ggplot(d, aes(x, y, colour = Condition, fill = Condition)) + geom_point(shape = 21) +
    scale_colour_manual(values = c(Ctrl = "#0072B2", Trt = "#D55E00")) +
    scale_fill_manual(values = c(Ctrl = "#0072B2", Trt = "#D55E00"))
  expect_equal(length(gregexpr("Condition:", legend_text_starter(p2))[[1]]), 1)
  expect_equal(legend_text_starter(p2), "Condition: Ctrl (blue); Trt (vermillion).")
  # one variable encoded by colour AND shape: each entry gets both words, once
  p3 <- ggplot(d, aes(x, y, colour = Replicate, shape = Replicate)) + geom_point() +
    scale_colour_manual(values = c("#CC79A7", "#56B4E9", "#999999")) + scale_shape_manual(values = c(21, 22, 24))
  expect_equal(legend_text_starter(p3),
               "Replicate: 1 (reddish purple, filled circles); 2 (sky blue, filled squares); 3 (grey, filled triangles).")
  # continuous fill reads as a colour scale, with its limits
  h <- ggplot(data.frame(a = 1:3, b = 1:3, v = c(0, 5, 10)), aes(a, b, fill = v)) + geom_tile() +
    scale_fill_gradient(low = "white", high = "#D55E00")
  expect_match(legend_text_starter(h), "v: colour scale from white \\(0\\) to vermillion \\(10\\)\\.")
  # nothing mapped -> empty
  expect_equal(legend_text_starter(ggplot(d, aes(x, y)) + geom_point()), "")
})

test_that("legend_key_plot re-enables hidden discrete legends and leaves continuous ones alone", {
  d <- data.frame(x = 1:4, y = 1:4, g = c("a", "a", "b", "b"))
  hidden <- ggplot(d, aes(x, y, fill = g)) + geom_col() + scale_fill_manual(values = c(a = "red", b = "blue"), guide = "none") +
    theme(legend.position = "none")
  expect_null(extract_legend_grob(hidden))
  shown <- legend_key_plot(hidden)
  g <- extract_legend_grob(shown)
  expect_false(is.null(g))
  expect_setequal(intersect(legend_grob_labels(g), c("a", "b")), c("a", "b"))
  # the embedded version carries the same labels as the standalone key
  g_in <- extract_legend_grob(legend_key_plot(hidden, inside = TRUE))
  expect_equal(sort(legend_grob_labels(g_in)), sort(legend_grob_labels(g)))
  # continuous: still a colourbar, not binned
  h <- ggplot(data.frame(a = 1:3, b = 1:3, v = c(0, 5, 10)), aes(a, b, fill = v)) + geom_tile()
  gh <- extract_legend_grob(legend_key_plot(h))
  expect_false(is.null(gh))
  expect_true(any(grepl("colourbar|colorbar", vapply(gh$grobs, function(x) x$name, character(1)), ignore.case = TRUE)) ||
              length(legend_grob_labels(gh)) >= 3)
  sz <- legend_key_size(g)
  expect_true(sz$w > 0.5 && sz$w < 5 && sz$h > 0.3 && sz$h < 5)
})
