# =============================================================================
# CellMorphR -- Single-Cell Morphometry Analysis Suite
# Version 3.1.2
# For time-resolved microscopy experiments (e.g., bacteriophage infection)
#
# Author: Michael Baffour Awuah <mbaffour890@gmail.com>
#
# Citation:
#   Awuah, M.B. CellMorphR: Single-cell morphometry analysis with
#   anti-pseudoreplication statistics. Version 3.1.2.
#   URL: https://github.com/mbaffour/CellMorphR
# =============================================================================
# Easiest way to run on Windows: double-click run_cellmorphr.bat (in this
# folder) -- it finds R, installs any missing packages, and opens the app in
# your browser.
#
# Manual install (one time):
# install.packages(c("shiny", "bslib", "ggplot2", "dplyr", "tidyr", "readxl",
#                     "lme4", "lmerTest", "emmeans", "ggridges", "patchwork",
#                     "colourpicker", "DT", "writexl", "scales", "ggbeeswarm",
#                     "cluster", "hexbin", "svglite", "officer", "rvg"))

library(shiny)
library(bslib)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)
library(DT)
library(scales)

# Null-coalescing operator (built-in from R 4.4.0; define for older R)
if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (!is.null(a)) a else b

# Single source of truth for the app version (navbar, hero, footer, session files)
APP_VERSION <- "3.1.2"

# Attempt to load optional packages
has_lme4     <- requireNamespace("lme4", quietly = TRUE)
has_lmerTest <- requireNamespace("lmerTest", quietly = TRUE)
has_emmeans  <- requireNamespace("emmeans", quietly = TRUE)
has_ridges   <- requireNamespace("ggridges", quietly = TRUE)
has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
has_writexl  <- requireNamespace("writexl", quietly = TRUE)
has_beeswarm <- requireNamespace("ggbeeswarm", quietly = TRUE)
has_colourpicker <- requireNamespace("colourpicker", quietly = TRUE)
has_cluster <- requireNamespace("cluster", quietly = TRUE)
has_hexbin <- requireNamespace("hexbin", quietly = TRUE)
has_officer <- requireNamespace("officer", quietly = TRUE)
has_rvg <- requireNamespace("rvg", quietly = TRUE)

# =============================================================================
# PUBLICATION GGPLOT THEME (parameterized)
# =============================================================================
theme_publication <- function(base_size = 12, legend_pos = "top",
                              show_grid = FALSE, show_border = FALSE,
                              font_family = "sans",
                              title_size = NULL, title_face = "bold",
                              axis_title_size = NULL, axis_title_face = "bold",
                              axis_text_size = NULL, axis_text_face = "plain",
                              legend_title_size = NULL, legend_title_face = "bold",
                              legend_text_size = NULL, legend_text_face = "plain",
                              strip_size = NULL, strip_face = "bold") {
  # Resolve sizes: use explicit values when provided, otherwise scale from base
  ts  <- title_size        %||% (base_size * 1.15)
  ats <- axis_title_size   %||% (base_size * 1.0)
  axs <- axis_text_size    %||% (base_size * 0.9)
  lts <- legend_title_size %||% (base_size * 0.9)
  lxs <- legend_text_size  %||% (base_size * 0.85)
  sts <- strip_size        %||% (base_size * 0.9)

  t <- theme_classic(base_size = base_size) %+replace%
    theme(
      text             = element_text(family = font_family, color = "black",
                                      size = base_size),
      plot.title       = element_text(family = font_family, face = title_face,
                                      size = ts, hjust = 0,
                                      margin = margin(b = 8), color = "black"),
      plot.subtitle    = element_text(family = font_family, size = base_size * 0.85,
                                      hjust = 0, color = "grey40",
                                      face = "plain", margin = margin(b = 10)),
      plot.caption     = element_text(family = font_family, size = base_size * 0.7,
                                      hjust = 1, color = "grey50",
                                      face = "plain", margin = margin(t = 8)),
      axis.title       = element_text(family = font_family, face = axis_title_face,
                                      size = ats, color = "black"),
      axis.text        = element_text(family = font_family, face = axis_text_face,
                                      size = axs, color = "black"),
      axis.line        = element_line(color = "black", linewidth = 0.5),
      axis.ticks       = element_line(color = "black", linewidth = 0.4),
      legend.title     = element_text(family = font_family, face = legend_title_face,
                                      size = lts, color = "black"),
      legend.text      = element_text(family = font_family, face = legend_text_face,
                                      size = lxs, color = "black"),
      legend.background = element_blank(),
      legend.key       = element_blank(),
      legend.key.spacing.y = unit(3, "pt"),
      # Stack multiple guides (e.g. Condition + Replicate, or colour + size)
      # instead of letting them compete for one row and run off the figure.
      legend.box       = "vertical",
      legend.position  = legend_pos,
      panel.background = element_blank(),
      panel.grid       = element_blank(),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text       = element_text(family = font_family, face = strip_face,
                                      size = sts, color = "black"),
      plot.margin      = margin(12, 12, 12, 12)
    )
  if (show_grid) {
    t <- t + theme(
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank()
    )
  }
  if (show_border) {
    t <- t + theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.line = element_blank()
    )
  }
  t
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# ---- Column mapping: guess which column serves each role ----
guess_col <- function(df, patterns, prefer_exact = NULL) {
  nms <- names(df)
  if (!is.null(prefer_exact) && prefer_exact %in% nms) return(prefer_exact)
  for (pat in patterns) {
    hit <- grep(pat, nms, ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  NULL
}

validate_data <- function(df, col_map = NULL) {
  # col_map: named list with keys condition, time, replicate -> actual column names
  # If NULL, auto-detect
  nms <- names(df)

  # Auto-detect if not provided
  if (is.null(col_map)) {
    cond_col <- guess_col(df, c("^condition$","^group$","^treatment$","^sample$","^cond$"))
    time_col <- guess_col(df, c("^time","^timepoint","^t_min","^min$","^hour","^day","^frame"))
    rep_col  <- guess_col(df, c("^replicate","^rep$","^rep_","biological","biol","^bio_rep","^experiment"))
  } else {
    cond_col <- col_map$condition
    time_col <- col_map$time
    rep_col  <- col_map$replicate
  }

  # Time is optional: many experiments are a single snapshot of several samples
  # with no time course at all. When absent we synthesise a single timepoint so
  # the rest of the app keeps working, and mark it so nothing is labelled "min".
  synthetic_time <- is.null(time_col)

  missing_roles <- c(
    if (is.null(cond_col)) "Condition/Group",
    if (is.null(rep_col))  "Replicate/Experiment"
  )

  if (length(missing_roles) > 0) {
    return(list(
      ok = FALSE, needs_mapping = TRUE,
      msg = paste("Could not auto-detect columns for:", paste(missing_roles, collapse=", "),
                  "\nPlease use the column mapping panel to assign each role."),
      cond_col = cond_col, time_col = time_col, rep_col = rep_col,
      all_cols = nms
    ))
  }

  # Rename to standard names internally
  df2 <- df
  if (cond_col != "Condition") { df2$Condition <- df2[[cond_col]] }
  if (synthetic_time) {
    df2$Time_min <- 0
  } else {
    # Always coerce, including when the column is already called Time_min. Files
    # with values like "30 min", "T30" or "0.5 h" used to stay character, which
    # made every downstream sort lexical: 0, 120, 30, 60 (audit CM-17).
    raw_time <- df2[[time_col]]
    coerced <- suppressWarnings(as.numeric(as.character(raw_time)))
    if (all(is.na(coerced)) && any(!is.na(raw_time) & nzchar(trimws(as.character(raw_time))))) {
      return(list(
        ok = FALSE, needs_mapping = TRUE,
        msg = paste0("The time column \"", time_col, "\" is not numeric -- values such as \"",
                     paste(utils::head(unique(stats::na.omit(as.character(raw_time))), 3), collapse = "\", \""),
                     "\" could not be read as minutes.\nStrip the units (use 30, not \"30 min\"), ",
                     "or map a different column to the time role below."),
        cond_col = cond_col, time_col = time_col, rep_col = rep_col, all_cols = nms))
    }
    df2$Time_min <- coerced
  }
  if (rep_col  != "Replicate") { df2$Replicate <- df2[[rep_col]] }

  # Exclude meta + ID-like from measurements
  meta_cols <- unique(c("Condition", "Time_min", "Replicate",
                         cond_col, time_col, rep_col,
                         "Image", "Cell_ID", "ID"))
  measure_cols <- setdiff(nms, meta_cols)
  # vapply, not sapply: sapply(character(0), ...) returns list(), and indexing a
  # character vector by a list is a hard error -- so the friendly "no numeric
  # measurement columns" message below used to be unreachable in exactly the
  # case it was written for (audit CM-5).
  numeric_flag <- vapply(measure_cols, function(m) is.numeric(df2[[m]]), logical(1))
  looks_numeric_but_is_not <- measure_cols[!numeric_flag & vapply(measure_cols, function(m) {
    v <- as.character(df2[[m]])
    v <- v[!is.na(v) & nzchar(trimws(v))]
    length(v) > 0 && mean(!is.na(suppressWarnings(as.numeric(v)))) > 0.5
  }, logical(1))]
  measure_cols <- measure_cols[numeric_flag]
  id_pattern <- grepl("^id$|_id$|^id_|cell.?id|object.?id|label$|^index$", measure_cols, ignore.case = TRUE)
  measure_cols <- measure_cols[!id_pattern]
  # Also exclude constant columns
  measure_cols <- measure_cols[vapply(measure_cols, function(m) {
    v <- df2[[m]]; length(unique(na.omit(v))) > 1
  }, logical(1))]

  parse_warning <- if (length(looks_numeric_but_is_not) > 0)
    paste0("These columns look numeric but were read as text, so they are not available as measurements: ",
           paste(looks_numeric_but_is_not, collapse = ", "),
           ". Values such as \"n/a\", \"NA (auto)\", \"#DIV/0!\" or a unit suffix make a whole column text.")
  else NULL

  if (length(measure_cols) == 0) {
    return(list(ok = FALSE, needs_mapping = FALSE,
                msg = paste(c("No numeric measurement columns found. All numeric columns are constant or excluded.",
                              parse_warning), collapse = "\n")))
  }

  list(ok = TRUE, needs_mapping = FALSE,
       msg = paste(c(if (synthetic_time)
                       "Data loaded successfully. No time column found - samples will be compared directly."
                     else "Data loaded successfully.",
                     parse_warning), collapse = "\n"),
       parse_warning = parse_warning,
       measures = measure_cols, df = df2, synthetic_time = synthetic_time,
       cond_col = cond_col, time_col = time_col, rep_col = rep_col)
}

compute_replicate_summaries <- function(df, measure) {
  df %>%
    group_by(Condition, Time_min, Replicate) %>%
    summarise(
      median_val = median(.data[[measure]], na.rm = TRUE),
      mean_val   = mean(.data[[measure]], na.rm = TRUE),
      sd_val     = sd(.data[[measure]], na.rm = TRUE),
      n_cells    = n(),
      q25        = quantile(.data[[measure]], 0.25, na.rm = TRUE),
      q75        = quantile(.data[[measure]], 0.75, na.rm = TRUE),
      .groups    = "drop"
    )
}

compute_condition_summaries <- function(rep_summary) {
  rep_summary %>%
    group_by(Condition, Time_min) %>%
    summarise(
      grand_median = mean(median_val),
      grand_mean   = mean(mean_val),
      sem_median   = sd(median_val) / sqrt(n()),
      sem_mean     = sd(mean_val) / sqrt(n()),
      n_reps       = n(),
      .groups      = "drop"
    )
}


# =============================================================================
# FLEXIBLE GROUPING + STATISTICS HELPERS (v2.6-flexstats)
# =============================================================================

# Detect likely experimental metadata columns while preserving the original
# Condition/Time_min/Replicate workflow.
detect_metadata_columns <- function(df, measure_cols = NULL) {
  nms <- names(df)
  if (is.null(measure_cols)) measure_cols <- character(0)
  id_like <- grepl("^id$|_id$|cell.?id|object.?id|label$|^index$", nms, ignore.case = TRUE)
  meta <- setdiff(nms, measure_cols)
  meta <- meta[!id_like[match(meta, nms)]]
  priority <- c("Condition", "Time_min", "Replicate", "Plasmid", "Vector", "Vector_Type",
                "Gene", "Infection", "Infection_Status", "Inducer", "Inducer_Concentration",
                "Treatment", "Strain")
  c(intersect(priority, meta), setdiff(meta, priority))
}

clean_group_piece <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "NA"
  x
}

build_group_label <- function(df, selected_cols, separator = " | ") {
  if (is.null(selected_cols) || length(selected_cols) == 0) return(as.character(df$Condition))
  selected_cols <- selected_cols[selected_cols %in% names(df)]
  if (length(selected_cols) == 0) return(as.character(df$Condition))
  pieces <- lapply(selected_cols, function(cl) {
    val <- clean_group_piece(df[[cl]])
    if (cl %in% c("Time_min", "Time", "Timepoint")) val <- paste0(val, " min")
    paste0(cl, "=", val)
  })
  do.call(paste, c(pieces, sep = separator))
}

# Wrap a single label string to a maximum line width so long titles, subtitles,
# and captions fold into multiple lines instead of running over the figure.
wrap_plot_text <- function(x, width = 100) {
  if (is.null(x) || !is.character(x) || length(x) != 1 || is.na(x) || !nzchar(x)) return(x)
  paste(strwrap(x, width = width), collapse = "\n")
}

# Fold long legend entries (e.g. long condition names) onto multiple lines so
# wide legends stop colliding with each other or the panel. Only touches
# discrete scales whose labels are still automatic, so explicit labels
# (wrapped x axes, named quartile bins, ...) are left alone.
wrap_legend_labels <- function(p, width = 24) {
  try({
    for (i in seq_along(p$scales$scales)) {
      sc <- p$scales$scales[[i]]
      if (inherits(sc, "ScaleDiscrete") && inherits(sc$labels, "waiver")) {
        p$scales$scales[[i]]$labels <- function(x) wrap_axis_labels(x, width)
      }
    }
  }, silent = TRUE)
  p
}

wrap_axis_labels <- function(labels, width = 18) {
  width <- suppressWarnings(as.integer(width))
  if (is.na(width) || width <= 0) width <- 18
  vapply(as.character(labels), function(z) paste(strwrap(z, width = width), collapse = "\n"), character(1))
}

make_significance_label <- function(p) {
  ifelse(is.na(p), "NA",
    ifelse(p < 0.001, "***",
      ifelse(p < 0.01, "**",
        ifelse(p < 0.05, "*", "ns"))))
}

format_p_value <- function(p) {
  ifelse(is.na(p), NA_character_, ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 4)))
}

# Point shapes offered wherever the user can choose one. Shapes 21-25 take a
# fill colour (outline + fill); the rest are drawn in a single colour.
SHAPE_CHOICES <- c(
  "Circle (filled)"        = 21,
  "Square (filled)"        = 22,
  "Triangle (filled)"      = 24,
  "Diamond (filled)"       = 23,
  "Triangle down (filled)" = 25,
  "Circle (solid)"         = 16,
  "Square (solid)"         = 15,
  "Triangle (solid)"       = 17,
  "Diamond (solid)"        = 18,
  "Circle (open)"          = 1,
  "Square (open)"          = 0,
  "Triangle (open)"        = 2,
  "Diamond (open)"         = 5,
  "Plus"                   = 3,
  "Cross"                  = 4,
  "Asterisk"               = 8
)
DEFAULT_REP_SHAPES <- c(21, 22, 24, 23, 25)

# Split Condition labels like "pSJExD_gp75_100nM" into a numeric dose column
# plus a base-label column (and generic part columns), so gradient designs can
# be analyzed as vector x dose instead of N unrelated groups.
split_condition_parts <- function(conds, delim = "_") {
  conds <- trimws(as.character(conds))
  # Dose = a trailing number with an optional unit suffix, matched at the END of
  # the label. Anchoring matters: naively stripping non-digits would read "75"
  # out of a gene name like "gp75" instead of the real dose.
  dose_re <- "([0-9]+\\.?[0-9]*)[[:space:]]*([A-Za-z]*)[[:space:]]*$"
  m <- regmatches(conds, regexec(dose_re, conds))
  num_val <- vapply(m, function(g) if (length(g) >= 2) suppressWarnings(as.numeric(g[2])) else NA_real_, numeric(1))
  unit <- vapply(m, function(g) if (length(g) >= 3) g[3] else "", character(1))
  unit_tab <- table(unit[nzchar(unit)])
  unit_name <- if (length(unit_tab) > 0) names(sort(unit_tab, decreasing = TRUE))[1] else ""
  # Base label = the condition with that trailing dose token removed
  base <- trimws(sub(dose_re, "", conds))
  base <- trimws(gsub(paste0("[", gsub("([][^\\-])", "\\\\\\1", delim), "[:space:]]+$"), "", base))
  base[!nzchar(base)] <- conds[!nzchar(base)]
  # Generic per-delimiter parts (useful when labels encode >2 factors)
  toks <- strsplit(conds, delim, fixed = TRUE)
  max_parts <- max(c(1L, lengths(toks)))
  parts <- lapply(seq_len(max_parts), function(i)
    vapply(toks, function(tt) if (length(tt) >= i) tt[i] else NA_character_, character(1)))
  names(parts) <- paste0("Cond_Part", seq_len(max_parts))
  list(num_val = num_val, unit_name = unit_name, base = base, parts = parts,
       dose_col = if (nzchar(unit_name)) paste0("Dose_", unit_name) else "Dose_Value")
}

# Minimum attainable two-sided p for a Wilcoxon rank-sum test with group sizes
# n1, n2 (exact test, no ties). With n=2 vs 2 this is 0.333 -- i.e., the test
# can NEVER reach p < 0.05. Used to warn users before they trust "ns" results.
min_wilcox_p <- function(n1, n2) {
  if (is.na(n1) || is.na(n2) || n1 < 1 || n2 < 1) return(NA_real_)
  min(2 / choose(n1 + n2, n1), 1)
}

# Shared significance computation for on-figure brackets (Distributions and
# Cell Metrics). Tests per-replicate summaries -- never raw cells -- using the
# test and correction chosen on the Statistics tab, so figures and tables agree.
# rs: replicate summary df with Condition, Time_min and the yvar column.
# pairs_list: list of c(cond1, cond2). Returns annotation df + caption text.
compute_sig_table <- function(rs, pairs_list, times, sig_yvar = "median_val",
                              test = "welch", correction = "holm",
                              label_style = "stars", show_ns = TRUE,
                              y_lookup = NULL, max_brackets = 12, max_per_time = 4) {
  sig_annot <- data.frame(Time_min = numeric(), y_pos = numeric(), p_raw = numeric(),
                          Cond1 = character(), Cond2 = character(), stringsAsFactors = FALSE)
  underpowered <- FALSE
  wilcox_impossible <- FALSE
  truncated_note <- NULL
  for (tm in times) {
    for (pair in pairs_list) {
      vals1 <- rs[[sig_yvar]][rs$Condition == pair[1] & rs$Time_min == tm]
      vals2 <- rs[[sig_yvar]][rs$Condition == pair[2] & rs$Time_min == tm]
      n1 <- length(vals1); n2 <- length(vals2)
      if (n1 >= 2 && n2 >= 2) {
        if (n1 < 3 || n2 < 3) underpowered <- TRUE
        if (test == "wilcox") {
          mp <- min_wilcox_p(n1, n2)
          if (!is.na(mp) && mp > 0.05) wilcox_impossible <- TRUE
        }
        tt <- tryCatch({
          if (test == "wilcox") wilcox.test(vals2, vals1, exact = FALSE)
          else t.test(vals2, vals1, var.equal = FALSE)
        }, error = function(e) NULL)
        if (!is.null(tt)) {
          y_top <- if (is.function(y_lookup)) y_lookup(tm, pair, vals1, vals2) else max(c(vals1, vals2), na.rm = TRUE)
          sig_annot <- rbind(sig_annot, data.frame(Time_min = tm, y_pos = y_top, p_raw = tt$p.value,
            Cond1 = pair[1], Cond2 = pair[2], stringsAsFactors = FALSE))
        }
      }
    }
  }
  n_family <- nrow(sig_annot)
  if (nrow(sig_annot) > 0) {
    corr <- if (is.null(correction) || correction == "none") "none" else correction
    # The family is every condition-within-timepoint comparison drawn on this
    # figure, adjusted together and BEFORE the display truncation below.
    sig_annot$p_adj <- if (corr == "none") sig_annot$p_raw else p.adjust(sig_annot$p_raw, method = corr)
    sig_annot$stars <- ifelse(sig_annot$p_adj < 0.001, "***",
                         ifelse(sig_annot$p_adj < 0.01, "**",
                                ifelse(sig_annot$p_adj < 0.05, "*", "ns")))
    sig_annot$label <- switch(label_style,
      "stars"   = sig_annot$stars,
      "stars_p" = ifelse(sig_annot$p_adj < 0.001, paste0(sig_annot$stars, "\n(p<0.001)"),
                         paste0(sig_annot$stars, "\n(p=", formatC(sig_annot$p_adj, format = "f", digits = 3), ")")),
      "pval"    = ifelse(sig_annot$p_adj < 0.001, "p<0.001",
                         paste0("p=", formatC(sig_annot$p_adj, format = "f", digits = 3))),
      sig_annot$stars)
    if (!isTRUE(show_ns)) sig_annot <- sig_annot[sig_annot$stars != "ns", , drop = FALSE]
    # Too many brackets bury the figure: keep only the most significant ones and
    # say so in the caption (full results live in the Statistics tab). Cap per
    # timepoint as well as overall, so one timepoint cannot become a ladder of
    # brackets that squashes the data into the bottom of the panel.
    n_total <- nrow(sig_annot)
    if (is.finite(max_per_time)) {
      keep_idx <- unlist(lapply(split(seq_len(nrow(sig_annot)), sig_annot$Time_min),
        function(ix) ix[order(sig_annot$p_adj[ix])][seq_len(min(length(ix), max_per_time))]))
      sig_annot <- sig_annot[sort(keep_idx), , drop = FALSE]
    }
    if (is.finite(max_brackets) && nrow(sig_annot) > max_brackets) {
      keep <- order(sig_annot$p_adj)[seq_len(max_brackets)]
      sig_annot <- sig_annot[sort(keep), , drop = FALSE]
    }
    if (nrow(sig_annot) < n_total) {
      truncated_note <- paste0("Showing the ", nrow(sig_annot), " most significant of ",
                               n_total, " comparisons; full results in the Statistics tab.")
    }
  }
  test_name <- if (test == "wilcox") "Wilcoxon rank-sum" else "Welch's t-test"
  corr_name <- switch(if (is.null(correction)) "holm" else correction,
    "holm" = "Holm", "BH" = "Benjamini-Hochberg", "bonferroni" = "Bonferroni", "none" = "uncorrected", correction)
  unit_label <- if (sig_yvar == "mean_val") "replicate means" else "replicate medians"
  caption <- paste0(test_name, " on ", unit_label,
    if (corr_name == "uncorrected") ", uncorrected"
    else paste0(", ", corr_name, "-corrected across all ", n_family,
                " condition-within-timepoint comparison", if (n_family == 1) "" else "s",
                " (one family)"),
    ". * p<0.05, ** p<0.01, *** p<0.001.")
  warn <- character(0)
  if (wilcox_impossible) warn <- c(warn,
    "Wilcoxon cannot reach p<0.05 with these replicate counts -- use Welch's t-test, the mixed model, or a trend test.")
  else if (underpowered) warn <- c(warn,
    "Fewer than 3 replicates in at least one group: tests are underpowered; consider the mixed model or a trend test.")
  if (!is.null(truncated_note)) warn <- c(warn, truncated_note)
  if (length(warn) > 0) caption <- paste0(caption, " Note: ", paste(warn, collapse = " "))
  list(annot = sig_annot, caption = caption, warnings = warn)
}

compute_flexible_replicate_summary <- function(df, measure, group_cols, yvar = "median_val") {
  stopifnot(measure %in% names(df))
  if (is.null(group_cols) || length(group_cols) == 0) group_cols <- "Condition"
  group_cols <- unique(group_cols[group_cols %in% names(df)])
  if (!"Replicate" %in% names(df)) stop("A biological replicate column named Replicate is required.")
  df$Statistics_Group <- build_group_label(df, group_cols)
  grp <- unique(c(group_cols, "Statistics_Group", "Replicate"))
  out <- df %>%
    group_by(across(all_of(grp))) %>%
    summarise(
      median_val = median(.data[[measure]], na.rm = TRUE),
      mean_val   = mean(.data[[measure]], na.rm = TRUE),
      sd_val     = sd(.data[[measure]], na.rm = TRUE),
      iqr_val    = IQR(.data[[measure]], na.rm = TRUE),
      n_cells    = sum(!is.na(.data[[measure]])),
      .groups = "drop"
    )
  out$Stat_Value <- out[[yvar]]
  out
}

# The stratified modes need their stratum column to be one of the grouping
# columns. Up to v3.1.2 a request for "within each time point" with Time_min
# unchecked fell through to vs_reference against an arbitrary first group, while
# the report still printed the mode the user asked for (audit CM-9). The mode
# actually built is now returned on the table as the "realised_mode" attribute,
# together with a "mode_message" explaining any downgrade, and the callers
# report that rather than the request.
build_comparison_table <- function(rep_df, mode = "vs_reference", reference_group = NULL,
                                   strict = FALSE) {
  requested_mode <- mode
  mode_message <- NULL
  needs <- c(within_time = "Time_min", within_condition = "Condition")[mode]
  if (!is.na(needs) && !needs %in% names(rep_df)) {
    msg <- paste0("Comparison mode \"", mode, "\" needs ", needs,
                  " to be one of the grouping columns, but the statistical groups are built from ",
                  paste(setdiff(names(rep_df), c("Statistics_Group", "Replicate", "Stat_Value",
                                                 "median_val", "mean_val", "sd_val", "iqr_val",
                                                 "n_cells")), collapse = " + "), ".")
    if (strict) stop(msg, call. = FALSE)
    mode <- "vs_reference"
    mode_message <- paste(msg, "Fell back to comparing every group against the reference.")
  }
  groups <- sort(unique(as.character(rep_df$Statistics_Group)))
  if (length(groups) < 2) return(data.frame())

  make_pairs <- function(gs) {
    gs <- sort(unique(as.character(gs)))
    if (length(gs) < 2) return(list())
    combn(gs, 2, simplify = FALSE)
  }

  pair_df <- data.frame(Group_1 = character(), Group_2 = character(), Stratum = character(), stringsAsFactors = FALSE)
  if (mode == "all_pairwise") {
    pairs <- make_pairs(groups)
    pair_df <- do.call(rbind, lapply(pairs, function(z) data.frame(Group_1 = z[1], Group_2 = z[2], Stratum = "All groups")))
  } else if (mode == "within_time" && "Time_min" %in% names(rep_df)) {
    for (tm in sort(unique(rep_df$Time_min))) {
      gs <- unique(rep_df$Statistics_Group[rep_df$Time_min == tm])
      pairs <- make_pairs(gs)
      if (length(pairs) > 0) pair_df <- rbind(pair_df, do.call(rbind, lapply(pairs, function(z)
        data.frame(Group_1 = z[1], Group_2 = z[2], Stratum = paste0("Time_min=", tm), stringsAsFactors = FALSE))))
    }
  } else if (mode == "within_condition" && "Condition" %in% names(rep_df)) {
    for (cd in sort(unique(as.character(rep_df$Condition)))) {
      gs <- unique(rep_df$Statistics_Group[as.character(rep_df$Condition) == cd])
      pairs <- make_pairs(gs)
      if (length(pairs) > 0) pair_df <- rbind(pair_df, do.call(rbind, lapply(pairs, function(z)
        data.frame(Group_1 = z[1], Group_2 = z[2], Stratum = paste0("Condition=", cd), stringsAsFactors = FALSE))))
    }
  } else {
    if (is.null(reference_group) || !reference_group %in% groups) reference_group <- groups[1]
    pairs <- lapply(setdiff(groups, reference_group), function(g) c(reference_group, g))
    pair_df <- do.call(rbind, lapply(pairs, function(z) data.frame(Group_1 = z[1], Group_2 = z[2], Stratum = "Selected reference", stringsAsFactors = FALSE)))
  }
  if (is.null(pair_df) || nrow(pair_df) == 0) return(data.frame())
  pair_df <- unique(pair_df)
  pair_df$Comparison_ID <- seq_len(nrow(pair_df))
  out <- pair_df[, c("Comparison_ID", "Stratum", "Group_1", "Group_2")]
  attr(out, "requested_mode") <- requested_mode
  attr(out, "realised_mode")  <- mode
  attr(out, "mode_message")   <- mode_message
  out
}

#' The comparison mode actually built, for reports. Never echoes the request.
realised_comparison_mode <- function(comparison_tbl, fallback = NA_character_) {
  attr(comparison_tbl, "realised_mode") %||% fallback
}

run_pairwise_tests <- function(rep_df, comparison_tbl, yvar = "Stat_Value",
                               test = "welch", correction = "holm", reference_group = NULL,
                               measurement = NA_character_, comparison_mode = NA_character_, stat_unit = NA_character_) {
  if (is.null(comparison_tbl) || nrow(comparison_tbl) == 0) return(data.frame())
  res <- lapply(seq_len(nrow(comparison_tbl)), function(i) {
    g1 <- comparison_tbl$Group_1[i]; g2 <- comparison_tbl$Group_2[i]
    d1 <- rep_df[rep_df$Statistics_Group == g1, , drop = FALSE]
    d2 <- rep_df[rep_df$Statistics_Group == g2, , drop = FALSE]
    v1 <- d1[[yvar]]; v2 <- d2[[yvar]]
    warn <- character()
    n1 <- length(unique(d1$Replicate)); n2 <- length(unique(d2$Replicate))
    if (n1 < 2 || n2 < 2) warn <- c(warn, "Too few biological replicates for reliable testing")
    if (max(n1, n2) >= 2 * max(1, min(n1, n2))) warn <- c(warn, "Unbalanced replicate counts")
    if (test == "wilcox") {
      mp <- min_wilcox_p(length(na.omit(v1)), length(na.omit(v2)))
      if (!is.na(mp) && mp > 0.05) warn <- c(warn,
        paste0("Wilcoxon cannot reach p<0.05 with these group sizes (minimum attainable p = ",
               formatC(mp, format = "f", digits = 3), ") -- use Welch's t-test, the mixed model, or a trend test"))
    }
    p_raw <- NA_real_; stat <- NA_real_; method_used <- test
    ci_low <- NA_real_; ci_high <- NA_real_; df_val <- NA_real_
    if (length(na.omit(v1)) >= 2 && length(na.omit(v2)) >= 2) {
      tt <- tryCatch({
        if (test == "wilcox") wilcox.test(v2, v1, exact = FALSE)
        else t.test(v2, v1, var.equal = FALSE)
      }, error = function(e) { warn <<- c(warn, paste("Test failed:", conditionMessage(e))); NULL })
      if (!is.null(tt)) {
        p_raw <- tt$p.value
        stat <- unname(tt$statistic[1])
        method_used <- tt$method
        if (!is.null(tt$conf.int)) { ci_low <- tt$conf.int[1]; ci_high <- tt$conf.int[2] }
        if (!is.null(tt$parameter)) df_val <- unname(tt$parameter[1])
      }
    }
    data.frame(
      Comparison_ID = comparison_tbl$Comparison_ID[i],
      Stratum = comparison_tbl$Stratum[i],
      Measurement = measurement,
      Comparison_Mode = comparison_mode,
      Statistical_Unit = stat_unit,
      Reference_Group = ifelse(is.null(reference_group), NA, reference_group),
      Group_1 = g1,
      Group_2 = g2,
      Test = method_used,
      Statistic = stat,
      df = round(df_val, 2),
      Mean_Group_1 = mean(v1, na.rm = TRUE),
      Mean_Group_2 = mean(v2, na.rm = TRUE),
      Difference_Group2_minus_Group1 = mean(v2, na.rm = TRUE) - mean(v1, na.rm = TRUE),
      CI95_low = ci_low,
      CI95_high = ci_high,
      N_BioRep_Group_1 = n1,
      N_BioRep_Group_2 = n2,
      N_Cells_Group_1 = sum(d1$n_cells, na.rm = TRUE),
      N_Cells_Group_2 = sum(d2$n_cells, na.rm = TRUE),
      p_raw = p_raw,
      Warning_Notes = paste(unique(warn), collapse = "; "),
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()
  if (nrow(res) > 0) {
    method <- if (is.null(correction) || correction == "none") "none" else correction
    res$p_adj <- if (method == "none") res$p_raw else p.adjust(res$p_raw, method = method)
    res$Significance <- make_significance_label(res$p_adj)
    res$p_raw_fmt <- format_p_value(res$p_raw)
    res$p_adj_fmt <- format_p_value(res$p_adj)
  }
  res
}

# =============================================================================
# STANDARDISED EFFECT SIZE
# =============================================================================
# Up to v3.1.2 the Effect Sizes tab reported a raw Cohen's d computed with
# sqrt((var1 + var2)/2) -- the unweighted average, which is only the pooled SD
# when the two groups have equal n -- with no small-sample correction and no
# confidence interval, while the figure drew dotted "large effect" guides at
# +/-0.8 (audit CM-3). At n = 3 vs 3 that d is 25% larger than the unbiased
# estimate. This computes Hedges' g instead: correctly pooled SD, Hedges'
# small-sample correction J, and a noncentral-t confidence interval, which at
# n = 3 is wide enough to make the point that three dishes cannot establish
# "large effect".

hedges_g <- function(ctrl, trt, conf_level = 0.95) {
  ctrl <- ctrl[is.finite(ctrl)]; trt <- trt[is.finite(trt)]
  n1 <- length(ctrl); n2 <- length(trt)
  blank <- function(note) list(g = NA_real_, d = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
                               pooled_sd = NA_real_, J = NA_real_, n_ctrl = n1, n_trt = n2,
                               conf_level = conf_level, note = note)
  if (n1 < 2 || n2 < 2) {
    return(blank(paste0("Not estimable: needs at least 2 biological replicates per group (got ",
                        n1, " and ", n2, ").")))
  }
  df <- n1 + n2 - 2
  # (n-1)-weighted pooled SD: valid for unequal n, unlike the average of variances.
  s_pool <- sqrt(((n1 - 1) * var(ctrl) + (n2 - 1) * var(trt)) / df)
  if (!is.finite(s_pool) || s_pool <= 0) {
    return(blank("Not estimable: the pooled standard deviation is zero (all replicates identical)."))
  }
  d <- (mean(trt) - mean(ctrl)) / s_pool
  J <- 1 - 3 / (4 * df - 1)            # Hedges' small-sample correction
  g <- J * d

  # Noncentral-t interval: invert the distribution of the observed t statistic.
  nfac <- sqrt((n1 * n2) / (n1 + n2))
  t_obs <- d * nfac
  alpha <- (1 - conf_level) / 2
  solve_ncp <- function(target) {
    # pt() warns about precision at the extremes of the bracket; that is noise
    # from the search, not from the answer uniroot converges on.
    f <- function(ncp) suppressWarnings(pt(t_obs, df, ncp)) - target
    lo <- t_obs - 10 * max(1, abs(t_obs)) - 10
    hi <- t_obs + 10 * max(1, abs(t_obs)) + 10
    tryCatch(uniroot(f, c(lo, hi), tol = 1e-6)$root, error = function(e) NA_real_)
  }
  ncp_hi <- solve_ncp(alpha)          # upper limit
  ncp_lo <- solve_ncp(1 - alpha)      # lower limit
  list(g = g, d = d,
       ci_low  = J * ncp_lo / nfac,
       ci_high = J * ncp_hi / nfac,
       pooled_sd = s_pool, J = J, n_ctrl = n1, n_trt = n2, conf_level = conf_level,
       note = if (min(n1, n2) < 4)
         paste0("Based on ", n1, " and ", n2,
                " biological replicates: the interval is very wide and a point estimate alone is not interpretable.")
       else NA_character_)
}

# =============================================================================
# ASSUMPTION DIAGNOSTICS (advisory only -- they never change the test that runs)
# =============================================================================
# The README has always documented a Levene test; until v3.2.0 none existed
# (audit CM-12), and Shapiro-Wilk was printed without any statement of what it
# did or did not do (audit CM-18). Both are computed here, on the same unit the
# tests use, and both are labelled advisory: nothing in the app switches test on
# the strength of them. Deciding to transform, or to move to the rank-based
# test, stays with the user.

#' Levene's test with median centering (Brown & Forsythe 1974), base R only.
levene_test <- function(y, group) {
  ok <- !is.na(y) & !is.na(group)
  y <- as.numeric(y[ok]); group <- droplevels(factor(group[ok]))
  if (nlevels(group) < 2) return(NULL)
  counts <- table(group)
  if (any(counts < 2)) return(NULL)
  centres <- tapply(y, group, median, na.rm = TRUE)
  z <- abs(y - centres[as.character(group)])
  if (all(abs(z - z[1]) < .Machine$double.eps^0.5)) return(NULL)  # no spread at all
  fit <- tryCatch(anova(lm(z ~ group)), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  list(statistic = fit[["F value"]][1], df1 = fit[["Df"]][1], df2 = fit[["Df"]][2],
       p.value = fit[["Pr(>F)"]][1], n_groups = nlevels(group),
       method = "Levene's test with median centering (Brown & Forsythe, 1974)")
}

#' Largest / smallest group SD. A ratio above ~3 is the usual rule of thumb for
#' variance heterogeneity that Welch's test handles and a pooled test does not.
variance_ratio <- function(y, group) {
  ok <- !is.na(y) & !is.na(group)
  y <- as.numeric(y[ok]); group <- droplevels(factor(group[ok]))
  if (nlevels(group) < 2) return(NULL)
  sds <- tapply(y, group, sd)
  sds <- sds[is.finite(sds) & sds > 0]
  if (length(sds) < 2) return(NULL)
  list(ratio = max(sds) / min(sds), min_sd = min(sds), max_sd = max(sds),
       min_group = names(sds)[which.min(sds)], max_group = names(sds)[which.max(sds)])
}

#' Run the advisory diagnostics on a set of replicate-level summaries.
assumption_diagnostics <- function(rep_df, yvar = "Stat_Value", group_col = "Statistics_Group",
                                   residuals = NULL) {
  if (is.null(rep_df) || nrow(rep_df) == 0 || !yvar %in% names(rep_df)) return(NULL)
  y <- rep_df[[yvar]]; g <- rep_df[[group_col]]
  sw_input <- if (!is.null(residuals)) residuals else
    unlist(lapply(split(y, g), function(v) v - mean(v, na.rm = TRUE)), use.names = FALSE)
  sw_input <- sw_input[is.finite(sw_input)]
  shapiro <- if (length(sw_input) >= 3 && length(unique(sw_input)) > 1)
    tryCatch(shapiro.test(sw_input), error = function(e) NULL) else NULL
  list(
    shapiro = shapiro,
    shapiro_on = if (!is.null(residuals)) "model residuals" else "within-group centred replicate summaries",
    levene = levene_test(y, g),
    variance_ratio = variance_ratio(y, g),
    advisory = paste(
      "These diagnostics are ADVISORY. CellMorphR does not switch test on them:",
      "the test you selected is the test that was run. Use them to decide whether to",
      "transform the data (Data tab) or to choose the Wilcoxon rank-sum test."
    )
  )
}

#' Printable block for the Statistics output and the downloaded report.
format_assumption_diagnostics <- function(d) {
  if (is.null(d)) return("Assumption diagnostics: not available for these groups.")
  out <- character(0)
  if (!is.null(d$shapiro)) {
    out <- c(out, sprintf("Normality  -- Shapiro-Wilk on %s: W = %.4f, p = %s",
                          d$shapiro_on, d$shapiro$statistic, format_p_value(d$shapiro$p.value)))
  } else out <- c(out, "Normality  -- Shapiro-Wilk: not computable (too few or constant values).")
  if (!is.null(d$levene)) {
    out <- c(out, sprintf("Variance   -- %s: F(%d, %d) = %.4f, p = %s",
                          d$levene$method, d$levene$df1, d$levene$df2,
                          d$levene$statistic, format_p_value(d$levene$p.value)))
  } else out <- c(out, "Variance   -- Levene's test: not computable (a group has fewer than 2 replicates).")
  if (!is.null(d$variance_ratio)) {
    out <- c(out, sprintf("Variance   -- largest/smallest group SD = %.2f (%s vs %s)%s",
                          d$variance_ratio$ratio, d$variance_ratio$max_group, d$variance_ratio$min_group,
                          if (d$variance_ratio$ratio > 3) "  <- heterogeneous; Welch's t-test already allows for this" else ""))
  }
  paste(c(out, d$advisory), collapse = "\n")
}

compute_skewness <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 3) return(NA)
  m <- mean(x)
  s <- sd(x)
  if (s == 0) return(0)
  (n / ((n - 1) * (n - 2))) * sum(((x - m) / s)^3)
}


# Suggest an appropriate statistical workflow from the selected comparison design.
# This intentionally suggests and explains; it does not hide the decision from the user.
suggest_statistical_workflow <- function(rep_df, comparison_tbl = NULL, comparison_mode = "vs_reference",
                                         stat_unit = "median_val", selected_group_cols = NULL) {
  if (is.null(rep_df) || nrow(rep_df) == 0) {
    return(list(test = "welch", correction = "holm", global_model = "skip",
                summary = "No replicate-summary data are available yet.",
                details = data.frame()))
  }
  groups <- sort(unique(as.character(rep_df$Statistics_Group)))
  n_groups <- length(groups)
  n_comp <- if (is.null(comparison_tbl)) 0 else nrow(comparison_tbl)

  group_info <- rep_df %>%
    group_by(Statistics_Group) %>%
    summarise(
      n_bio_reps = n_distinct(Replicate),
      n_cells = sum(n_cells, na.rm = TRUE),
      mean_stat = mean(Stat_Value, na.rm = TRUE),
      sd_stat = sd(Stat_Value, na.rm = TRUE),
      .groups = "drop"
    )

  min_reps <- suppressWarnings(min(group_info$n_bio_reps, na.rm = TRUE))
  max_reps <- suppressWarnings(max(group_info$n_bio_reps, na.rm = TRUE))
  if (!is.finite(min_reps)) min_reps <- NA_integer_
  if (!is.finite(max_reps)) max_reps <- NA_integer_

  # Pairwise test suggestion: Welch is the default because it is robust to unequal variances.
  suggested_test <- "welch"
  suggested_correction <- if (n_comp <= 1) "none" else "holm"

  # Optional global model suggestion. Pairwise tests are still run from the explicit comparison table.
  has_condition <- "Condition" %in% names(rep_df)
  has_time <- "Time_min" %in% names(rep_df)
  n_cond <- if (has_condition) length(unique(rep_df$Condition)) else 0
  n_time <- if (has_time) length(unique(rep_df$Time_min)) else 0
  same_reps_across_time <- FALSE
  if (has_time && "Replicate" %in% names(rep_df) && n_time > 1) {
    rep_time_counts <- rep_df %>% distinct(Replicate, Time_min) %>% count(Replicate, name = "n_times")
    same_reps_across_time <- any(rep_time_counts$n_times > 1)
  }
  suggested_global <- "skip"
  if (has_condition && has_time && n_cond > 1 && n_time > 1) suggested_global <- "anova"
  if (same_reps_across_time && has_condition && has_time && n_cond > 1 && n_time > 1) suggested_global <- "lmm"

  unit_label <- ifelse(stat_unit == "mean_val", "replicate means", "replicate medians")
  comparison_label <- switch(comparison_mode,
    vs_reference = "each group against the selected reference",
    all_pairwise = "all pairwise group comparisons",
    within_time = "pairwise comparisons within each time point",
    within_condition = "pairwise comparisons within each condition",
    "the selected comparison table"
  )

  warnings <- character()
  if (!is.na(min_reps) && min_reps < 2) warnings <- c(warnings, "At least one group has fewer than 2 biological replicates. The app can summarize it, but formal pairwise testing is not reliable.")
  if (!is.na(min_reps) && min_reps < 3) warnings <- c(warnings, "At least one group has fewer than 3 biological replicates. Interpret p-values cautiously and emphasize effect sizes/plots.")
  if (!is.na(min_reps) && !is.na(max_reps) && min_reps > 0 && max_reps >= 2 * min_reps) warnings <- c(warnings, "Biological replicate counts are unbalanced across groups.")
  if (n_comp > 15) warnings <- c(warnings, "Many comparisons are selected. Consider filtering brackets/results to the biologically most important comparisons.")

  rationale <- c(
    paste0("Recommended pairwise test: ", ifelse(suggested_test == "welch", "Welch's t-test", "Wilcoxon rank-sum test"), "."),
    paste0("Why: your current setup defines ", n_groups, " statistical group(s) and ", n_comp, " planned comparison(s), using ", unit_label, " as the statistical unit."),
    paste0("The pairwise workflow is: ", comparison_label, "."),
    "Welch's t-test is the safest default for two-group replicate-summary comparisons because it does not assume equal variance.",
    "Raw single-cell values should be visualized, but p-values should be based on biological replicate summaries to avoid pseudoreplication.",
    paste0("Recommended p-value correction: ", ifelse(suggested_correction == "none", "none, because only one comparison is planned", suggested_correction), "."),
    paste0("Optional global model suggestion: ", switch(suggested_global,
      anova = "repeated-measures ANOVA with the dish as the error stratum, because Condition and Time_min appear available with multiple levels",
      lmm = "linear mixed model, because replicate IDs appear across time and a time-course structure may be present",
      skip = "skip global model unless you have a specific factorial question"), ".")
  )
  if (length(warnings) > 0) rationale <- c(rationale, "", "Warnings:", paste0("- ", warnings))

  list(
    test = suggested_test,
    correction = suggested_correction,
    global_model = suggested_global,
    summary = paste(rationale, collapse = "\n"),
    details = group_info
  )
}

# =============================================================================
# DATA TRANSFORMS (one pipeline, always reported)
# =============================================================================
# Every analysis in the app -- the pairwise tests, the global models, the figure
# axes -- must agree about what scale the numbers are on, and the report must
# say so. In particular the log transform applies a DATA-DEPENDENT shift when a
# column contains non-positive values, so the shift constant is part of the
# answer and is carried around with the values (audit CM-13, CM-14).

apply_measure_transform <- function(x, transform = "none") {
  transform <- transform %||% "none"
  if (transform == "log") {
    shift <- 0
    if (any(x <= 0, na.rm = TRUE)) shift <- -min(x, na.rm = TRUE) + 1
    list(values = log(x + shift), transform = "log", shift = shift,
         label = if (shift != 0)
                   paste0("log(x + ", formatC(shift, format = "g", digits = 6, width = 1),
                          ")  [shift chosen so the smallest value maps to 1]")
                 else "log(x)")
  } else if (transform == "sqrt") {
    list(values = sqrt(pmax(x, 0)), transform = "sqrt", shift = 0,
         label = "sqrt(max(x, 0))")
  } else {
    list(values = x, transform = "none", shift = 0, label = "none (raw values)")
  }
}

#' Apply one transform to several measurement columns, keeping a per-column
#' record of exactly what was done (attached as the "transform_info" attribute).
apply_transform_to_columns <- function(df, cols, transform = "none") {
  transform <- transform %||% "none"
  info <- list()
  for (m in cols) {
    if (!m %in% names(df) || !is.numeric(df[[m]])) next
    tr <- tryCatch(apply_measure_transform(df[[m]], transform), error = function(e) NULL)
    if (is.null(tr)) next
    df[[m]] <- tr$values
    info[[m]] <- tr[c("transform", "shift", "label")]
  }
  attr(df, "transform_info") <- info
  attr(df, "transform") <- transform
  df
}

#' One-line description of the transform applied to a single measurement, for
#' reports, captions and axis labels.
describe_transform <- function(df_or_info, measure, transform = NULL) {
  info <- if (is.list(df_or_info) && !is.null(attr(df_or_info, "transform_info")))
            attr(df_or_info, "transform_info") else df_or_info
  entry <- if (is.list(info) && !is.null(info[[measure]])) info[[measure]] else NULL
  if (!is.null(entry)) return(entry$label)
  apply_measure_transform(numeric(0), transform %||% "none")$label
}

# =============================================================================
# ROBUST LMM FITTING
# =============================================================================

# The LMM is fitted on cells, with the experimental design carried by the random
# effects. A random INTERCEPT per dish, (1|Rep_unique), only absorbs variation
# that is constant over the whole time course. Condition varies between dishes,
# so its test is protected by that term -- but Time and Condition:Time vary
# WITHIN a dish, and up to v3.1.2 they were therefore tested against the
# cell-to-cell residual: 99-100% false positives once dishes drift over time
# (audit CM-2). Adding (1|Rep_unique:Time_fac) gives each dish its own offset at
# each timepoint -- focus drift, medium evaporation, per-session segmentation
# quality -- which is the error term those two effects belong to.

LMM_LABEL <- "Linear mixed model on cells (dish and dish x timepoint as random effects)"

fit_lmm_safe <- function(df, measure, transform = "none") {
  tr <- apply_measure_transform(df[[measure]], transform)
  df[[measure]] <- tr$values
  df <- df[!is.na(df[[measure]]), ]
  df$Condition <- droplevels(factor(df$Condition))
  df$Time_fac  <- droplevels(factor(df$Time_min))
  df$Replicate <- factor(df$Replicate)
  df$Rep_unique <- droplevels(interaction(df$Condition, df$Replicate, drop = TRUE))

  no_model <- function(msg) list(
    model = NULL, method = "No global model", warnings = msg,
    fallback_used = FALSE, fitted = FALSE, singular = FALSE,
    transform_used = tr$transform, transform_label = tr$label, transform_shift = tr$shift)

  n_rep_levels <- nlevels(df$Rep_unique)
  if (n_rep_levels < 3) {
    return(no_model(paste0("Only ", n_rep_levels,
      " unique replicate groups (need at least 3). No global model was fitted -- ",
      "the Pairwise results above are still valid.")))
  }
  if (any(table(df$Condition, df$Time_fac) == 0)) {
    return(no_model(paste0("At least one Condition x Timepoint combination has no cells. ",
      "No global model was fitted -- the Pairwise results above are still valid.")))
  }

  # The dish x timepoint term is only estimable when dishes are measured at more
  # than one timepoint with more than one cell per dish per timepoint; otherwise
  # it is confounded with the residual.
  n_times <- nlevels(df$Time_fac)
  cells_per_cell <- table(df$Rep_unique, df$Time_fac)
  can_fit_reptime <- n_times > 1 && min(cells_per_cell[cells_per_cell > 0]) >= 2
  rand <- if (can_fit_reptime) "(1|Rep_unique) + (1|Rep_unique:Time_fac)" else "(1|Rep_unique)"
  formula_lmm <- as.formula(paste0(measure, " ~ Condition * Time_fac + ", rand))

  warn_msgs <- character()
  model <- tryCatch(
    withCallingHandlers({
      if (has_lmerTest) lmerTest::lmer(formula_lmm, data = df, na.action = na.omit,
                       control = lme4::lmerControl(optimizer = "bobyqa",
                                                   optCtrl = list(maxfun = 50000)))
      else lme4::lmer(formula_lmm, data = df, na.action = na.omit,
                       control = lme4::lmerControl(optimizer = "bobyqa",
                                                   optCtrl = list(maxfun = 50000)))
    }, warning = function(w) { warn_msgs <<- c(warn_msgs, conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = function(e) NULL)

  if (is.null(model)) {
    return(no_model(paste0("The mixed model did not converge",
      if (length(warn_msgs)) paste0(" (", paste(unique(warn_msgs), collapse = "; "), ")") else "",
      ". No global model was fitted -- the Pairwise results above are still valid.")))
  }

  notes <- character(0)
  if (!can_fit_reptime) notes <- c(notes, paste0(
    "Only one timepoint, or only one cell per dish per timepoint: the dish x timepoint ",
    "random term is not estimable and was dropped. Time and Condition:Time are NOT ",
    "protected against pseudoreplication in this fit -- prefer the repeated-measures ANOVA."))

  # A singular fit means a variance component was estimated at exactly zero. That
  # is reported, not silently swapped for another model: the fit is still usable
  # (its tests become conservative) and the user needs to know which term
  # collapsed. Up to v3.1.2 this branch returned NULL while promising a fallback
  # ANOVA that was never computed (audit CM-7).
  is_singular <- tryCatch(isTRUE(lme4::isSingular(model)), error = function(e) FALSE)
  if (is_singular) {
    vc <- tryCatch(as.data.frame(lme4::VarCorr(model)), error = function(e) NULL)
    zero_terms <- if (!is.null(vc)) vc$grp[!is.na(vc$vcov) & vc$vcov < 1e-10] else character(0)
    notes <- c(notes, paste0(
      "Singular fit: ",
      if (length(zero_terms)) paste0("the random term(s) ", paste(zero_terms, collapse = ", "),
                                     " were estimated at zero variance. ")
      else "at least one random-effect variance was estimated at zero. ",
      "The model below is still reported, and its tests are conservative rather than wrong, ",
      "but with this few replicates the variance components are barely identified. ",
      "Treat the repeated-measures ANOVA on replicate summaries as the primary global model."))
  }
  if (!has_lmerTest) notes <- c(notes,
    "lmerTest is not installed, so the mixed-model table has no p-value column. Install lmerTest for Satterthwaite degrees of freedom.")

  list(model = model, method = LMM_LABEL,
       formula = formula_lmm,
       random_structure = rand,
       has_reptime_term = can_fit_reptime,
       warnings = c(warn_msgs, notes), notes = notes,
       fallback_used = FALSE, fitted = TRUE, singular = is_singular,
       df_method = if (has_lmerTest) "Satterthwaite (lmerTest)" else "none (lmerTest not installed)",
       transform_used = tr$transform, transform_label = tr$label, transform_shift = tr$shift,
       data_used = df)
}

# =============================================================================
# PCA / MORPHOSPACE HELPER
# =============================================================================

run_pca <- function(df, measure_cols, scale = TRUE) {
  complete_rows <- complete.cases(df[, measure_cols, drop = FALSE])
  df_clean <- df[complete_rows, ]
  if (nrow(df_clean) < 10) return(list(ok = FALSE, msg = "Fewer than 10 complete observations."))
  if (length(measure_cols) < 2) return(list(ok = FALSE, msg = "Need >= 2 numeric features."))
  sds <- sapply(df_clean[measure_cols], sd, na.rm = TRUE)
  valid_cols <- measure_cols[sds > 0]
  if (length(valid_cols) < 2) return(list(ok = FALSE, msg = "Fewer than 2 features with nonzero variance."))
  pca_result <- tryCatch(prcomp(df_clean[valid_cols], center = TRUE, scale. = scale), error = function(e) NULL)
  if (is.null(pca_result)) return(list(ok = FALSE, msg = "PCA computation failed."))
  list(ok = TRUE, pca = pca_result, df = df_clean, features_used = valid_cols,
       var_explained = summary(pca_result)$importance)
}

# =============================================================================
# CLUSTERING HELPER
# =============================================================================

find_optimal_k <- function(data_matrix, max_k = 10) {
  max_k <- min(max_k, nrow(data_matrix) - 1, 15)
  if (max_k < 2) return(list(optimal_k = 2, wss = NULL, sil = NULL))
  wss <- sil <- numeric(max_k); wss[1] <- sil[1] <- NA
  for (k in 2:max_k) {
    km <- tryCatch(kmeans(data_matrix, centers = k, nstart = 25, iter.max = 100), error = function(e) NULL)
    if (is.null(km)) { wss[k] <- NA; sil[k] <- NA; next }
    wss[k] <- km$tot.withinss
    if (has_cluster) { ss <- tryCatch(cluster::silhouette(km$cluster, dist(data_matrix)), error = function(e) NULL)
      sil[k] <- if (!is.null(ss)) mean(ss[, 3]) else NA }
  }
  valid_sil <- sil[2:max_k]
  optimal_k <- if (any(!is.na(valid_sil))) which.max(valid_sil) + 1 else 3
  list(optimal_k = optimal_k, wss = wss[2:max_k], sil = sil[2:max_k], k_range = 2:max_k)
}

# =============================================================================
# MULTI-FEATURE ANALYSIS HELPERS
# =============================================================================

# -----------------------------------------------------------------------------
# Repeated-measures ANOVA on replicate-level summaries.
#
# In a time course the SAME dish is imaged at every timepoint, so a dish is one
# experimental unit, not k of them. Up to v3.1.2 this function fitted a plain
# `aov(y ~ Condition * Time_fac)` on the replicate summaries, which counted each
# dish once per timepoint and gave a measured 41% false-positive rate for the
# Condition effect on null data (audit CM-1).
#
# The model below puts the dish in its own error stratum, so each term is tested
# against the variation it actually varies over:
#   * Condition varies BETWEEN dishes -> tested on dish-to-dish variation
#     (Error: Rep_unique), df = n_dishes - n_conditions
#   * Time and Condition:Time vary WITHIN a dish -> tested on the dish x time
#     variation (Error: Within), never on the cell-to-cell residual
#
# Note there is exactly one summary value per (dish, timepoint), so the within
# stratum residual IS the dish x time interaction. Nesting time explicitly --
# Error(Rep_unique/Time_fac) -- would leave that stratum with zero residual df
# and no p-values at all, which is why Error(Rep_unique) is used.
# -----------------------------------------------------------------------------

RM_ANOVA_LABEL <- "Repeated-measures ANOVA on replicate summaries (dish as error stratum)"

fit_rm_anova <- function(rs, yvar = "median_val") {
  if (is.null(rs) || nrow(rs) == 0 || !yvar %in% names(rs)) {
    return(list(ok = FALSE, msg = "No replicate-level summaries available."))
  }
  if (!all(c("Condition", "Replicate") %in% names(rs))) {
    return(list(ok = FALSE, msg = "The global model needs Condition and Replicate columns."))
  }
  rs <- rs[!is.na(rs[[yvar]]), , drop = FALSE]
  rs$Condition  <- droplevels(factor(rs$Condition))
  rs$Replicate  <- factor(rs$Replicate)
  rs$Rep_unique <- droplevels(interaction(rs$Condition, rs$Replicate, drop = TRUE))
  has_time <- "Time_min" %in% names(rs)
  rs$Time_fac <- if (has_time) droplevels(factor(rs$Time_min)) else factor(rep("all", nrow(rs)))

  n_cond <- nlevels(rs$Condition)
  n_time <- nlevels(rs$Time_fac)
  n_dish <- nlevels(rs$Rep_unique)
  if (n_cond < 2) return(list(ok = FALSE, msg = "Need at least two conditions for a global model."))
  if (n_dish <= n_cond) {
    return(list(ok = FALSE, msg = paste0(
      "Only ", n_dish, " biological replicates across ", n_cond,
      " conditions: there is no dish-to-dish variation left to test Condition against. ",
      "A global model needs at least two replicates in at least one condition.")))
  }

  repeated <- n_time > 1
  form <- if (repeated) {
    as.formula(paste0(yvar, " ~ Condition * Time_fac + Error(Rep_unique)"))
  } else {
    # One timepoint: the design collapses to a one-way ANOVA over dishes, which
    # already has the dish as its unit. No error stratum is needed or estimable.
    as.formula(paste0(yvar, " ~ Condition"))
  }
  model <- tryCatch(aov(form, data = rs), error = function(e) NULL)
  if (is.null(model)) return(list(ok = FALSE, msg = "The global model could not be fitted."))

  tbl <- tidy_rm_anova(model)
  if (is.null(tbl) || nrow(tbl) == 0) {
    return(list(ok = FALSE, msg = "The global model produced no testable terms."))
  }
  notes <- character(0)
  if (!repeated) notes <- c(notes,
    "Only one timepoint: fitted as a one-way ANOVA on replicate summaries (no Time term).")
  if (any(is.na(tbl$p))) notes <- c(notes,
    "Some terms have no residual degrees of freedom left and cannot be tested; add replicates or timepoints.")

  list(ok = TRUE, model = model, table = tbl, data = rs,
       method = if (repeated) RM_ANOVA_LABEL else
         "One-way ANOVA on replicate summaries (one timepoint; dish is the unit)",
       repeated = repeated, n_dishes = n_dish, n_conditions = n_cond, n_times = n_time,
       notes = notes, yvar = yvar)
}

#' Flatten an aov / aovlist summary into one tidy table, recording which error
#' stratum -- and therefore which denominator df -- tested each term.
tidy_rm_anova <- function(model) {
  s <- summary(model)
  strata <- if (inherits(model, "aovlist")) s else list("Error: Within" = s)
  out <- lapply(names(strata), function(nm) {
    tab <- strata[[nm]]
    if (is.list(tab) && !is.data.frame(tab)) tab <- tab[[1]]
    if (is.null(tab)) return(NULL)
    terms <- trimws(rownames(tab))
    ss <- tab[["Sum Sq"]]
    resid_i <- which(terms == "Residuals")
    ss_resid <- if (length(resid_i) == 1) ss[resid_i] else NA_real_
    df_resid <- if (length(resid_i) == 1) tab[["Df"]][resid_i] else NA_real_
    keep <- setdiff(seq_along(terms), resid_i)
    if (length(keep) == 0) return(NULL)
    data.frame(
      Term = terms[keep],
      Error_Stratum = sub("^Error: ", "", nm),
      Df = tab[["Df"]][keep],
      DenDF = df_resid,
      Sum_Sq = ss[keep],
      Mean_Sq = tab[["Mean Sq"]][keep],
      F_value = tab[["F value"]][keep],
      p = tab[["Pr(>F)"]][keep],
      partial_eta2 = if (is.na(ss_resid)) NA_real_ else ss[keep] / (ss[keep] + ss_resid),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out[!vapply(out, is.null, logical(1))])
  if (is.null(out)) return(NULL)
  rownames(out) <- NULL
  out
}

#' Human-readable print of a fit_rm_anova() result, used by the Statistics tab
#' and the downloaded report.
format_rm_anova <- function(fit, digits = 4) {
  if (is.null(fit) || !isTRUE(fit$ok)) {
    return(paste("No global model:", fit$msg %||% "not available."))
  }
  tb <- fit$table
  df <- data.frame(
    Term = tb$Term,
    `Tested against` = tb$Error_Stratum,
    Df = tb$Df, DenDF = tb$DenDF,
    `Sum Sq` = signif(tb$Sum_Sq, digits), `Mean Sq` = signif(tb$Mean_Sq, digits),
    `F` = signif(tb$F_value, digits),
    `p` = format_p_value(tb$p),
    `partial eta2` = round(tb$partial_eta2, 3),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  lines <- c(
    fit$method,
    paste0("Unit of replication: biological replicate (", fit$n_dishes,
           " replicate groups across ", fit$n_conditions, " conditions",
           if (fit$repeated) paste0(", ", fit$n_times, " timepoints") else "", ")"),
    paste(capture.output(print(df, row.names = FALSE)), collapse = "\n")
  )
  if (length(fit$notes) > 0) lines <- c(lines, paste("Note:", fit$notes))
  paste(lines, collapse = "\n")
}

run_anova_for_feature <- function(df, measure, yvar = "median_val") {
  rs <- compute_replicate_summaries(df, measure)
  fit <- fit_rm_anova(rs, yvar = yvar)
  if (!isTRUE(fit$ok)) return(NULL)
  tb <- fit$table
  pick <- function(term, col) {
    i <- match(term, tb$Term)
    if (is.na(i)) NA_real_ else tb[[col]][i]
  }
  data.frame(
    Feature = measure,
    Condition_p = pick("Condition", "p"),
    Time_p = pick("Time_fac", "p"),
    Interaction_p = pick("Condition:Time_fac", "p"),
    Condition_partial_eta2 = round(pick("Condition", "partial_eta2"), 4),
    Time_partial_eta2 = round(pick("Time_fac", "partial_eta2"), 4),
    Interaction_partial_eta2 = round(pick("Condition:Time_fac", "partial_eta2"), 4),
    stringsAsFactors = FALSE
  )
}

detect_earliest_divergence <- function(df, measure, ref_cond, alpha = 0.05) {
  rs <- compute_replicate_summaries(df, measure)
  times <- sort(unique(rs$Time_min))
  trt_conds <- setdiff(unique(as.character(rs$Condition)), ref_cond)
  results <- list()
  for (trt in trt_conds) {
    p_vals <- sapply(times, function(tm) {
      ctrl <- rs$median_val[rs$Condition == ref_cond & rs$Time_min == tm]
      treat <- rs$median_val[rs$Condition == trt & rs$Time_min == tm]
      if (length(ctrl) < 2 || length(treat) < 2) return(NA)
      tryCatch(t.test(treat, ctrl, var.equal = FALSE)$p.value, error = function(e) NA) })
    p_adj <- p.adjust(p_vals, method = "holm")
    first_sig <- which(p_adj < alpha)
    earliest <- if (length(first_sig) > 0) times[first_sig[1]] else NA
    earliest_p <- if (length(first_sig) > 0) p_adj[first_sig[1]] else NA
    results[[length(results) + 1]] <- data.frame(Feature = measure, Treatment = trt,
      Earliest_Significant_Time = earliest, p_adj_at_earliest = round(earliest_p, 4), stringsAsFactors = FALSE)
  }
  do.call(rbind, results)
}

# =============================================================================
# CELL METRICS HELPER
# =============================================================================
compute_cell_metrics <- function(df, feature) {
  df %>%
    group_by(Condition, Time_min) %>%
    summarise(
      n_cells    = n(),
      mean_val   = mean(.data[[feature]], na.rm = TRUE),
      median_val = median(.data[[feature]], na.rm = TRUE),
      sd_val     = sd(.data[[feature]], na.rm = TRUE),
      iqr_val    = IQR(.data[[feature]], na.rm = TRUE),
      cv_val     = if (abs(mean(.data[[feature]], na.rm = TRUE)) > 1e-10)
                     sd(.data[[feature]], na.rm = TRUE) / abs(mean(.data[[feature]], na.rm = TRUE)) * 100 else NA_real_,
      skewness   = compute_skewness(.data[[feature]]),
      q25        = quantile(.data[[feature]], 0.25, na.rm = TRUE),
      q75        = quantile(.data[[feature]], 0.75, na.rm = TRUE),
      .groups    = "drop"
    )
}

# =============================================================================
# RAW HEATMAP DATA HELPER
# =============================================================================
compute_raw_heatmap_data <- function(df, features, stat_fn = "mean",
                                     time_fmt = function(v) paste0(v, " min")) {
  fn_val <- if (stat_fn == "median") "median_val" else "mean_val"
  rs_list <- lapply(features, function(m) {
    compute_replicate_summaries(df, m) %>%
      group_by(Condition, Time_min) %>%
      summarise(value = mean(.data[[fn_val]], na.rm = TRUE), .groups = "drop") %>%
      mutate(Feature = m)
  })
  res <- bind_rows(rs_list)
  # Single-line label: these are drawn on an angled axis, where an embedded
  # line break makes the two lines cross each other.
  res$Group <- paste0(res$Condition, " | ", time_fmt(res$Time_min))
  res
}

# =============================================================================
# SPIDER/RADAR PLOT HELPERS
# =============================================================================
prep_radar_data <- function(df, features, stat_fn = "mean") {
  fn_val <- if (stat_fn == "median") "median_val" else "mean_val"
  rs_list <- lapply(features, function(m) {
    compute_replicate_summaries(df, m) %>%
      group_by(Condition) %>%
      summarise(value = mean(.data[[fn_val]], na.rm = TRUE), .groups = "drop") %>%
      mutate(Feature = m)
  })
  raw <- bind_rows(rs_list)
  raw <- raw %>% group_by(Feature) %>%
    mutate(scaled = (value - min(value, na.rm = TRUE)) /
                    (max(value, na.rm = TRUE) - min(value, na.rm = TRUE) + 1e-10)) %>%
    ungroup()
  raw
}

plot_radar <- function(radar_df, cols, base_sz = 12, show_points = TRUE,
                        show_grid = TRUE, alpha_fill = 0.2) {
  features <- unique(radar_df$Feature)
  n_feat <- length(features)
  if (n_feat < 3) return(ggplot() + annotate("text", x = 0.5, y = 0.5,
    label = "Need >= 3 features for radar chart", size = 5) + theme_void())

  angle_df <- data.frame(Feature = features,
    angle = seq(0, 2 * pi * (1 - 1 / n_feat), length.out = n_feat),
    stringsAsFactors = FALSE)
  radar_df <- merge(radar_df, angle_df, by = "Feature")
  radar_df$x <- radar_df$scaled * cos(radar_df$angle - pi / 2)
  radar_df$y <- radar_df$scaled * sin(radar_df$angle - pi / 2)

  # Close each polygon
  close_rows <- radar_df %>% group_by(Condition) %>%
    arrange(angle) %>% dplyr::slice(1) %>% ungroup()
  radar_df <- bind_rows(radar_df %>% arrange(Condition, angle), close_rows) %>%
    arrange(Condition, angle)

  # Axis labels
  axis_df <- angle_df %>% mutate(
    x = 1.18 * cos(angle - pi / 2), y = 1.18 * sin(angle - pi / 2))

  # Grid rings
  grid_rings <- do.call(rbind, lapply(c(0.25, 0.5, 0.75, 1.0), function(r) {
    a <- seq(0, 2 * pi, length.out = 100)
    data.frame(x = r * cos(a - pi / 2), y = r * sin(a - pi / 2), r = r)
  }))

  # Spoke lines
  spoke_df <- data.frame(x0 = 0, y0 = 0,
    x1 = cos(angle_df$angle - pi / 2), y1 = sin(angle_df$angle - pi / 2))

  p <- ggplot()
  if (show_grid) {
    p <- p + geom_path(data = grid_rings, aes(x = x, y = y, group = r),
                       color = "grey80", linewidth = 0.3) +
      geom_segment(data = spoke_df, aes(x = x0, y = y0, xend = x1, yend = y1),
                   color = "grey70", linewidth = 0.3)
  }
  p <- p + geom_polygon(data = radar_df,
    aes(x = x, y = y, fill = Condition, color = Condition, group = Condition),
    alpha = alpha_fill, linewidth = 0.8)
  if (show_points) {
    p <- p + geom_point(data = radar_df,
      aes(x = x, y = y, color = Condition), size = 3, alpha = 0.9)
  }
  p <- p + geom_text(data = axis_df, aes(x = x, y = y, label = Feature),
    size = base_sz / 3.5, fontface = "bold", color = "grey20") +
    scale_fill_manual(values = cols) + scale_color_manual(values = cols) +
    # Extra room around the wheel so the outermost feature labels (drawn at
    # radius 1.18) are not clipped at the left/right edges of the panel.
    coord_fixed(xlim = c(-1.6, 1.6), ylim = c(-1.3, 1.3)) +
    theme_void(base_size = base_sz) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold", hjust = 0.5))
  p
}

# ---- END OF STATISTICAL CORE (tests/helper-cellmorphr.R sources up to here) ----
# Everything above this marker is plain R with no Shiny reactivity, so the test
# suite can source it standalone. Keep new statistical helpers above the marker.

# =============================================================================
# SHARED FIGURE CONTROLS (reusable UI component)
# =============================================================================

figure_controls_ui <- function(prefix, def_width = 8, def_height = 5, show_legend = TRUE,
                                show_grid = TRUE, show_theme = TRUE) {
  face_choices <- c("Bold" = "bold", "Plain" = "plain",
                    "Italic" = "italic", "Bold Italic" = "bold.italic")
  controls <- tagList(
    tags$p(tags$strong("Figure Appearance"), class = "mb-1 mt-2"),
    sliderInput(paste0(prefix, "_base_size"), "Base font size", 8, 26, 12, 1),
    selectInput(paste0(prefix, "_font_family"), "Font family",
      choices = c("Sans-serif (Arial)" = "sans", "Serif (Times)" = "serif",
                  "Mono (Courier)" = "mono"),
      selected = "sans"),
    if (show_legend) selectInput(paste0(prefix, "_legend_pos"), "Legend position",
      choices = c("Top" = "top", "Right" = "right", "Bottom" = "bottom",
                  "Left" = "left", "None" = "none"), selected = "top") else NULL,
    if (show_grid) checkboxInput(paste0(prefix, "_grid"), "Show gridlines", FALSE) else NULL,
    if (show_theme) checkboxInput(paste0(prefix, "_border"), "Panel border (instead of axes)", FALSE) else NULL,
    # --- Per-element typography (collapsible) ---
    tags$details(class = "mt-2",
      tags$summary(tags$strong("Typography details", style = "cursor:pointer; font-size:0.88rem;")),
      tags$div(class = "ps-1 pt-2",
        tags$p(class = "small text-muted mb-1", "Override individual sizes and styles."),
        # Title
        tags$label("Title", class = "small fw-bold"),
        layout_columns(col_widths = c(6, 6), gap = "4px",
          numericInput(paste0(prefix, "_title_size"), NULL, value = NA, min = 6, max = 36, step = 1),
          selectInput(paste0(prefix, "_title_face"), NULL, choices = face_choices, selected = "bold")
        ),
        # Axis labels
        tags$label("Axis labels", class = "small fw-bold"),
        layout_columns(col_widths = c(6, 6), gap = "4px",
          numericInput(paste0(prefix, "_axis_title_size"), NULL, value = NA, min = 6, max = 36, step = 1),
          selectInput(paste0(prefix, "_axis_title_face"), NULL, choices = face_choices, selected = "bold")
        ),
        # Axis tick text
        tags$label("Axis tick text", class = "small fw-bold"),
        layout_columns(col_widths = c(6, 6), gap = "4px",
          numericInput(paste0(prefix, "_axis_text_size"), NULL, value = NA, min = 6, max = 36, step = 1),
          selectInput(paste0(prefix, "_axis_text_face"), NULL, choices = face_choices, selected = "plain")
        ),
        # Legend title
        tags$label("Legend title", class = "small fw-bold"),
        layout_columns(col_widths = c(6, 6), gap = "4px",
          numericInput(paste0(prefix, "_legend_title_size"), NULL, value = NA, min = 6, max = 36, step = 1),
          selectInput(paste0(prefix, "_legend_title_face"), NULL, choices = face_choices, selected = "bold")
        ),
        # Legend text
        tags$label("Legend text", class = "small fw-bold"),
        layout_columns(col_widths = c(6, 6), gap = "4px",
          numericInput(paste0(prefix, "_legend_text_size"), NULL, value = NA, min = 6, max = 36, step = 1),
          selectInput(paste0(prefix, "_legend_text_face"), NULL, choices = face_choices, selected = "plain")
        ),
        # Strip (facet) labels
        tags$label("Facet labels", class = "small fw-bold"),
        layout_columns(col_widths = c(6, 6), gap = "4px",
          numericInput(paste0(prefix, "_strip_size"), NULL, value = NA, min = 6, max = 36, step = 1),
          selectInput(paste0(prefix, "_strip_face"), NULL, choices = face_choices, selected = "bold")
        ),
        tags$p(class = "small text-muted mt-1", "Leave size blank to auto-scale from base."),
        actionButton(paste0(prefix, "_apply_global"), "Apply to all tabs",
          class = "btn-sm btn-outline-secondary w-100 mt-2", icon = icon("globe"))
      )
    ),
    # --- Custom labels (collapsible) ---
    tags$details(class = "mt-2",
      tags$summary(tags$strong("Labels & caption", style = "cursor:pointer; font-size:0.88rem;")),
      tags$div(class = "ps-1 pt-2",
        tags$p(class = "small text-muted mb-1", "Leave blank to keep automatic text. Applies to this tab's main plot."),
        textInput(paste0(prefix, "_custom_title"),    "Title", ""),
        textInput(paste0(prefix, "_custom_subtitle"), "Subtitle", ""),
        textInput(paste0(prefix, "_custom_x_label"),  "X axis label", ""),
        textInput(paste0(prefix, "_custom_y_label"),  "Y axis label", ""),
        textInput(paste0(prefix, "_custom_legend"),   "Legend title", ""),
        textInput(paste0(prefix, "_custom_caption"),  "Caption", "")
      )
    ),
    tags$hr(),
    tags$p(tags$strong("Export Settings"), class = "mb-1"),
    numericInput(paste0(prefix, "_width"), "Width (in)", def_width, 3, 20, 0.5),
    numericInput(paste0(prefix, "_height"), "Height (in)", def_height, 2, 16, 0.5),
    numericInput(paste0(prefix, "_dpi"), "DPI", 300, 72, 600, 50),
    selectInput(paste0(prefix, "_format"), "File format",
      choices = c("PDF" = "pdf", "PNG" = "png", "SVG" = "svg", "TIFF" = "tiff", "EPS" = "eps",
                  "PPTX (editable)" = "pptx"), selected = "pdf")
  )
  controls
}

# Helper: get theme from shared controls
# Reads base + per-element typography overrides from input[[prefix_*]] inputs
get_theme <- function(input, prefix) {
  base <- input[[paste0(prefix, "_base_size")]] %||% 12
  legend <- input[[paste0(prefix, "_legend_pos")]] %||% "top"
  grid <- isTRUE(input[[paste0(prefix, "_grid")]])
  border <- isTRUE(input[[paste0(prefix, "_border")]])
  font_fam <- input[[paste0(prefix, "_font_family")]] %||% "sans"
  # Per-element overrides: NA means "use default from base_size"
  na_to_null <- function(x) { if (is.null(x) || is.na(x)) NULL else x }
  theme_publication(
    base_size = base, legend_pos = legend, show_grid = grid,
    show_border = border, font_family = font_fam,
    title_size        = na_to_null(input[[paste0(prefix, "_title_size")]]),
    title_face        = input[[paste0(prefix, "_title_face")]] %||% "bold",
    axis_title_size   = na_to_null(input[[paste0(prefix, "_axis_title_size")]]),
    axis_title_face   = input[[paste0(prefix, "_axis_title_face")]] %||% "bold",
    axis_text_size    = na_to_null(input[[paste0(prefix, "_axis_text_size")]]),
    axis_text_face    = input[[paste0(prefix, "_axis_text_face")]] %||% "plain",
    legend_title_size = na_to_null(input[[paste0(prefix, "_legend_title_size")]]),
    legend_title_face = input[[paste0(prefix, "_legend_title_face")]] %||% "bold",
    legend_text_size  = na_to_null(input[[paste0(prefix, "_legend_text_size")]]),
    legend_text_face  = input[[paste0(prefix, "_legend_text_face")]] %||% "plain",
    strip_size        = na_to_null(input[[paste0(prefix, "_strip_size")]]),
    strip_face        = input[[paste0(prefix, "_strip_face")]] %||% "bold"
  )
}

# Helper: apply user label overrides from the shared "Labels & caption" controls.
# Must run as the LAST step of each plot reactive so downloads (which reuse the
# reactives) match the on-screen preview exactly.
apply_custom_labs <- function(p, input, prefix) {
  get_txt <- function(id) {
    v <- input[[paste0(prefix, id)]]
    if (!is.null(v) && nzchar(trimws(v))) v else NULL
  }
  ttl <- get_txt("_custom_title");    if (!is.null(ttl)) p <- p + labs(title = ttl)
  sub <- get_txt("_custom_subtitle"); if (!is.null(sub)) p <- p + labs(subtitle = sub)
  xl  <- get_txt("_custom_x_label");  if (!is.null(xl))  p <- p + labs(x = xl)
  yl  <- get_txt("_custom_y_label");  if (!is.null(yl))  p <- p + labs(y = yl)
  lg  <- get_txt("_custom_legend")
  if (!is.null(lg)) {
    if (!is.null(p$labels$fill)) p <- p + labs(fill = lg)
    if (!is.null(p$labels$colour)) p <- p + labs(color = lg)
  }
  cap <- get_txt("_custom_caption");  if (!is.null(cap)) p <- p + labs(caption = cap)
  # Fold long text so it never runs past the figure edge (auto or custom alike)
  p$labels$title    <- wrap_plot_text(p$labels$title,    width = 70)
  p$labels$subtitle <- wrap_plot_text(p$labels$subtitle, width = 105)
  p$labels$caption  <- wrap_plot_text(p$labels$caption,  width = 110)
  wrap_legend_labels(p)
}


# =============================================================================
# UI
# =============================================================================
## ---- Inline SVG Strings for Decorative Cells ----
svg_bacterium_rod <- '<svg viewBox="0 0 120 40" xmlns="http://www.w3.org/2000/svg"><rect x="10" y="5" width="100" height="30" rx="15" ry="15" fill="#4FC3F7" opacity="0.18" stroke="#4FC3F7" stroke-width="1" stroke-opacity="0.25"/><ellipse cx="35" cy="20" rx="4" ry="3" fill="#81C784" opacity="0.22"/><ellipse cx="75" cy="20" rx="4" ry="3" fill="#81C784" opacity="0.22"/><path d="M5 20 Q-8 8 -2 0" fill="none" stroke="#4FC3F7" stroke-width="1.2" opacity="0.2"/><path d="M5 22 Q-10 30 -4 38" fill="none" stroke="#4FC3F7" stroke-width="1.2" opacity="0.2"/></svg>'

svg_phage_t4 <- '<svg viewBox="0 0 80 140" xmlns="http://www.w3.org/2000/svg"><polygon points="40,5 65,20 65,45 40,60 15,45 15,20" fill="#81C784" opacity="0.15" stroke="#81C784" stroke-width="1" stroke-opacity="0.25"/><line x1="40" y1="60" x2="40" y2="105" stroke="#81C784" stroke-width="2" opacity="0.2"/><rect x="30" y="100" width="20" height="6" rx="2" fill="#81C784" opacity="0.18"/><line x1="30" y1="106" x2="15" y2="130" stroke="#81C784" stroke-width="1.2" opacity="0.18"/><line x1="50" y1="106" x2="65" y2="130" stroke="#81C784" stroke-width="1.2" opacity="0.18"/><line x1="35" y1="106" x2="25" y2="135" stroke="#81C784" stroke-width="1" opacity="0.15"/><line x1="45" y1="106" x2="55" y2="135" stroke="#81C784" stroke-width="1" opacity="0.15"/><line x1="40" y1="106" x2="40" y2="132" stroke="#81C784" stroke-width="1" opacity="0.15"/></svg>'

svg_coccus <- '<svg viewBox="0 0 50 50" xmlns="http://www.w3.org/2000/svg"><circle cx="25" cy="25" r="20" fill="#FFB74D" opacity="0.14" stroke="#FFB74D" stroke-width="1" stroke-opacity="0.22"/><line x1="25" y1="6" x2="25" y2="44" stroke="#FFB74D" stroke-width="0.8" opacity="0.15" stroke-dasharray="3,2"/></svg>'

svg_phage_small <- '<svg viewBox="0 0 50 90" xmlns="http://www.w3.org/2000/svg"><polygon points="25,3 42,14 42,33 25,44 8,33 8,14" fill="#CE93D8" opacity="0.14" stroke="#CE93D8" stroke-width="0.8" stroke-opacity="0.2"/><line x1="25" y1="44" x2="25" y2="72" stroke="#CE93D8" stroke-width="1.5" opacity="0.16"/><line x1="18" y1="72" x2="10" y2="88" stroke="#CE93D8" stroke-width="1" opacity="0.14"/><line x1="32" y1="72" x2="40" y2="88" stroke="#CE93D8" stroke-width="1" opacity="0.14"/><line x1="25" y1="72" x2="25" y2="86" stroke="#CE93D8" stroke-width="1" opacity="0.14"/></svg>'

# =============================================================================
# CUSTOM CSS FOR DARK THEME
# =============================================================================
dark_css <- "
/* ---------- Soft Dark Base ---------- */
body { background-color: #1e1e2e !important; color: #e0e0e8 !important; }
.navbar { background-color: #16162a !important; border-bottom: 1px solid #3a3a55; }
.navbar-brand, .navbar .nav-link { color: #e0e0e8 !important; }
.navbar .nav-link:hover, .navbar .nav-link.active { color: #4FC3F7 !important; }
.nav-underline .nav-link.active { border-bottom-color: #4FC3F7 !important; }

/* ---------- Sidebar ---------- */
.bslib-sidebar-layout > .sidebar { background-color: #242438 !important; border-right: 1px solid #3a3a55 !important; }
.sidebar .sidebar-title { color: #b0b0c8 !important; font-weight: 600; }
.sidebar label, .sidebar .control-label { color: #c8c8da !important; }
.sidebar .form-text, .sidebar .text-muted { color: #8888a8 !important; }
.sidebar p, .sidebar li, .sidebar span { color: #b8b8d0 !important; }
.sidebar code { color: #4FC3F7 !important; background: #2a2a42 !important; }
.sidebar hr { border-color: #3a3a55 !important; }

/* ---------- Cards ---------- */
.card { background-color: #252540 !important; border: 1px solid #3a3a55 !important; border-radius: 10px !important; }
.card-header { background-color: #2a2a48 !important; color: #c8c8e0 !important; border-bottom: 1px solid #3a3a55 !important; }
.card-body { color: #d0d0e0 !important; }

/* ---------- Inputs ---------- */
.form-control, .form-select, .shiny-input-container input[type='text'],
.shiny-input-container input[type='number'] {
  background-color: #2a2a42 !important; color: #e0e0e8 !important;
  border: 1px solid #4a4a65 !important; border-radius: 6px !important;
}
.form-control:focus, .form-select:focus { border-color: #4FC3F7 !important;
  box-shadow: 0 0 0 0.15rem rgba(79,195,247,0.25) !important; }
.btn-primary { background-color: #4FC3F7 !important; border-color: #4FC3F7 !important; color: #111 !important; font-weight: 600; }
.btn-primary:hover { background-color: #39b0e0 !important; }
.btn-outline-primary { color: #4FC3F7 !important; border-color: #4FC3F7 !important; }
.btn-outline-primary:hover { background-color: #4FC3F7 !important; color: #111 !important; }
.btn-outline-secondary { color: #81C784 !important; border-color: #81C784 !important; }
.btn-outline-secondary:hover { background-color: #81C784 !important; color: #111 !important; }

/* ---------- Radio / Checkbox ---------- */
.form-check-label { color: #c8c8da !important; }
.form-check-input { background-color: #2a2a42 !important; border-color: #5a5a75 !important; }
.form-check-input:checked { background-color: #4FC3F7 !important; border-color: #4FC3F7 !important; }
.shiny-options-group label { color: #c8c8da !important; }

/* ---------- Tabs inside cards ---------- */
.nav-tabs { border-bottom-color: #3a3a55 !important; }
.nav-tabs .nav-link { color: #9898b8 !important; }
.nav-tabs .nav-link:hover { color: #c8c8e0 !important; border-color: transparent !important; }
.nav-tabs .nav-link.active { background-color: #252540 !important; color: #4FC3F7 !important;
  border-color: #3a3a55 #3a3a55 #252540 !important; }

/* ---------- DataTables ---------- */
.dataTables_wrapper { color: #c8c8da !important; }
table.dataTable { background-color: #252540 !important; color: #d0d0e0 !important; }
table.dataTable thead th { background-color: #2a2a48 !important; color: #b8b8d0 !important;
  border-bottom-color: #4a4a65 !important; }
table.dataTable tbody tr { background-color: #252540 !important; }
table.dataTable tbody tr:hover { background-color: #30305a !important; }
table.dataTable tbody td { border-top-color: #3a3a55 !important; }
.dataTables_info, .dataTables_paginate .paginate_button { color: #9898b8 !important; }
.dataTables_paginate .paginate_button.current { background: #4FC3F7 !important; color: #111 !important;
  border-color: #4FC3F7 !important; }
div.dataTables_filter input { background-color: #2a2a42 !important; color: #e0e0e8 !important;
  border: 1px solid #4a4a65 !important; }

/* ---------- Value Boxes ---------- */
.bslib-value-box { border-radius: 10px !important; }
.bslib-value-box .value-box-title { color: #b0b0c8 !important; font-size: 0.85rem; }
.bslib-value-box .value-box-value { color: #e8e8f4 !important; font-weight: 700; }

/* ---------- Alerts ---------- */
.alert-info { background-color: #1e2e3e !important; color: #90c8e8 !important; border-color: #2a4a5e !important; }

/* ---------- Scrollbar ---------- */
::-webkit-scrollbar { width: 8px; height: 8px; }
::-webkit-scrollbar-track { background: #1e1e2e; }
::-webkit-scrollbar-thumb { background: #4a4a65; border-radius: 4px; }
::-webkit-scrollbar-thumb:hover { background: #5a5a75; }

/* ---------- Conditional Panel ---------- */
.shiny-html-output p, .shiny-html-output li, .shiny-html-output span { color: #c8c8da; }
.shiny-html-output h4, .shiny-html-output h5, .shiny-html-output h6 { color: #e0e0e8; }
.shiny-html-output strong { color: #e8e8f4; }

/* ---------- Hero Landing Page ---------- */
.hero-landing { text-align: center; padding: 2.5rem 1rem 2rem; position: relative; overflow: hidden; }
.hero-landing h1 { font-size: 2.8rem; font-weight: 800; color: #4FC3F7;
  letter-spacing: -0.5px; margin-bottom: 0.2rem; }
.hero-landing h1 .morph-highlight { color: #81C784; }
.hero-landing .version-tag { color: #8888a8; font-size: 1rem; font-weight: 300; }
.hero-landing .mantra { color: #b0b0c8; font-size: 1.05rem; font-style: italic;
  margin: 0.8rem 0 0.4rem; }
.hero-landing .motto { color: #FFB74D; font-size: 0.92rem; font-weight: 700;
  letter-spacing: 1.5px; text-transform: uppercase; margin-top: 0.3rem; }
.hero-landing .trademark { color: #6a6a88; font-size: 0.78rem; margin-top: 1.2rem; }
.hero-card { background: #2a2a48 !important; border: 1px solid #3a3a55 !important;
  border-radius: 12px !important; transition: border-color 0.2s; }
.hero-card:hover { border-color: #4FC3F7 !important; }
.hero-card .card-header { font-weight: 600; }
.hero-card h6 { color: #c8c8e0; margin-bottom: 0.5rem; }
.feature-icon-card { text-align: center; padding: 1rem 0.5rem; }
.feature-icon-card i { font-size: 1.6rem; color: #4FC3F7; margin-bottom: 0.4rem; }
.feature-icon-card .feat-label { font-size: 0.8rem; color: #9898b8; line-height: 1.3; }

/* ---------- Decorative cells (background) ---------- */
.deco-cell { position: absolute; pointer-events: none; z-index: 0; }
.hero-content { position: relative; z-index: 1; }

/* ---------- Misc fixes ---------- */
.progress-bar { background-color: #4FC3F7 !important; }
.modal-content { background-color: #252540 !important; color: #e0e0e8 !important; }
.small, small { color: #9898b8 !important; }
.text-muted { color: #8888a8 !important; }
p strong, li strong { color: #d8d8ec !important; }
"

light_css <- "
/* ---------- Light Mode Overrides ---------- */
body { background-color: #f8f9fa !important; color: #2c3e50 !important; }
.navbar { background-color: #2C3E50 !important; border-bottom: 1px solid #ddd; }
.navbar-brand, .navbar .nav-link { color: #ecf0f1 !important; }
.navbar .nav-link:hover, .navbar .nav-link.active { color: #18BC9C !important; }

.bslib-sidebar-layout > .sidebar { background-color: #ffffff !important; border-right: 1px solid #dee2e6 !important; }
.sidebar .sidebar-title { color: #2c3e50 !important; }
.sidebar label, .sidebar .control-label { color: #495057 !important; }
.sidebar .form-text, .sidebar .text-muted { color: #6c757d !important; }
.sidebar p, .sidebar li, .sidebar span { color: #495057 !important; }
.sidebar code { color: #e74c3c !important; background: #f1f3f5 !important; }
.sidebar hr { border-color: #dee2e6 !important; }

.card { background-color: #ffffff !important; border: 1px solid #dee2e6 !important; border-radius: 10px !important; }
.card-header { background-color: #f8f9fa !important; color: #2c3e50 !important; border-bottom: 1px solid #dee2e6 !important; }
.card-body { color: #2c3e50 !important; }

.form-control, .form-select { background-color: #ffffff !important; color: #2c3e50 !important;
  border: 1px solid #ced4da !important; }
.form-control:focus, .form-select:focus { border-color: #18BC9C !important;
  box-shadow: 0 0 0 0.15rem rgba(24,188,156,0.25) !important; }
.btn-primary { background-color: #2C3E50 !important; border-color: #2C3E50 !important; color: #fff !important; }
.btn-primary:hover { background-color: #1a252f !important; }
.btn-outline-primary { color: #2C3E50 !important; border-color: #2C3E50 !important; }
.btn-outline-secondary { color: #18BC9C !important; border-color: #18BC9C !important; }

.form-check-label { color: #495057 !important; }
.form-check-input { background-color: #fff !important; border-color: #ced4da !important; }
.form-check-input:checked { background-color: #18BC9C !important; border-color: #18BC9C !important; }
.shiny-options-group label { color: #495057 !important; }

.nav-tabs { border-bottom-color: #dee2e6 !important; }
.nav-tabs .nav-link { color: #6c757d !important; }
.nav-tabs .nav-link.active { background-color: #fff !important; color: #2C3E50 !important;
  border-color: #dee2e6 #dee2e6 #fff !important; }

table.dataTable { background-color: #fff !important; color: #2c3e50 !important; }
table.dataTable thead th { background-color: #f8f9fa !important; color: #495057 !important; }
table.dataTable tbody tr:hover { background-color: #e9ecef !important; }
.dataTables_info, .dataTables_paginate .paginate_button { color: #6c757d !important; }
.dataTables_paginate .paginate_button.current { background: #2C3E50 !important; color: #fff !important; }
div.dataTables_filter input { background-color: #fff !important; color: #2c3e50 !important; }

.bslib-value-box .value-box-title { color: #6c757d !important; }
.bslib-value-box .value-box-value { color: #2c3e50 !important; }

.alert-info { background-color: #d1ecf1 !important; color: #0c5460 !important; border-color: #bee5eb !important; }
.shiny-html-output p, .shiny-html-output li, .shiny-html-output span { color: #495057; }
.shiny-html-output h4, .shiny-html-output h5, .shiny-html-output h6 { color: #2c3e50; }
.shiny-html-output strong { color: #2c3e50; }

.hero-landing h1 { color: #2C3E50; }
.hero-landing h1 .morph-highlight { color: #18BC9C; }
.hero-landing .version-tag { color: #6c757d; }
.hero-landing .mantra { color: #495057; }
.hero-landing .motto { color: #e67e22; }
.hero-landing .trademark { color: #adb5bd; }
.hero-card { background: #ffffff !important; border: 1px solid #dee2e6 !important; }
.hero-card:hover { border-color: #18BC9C !important; }
.feature-icon-card i { color: #2C3E50; }
.feature-icon-card .feat-label { color: #6c757d; }

.modal-content { background-color: #fff !important; color: #2c3e50 !important; }
.small, small { color: #6c757d !important; }
.text-muted { color: #6c757d !important; }
p strong, li strong { color: #2c3e50 !important; }
.progress-bar { background-color: #18BC9C !important; }
"

# Store themes for toggling
dark_theme <- bs_theme(version = 5, bootswatch = "darkly",
                       primary = "#4FC3F7", secondary = "#81C784",
                       success = "#18BC9C", info = "#4FC3F7",
                       font_scale = 0.92, "navbar-bg" = "#16162a")
light_theme <- bs_theme(version = 5, bootswatch = "flatly",
                        primary = "#2C3E50", secondary = "#18BC9C",
                        success = "#18BC9C", info = "#3498DB",
                        font_scale = 0.92, "navbar-bg" = "#2C3E50")

# Small reusable hint shown under every editable figure
edit_hint_ui <- function() {
  tags$p(class = "small text-muted px-2 pt-1 mb-1",
    icon("pen"), " Double-click the figure to edit its title, axis labels, legend title, and group names.")
}

ui <- page_navbar(
  title = tags$span(
    icon("microscope", style = "color:#4FC3F7; margin-right:6px;"),
    tags$strong("CellMorphR", style = "color:#e0e0e8;"),
    tags$span(paste0(" v", APP_VERSION), style = "font-weight:300; opacity:0.7; color:#8888a8;")
  ),
  theme = dark_theme,
  header = tags$head(
    tags$style(id = "theme-css", HTML(dark_css)),
    tags$script(HTML("
      Shiny.addCustomMessageHandler('swap-css', function(css) {
        var el = document.getElementById('theme-css');
        if (el) el.innerHTML = css;
      });
    "))
  ),
  fillable = TRUE,

  # === TAB 1: DATA ============================================================
  nav_panel(
    title = "Data",
    icon  = icon("upload"),

    # --- Hero Landing Page (shown when no data loaded) ---
    conditionalPanel(
      condition = "!output.has_data",
      tags$div(class = "hero-landing",
        # Decorative background cells
        tags$div(class = "deco-cell", style = "top:8%; left:5%; width:110px;", HTML(svg_bacterium_rod)),
        tags$div(class = "deco-cell", style = "top:15%; right:8%; width:65px;", HTML(svg_phage_t4)),
        tags$div(class = "deco-cell", style = "bottom:18%; left:12%; width:45px;", HTML(svg_coccus)),
        tags$div(class = "deco-cell", style = "top:35%; right:15%; width:95px;", HTML(svg_bacterium_rod)),
        tags$div(class = "deco-cell", style = "bottom:10%; right:5%; width:42px;", HTML(svg_phage_small)),
        tags$div(class = "deco-cell", style = "top:55%; left:3%; width:50px;", HTML(svg_phage_small)),
        tags$div(class = "deco-cell", style = "bottom:30%; right:22%; width:40px;", HTML(svg_coccus)),
        tags$div(class = "deco-cell", style = "top:5%; left:45%; width:80px;", HTML(svg_phage_t4)),

        # Main hero content
        tags$div(class = "hero-content",
          tags$h1(
            tags$span("Cell"),
            tags$span("Morph", class = "morph-highlight"),
            tags$span("R")
          ),
          tags$p(class = "version-tag", paste0("v", APP_VERSION, " | Single-Cell Morphometry Analysis Suite")),
          tags$p(class = "mantra", icon("heart", style = "color:#e57373; font-size:0.9rem;"),
            " Built with passion for science and discovery"),
          tags$p(class = "motto", "YOU CAN JUST BUILD THINGS"),
          tags$hr(style = "border-color:#3a3a55; max-width:500px; margin:1.5rem auto;"),

          # --- Quick-start cards ---
          tags$div(style = "max-width:900px; margin:0 auto;",
            layout_columns(col_widths = c(4, 4, 4),
              tags$div(class = "hero-card card p-3",
                tags$h6(icon("upload", style = "color:#4FC3F7;"), " Upload Data"),
                fileInput("file", NULL, accept = c(".xlsx", ".xls", ".csv", ".tsv"),
                          placeholder = "Excel or CSV", width = "100%"),
                tags$p(class = "small", style = "color:#8888a8;",
                  "Upload .xlsx, .xls, .csv, or .tsv files with columns: ",
                  tags$code("Condition"), ", ", tags$code("Time_min"), ", ",
                  tags$code("Replicate"), " + measurement columns")
              ),
              tags$div(class = "hero-card card p-3",
                tags$h6(icon("play-circle", style = "color:#81C784;"), " Demo Data"),
                tags$p(class = "small", style = "color:#9898b8;",
                  "Try CellMorphR instantly with simulated phage infection data:
                   2 conditions, 4 timepoints, 3 replicates."),
                actionButton("load_demo", "Load Demo Data",
                  class = "btn-outline-secondary btn-sm w-100 mt-2", icon = icon("flask"))
              ),
              tags$div(class = "hero-card card p-3",
                tags$h6(icon("folder-open", style = "color:#FFB74D;"), " Resume Session"),
                tags$p(class = "small", style = "color:#9898b8;",
                  "Resume a previous analysis from a saved .rds session file."),
                fileInput("load_session", NULL, accept = ".rds",
                          placeholder = "Load .rds session", width = "100%")
              )
            )
          ),

          # --- Feature highlights ---
          tags$div(style = "max-width:800px; margin:1.5rem auto 0;",
            layout_columns(col_widths = c(2, 2, 2, 2, 2, 2),
              tags$div(class = "feature-icon-card", icon("chart-area"),
                tags$div(class = "feat-label", "10 Analysis", tags$br(), "Tabs")),
              tags$div(class = "feature-icon-card", icon("image"),
                tags$div(class = "feat-label", "Publication", tags$br(), "Figures")),
              tags$div(class = "feature-icon-card", icon("shield-halved"),
                tags$div(class = "feat-label", "Anti-Pseudo", tags$br(), "replication")),
              tags$div(class = "feature-icon-card", icon("chart-line"),
                tags$div(class = "feat-label", "Spider &", tags$br(), "Radar Plots")),
              tags$div(class = "feature-icon-card", icon("shapes"),
                tags$div(class = "feat-label", "Phenotype", tags$br(), "Discovery")),
              tags$div(class = "feature-icon-card", icon("download"),
                tags$div(class = "feat-label", "Export", tags$br(), "Everything"))
            )
          ),
          tags$p(class = "trademark",
            "\u00a9 CellMorphR by Michael Baffour Awuah | github.com/mbaffour/CellMorphR")
        )
      )
    ),

    # --- Data Panel (shown after data is loaded) ---
    conditionalPanel(
      condition = "output.has_data",
      layout_sidebar(
        sidebar = sidebar(
          width = 320, title = "Configure",
          tags$p(tags$strong("Required columns:"), class = "mb-1"),
          tags$ul(class = "small",
            tags$li(tags$code("Condition"), " -- e.g., Uninfected, +N4"),
            tags$li(tags$code("Time_min"), " -- numeric (30, 60, 90, 120)"),
            tags$li(tags$code("Replicate"), " -- biological replicate ID (1, 2, 3)")
          ),
          tags$p(class = "small", style = "color:#8888a8;", "Plus numeric measurement columns."),
          tags$hr(),
          uiOutput("measure_selector"),
          uiOutput("condition_order_ui"),
          tags$hr(),
          tags$details(
            tags$summary(tags$strong("Rename samples (labels on graphs)", style = "cursor:pointer;")),
            tags$div(class = "pt-2",
              tags$p(class = "small text-muted",
                "Edit how each sample is written everywhere -- legends, axes, tables, exports.
                 Names must stay unique."),
              tags$p(tags$strong("Conditions"), class = "small mb-1"),
              actionButton("rename_remove_underscores", "Replace all underscores with spaces",
                class = "btn-outline-secondary btn-sm w-100 mb-2", icon = icon("broom")),
              uiOutput("condition_rename_ui"),
              actionButton("rename_conditions_apply", "Apply new names",
                class = "btn-outline-primary btn-sm w-100", icon = icon("pen")),
              tags$hr(),
              tags$p(tags$strong("Timepoints"), class = "small mb-1"),
              tags$p(class = "small text-muted",
                "Rename how timepoints appear (e.g. \"Uninduced\", \"1 h\"). The underlying
                 numeric time is unchanged, so ordering, trajectories and the trend test still work."),
              uiOutput("time_rename_ui"),
              actionButton("rename_times_apply", "Apply timepoint labels",
                class = "btn-outline-primary btn-sm w-100", icon = icon("clock"))
            )
          ),
          tags$hr(),
          tags$details(
            tags$summary(tags$strong("Split Condition into factors", style = "cursor:pointer;")),
            tags$div(class = "pt-2",
              tags$p(class = "small text-muted",
                "If your Condition labels encode several factors (e.g., 'pSJExD_gp75_100nM' = vector + dose),
                 split them into separate columns. The new columns become available for statistical grouping
                 and the dose-response Trend test."),
              textInput("split_delim", "Delimiter", "_"),
              verbatimTextOutput("split_preview"),
              actionButton("split_condition_apply", "Create columns",
                class = "btn-outline-primary btn-sm w-100", icon = icon("scissors"))
            )
          ),
          tags$hr(),
          selectInput("data_transform", "Transform Numeric Columns",
                      choices = c("None" = "none", "Log (x+1)" = "log", "Square Root" = "sqrt"),
                      selected = "none"),
          tags$p(class = "small", style = "color:#8888a8;", "Applied globally. Useful for right-skewed data."),
          tags$hr(),
          uiOutput("color_pickers"),
          tags$hr(),
          tags$p(tags$strong("Session"), class = "mb-1"),
          downloadButton("save_session", "Save Session (.rds)", class = "btn-outline-primary btn-sm w-100"),
          tags$p(class = "small mt-1", style = "color:#8888a8;", "Save data + settings to resume later."),
          tags$hr(),
          fileInput("file_replace", "Upload New Data",
            accept = c(".xlsx", ".xls", ".csv", ".tsv"), placeholder = "Replace current data"),
          actionButton("load_demo_again", "Load Demo Data",
            class = "btn-outline-secondary btn-sm w-100", icon = icon("flask"))
        ),
        card(card_header("Data Preview"), DTOutput("data_table")),
        layout_columns(col_widths = c(4, 4, 4),
          value_box(title = "Total Cells", value = textOutput("n_cells_txt"),
                    showcase = icon("microscope"), theme = "primary"),
          value_box(title = "Conditions x Time Points", value = textOutput("n_groups_txt"),
                    showcase = icon("layer-group"), theme = "secondary"),
          value_box(title = "Biological Replicates", value = textOutput("n_reps_txt"),
                    showcase = icon("flask"), theme = "info")
        ),
        card(card_header("Feature Diagnostics"), DTOutput("feature_diagnostics_table"))
      )
    )
  ),

  # === TAB 2: DISTRIBUTIONS (with sub-tabs) ===================================
  nav_panel(
    title = "Distributions",
    icon  = icon("chart-area"),
    layout_sidebar(
      sidebar = sidebar(
        width = 300, title = "Distribution Options",
        # ── Subset / Filter ──────────────────────────────────────────
        tags$p(tags$strong("Show Samples"), class = "mb-1"),
        uiOutput("dist_sample_filter_ui"),
        tags$hr(),
        conditionalPanel("output.has_timecourse",
          tags$p(tags$strong("Show Timepoints"), class = "mb-1"),
          uiOutput("dist_time_filter_ui"),
          tags$hr()),
        # ── X-axis grouping ──────────────────────────────────────────
        conditionalPanel("output.has_timecourse",
        radioButtons("dist_grouping_mode", "Arrange X axis by",
          choices = c(
            "Timepoint — conditions dodge"  = "by_time",
            "Condition — timepoints dodge"  = "by_condition"
          ),
          selected = "by_time", inline = FALSE),
        tags$p(class="small text-muted", style="margin-top:-6px;margin-bottom:6px;font-size:10px;",
          "Pick which variable sits on the X axis. The other variable becomes the fill/dodge colour.")),
        tags$hr(),
        radioButtons("dist_type", "Plot Type",
          choices = c("Violin + Medians" = "violin", "Strip (Jitter)" = "strip",
                      "Sina (Density-Jitter)" = "sina", "Box + Strip" = "boxstrip",
                      "Dot Plot (Replicates Only)" = "dotplot",
                      "Ridge Plot" = "ridge", "ECDF" = "ecdf",
                      "Histogram" = "histogram",
                      "Stacked Bar (quartiles)" = "stacked_bar",
                      "Pie Chart" = "pie", "Donut Chart" = "donut"),
          selected = "violin"),
        conditionalPanel("input.dist_type == 'histogram'",
          numericInput("dist_hist_bins", "Histogram bins", 40, 5, 200, 5)),
        tags$hr(),
        checkboxInput("show_rep_medians", "Show replicate medians", TRUE),
        checkboxInput("show_cell_points", "Show individual cells", FALSE),
        conditionalPanel("input.show_cell_points",
          sliderInput("point_alpha", "Point opacity", 0.02, 0.5, 0.08, 0.02)),
        tags$hr(),
        tags$p(tags$strong("SuperPlot Options"), class = "mb-1"),
        radioButtons("dist_rep_encode", "Distinguish replicates by",
          choices = c("Nothing (identical dots)" = "none",
                      "Point shape per replicate" = "shape"),
          selected = "none"),
        conditionalPanel("input.dist_rep_encode == 'shape'",
          tags$p(class = "small text-muted mb-1", "Shape for each replicate (also shown in the legend):"),
          uiOutput("dist_rep_shape_ui")),
        conditionalPanel("input.show_cell_points",
          selectInput("dist_cell_shape", "Individual cell point shape",
                      choices = SHAPE_CHOICES, selected = 19)),
        checkboxInput("dist_cells_by_rep", "Color individual cells by replicate", FALSE),
        checkboxInput("dist_show_n", "Show n (cells / replicates) under each group", FALSE),
        tags$p(class = "small text-muted", style = "font-size:10px; margin-top:-4px;",
          "SuperPlot style (Lord et al. 2020): show cell-level spread AND replicate-level agreement without needing p-values."),
        tags$hr(),
        tags$p(tags$strong("Mean/Median Markers"), class = "mb-1"),
        checkboxInput("show_group_mean", "Show group mean (horizontal bar)", TRUE),
        radioButtons("dist_stat_type", "Statistic for annotations",
          choices = c("Mean" = "mean", "Median" = "median"), selected = "mean", inline = TRUE),
        tags$hr(),
        tags$p(tags$strong("Significance Brackets"), class = "mb-1"),
        checkboxInput("show_sig", "Show significance brackets", FALSE),
        conditionalPanel("input.show_sig",
          radioButtons("sig_compare_mode", "Comparison type",
            choices = c("Each treatment vs control" = "vs_control",
                        "All pairwise (every pair)" = "all_pairwise"),
            selected = "vs_control"),
          radioButtons("sig_label_style", "Label style",
            choices = c("Stars (*, **, ***)" = "stars", "Stars + p-value" = "stars_p",
                        "p-value only" = "pval"), selected = "stars"),
          checkboxInput("show_ns", "Show ns brackets", FALSE),
          radioButtons("sig_stat_unit", "Test using replicate-level:",
            choices = c("Medians" = "median_val", "Means" = "mean_val"),
            selected = "median_val", inline = TRUE),
          tags$div(class = "alert alert-info py-1 px-2 mt-2 small",
            tags$strong("vs control:"), "Compares each treatment to the reference condition at each timepoint.",
            tags$br(), tags$strong("All pairwise:"), "Compares every pair of conditions at each timepoint.",
            tags$br(), tags$strong("Test & correction:"), "Follow the Statistics tab settings (default: Welch's t-test, Holm). Always computed on per-replicate summaries (N = biological replicates), never raw cells.",
            tags$br(), tags$strong("Thresholds:"), "* p<0.05, ** p<0.01, *** p<0.001.")
        ),
        tags$hr(),
        tags$p(tags$strong("Violin Appearance"), class = "mb-1"),
        sliderInput("violin_alpha",    "Violin fill opacity",   0.1, 1.0, 0.70, 0.05),
        sliderInput("violin_width_scale","Violin width scale",  0.3, 2.0, 1.0, 0.05),
        checkboxInput("violin_trim",   "Trim violin tails",     FALSE),
        sliderInput("box_notch_w",     "Box width (violin+box)",0.05, 0.5, 0.18, 0.01),
        checkboxInput("violin_outline","Show violin outline",   TRUE),
        sliderInput("violin_line_w",   "Outline thickness",     0.1, 2.0, 0.3, 0.1),
        tags$hr(),
        tags$p(tags$strong("Point Appearance"), class = "mb-1"),
        sliderInput("rep_dot_size",    "Replicate dot size",    0.5, 8.0, 3.5, 0.25),
        colourpicker::colourInput("rep_dot_color", "Replicate dot border", "#000000", showColour="background"),
        sliderInput("rep_dot_alpha",   "Replicate dot opacity", 0.1, 1.0, 0.9, 0.05),
        tags$hr(),
        tags$p(tags$strong("Axis Customization"), class = "mb-1"),
        checkboxInput("dist_log_y",    "Log Y axis",            FALSE),
        selectInput("dist_x_angle", "X label angle", choices = c("0" = 0, "30" = 30, "45" = 45, "60" = 60, "90" = 90), selected = 45),
        numericInput("dist_label_wrap", "Wrap X labels after characters", value = 18, min = 6, max = 60, step = 2),
        numericInput("dist_bottom_margin", "Bottom margin", value = 18, min = 5, max = 80, step = 2),
        checkboxInput("dist_sample_x_mode", "Show combined sample names on X axis and hide legend", FALSE),
        numericInput("dist_y_min",     "Y axis min (blank=auto)", value=NA, step=0.1),
        numericInput("dist_y_max",     "Y axis max (blank=auto)", value=NA, step=0.1),
        tags$hr(),
        figure_controls_ui("dist", 8, 5.5),
        tags$hr(),
        downloadButton("dl_dist", "Download Figure", class = "btn-primary w-100")
      ),
      navset_card_tab(
        nav_panel("Plot",
          plotOutput("dist_plot", height = "600px", dblclick = "dist_plot_dblclick"),
          edit_hint_ui()
        ),
        nav_panel("About",
          tags$div(class = "p-3 small",
            tags$h6("Understanding Distribution Plots"),
            tags$p(tags$strong("SuperPlot options"), " let you encode biological replicates in the figure itself:
                    per-replicate point shapes on the summary dots, or coloring individual cells by replicate
                    (Lord et al. 2020, J Cell Biol). This shows sample-to-sample variation honestly without
                    requiring any hypothesis test."),
            tags$p(tags$strong("Violin plots"), " show kernel density estimates of single-cell
                    measurements. Width = density. Best for large datasets (200+ cells)."),
            tags$p(tags$strong("Strip/Sina"), " show every cell as a point. Sina spreads points
                    proportional to density."),
            tags$p(tags$strong("Overlaid replicate medians"), " (colored dots) are your true N.
                    Each dot = median of one biological replicate."),
            tags$p(tags$strong("Group mean bar"), " shows the mean of replicate medians --
                    the honest summary of your data."),
            tags$p(tags$strong("Significance brackets"), " use Welch's t-test on replicate-level
                    summaries (not raw cells). This avoids pseudoreplication. Holm correction
                    applied across all comparisons."),
            tags$p(tags$strong("vs control:"), " compares each treatment to the reference condition
                    at each timepoint. Best when you have a clear control group."),
            tags$p(tags$strong("All pairwise:"), " compares every pair of conditions at each timepoint.
                    Best when you have multiple treatments and want to compare them all."),
            tags$p(tags$strong("Mean vs Median:"), " Median is more robust to outliers and skewed data
                    (common in cell biology). Mean is appropriate for symmetric distributions.")
          )
        )
      )
    )
  ),

  # === TAB 3: TRAJECTORIES ====================================================
  nav_panel(
    title = "Trajectories",
    icon  = icon("chart-line"),
    layout_sidebar(
      sidebar = sidebar(
        width = 300, title = "Trajectory Options",
        radioButtons("traj_summary", "Summary Statistic",
                     choices = c("Median" = "median_val", "Mean" = "mean_val"), selected = "median_val"),
        checkboxInput("traj_smooth", "Add trend line (loess)", TRUE),
        checkboxInput("traj_error", "Show SEM ribbon", TRUE),
        checkboxInput("traj_connect_reps", "Connect individual replicates", TRUE),
        radioButtons("traj_style", "Style",
                     choices = c("Points + Lines" = "pointline", "Points + Ribbon" = "ribbon"),
                     selected = "pointline"),
        sliderInput("traj_point_size", "Point size", 0.5, 8, 3, 0.25),
        sliderInput("traj_line_width", "Trend line width", 0.3, 3, 1.2, 0.1),
        tags$hr(),
        figure_controls_ui("traj", 7, 5),
        tags$hr(),
        downloadButton("dl_traj", "Download Figure", class = "btn-primary w-100")
      ),
      navset_card_tab(
        nav_panel("Plot",
          plotOutput("traj_plot", height = "550px", dblclick = "traj_plot_dblclick"),
          edit_hint_ui()),
        nav_panel("About", tags$div(class = "p-3 small",
          tags$h6("Understanding Trajectory Plots"),
          tags$p(tags$strong("Each point"), " is the median (or mean) of one biological replicate."),
          tags$p(tags$strong("Diverging lines"), " between conditions = Condition x Time interaction."),
          tags$p(tags$strong("SEM ribbons"), " computed from replicate-level summaries (N=3).")
        ))
      )
    )
  ),

  # === TAB 4: MORPHOSPACE (PCA) ================================================
  nav_panel(
    title = "Morphospace",
    icon  = icon("cube"),
    layout_sidebar(
      sidebar = sidebar(
        width = 300, title = "PCA Options",
        checkboxInput("pca_scale", "Scale features (recommended)", TRUE),
        selectInput("pca_color_by", "Color by",
          choices = c("Condition" = "Condition", "Time" = "Time_min", "Condition x Time" = "cond_time"),
          selected = "Condition"),
        uiOutput("pca_exclude_features_ui"),
        tags$hr(),
        checkboxInput("pca_show_density", "Show density contours", FALSE),
        checkboxInput("pca_show_centroids", "Show group centroids", TRUE),
        checkboxInput("pca_show_trajectories", "Show centroid trajectories over time", FALSE),
        checkboxInput("pca_show_loadings", "Show loading vectors", FALSE),
        conditionalPanel("input.pca_show_loadings",
          sliderInput("pca_loading_scale", "Loading vector scale", 1, 10, 3, 0.5)),
        sliderInput("pca_point_alpha", "Point opacity", 0.05, 1.0, 0.2, 0.05),
        sliderInput("pca_point_size", "Point size", 0.3, 6, 1, 0.1),
        tags$hr(),
        figure_controls_ui("pca", 8, 7),
        tags$hr(),
        downloadButton("dl_pca", "Download Morphospace", class = "btn-primary w-100"),
        downloadButton("dl_pca_traj", "Download Trajectory", class = "btn-outline-primary w-100 mt-2")
      ),
      navset_card_tab(
        nav_panel("Morphospace",
          plotOutput("pca_scatter", height = "550px", dblclick = "pca_scatter_dblclick"),
          edit_hint_ui()),
        nav_panel("Loadings",
          plotOutput("pca_loadings_bar", height = "500px"),
          downloadButton("dl_pca_loadings", "Download Loadings Figure", class = "btn-sm btn-outline-primary m-2")),
        nav_panel("Trajectories", plotOutput("pca_trajectory", height = "500px")),
        nav_panel("Diversity",
          layout_columns(col_widths = c(12),
            card(card_header("Morphological Diversity Metrics"), DTOutput("pca_diversity_table")),
            card(card_header("PC Variance by Group"), plotOutput("pca_diversity_plot", height = "350px"),
                 downloadButton("dl_pca_diversity", "Download Diversity Figure", class = "btn-sm btn-outline-primary m-2"))
          )
        ),
        nav_panel("Scree",
          layout_columns(col_widths = c(6, 6),
            card(card_header("Variance Explained"), plotOutput("pca_scree", height = "300px"),
                 downloadButton("dl_pca_scree", "Download Scree Figure", class = "btn-sm btn-outline-primary m-2")),
            card(card_header("Loadings Table"), DTOutput("pca_loadings_table"))
          )
        ),
        nav_panel("About", tags$div(class = "p-3 small",
          tags$h6("About Morphospace Analysis"),
          tags$p(tags$strong("PCA"), " reduces all features into 2 axes. Each point = one cell."),
          tags$p(tags$strong("Density contours"), " show where cells concentrate in morphospace."),
          tags$p(tags$strong("Loadings bar plot"), " shows which features drive PC1 and PC2 most."),
          tags$p(tags$strong("Trajectories"), " connect group centroids over time, showing morphological trajectory of infection."),
          tags$p(tags$strong("Diversity metrics"), " quantify morphological heterogeneity:"),
          tags$ul(
            tags$li(tags$strong("Convex hull area:"), " Area occupied in PC1-PC2 space. Large = diverse morphology."),
            tags$li(tags$strong("PC1/PC2 variance:"), " Spread along each principal component."),
            tags$li(tags$strong("Mean distance to centroid:"), " Average morphological deviation from group center.")
          ),
          tags$p(tags$strong("Feature exclusion:"), " Use the dropdown above to exclude columns (e.g., ID) from PCA.")
        ))
      )
    )
  ),

  # === TAB 5: PHENOTYPES (CLUSTERING) ==========================================
  nav_panel(
    title = "Phenotypes",
    icon  = icon("shapes"),
    layout_sidebar(
      sidebar = sidebar(
        width = 300, title = "Clustering Options",
        radioButtons("clust_method", "Method",
                     choices = c("K-means" = "kmeans", "Hierarchical (Ward)" = "hclust"), selected = "kmeans"),
        checkboxInput("clust_on_pca", "Cluster on PCA (recommended)", TRUE),
        conditionalPanel("input.clust_on_pca",
          sliderInput("clust_n_pcs", "Number of PCs", 2, 10, 3, 1)),
        tags$hr(),
        radioButtons("clust_k_method", "Number of clusters",
                     choices = c("Auto-detect (silhouette)" = "auto", "Manual" = "manual"), selected = "auto"),
        conditionalPanel("input.clust_k_method == 'manual'",
          sliderInput("clust_k", "k", 2, 10, 3, 1)),
        actionButton("run_clustering", "Run Clustering", class = "btn-primary w-100 mt-2", icon = icon("play")),
        tags$hr(),
        figure_controls_ui("clust", 10, 5),
        tags$hr(),
        downloadButton("dl_clust", "Download Figure", class = "btn-primary w-100")
      ),
      navset_card_tab(
        nav_panel("Clusters",
          layout_columns(col_widths = c(6, 6),
            card(card_header("Optimal k"), plotOutput("clust_k_plot", height = "300px"),
                 downloadButton("dl_clust_k", "Download", class = "btn-sm btn-outline-primary m-2")),
            card(card_header("Morphospace by Phenotype"), plotOutput("clust_pca_plot", height = "300px"),
                 downloadButton("dl_clust_pca", "Download", class = "btn-sm btn-outline-primary m-2"))
          )
        ),
        nav_panel("Composition",
          plotOutput("clust_composition", height = "450px", dblclick = "clust_composition_dblclick"),
          edit_hint_ui()),
        nav_panel("Composition Table", DTOutput("clust_comp_table")),
        nav_panel("Characterization", DTOutput("clust_char_table")),
        nav_panel("About", tags$div(class = "p-3",
          tags$h5("About Phenotype Discovery"),
          tags$p("Phenotype discovery uses unsupervised clustering to identify distinct morphological
                  subpopulations ('phenotypes') in your data. Cells are grouped purely on their measured
                  features, without prior knowledge of what morphological states might exist."),
          tags$h6("Clustering Methods"),
          tags$ul(
            tags$li(tags$strong("K-means:"), " Partitions cells into k groups by minimizing within-cluster
                     variance. Fast and widely used. Best when clusters are roughly spherical and equally sized."),
            tags$li(tags$strong("Hierarchical (Ward):"), " Agglomerative approach that merges nearest clusters
                     using Ward's minimum variance criterion. Produces a dendrogram and tends to find compact,
                     roughly equal-sized clusters.")
          ),
          tags$h6("Selecting k (Number of Clusters)"),
          tags$p("In 'Auto-detect' mode, silhouette analysis evaluates k = 2 through 10. For each k,
                  the mean silhouette width measures how well each cell fits its assigned cluster vs. the
                  nearest neighboring cluster. Values near 1 indicate strong cluster structure; near 0
                  suggests overlapping clusters. The k with the highest silhouette score is selected."),
          tags$h6("PCA Space"),
          tags$p("By default, clustering operates on PCA-transformed data rather than raw features.
                  This reduces noise, removes correlated features, and often improves cluster separation.
                  The 'Morphospace by Phenotype' scatter shows clusters projected onto PC1 vs PC2."),
          tags$h6("Interpreting Results"),
          tags$ul(
            tags$li(tags$strong("Composition Plot:"), " Stacked bars showing the % of cells in each
                     phenotype per Condition × Time. Shifts in composition over time indicate
                     treatment-specific morphological responses."),
            tags$li(tags$strong("Composition Table:"), " Exact counts and percentages per
                     Condition × Time × Phenotype for precise comparisons."),
            tags$li(tags$strong("Characterization:"), " Summary statistics (mean, median, SD, IQR, CV%)
                     of each measured feature per phenotype. Helps interpret what each phenotype
                     represents biologically (e.g., P1 = small round cells, P2 = large elongated cells).")
          )
        ))
      )
    )
  ),

  # === TAB 6: STATISTICS ======================================================
  nav_panel(
    title = "Statistics",
    icon  = icon("square-root-variable"),
    layout_sidebar(
      sidebar = sidebar(
        width = 340, title = "Flexible Statistics Setup",
        uiOutput("stat_measure_ui"),
        tags$hr(),
        tags$p(tags$strong("Build statistical groups from"), class = "mb-1"),
        uiOutput("stat_group_cols_ui"),
        tags$p(class = "small text-muted", "Choose the metadata columns that define the groups you want to compare. This can be Condition only, Condition + Time_min, Plasmid + Inducer, etc."),
        tags$hr(),
        radioButtons("stat_summary", "Statistical unit",
                     choices = c("Replicate median" = "median_val", "Replicate mean" = "mean_val"), selected = "median_val"),
        radioButtons("stat_compare_mode", "Comparison mode",
          choices = c("Each group vs selected reference" = "vs_reference",
                      "All pairwise groups" = "all_pairwise",
                      "Pairwise within each time point" = "within_time",
                      "Pairwise within each condition" = "within_condition"), selected = "within_time"),
        uiOutput("stat_reference_ui"),
        selectInput("stat_test", "Pairwise test",
          choices = c("Welch's t-test" = "welch", "Wilcoxon rank-sum" = "wilcox"), selected = "welch"),
        selectInput("stat_p_adjust", "P-value correction",
          choices = c("Holm" = "holm", "Benjamini-Hochberg" = "BH", "Bonferroni" = "bonferroni", "None" = "none"), selected = "holm"),
        actionButton("suggest_stats", "Suggest Test & Settings",
          class = "btn-outline-secondary w-100 mt-2", icon = icon("magic")),
        tags$p(class = "small text-muted mt-1", "This inspects your selected groups, comparisons, and biological replicate counts, then explains the recommended test before you run it."),
        tags$hr(),
        tags$details(
          tags$summary(tags$strong("Optional global model"), style = "cursor:pointer;"),
          radioButtons("stat_method", "Global model",
            choices = c("Repeated-measures ANOVA (dish as error stratum)" = "anova",
                        "Linear mixed model on cells" = "lmm",
                        "Skip global model" = "skip"), selected = "anova"),
          tags$p(class = "small text-muted",
            "Both models treat one dish measured at k timepoints as one experimental unit. ",
            "The repeated-measures ANOVA runs on the replicate summaries; the mixed model runs on cells, ",
            "with a random effect per dish and per dish-and-timepoint. Both use the Data tab's transform, ",
            "the same one the pairwise tests use.")
        ),
        tags$hr(),
        actionButton("run_stats", "Run Flexible Statistics", class = "btn-primary w-100", icon = icon("play")),
        tags$hr(),
        downloadButton("dl_stats", "Download Results", class = "btn-outline-primary w-100"),
        downloadButton("dl_stats_csv", "Download Pairwise CSV", class = "btn-outline-secondary w-100 mt-2")
      ),
      navset_card_tab(
        nav_panel("Group Preview", DTOutput("stat_group_preview_table")),
        nav_panel("Comparison Preview", DTOutput("comparison_preview_table")),
        nav_panel("Test Suggestion", verbatimTextOutput("test_suggestion_output"), DTOutput("test_suggestion_group_table")),
        nav_panel("Descriptive", DTOutput("stat_descriptive_table")),
        nav_panel("Results", verbatimTextOutput("stat_output")),
        nav_panel("Pairwise", DTOutput("pairwise_table")),
        nav_panel("Trend",
          tags$div(class = "p-2",
            layout_columns(col_widths = c(4, 4, 4),
              uiOutput("trend_predictor_ui"),
              uiOutput("trend_strata_ui"),
              tags$div(
                actionButton("run_trend", "Run Trend Test", class = "btn-primary w-100 mt-4", icon = icon("chart-line"))
              )
            ),
            tags$p(class = "small text-muted",
              "Fits a linear model of replicate summaries against a numeric predictor (e.g., inducer dose) within each stratum. ",
              "This asks 'does the measurement change along the gradient?' -- usually far more powerful than dozens of pairwise tests, ",
              "especially with few biological replicates. Use 'Split Condition' on the Data tab to extract a numeric dose column first if needed."),
            plotOutput("trend_plot", height = "380px"),
            layout_columns(col_widths = c(6, 6),
              downloadButton("dl_trend_plot", "Download Trend Figure", class = "btn-outline-primary w-100"),
              downloadButton("dl_trend_csv", "Download Trend Table (CSV)", class = "btn-outline-secondary w-100")
            ),
            DTOutput("trend_table")
          )
        ),
        nav_panel("EMM Plot",
          plotOutput("emm_plot", height = "450px"),
          tags$div(class = "p-2",
            layout_columns(col_widths = c(4, 8),
              selectInput("emm_format", "Format",
                choices = c("PDF" = "pdf", "PNG" = "png", "SVG" = "svg", "PPTX (editable)" = "pptx"), selected = "pdf"),
              downloadButton("dl_emm_plot", "Download EMM Figure", class = "btn-outline-primary w-100 mt-4")
            )
          )
        ),
        nav_panel("About", tags$div(class = "p-3 small",
          tags$h6("Understanding the Statistical Analysis"),
          tags$p(tags$strong("Step 1:"), " Collapse cells to replicate-level summaries (median or mean).
                 P-values are computed on these summaries -- never on raw single cells -- so thousands of cells
                 cannot masquerade as thousands of independent observations (anti-pseudoreplication;
                 see Lord et al. 2020, J Cell Biol, 'SuperPlots')."),
          tags$p(tags$strong("Step 2:"), " Choose the metadata columns that define the comparison groups."),
          tags$p(tags$strong("Reference group:"), " Pick it manually. The app no longer silently uses the earliest time point or lowest inducer as the reference."),
          tags$p(tags$strong("Step 3:"), " Pairwise tests compare replicate summaries, with adjustable p-value correction.
                 The significance brackets on the Distributions and Cell Metrics figures use the test and correction selected here,
                 so figures and tables always agree."),
          tags$p(tags$strong("Few replicates?"), " With only 2 biological replicates per group, the Wilcoxon rank-sum test
                 mathematically cannot reach p < 0.05 (minimum attainable p = 0.33), and Welch's t-test has ~1 degree of freedom.
                 The app warns you when this applies. In that situation prefer the linear mixed model (global model below),
                 a dose-response trend test (Trend tab), or simply present the SuperPlot with replicate-level points and no p-values."),
          tags$p(tags$strong("Trend test:"), " For gradient/dose designs, the Trend tab fits replicate summaries against a
                 numeric predictor (slope, 95% CI, p, R²) within each stratum -- one focused question instead of many pairwise ones."),
          tags$p(tags$strong("Optional global model:"), " ANOVA/LMM is retained for simple Condition × Time_min workflows.")
        ))
      )
    )
  ),

  # === TAB 7: EFFECT SIZES ====================================================
  nav_panel(
    title = "Effect Sizes",
    icon  = icon("arrows-left-right"),
    layout_sidebar(
      sidebar = sidebar(
        width = 300, title = "Options",
        radioButtons("es_summary", "Summary Statistic",
          choices = c("Median" = "median_val", "Mean" = "mean_val"), selected = "median_val"),
        radioButtons("es_type", "Metric",
          choices = c("Absolute difference" = "abs_diff", "Percent change" = "pct_change",
                      "Hedges' g (95% CI)" = "hedges_g"), selected = "pct_change"),
        tags$hr(),
        checkboxInput("es_show_control", "Show control trajectory (baseline)", FALSE),
        radioButtons("es_plot_style", "Plot style",
          choices = c("Lines + Points" = "linepoint", "Points only" = "points",
                      "Bubble (size = magnitude)" = "bubble"),
          selected = "linepoint"),
        tags$hr(),
        figure_controls_ui("es", 7, 4.5),
        tags$hr(),
        downloadButton("dl_es", "Download Figure", class = "btn-primary w-100")
      ),
      navset_card_tab(
        nav_panel("Plot",
          plotOutput("es_plot", height = "450px", dblclick = "es_plot_dblclick"),
          edit_hint_ui()),
        nav_panel("Table", DTOutput("es_table")),
        nav_panel("About", tags$div(class = "p-3 small",
          tags$h6("Understanding Effect Sizes"),
          tags$p("P-values show if an effect exists; effect sizes show ", tags$strong("how big"), "."),
          tags$p(tags$strong("Percent change"), " is most intuitive: 'Cells are 25% larger at 120 min.'"),
          tags$p(tags$strong("Hedges' g"), " standardizes by variability: the difference in means divided by the pooled SD of the biological replicates, with Hedges' small-sample correction. Cohen's conventional benchmarks put |g| > 0.8 at \"large\", but with three dishes per group the 95% interval shown alongside each point is usually several units wide -- read the interval, not the point.")
        ))
      )
    )
  ),

  # === TAB 8: MULTI-FEATURE SUMMARY ===========================================
  nav_panel(
    title = "Multi-Feature",
    icon  = icon("table-cells"),
    layout_sidebar(
      sidebar = sidebar(
        width = 300, title = "Multi-Feature Analysis",
        radioButtons("mf_summary", "Summarize by",
          choices = c("Median" = "median_val", "Mean" = "mean_val"), selected = "median_val"),
        numericInput("mf_alpha", "Significance threshold", 0.05, 0.001, 0.10, 0.01),
        tags$hr(),
        tags$p(tags$strong("Stat Heatmap Appearance"), class = "mb-1"),
        selectInput("mf_palette", "Color palette",
          choices = c("White-Yellow-Red" = "wyr", "White-Blue" = "wb",
                      "Viridis" = "viridis", "White-Green-Black" = "wgb"),
          selected = "wyr"),
        checkboxInput("mf_show_stars", "Show significance stars", TRUE),
        checkboxInput("mf_show_pvals", "Show p-values in cells", FALSE),
        uiOutput("mf_feature_selector"),
        tags$hr(),
        actionButton("run_multifeature", "Run All Features", class = "btn-primary w-100", icon = icon("play")),
        tags$p(class = "small text-muted mt-2", "Repeated-measures ANOVA (dish as error stratum) across all features simultaneously."),
        tags$hr(),
        tags$p(tags$strong("Raw Heatmap Settings"), class = "mb-1"),
        radioButtons("mf_raw_stat", "Statistic",
          choices = c("Mean" = "mean", "Median" = "median"), selected = "mean"),
        checkboxInput("mf_raw_zscore", "Z-score scale (per feature)", FALSE),
        selectInput("mf_raw_palette", "Color palette",
          choices = c("Blue-White-Red" = "bwr", "Green-White-Red" = "gwr",
                      "Viridis" = "viridis", "Purple-White-Orange" = "pwo"),
          selected = "bwr"),
        checkboxInput("mf_raw_show_values", "Show values in cells", TRUE),
        tags$hr(),
        tags$p(tags$strong("Spider/Radar Settings"), class = "mb-1"),
        radioButtons("mf_radar_stat", "Statistic",
          choices = c("Mean" = "mean", "Median" = "median"), selected = "mean"),
        uiOutput("mf_radar_time_ui"),
        checkboxInput("mf_radar_points", "Show data points", TRUE),
        checkboxInput("mf_radar_grid", "Show grid rings", TRUE),
        sliderInput("mf_radar_alpha", "Fill opacity", 0, 1, 0.25, 0.05),
        tags$hr(),
        figure_controls_ui("mf", 8, 6),
        tags$hr(),
        downloadButton("dl_mf_heatmap", "Download Heatmap", class = "btn-primary w-100"),
        downloadButton("dl_mf_raw_heatmap", "Download Raw Heatmap", class = "btn-outline-primary w-100 mt-2"),
        downloadButton("dl_mf_radar", "Download Radar Plot", class = "btn-outline-primary w-100 mt-2"),
        downloadButton("dl_mf_table", "Download Results (CSV)", class = "btn-outline-primary w-100 mt-2")
      ),
      navset_card_tab(
        nav_panel("Stat Heatmap",
          plotOutput("mf_heatmap", height = "500px", dblclick = "mf_heatmap_dblclick"),
          edit_hint_ui()),
        nav_panel("Raw Heatmap", plotOutput("mf_raw_heatmap", height = "500px")),
        nav_panel("Spider/Radar", plotOutput("mf_radar_plot", height = "550px")),
        nav_panel("Importance",
          plotOutput("mf_importance", height = "400px"),
          downloadButton("dl_mf_importance", "Download Importance Figure", class = "btn-sm btn-outline-primary m-2")),
        nav_panel("Divergence", DTOutput("mf_divergence_table")),
        nav_panel("Full Table", DTOutput("mf_results_table"))
      )
    )
  ),

  # === TAB 9: CELL METRICS ====================================================
  nav_panel(
    title = "Cell Metrics",
    icon  = icon("chart-line"),
    layout_sidebar(
      sidebar = sidebar(
        width = 320, title = "Cell Metrics Explorer", open = TRUE,
        # --- Measurement selector ---
        uiOutput("cm_measure_ui"),
        tags$hr(),
        tags$p(tags$strong("Chart Type"), class = "mb-1"),
        radioButtons("cm_chart_type", NULL,
          choices = c(
            "Bar chart" = "bar", "Line chart" = "line",
            "Violin plot" = "violin", "Boxplot" = "boxplot",
            "Ridgeline" = "ridgeline", "ECDF" = "ecdf",
            "Scatter (2 features)" = "scatter", "Hexbin (2 features)" = "hexbin",
            "Heatmap" = "heatmap", "Bubble plot" = "bubble",
            "Spider / Radar" = "radar"),
          selected = "violin"),
        tags$hr(),
        # --- Metric selector (for summary-based charts) ---
        conditionalPanel("['bar','line','heatmap','bubble'].indexOf(input.cm_chart_type) >= 0",
          selectInput("cm_metric", "Summary metric",
            choices = c("Mean" = "mean_val", "Median" = "median_val", "SD" = "sd_val",
                        "IQR" = "iqr_val", "CV (%)" = "cv_val", "Skewness" = "skewness",
                        "N Cells" = "n_cells", "Q25" = "q25", "Q75" = "q75"),
            selected = "mean_val")
        ),
        # --- Bubble size metric ---
        conditionalPanel("input.cm_chart_type == 'bubble'",
          selectInput("cm_bubble_size", "Bubble size by",
            choices = c("N Cells" = "n_cells", "SD" = "sd_val", "CV (%)" = "cv_val"),
            selected = "n_cells")
        ),
        # --- Scatter / Hexbin: two feature selectors ---
        conditionalPanel("input.cm_chart_type == 'scatter' || input.cm_chart_type == 'hexbin'",
          tags$div(class = "alert alert-secondary px-2 py-1 small mb-2",
            "Compare two measurements against each other (e.g., Area vs Perimeter).
             Groups are shown via color/facets."),
          uiOutput("cm_scatter_features")
        ),
        # --- Radar: timepoint selector ---
        conditionalPanel("input.cm_chart_type == 'radar'",
          uiOutput("cm_radar_time_ui")
        ),
        tags$hr(),
        # --- Grouping control (hidden for scatter/hexbin/radar) ---
        conditionalPanel("['scatter','hexbin','radar'].indexOf(input.cm_chart_type) < 0",
          tags$p(tags$strong("Layout"), class = "mb-1"),
          radioButtons("cm_grouping", "Group axis by",
            choices = c("By timepoint (compare conditions)" = "by_time",
                        "By condition (compare timepoints)" = "by_cond"),
            selected = "by_time")
        ),
        checkboxInput("cm_show_values", "Show values on chart", FALSE),
        tags$hr(),
        # --- Significance testing ---
        tags$p(tags$strong("Significance Testing"), class = "mb-1"),
        checkboxInput("cm_show_sig", "Show significance brackets", FALSE),
        conditionalPanel("input.cm_show_sig",
          radioButtons("cm_sig_mode", "Comparison type",
            choices = c("Each treatment vs control" = "vs_control",
                        "All pairwise (every pair)" = "all_pairwise"),
            selected = "vs_control"),
          radioButtons("cm_sig_stat", "Test using replicate-level:",
            choices = c("Medians" = "median_val", "Means" = "mean_val"),
            selected = "median_val", inline = TRUE),
          radioButtons("cm_sig_label", "Label style",
            choices = c("Stars (*, **, ***)" = "stars",
                        "Stars + p-value" = "stars_p",
                        "p-value only" = "pval"), selected = "stars"),
          checkboxInput("cm_show_ns", "Show ns brackets", FALSE),
          tags$div(class = "alert alert-info px-2 py-1 small",
            "Test and correction follow the Statistics tab settings (default: Welch's t-test, Holm). Always on biological replicate summaries.")
        ),
        tags$hr(),
        figure_controls_ui("cm", 8, 5),
        tags$hr(),
        downloadButton("dl_cm_plot", "Download Figure", class = "btn-primary w-100"),
        downloadButton("dl_cm_table", "Download Table (CSV)", class = "btn-outline-primary w-100 mt-2")
      ),
      navset_card_tab(
        nav_panel("Chart",
          plotOutput("cm_chart", height = "550px", dblclick = "cm_chart_dblclick"),
          edit_hint_ui()),
        nav_panel("Table", DTOutput("cm_table")),
        nav_panel("About", uiOutput("cm_about"))
      )
    )
  ),

  # === TAB 10: PUBLICATION FIGURE ==============================================
  nav_panel(
    title = "Pub Figure",
    icon  = icon("image"),
    layout_sidebar(
      sidebar = sidebar(
        width = 300, title = "Composite Figure",
        checkboxGroupInput("pub_panels", "Include panels:",
          choices = c("A: Distributions" = "dist", "B: Trajectories" = "traj",
                      "C: Effect sizes" = "effect", "D: Morphospace" = "morphospace",
                      "E: Phenotype composition" = "phenotype",
                      "F: Cell Metrics (line)" = "cellmetrics",
                      "G: Spider/Radar" = "radar",
                      "H: Morph. Trajectory" = "morph_traj",
                      "I: Morphospace Diversity" = "morph_div"),
          selected = c("dist", "traj")),
        tags$hr(),
        figure_controls_ui("pub", 12, 5, show_legend = FALSE),
        tags$hr(),
        downloadButton("dl_pub", "Download Figure", class = "btn-primary w-100"),
        downloadButton("dl_pub_deck", "Export Panels to PPTX Deck", class = "btn-outline-primary w-100 mt-2"),
        tags$p(class = "small text-muted mt-1",
          "PPTX deck: one editable vector slide per panel (requires officer + rvg).")
      ),
      card(card_header("Composite Publication Figure"),
           plotOutput("pub_plot", height = "550px", dblclick = "pub_plot_dblclick"),
           edit_hint_ui())
    )
  ),

  # === TAB 11: GUIDE & SCRIPT BUILDER =========================================
  nav_panel(
    title = "Guide",
    icon  = icon("book-open"),
    navset_card_tab(
      nav_panel("Quick Start",
        tags$div(class = "p-3", style = "max-width:900px;",
          uiOutput("guide_status"),
          tags$h5("1. Load your data"),
          tags$p("Go to the ", tags$strong("Data"), " tab and upload a .csv, .tsv or Excel file with ",
                 tags$strong("one row per cell"), ". Two example files ship with the app in the ",
                 tags$code("examples/"), " folder. You need:"),
          tags$ul(
            tags$li(tags$code("Condition"), " (or Sample / Group / Treatment) - the experimental group"),
            tags$li(tags$code("Replicate"), " - the biological replicate. This is your true n."),
            tags$li("at least one numeric measurement column (Length, Area, Circularity, ...)"),
            tags$li(tags$code("Time_min"), " is ", tags$strong("optional"),
                    " - leave it out for a single snapshot and samples are compared directly.")
          ),
          tags$h5("2. Pick a measurement"),
          tags$p("Choose the Measurement Variable and reference (control) condition in the Data sidebar. ",
                 "Everything else follows from that choice."),
          tags$h5("3. Make the figure"),
          tags$p("The ", tags$strong("Distributions"), " tab is the one you will use most. Turn on ",
                 tags$em("Distinguish replicates by shape"), " to show sample-to-sample variation ",
                 "honestly without any p-values (a SuperPlot; Lord et al. 2020, J Cell Biol)."),
          tags$h5("4. Add statistics only if you need them"),
          tags$p("On the ", tags$strong("Statistics"), " tab choose the grouping, test and correction, ",
                 "then switch brackets on in the Distributions sidebar. Figures always use the same ",
                 "settings as the tables. Every test runs on biological-replicate summaries, never on ",
                 "individual cells."),
          tags$h5("5. Edit and export"),
          tags$p(icon("pen"), tags$strong(" Double-click any figure"),
                 " to rename samples and timepoints or change any label. Then export as PDF, SVG, ",
                 "TIFF, PNG, EPS, or an editable PowerPoint."),
          tags$hr(),
          tags$p(class = "text-muted small",
            "Full documentation: README.md - data layout: DATA_FORMAT.md - citable methods: METHODS.md")
        )
      ),
      nav_panel("Check My Data",
        tags$div(class = "p-3", style = "max-width:900px;",
          tags$p("A plain-language report of what the app detected in your file, and whether the ",
                 "design supports the statistics you might want."),
          actionButton("guide_check", "Check my data", class = "btn-primary", icon = icon("stethoscope")),
          tags$hr(),
          uiOutput("guide_report")
        )
      ),
      nav_panel("Script Builder",
        tags$div(class = "p-3",
          tags$p("Build a ready-to-run R script that reproduces this analysis outside the app - ",
                 "same helper functions, same statistics, same theme. Use it to re-run when data is ",
                 "updated, batch over many measurements, or attach a reproducible script to a paper."),
          layout_columns(col_widths = c(5, 7),
            tags$div(
              uiOutput("script_measure_ui"),
              uiOutput("script_group_ui"),
              selectInput("script_plot", "Figure",
                choices = c("Violin + replicate points (SuperPlot)" = "violin",
                            "Sina + replicate points" = "sina",
                            "Box + replicate points" = "box",
                            "Replicate dot plot only" = "dot"),
                selected = "violin"),
              selectInput("script_stat_unit", "Statistical unit",
                choices = c("Replicate median" = "median_val", "Replicate mean" = "mean_val"),
                selected = "median_val"),
              selectInput("script_test", "Pairwise test",
                choices = c("Welch's t-test" = "welch", "Wilcoxon rank-sum" = "wilcox"),
                selected = "welch"),
              selectInput("script_corr", "P-value correction",
                choices = c("Holm" = "holm", "Benjamini-Hochberg" = "BH",
                            "Bonferroni" = "bonferroni", "None" = "none"),
                selected = "holm"),
              checkboxInput("script_all_measures", "Loop over every measurement column", FALSE),
              checkboxInput("script_trend", "Include dose-response trend test", FALSE),
              conditionalPanel("input.script_trend",
                uiOutput("script_trend_pred_ui")),
              selectInput("script_fmt", "Figure format",
                choices = c("PDF" = "pdf", "PNG" = "png", "SVG" = "svg", "TIFF" = "tiff"),
                selected = "pdf"),
              actionButton("script_build", "Generate script",
                class = "btn-primary w-100 mt-2", icon = icon("wand-magic-sparkles")),
              downloadButton("script_download", "Download .R file",
                class = "btn-outline-primary w-100 mt-2")
            ),
            tags$div(
              tags$p(tags$strong("Generated script"), class = "mb-1"),
              tags$p(class = "small text-muted",
                "Save this next to CellMorphR_v3.R and run:  Rscript my_analysis.R"),
              verbatimTextOutput("script_preview")
            )
          )
        )
      ),
      nav_panel("Statistics Help",
        tags$div(class = "p-3", style = "max-width:900px;",
          tags$h5("Your n is the number of biological replicates"),
          tags$p("Thousands of cells from one experiment are not thousands of independent ",
                 "observations. CellMorphR collapses each replicate to a single value before ",
                 "any test, so cells cannot inflate significance (pseudoreplication)."),
          tags$h5("Which test?"),
          tags$ul(
            tags$li(tags$strong("Welch's t-test"), " - the safe default for comparing two groups of ",
                    "replicate summaries. Does not assume equal variance."),
            tags$li(tags$strong("Wilcoxon rank-sum"), " - rank based, but needs a reasonable number ",
                    "of replicates. With 2 vs 2 it can never reach p < 0.05."),
            tags$li(tags$strong("Linear mixed model"), " - uses cell-level data with replicate as a ",
                    "random effect. The best option when you only have 2 replicates."),
            tags$li(tags$strong("Trend test"), " - for dose gradients, asks whether the measurement ",
                    "changes with dose. One focused question instead of many pairwise ones.")
          ),
          tags$h5("Which correction?"),
          tags$p(tags$strong("Holm"), " is a safe default. ", tags$strong("Benjamini-Hochberg"),
                 " is appropriate when you are screening many features and can tolerate a ",
                 "controlled false-discovery rate. Use ", tags$strong("none"),
                 " only for a single pre-planned comparison."),
          tags$h5("How many replicates do I need?"),
          tags$p("Three or more independent experiments is the practical minimum for pairwise ",
                 "p-values. With two, prefer the mixed model or a trend test, and present a ",
                 "SuperPlot showing the replicate-level agreement."),
          tags$h5("Reporting"),
          tags$p("Report the test, the correction, the number of biological replicates, and an ",
                 "effect size with its confidence interval - not a p-value alone. The Pairwise ",
                 "table gives you all of these. METHODS.md contains a citable description.")
        )
      )
    )
  ),

  # === FOOTER ==================================================================
  nav_spacer(),
  nav_item(actionLink("theme_toggle", label = NULL,
    icon = icon("moon"), class = "btn btn-sm",
    style = "color:#FFB74D; font-size:1.1rem; padding:4px 8px;",
    title = "Toggle Light/Dark Mode")),
  nav_item(tags$span(class = "navbar-text small",
    style = "color:#8888a8;",
    paste0("CellMorphR v", APP_VERSION, " | "),
    tags$em("Built with passion for science and discovery"),
    tags$span(" | ", style = "opacity:0.5;"),
    tags$strong("YOU CAN JUST BUILD THINGS", style = "color:#FFB74D; letter-spacing:0.5px;")
  ))
)


# =============================================================================
# SERVER
# =============================================================================
server <- function(input, output, session) {

  # --- Theme toggle (light/dark) ---
  theme_mode <- reactiveVal("dark")
  observeEvent(input$theme_toggle, {
    if (theme_mode() == "dark") {
      theme_mode("light")
      session$setCurrentTheme(light_theme)
      session$sendCustomMessage("swap-css", light_css)
    } else {
      theme_mode("dark")
      session$setCurrentTheme(dark_theme)
      session$sendCustomMessage("swap-css", dark_css)
    }
  })

  # --- Apply typography to all tabs ---
  all_fig_prefixes <- c("dist", "traj", "pca", "clust", "es", "mf", "cm", "pub")
  typo_fields <- c("base_size", "font_family", "legend_pos",
    "title_size", "title_face", "axis_title_size", "axis_title_face",
    "axis_text_size", "axis_text_face", "legend_title_size", "legend_title_face",
    "legend_text_size", "legend_text_face", "strip_size", "strip_face")

  lapply(all_fig_prefixes, function(src) {
    observeEvent(input[[paste0(src, "_apply_global")]], {
      for (tgt in setdiff(all_fig_prefixes, src)) {
        for (fld in typo_fields) {
          val <- input[[paste0(src, "_", fld)]]
          if (!is.null(val)) {
            tgt_id <- paste0(tgt, "_", fld)
            if (fld == "base_size") {
              updateSliderInput(session, tgt_id, value = val)
            } else if (grepl("_size$", fld) && fld != "base_size") {
              updateNumericInput(session, tgt_id, value = val)
            } else {
              updateSelectInput(session, tgt_id, selected = val)
            }
          }
        }
      }
      showNotification("Typography applied to all tabs.", type = "message")
    })
  })

  # --- Helper: Save figures ---
  save_figure <- function(file, plot_obj, width, height, dpi, format) {
    if (format == "tiff") {
      ggsave(file, plot_obj, width = width, height = height, dpi = dpi, device = "tiff",
             compression = "lzw", bg = "white")
    } else if (format == "eps") {
      if (isTRUE(capabilities("cairo"))) {
        ggsave(file, plot_obj, width = width, height = height, dpi = dpi, device = cairo_ps, bg = "white")
      } else {
        showNotification("Cairo is unavailable on this R build; using the standard postscript device (fonts are not embedded).",
                         type = "warning", duration = 8)
        ggsave(file, plot_obj, width = width, height = height, dpi = dpi, device = "eps", bg = "white")
      }
    } else if (format == "pptx") {
      if (has_officer && has_rvg) {
        tryCatch({
          # Scale the figure into a standard 10 x 7.5 in slide, preserving aspect
          s <- min(10 / width, 7.5 / height, 1)
          doc <- officer::read_pptx()
          doc <- officer::add_slide(doc, layout = "Blank", master = "Office Theme")
          doc <- officer::ph_with(doc, rvg::dml(ggobj = plot_obj),
                   location = officer::ph_location(
                     left = (10 - width * s) / 2, top = (7.5 - height * s) / 2,
                     width = width * s, height = height * s))
          print(doc, target = file)
        }, error = function(e) {
          showNotification(paste("Editable PPTX export failed for this plot:", conditionMessage(e),
                                 "- saving a PDF-in-pptx fallback instead."), type = "warning", duration = 10)
          ggsave(file, plot_obj, width = width, height = height, dpi = dpi, device = "pdf")
        })
      } else {
        showNotification("Editable PPTX export needs the 'officer' and 'rvg' packages. Run run_cellmorphr.bat to install them, or install.packages(c('officer','rvg')). Saving PDF instead.",
                         type = "warning", duration = 10)
        ggsave(file, plot_obj, width = width, height = height, dpi = dpi, device = "pdf")
      }
    } else if (format == "svg") {
      if (requireNamespace("svglite", quietly = TRUE))
        ggsave(file, plot_obj, width = width, height = height, dpi = dpi, device = svglite::svglite)
      else ggsave(file, plot_obj, width = width, height = height, dpi = dpi, device = "svg")
    } else if (format == "png") {
      ggsave(file, plot_obj, width = width, height = height, dpi = dpi, device = "png", bg = "white")
    } else {
      ggsave(file, plot_obj, width = width, height = height, dpi = dpi, device = "pdf")
    }
  }
  fmt_ext <- function(format) format

  # === REACTIVE: loaded data ===
  raw_data <- reactiveVal(NULL)
  measure_cols <- reactiveVal(NULL)

  # --- Conditional landing page: tells UI whether data exists ---
  output$has_data <- reactive({ !is.null(raw_data()) })
  outputOptions(output, "has_data", suspendWhenHidden = FALSE)

  # Lets the UI hide time-specific controls for single-snapshot experiments
  output$has_timecourse <- reactive({ !is.null(raw_data()) && has_time() })
  outputOptions(output, "has_timecourse", suspendWhenHidden = FALSE)

  # --- Column mapping state ---
  pending_raw <- reactiveVal(NULL)   # holds raw df while user maps columns
  col_map_rv  <- reactiveVal(NULL)   # holds confirmed col_map

  # --- Timepoint display labels -------------------------------------------
  # Time_min stays numeric (ordering, trend tests, trajectories depend on it);
  # only how each timepoint is *written* on figures is user-editable.
  time_label_map <- reactiveVal(NULL)   # named chr: names = as.character(Time_min)
  time_synthetic <- reactiveVal(FALSE)  # TRUE when the file had no time column

  # TRUE when the experiment really has a time dimension worth plotting
  has_time <- reactive({
    if (isTRUE(time_synthetic())) return(FALSE)
    length(time_values()) > 1 || !identical(as.numeric(time_values()), 0)
  })

  fmt_time <- function(v) {
    key <- as.character(v)
    m <- time_label_map()
    out <- if (!is.null(m)) unname(m[key]) else rep(NA_character_, length(key))
    dflt <- if (isTRUE(time_synthetic())) rep("All samples", length(key)) else paste0(key, " min")
    ifelse(is.na(out) | !nzchar(out), dflt, out)
  }

  # Timepoints present in the data, in numeric order
  time_values <- reactive({
    req(raw_data())
    tv <- suppressWarnings(as.numeric(as.character(raw_data()$Time_min)))
    sort(unique(tv[!is.na(tv)]))
  })

  # Apply a set of timepoint labels; blanks fall back to the "<n> min" default
  apply_time_labels <- function(labels) {
    tv <- isolate(time_values())
    if (length(tv) == 0) return(invisible(FALSE))
    labels <- trimws(as.character(labels))
    labels[!nzchar(labels)] <- paste0(tv[!nzchar(labels)], " min")
    if (anyDuplicated(labels)) {
      showNotification("Timepoint labels must be unique.", type = "warning"); return(invisible(FALSE))
    }
    time_label_map(setNames(labels, as.character(tv)))
    invisible(TRUE)
  }

  # --- Rename conditions: change display names everywhere (legends, axes,
  # tables, stats, exports) by rewriting the values in the data itself. ---
  apply_condition_rename <- function(mapping) {
    # mapping: named chr vector, names = old values, values = new values
    new_names <- trimws(mapping)
    if (any(!nzchar(new_names))) {
      showNotification("Names cannot be empty.", type = "warning"); return(invisible(FALSE))
    }
    if (anyDuplicated(new_names)) {
      showNotification("Names must be unique -- two conditions would merge.", type = "warning"); return(invisible(FALSE))
    }
    rd <- raw_data()
    old_ref <- isolate(input$ref_condition)
    rd$Condition <- unname(new_names[as.character(rd$Condition)])
    raw_data(rd)
    # Keep the reference-condition selection pointing at the renamed value
    if (!is.null(old_ref) && old_ref %in% names(new_names)) {
      updateSelectInput(session, "ref_condition", selected = unname(new_names[old_ref]))
    }
    showNotification("Condition names updated everywhere.", type = "message")
    invisible(TRUE)
  }

  output$condition_rename_ui <- renderUI({
    req(raw_data())
    conds <- unique(as.character(raw_data()$Condition))
    tagList(lapply(seq_along(conds), function(i) {
      textInput(paste0("cond_rename_", i), NULL, value = conds[i], width = "100%")
    }))
  })

  observeEvent(input$rename_conditions_apply, {
    req(raw_data())
    conds <- unique(as.character(raw_data()$Condition))
    new_names <- vapply(seq_along(conds), function(i) input[[paste0("cond_rename_", i)]] %||% conds[i], character(1))
    mapping <- setNames(new_names, conds)
    if (identical(unname(mapping), conds)) {
      showNotification("No changes to apply.", type = "message"); return()
    }
    apply_condition_rename(mapping)
  })

  output$dist_rep_shape_ui <- renderUI({
    req(raw_data())
    reps <- levels(factor(raw_data()$Replicate))
    tagList(lapply(seq_along(reps), function(i) {
      selectInput(paste0("dist_rep_shape_", i), paste("Replicate", reps[i]),
                  choices = SHAPE_CHOICES,
                  selected = DEFAULT_REP_SHAPES[((i - 1) %% length(DEFAULT_REP_SHAPES)) + 1])
    }))
  })

  output$time_rename_ui <- renderUI({
    req(time_values())
    if (!has_time()) return(tags$p(class = "small text-muted",
      "This dataset has no time column, so there are no timepoint labels to edit."))
    tv <- time_values(); cur <- fmt_time(tv)
    tagList(lapply(seq_along(tv), function(i)
      textInput(paste0("time_rename_", i), NULL, value = cur[i], width = "100%")))
  })

  observeEvent(input$rename_times_apply, {
    tv <- time_values(); req(length(tv) > 0)
    cur <- fmt_time(tv)
    new_t <- vapply(seq_along(tv), function(i) input[[paste0("time_rename_", i)]] %||% cur[i], character(1))
    if (identical(trimws(new_t), cur)) {
      showNotification("No changes to apply.", type = "message"); return()
    }
    if (isTRUE(apply_time_labels(new_t)))
      showNotification("Timepoint labels updated on all figures.", type = "message")
  })

  observeEvent(input$rename_remove_underscores, {
    req(raw_data())
    conds <- unique(as.character(raw_data()$Condition))
    mapping <- setNames(gsub("_", " ", conds, fixed = TRUE), conds)
    if (identical(unname(mapping), conds)) {
      showNotification("No underscores found in condition names.", type = "message"); return()
    }
    apply_condition_rename(mapping)
  })

  # --- Double-click-to-edit: modal editor for the displayed figure ------------
  # Double-clicking a main figure opens this editor; applied changes feed the
  # tab's custom-label inputs (so preview and downloads stay identical) and the
  # condition renames rewrite the data itself (so every tab and table follows).
  label_editor_prefix <- reactiveVal(NULL)

  open_label_editor <- function(prefix, tab_title) {
    label_editor_prefix(prefix)
    cur <- function(id) isolate(input[[paste0(prefix, id)]]) %||% ""
    conds <- if (!is.null(isolate(raw_data()))) unique(as.character(isolate(raw_data())$Condition)) else character(0)
    tvals <- if (!is.null(isolate(raw_data())) && isolate(has_time())) isolate(time_values()) else numeric(0)
    tlabs <- if (length(tvals) > 0) fmt_time(tvals) else character(0)
    showModal(modalDialog(
      title = tagList(icon("pen"), paste0(" Edit figure text -- ", tab_title)),
      size = "l", easyClose = TRUE,
      layout_columns(col_widths = c(6, 6),
        tags$div(
          tags$p(tags$strong("Figure text"), class = "mb-1"),
          tags$p(class = "small text-muted", "Blank fields keep the automatic text."),
          textInput("led_title",    "Title",        cur("_custom_title")),
          textInput("led_subtitle", "Subtitle",     cur("_custom_subtitle")),
          textInput("led_x",        "X axis label", cur("_custom_x_label")),
          textInput("led_y",        "Y axis label", cur("_custom_y_label")),
          textInput("led_legend",   "Legend title", cur("_custom_legend")),
          textInput("led_caption",  "Caption",      cur("_custom_caption"))
        ),
        tags$div(
          tags$p(tags$strong("Sample names"), class = "mb-1"),
          tags$p(class = "small text-muted",
            "Renaming applies everywhere: all figures, statistics tables, and exports."),
          actionButton("led_unders", "Replace underscores with spaces",
            class = "btn-outline-secondary btn-sm w-100 mb-2", icon = icon("broom")),
          tags$p(tags$strong("Conditions"), class = "small mb-1"),
          lapply(seq_along(conds), function(i)
            textInput(paste0("led_cond_", i), NULL, value = conds[i], width = "100%")),
          if (length(tvals) > 0) tagList(
            tags$p(tags$strong("Timepoints"), class = "small mb-1 mt-2"),
            lapply(seq_along(tvals), function(i)
              textInput(paste0("led_time_", i), NULL, value = tlabs[i], width = "100%"))
          )
        )
      ),
      footer = tagList(modalButton("Cancel"),
                       actionButton("led_apply", "Apply", class = "btn-primary", icon = icon("check")))
    ))
  }

  # One observer per editable figure
  local({
    editor_map <- list(
      dist_plot_dblclick         = c("dist",  "Distributions"),
      traj_plot_dblclick         = c("traj",  "Trajectories"),
      pca_scatter_dblclick       = c("pca",   "Morphospace"),
      clust_composition_dblclick = c("clust", "Phenotype Composition"),
      es_plot_dblclick           = c("es",    "Effect Sizes"),
      mf_heatmap_dblclick        = c("mf",    "Multi-Feature Heatmap"),
      cm_chart_dblclick          = c("cm",    "Cell Metrics"),
      pub_plot_dblclick          = c("pub",   "Publication Figure")
    )
    for (id in names(editor_map)) {
      local({
        id2 <- id; m <- editor_map[[id]]
        observeEvent(input[[id2]], open_label_editor(m[1], m[2]), ignoreInit = TRUE)
      })
    }
  })

  observeEvent(input$led_unders, {
    req(raw_data())
    conds <- unique(as.character(raw_data()$Condition))
    for (i in seq_along(conds)) {
      curv <- input[[paste0("led_cond_", i)]] %||% conds[i]
      updateTextInput(session, paste0("led_cond_", i), value = gsub("_", " ", curv, fixed = TRUE))
    }
    tv <- time_values(); cur <- fmt_time(tv)
    for (i in seq_along(tv)) {
      curv <- input[[paste0("led_time_", i)]] %||% cur[i]
      updateTextInput(session, paste0("led_time_", i), value = gsub("_", " ", curv, fixed = TRUE))
    }
  })

  observeEvent(input$led_apply, {
    prefix <- label_editor_prefix(); req(prefix)
    updateTextInput(session, paste0(prefix, "_custom_title"),    value = input$led_title %||% "")
    updateTextInput(session, paste0(prefix, "_custom_subtitle"), value = input$led_subtitle %||% "")
    updateTextInput(session, paste0(prefix, "_custom_x_label"),  value = input$led_x %||% "")
    updateTextInput(session, paste0(prefix, "_custom_y_label"),  value = input$led_y %||% "")
    updateTextInput(session, paste0(prefix, "_custom_legend"),   value = input$led_legend %||% "")
    updateTextInput(session, paste0(prefix, "_custom_caption"),  value = input$led_caption %||% "")
    if (!is.null(raw_data())) {
      conds <- unique(as.character(raw_data()$Condition))
      new_names <- vapply(seq_along(conds), function(i) input[[paste0("led_cond_", i)]] %||% conds[i], character(1))
      if (!identical(trimws(new_names), conds)) apply_condition_rename(setNames(new_names, conds))
      tv <- time_values()
      if (length(tv) > 0) {
        cur <- fmt_time(tv)
        new_t <- vapply(seq_along(tv), function(i) input[[paste0("led_time_", i)]] %||% cur[i], character(1))
        if (!identical(trimws(new_t), cur)) apply_time_labels(new_t)
      }
    }
    removeModal()
  })

  # --- Split Condition into factor columns (e.g., vector + numeric dose) ---
  output$split_preview <- renderPrint({
    req(raw_data())
    uc <- unique(as.character(raw_data()$Condition))
    sp <- split_condition_parts(uc, input$split_delim %||% "_")
    prev <- data.frame(Condition = uc, Base = sp$base, check.names = FALSE)
    if (any(!is.na(sp$num_val))) prev[[sp$dose_col]] <- sp$num_val
    print(utils::head(prev, 8), row.names = FALSE)
  })

  observeEvent(input$split_condition_apply, {
    req(raw_data())
    rd <- raw_data()
    sp <- split_condition_parts(rd$Condition, input$split_delim %||% "_")
    added <- character(0)
    if (any(!is.na(sp$num_val))) {
      rd[[sp$dose_col]] <- sp$num_val
      rd[["Condition_Base"]] <- sp$base
      added <- c(added, sp$dose_col, "Condition_Base")
    }
    if (length(sp$parts) > 1) {
      for (nm in names(sp$parts)) rd[[nm]] <- sp$parts[[nm]]
      added <- c(added, names(sp$parts))
    }
    if (length(added) == 0) {
      showNotification("Nothing to split: Condition has no delimiter or numeric token.", type = "warning")
    } else {
      raw_data(rd)
      showNotification(paste("Added columns:", paste(added, collapse = ", "),
                             "- now available for statistical grouping and the Trend test."), type = "message", duration = 8)
    }
  })

  # --- Helper: load and validate a data file ---
  load_data_file <- function(filepath, filename) {
    ext <- tolower(tools::file_ext(filename))
    df <- tryCatch({
      if (ext %in% c("xlsx", "xls")) readxl::read_excel(filepath)
      else if (ext == "csv")  read.csv(filepath, stringsAsFactors = FALSE, check.names = FALSE)
      else if (ext == "tsv")  read.delim(filepath, stringsAsFactors = FALSE, check.names = FALSE)
      else {
        # Try CSV as fallback for any text file
        tryCatch(read.csv(filepath, stringsAsFactors = FALSE, check.names = FALSE),
                 error = function(e) NULL)
      }
    }, error = function(e) NULL)

    if (is.null(df) || nrow(df) == 0) {
      showNotification("Could not read the file. Check that it is a valid CSV, TSV, or Excel file.",
                       type = "error", duration = 8)
      return()
    }

    val <- tryCatch(validate_data(df), error = function(e) {
      list(ok = FALSE, needs_mapping = FALSE,
           msg = paste("Validation error:", conditionMessage(e)))
    })

    if (val$ok) {
      # Perfect: load directly
      finalize_load(val$df, val$measures, synthetic_time = isTRUE(val$synthetic_time))
    } else if (isTRUE(val$needs_mapping)) {
      # Show column mapping modal
      pending_raw(df)
      showModal(col_map_modal(names(df), val$cond_col, val$time_col, val$rep_col))
    } else {
      showNotification(val$msg, type = "error", duration = 12)
    }
  }

  finalize_load <- function(df, measures, synthetic_time = FALSE) {
    time_synthetic(isTRUE(synthetic_time))
    time_label_map(NULL)
    df$Condition <- as.character(df$Condition)
    df$Time_min  <- suppressWarnings(as.numeric(df$Time_min))
    if (any(is.na(df$Time_min))) {
      # Time might be non-numeric (e.g. "0h","24h") — convert by stripping non-numeric suffix
      raw_time <- as.character(df$Time_min)
      df$Time_min <- suppressWarnings(as.numeric(gsub("[^0-9.-]", "", raw_time)))
    }
    df$Replicate <- as.factor(df$Replicate)
    raw_data(df)
    measure_cols(measures)
    pending_raw(NULL)
    removeModal()
    showNotification(
      paste0("Loaded ", format(nrow(df), big.mark=","), " rows, ",
             length(measures), " measurement columns.",
             if (isTRUE(synthetic_time)) " No time column found - samples are compared directly." else ""),
      type = "message", duration = if (isTRUE(synthetic_time)) 9 else 5
    )
  }

  # --- Column mapping modal ---
  col_map_modal <- function(all_cols, guess_cond, guess_time, guess_rep) {
    make_sel <- function(id, label, guess) {
      selectInput(id, label, choices = all_cols,
                  selected = if (!is.null(guess) && guess %in% all_cols) guess else all_cols[1])
    }
    modalDialog(
      title = "Map Column Roles",
      tags$p(class="text-muted small",
        "The app couldn't auto-detect required column roles. ",
        "Please assign which column represents each role below."),
      make_sel("map_cond_col", "Condition / Group column", guess_cond),
      selectInput("map_time_col", "Time / Timepoint column",
                  choices = c("(no timepoints - compare samples directly)" = "__none__", all_cols),
                  selected = if (!is.null(guess_time) && guess_time %in% all_cols) guess_time else "__none__"),
      make_sel("map_rep_col",  "Replicate / Experiment column", guess_rep),
      tags$hr(),
      tags$p(class="text-muted small",
        "All remaining numeric columns will be treated as measurements."),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_col_map", "Load Data", class="btn-primary")
      ),
      easyClose = FALSE, size = "m"
    )
  }

  # --- Confirm column mapping ---
  observeEvent(input$confirm_col_map, {
    req(pending_raw(), input$map_cond_col, input$map_time_col, input$map_rep_col)
    df <- pending_raw()
    cmap <- list(condition = input$map_cond_col,
                 time      = if (identical(input$map_time_col, "__none__")) NULL else input$map_time_col,
                 replicate = input$map_rep_col)
    val <- tryCatch(validate_data(df, col_map = cmap), error = function(e) {
      list(ok=FALSE, needs_mapping=FALSE, msg=paste("Error:", conditionMessage(e)))
    })
    if (val$ok) {
      finalize_load(val$df, val$measures, synthetic_time = isTRUE(val$synthetic_time))
    } else {
      showNotification(val$msg, type="error", duration=10)
    }
  })

  # File upload (landing page)
  observeEvent(input$file, { req(input$file); load_data_file(input$file$datapath, input$file$name) })

  # File replace (sidebar, after data loaded)
  observeEvent(input$file_replace, { req(input$file_replace); load_data_file(input$file_replace$datapath, input$file_replace$name) })

  # --- Generate demo data ---
  generate_demo <- function() {
    set.seed(42)
    conditions <- c("Uninfected", "+N4"); times <- c(30, 60, 90, 120); reps <- 1:3
    rows <- list()
    for (cond in conditions) for (tm in times) for (rp in reps) {
      n_cells <- sample(180:220, 1)
      if (cond == "Uninfected") { ba <- 3.0 + (tm - 30) * 0.002; areas <- rlnorm(n_cells, log(ba), 0.25)
      } else { ba <- 3.0 + (tm - 30) * 0.025 + (tm / 120)^2 * 2; areas <- rlnorm(n_cells, log(ba), 0.3 + tm * 0.001) }
      areas <- areas * (1 + (rp - 2) * 0.03)
      perimeters <- 2 * sqrt(pi * areas) * runif(n_cells, 0.95, 1.15)
      rows[[length(rows) + 1]] <- data.frame(Condition = cond, Time_min = tm, Replicate = rp,
        Area = round(areas, 3), Perimeter = round(perimeters, 3),
        Circularity = round(runif(n_cells, 0.6, 1.0), 3))
    }
    df <- bind_rows(rows); df$Replicate <- as.factor(df$Replicate)
    raw_data(df); measure_cols(c("Area", "Perimeter", "Circularity"))
    showNotification("Demo data loaded.", type = "message")
  }

  # Demo data buttons (landing page + sidebar)
  observeEvent(input$load_demo, { generate_demo() })
  observeEvent(input$load_demo_again, { generate_demo() })

  # --- Save session (.rds) ---
  output$save_session <- downloadHandler(
    filename = function() paste0("CellMorphR_session_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds"),
    content = function(file) {
      session_data <- list(
        data = raw_data(),
        measure_cols = measure_cols(),
        settings = list(
          measure = input$measure,
          ref_condition = input$ref_condition,
          data_transform = input$data_transform
        ),
        app_version = APP_VERSION,
        saved_at = Sys.time()
      )
      saveRDS(session_data, file)
      showNotification("Session saved successfully.", type = "message")
    }
  )

  # --- Load session (.rds) ---
  # Settings restore uses ONE permanent observer + a pending reactiveVal so that
  # repeated session loads cannot stack up stale one-shot observers.
  pending_session_settings <- reactiveVal(NULL)
  observeEvent(measure_cols(), {
    s <- pending_session_settings()
    if (is.null(s)) return()
    if (!is.null(s$measure)) updateSelectInput(session, "measure", selected = s$measure)
    if (!is.null(s$ref_condition)) updateSelectInput(session, "ref_condition", selected = s$ref_condition)
    if (!is.null(s$data_transform)) updateSelectInput(session, "data_transform", selected = s$data_transform)
    pending_session_settings(NULL)
  }, ignoreInit = TRUE)

  observeEvent(input$load_session, {
    req(input$load_session)
    sess <- tryCatch(readRDS(input$load_session$datapath), error = function(e) NULL)
    if (is.null(sess) || is.null(sess$data)) {
      showNotification("Failed to load session file. Invalid format.", type = "error"); return()
    }
    # Restore data
    df <- sess$data
    df$Condition <- as.character(df$Condition)
    df$Time_min  <- as.numeric(df$Time_min)
    df$Replicate <- as.factor(df$Replicate)
    pending_session_settings(sess$settings)
    raw_data(df); measure_cols(sess$measure_cols)

    saved_time <- if (!is.null(sess$saved_at)) format(sess$saved_at, "%Y-%m-%d %H:%M") else "unknown"
    showNotification(
      paste0("Session restored (saved ", saved_time, "). ",
             nrow(df), " cells, ", length(sess$measure_cols), " features."),
      type = "message", duration = 8
    )
  })

  output$measure_selector <- renderUI({
    req(measure_cols())
    selectInput("measure", "Measurement Variable", choices = measure_cols(), selected = measure_cols()[1])
  })

  output$condition_order_ui <- renderUI({
    req(raw_data())
    conds <- unique(raw_data()$Condition)
    cur <- isolate(input$ref_condition)
    selectInput("ref_condition", "Reference (control) condition", choices = conds,
      selected = if (!is.null(cur) && cur %in% conds) cur else {
        idx <- grep("uninf|control|wt|mock", conds, ignore.case = TRUE)
        if (length(idx) > 0) conds[idx[1]] else conds[1] })
  })


  # --- Flexible statistics UI helpers ---
  output$stat_measure_ui <- renderUI({
    req(measure_cols())
    selectInput("stat_measure", "Measurement variable", choices = measure_cols(),
                selected = input$measure %||% measure_cols()[1])
  })

  output$stat_group_cols_ui <- renderUI({
    req(raw_data())
    meta <- detect_metadata_columns(raw_data(), measure_cols())
    default <- intersect(c("Condition", "Time_min"), meta)
    if (length(default) == 0 && length(meta) > 0) default <- meta[1]
    checkboxGroupInput("stat_group_cols", NULL,
      choices = meta, selected = default, inline = FALSE)
  })

  flexible_rep_summ <- reactive({
    req(transformed_data())
    m <- input$stat_measure %||% input$measure
    req(m)
    cols <- input$stat_group_cols
    if (is.null(cols) || length(cols) == 0) cols <- "Condition"
    yvar <- input$stat_summary %||% "median_val"
    compute_flexible_replicate_summary(transformed_data(), m, cols, yvar = yvar)
  })

  comparison_preview <- reactive({
    req(flexible_rep_summ())
    ref <- input$stat_reference_group
    build_comparison_table(flexible_rep_summ(), input$stat_compare_mode %||% "within_time", ref)
  })

  output$stat_reference_ui <- renderUI({
    req(flexible_rep_summ())
    groups <- sort(unique(as.character(flexible_rep_summ()$Statistics_Group)))
    idx <- grep("uninf|control|wt|mock|empty|vector|0nm|0 nM|0_nM", groups, ignore.case = TRUE)
    default <- if (length(idx) > 0) groups[idx[1]] else groups[1]
    if ((input$stat_compare_mode %||% "within_time") == "all_pairwise") {
      tags$p(class = "small text-muted", "Reference group is not used for all-pairwise mode.")
    } else {
      selectInput("stat_reference_group", "Reference group", choices = groups, selected = default)
    }
  })

  output$stat_group_preview_table <- renderDT({
    req(flexible_rep_summ())
    rs <- flexible_rep_summ()
    preview <- rs %>% group_by(Statistics_Group) %>%
      summarise(n_bio_reps = n_distinct(Replicate), n_cells = sum(n_cells, na.rm = TRUE),
                mean_of_replicate_medians = mean(median_val, na.rm = TRUE),
                mean_of_replicate_means = mean(mean_val, na.rm = TRUE), .groups = "drop")
    preview[] <- lapply(preview, function(x) if (is.numeric(x)) round(x, 4) else x)
    datatable(preview, options = list(pageLength = 15, scrollX = TRUE),
              caption = "Preview of user-defined statistical groups. Each row is tested using biological replicate summaries.")
  })

  output$comparison_preview_table <- renderDT({
    req(comparison_preview())
    datatable(comparison_preview(), options = list(pageLength = 20, scrollX = TRUE),
              caption = "Planned comparisons. These are the only pairwise comparisons used for the flexible statistics output.")
  })


  stat_suggestion <- reactiveVal(NULL)

  observeEvent(input$suggest_stats, {
    req(flexible_rep_summ())
    suggestion <- suggest_statistical_workflow(
      rep_df = flexible_rep_summ(),
      comparison_tbl = comparison_preview(),
      comparison_mode = input$stat_compare_mode %||% "within_time",
      stat_unit = input$stat_summary %||% "median_val",
      selected_group_cols = input$stat_group_cols
    )
    stat_suggestion(suggestion)
    updateSelectInput(session, "stat_test", selected = suggestion$test)
    updateSelectInput(session, "stat_p_adjust", selected = suggestion$correction)
    if (!is.null(suggestion$global_model) && suggestion$global_model %in% c("anova", "lmm", "skip")) {
      updateRadioButtons(session, "stat_method", selected = suggestion$global_model)
    }
    showNotification("Suggested statistics settings were applied. Review the Test Suggestion tab before running.", type = "message")
  })

  current_stat_suggestion <- reactive({
    req(flexible_rep_summ())
    suggest_statistical_workflow(
      rep_df = flexible_rep_summ(),
      comparison_tbl = comparison_preview(),
      comparison_mode = input$stat_compare_mode %||% "within_time",
      stat_unit = input$stat_summary %||% "median_val",
      selected_group_cols = input$stat_group_cols
    )
  })

  output$test_suggestion_output <- renderText({
    req(current_stat_suggestion())
    current_stat_suggestion()$summary
  })

  output$test_suggestion_group_table <- renderDT({
    req(current_stat_suggestion())
    details <- current_stat_suggestion()$details
    if (is.null(details) || nrow(details) == 0) return(NULL)
    details[] <- lapply(details, function(x) if (is.numeric(x)) round(x, 4) else x)
    datatable(details, options = list(pageLength = 15, scrollX = TRUE),
              caption = "Biological replicate counts used to generate the recommendation.")
  })

  # Multi-feature feature selector
  output$mf_feature_selector <- renderUI({
    req(measure_cols())
    checkboxGroupInput("mf_features", "Features to include", choices = measure_cols(),
                       selected = measure_cols())
  })

  # ── Distributions filter UI: one checkbox per unique Condition ──
  output$dist_sample_filter_ui <- renderUI({
    req(raw_data())
    conds <- as.character(unique(raw_data()$Condition))
    # Detect "vector" group (empty vector) vs experimental
    ev_pat <- "uninf|control|wt|mock|empty|ev|psj|vector|0nm"
    is_ev  <- grepl(ev_pat, conds, ignore.case = TRUE)
    tagList(
      # Quick-select buttons
      div(style = "display:flex; gap:4px; margin-bottom:6px; flex-wrap:wrap;",
        actionButton("dist_sel_all",  "All",          class = "btn-xs btn-outline-primary",  style = "font-size:10px; padding:2px 7px;"),
        actionButton("dist_sel_none", "None",         class = "btn-xs btn-outline-secondary", style = "font-size:10px; padding:2px 7px;"),
        actionButton("dist_sel_ev",   "Empty Vector", class = "btn-xs btn-outline-secondary", style = "font-size:10px; padding:2px 7px;"),
        actionButton("dist_sel_exp",  "Experimental", class = "btn-xs btn-outline-secondary", style = "font-size:10px; padding:2px 7px;")
      ),
      checkboxGroupInput("dist_cond_sel", NULL,
        choices  = setNames(conds, conds),
        selected = conds,
        inline   = FALSE)
    )
  })

  # ── Distributions filter UI: one checkbox per unique timepoint ──
  output$dist_time_filter_ui <- renderUI({
    req(raw_data())
    times <- sort(unique(suppressWarnings(as.numeric(raw_data()$Time_min))))
    times <- times[!is.na(times)]
    labels <- paste0(times, " min")
    div(style = "display:flex; flex-wrap:wrap; gap:4px;",
      lapply(seq_along(times), function(i) {
        div(style = "min-width:80px;",
          checkboxInput(
            paste0("dist_time_", times[i]),
            labels[i],
            value = TRUE
          )
        )
      })
    )
  })

  # ── Quick-select observers ──
  observeEvent(input$dist_sel_all, {
    req(raw_data())
    conds <- as.character(unique(raw_data()$Condition))
    updateCheckboxGroupInput(session, "dist_cond_sel", selected = conds)
  })
  observeEvent(input$dist_sel_none, {
    updateCheckboxGroupInput(session, "dist_cond_sel", selected = character(0))
  })
  observeEvent(input$dist_sel_ev, {
    req(raw_data())
    conds <- as.character(unique(raw_data()$Condition))
    ev_pat <- "uninf|control|wt|mock|empty|ev|psj|vector|0nm"
    ev_conds <- conds[grepl(ev_pat, conds, ignore.case = TRUE)]
    if (length(ev_conds) == 0) ev_conds <- conds[1]
    updateCheckboxGroupInput(session, "dist_cond_sel", selected = ev_conds)
  })
  observeEvent(input$dist_sel_exp, {
    req(raw_data())
    conds <- as.character(unique(raw_data()$Condition))
    ev_pat <- "uninf|control|wt|mock|empty|ev|psj|vector|0nm"
    exp_conds <- conds[!grepl(ev_pat, conds, ignore.case = TRUE)]
    if (length(exp_conds) == 0) exp_conds <- conds
    updateCheckboxGroupInput(session, "dist_cond_sel", selected = exp_conds)
  })

  # Ordered data
  plot_data <- reactive({
    req(raw_data())
    df <- raw_data()
    conds <- unique(as.character(df$Condition))
    # Use ref_condition if set, otherwise default to first condition
    ref <- if (isTruthy(input$ref_condition) && input$ref_condition %in% conds)
             input$ref_condition
           else
             conds[1]
    conds <- c(ref, setdiff(conds, ref))
    df$Condition <- factor(df$Condition, levels = conds)
    df$Time_min  <- suppressWarnings(as.numeric(df$Time_min))
    df
  })

  # Transformed data
  transformed_data <- reactive({
    req(plot_data())
    df <- plot_data()
    mc <- tryCatch(measure_cols(), error = function(e) NULL)
    if (is.null(mc) || length(mc) == 0) return(df)
    transform <- input$data_transform %||% "none"
    # One pipeline for every analysis in the app, and it records the transform
    # (including the log shift) so the reports can state it -- see CM-13/CM-14.
    apply_transform_to_columns(df, mc, transform)
  })

  measure_label <- reactive({
    m <- input$measure; tr <- input$data_transform %||% "none"
    if (tr == "log") paste0("log(", m, ")") else if (tr == "sqrt") paste0("sqrt(", m, ")") else m
  })

  # Color pickers
  output$color_pickers <- renderUI({
    req(plot_data())
    conds <- levels(plot_data()$Condition)
    n <- length(conds)

    # Named palette presets
    PALETTE_PRESETS <- list(
      "Custom"       = NULL,
      "Blue/Red"     = c("#4A90D9","#D94A4A","#F5A623","#7ED321","#9B59B6","#E74C3C","#1ABC9C","#F39C12"),
      "Teal/Orange"  = c("#2a9d8f","#e76f51","#e9c46a","#264653","#f4a261","#457b9d","#a8dadc","#1d3557"),
      "Pastel"       = c("#aec6cf","#ffb7b2","#b5ead7","#ffdac1","#c7ceea","#ff9aa2","#e2f0cb","#f8c8d4"),
      "Bold"         = c("#e63946","#2a9d8f","#e9c46a","#264653","#f4a261","#457b9d","#a8dadc","#1d3557"),
      "Publication"  = c("#0077b6","#e63946","#2a9d8f","#f4a261","#6a0dad","#e9c46a","#264653","#f77f00"),
      "Colorblind"   = c("#0077bb","#ee7733","#009988","#cc3311","#33bbee","#ee3377","#bbbbbb","#000000"),
      "Viridis"      = c("#440154","#31688e","#35b779","#fde725","#21918c","#5ec962","#90d743","#c8e020"),
      "Grayscale"    = paste0("grey", round(seq(20, 80, length.out = 8)))
    )

    defaults_by_n <- list(
      "2" = c("#4A90D9","#D94A4A"),
      "3" = c("#4A90D9","#D94A4A","#F5A623"),
      "4" = c("#4A90D9","#D94A4A","#F5A623","#7ED321"),
      "5" = c("#4A90D9","#D94A4A","#F5A623","#7ED321","#9B59B6"),
      "6" = c("#4A90D9","#D94A4A","#F5A623","#7ED321","#9B59B6","#E74C3C")
    )
    dcols <- if (n <= 6) defaults_by_n[[as.character(n)]] else scales::hue_pal()(n)

    sel_preset <- input$palette_preset_sel %||% "Custom"
    if (!is.null(sel_preset) && sel_preset != "Custom" && sel_preset %in% names(PALETTE_PRESETS)) {
      preset_cols <- rep_len(PALETTE_PRESETS[[sel_preset]], n)
    } else {
      preset_cols <- dcols
    }

    tagList(
      tags$p(tags$strong("Color Palette"), class = "mb-1"),
      selectInput("palette_preset_sel", "Quick palette",
        choices = names(PALETTE_PRESETS), selected = sel_preset),
      if (has_colourpicker)
        tagList(lapply(seq_along(conds), function(i) {
          cur <- input[[paste0("color_", i)]]
          init <- if (!is.null(cur) && nchar(cur) > 0) cur else preset_cols[i]
          colourpicker::colourInput(paste0("color_", i), conds[i], init, showColour = "both")
        }))
      else
        tags$p(class = "small text-muted", "Install 'colourpicker' for custom per-condition colors.")
    )
  })

  # When palette preset changes, update color pickers
  observeEvent(input$palette_preset_sel, {
    req(plot_data())
    conds <- levels(plot_data()$Condition); n <- length(conds)
    PALETTE_PRESETS <- list(
      "Blue/Red"    = c("#4A90D9","#D94A4A","#F5A623","#7ED321","#9B59B6","#E74C3C","#1ABC9C","#F39C12"),
      "Teal/Orange" = c("#2a9d8f","#e76f51","#e9c46a","#264653","#f4a261","#457b9d","#a8dadc","#1d3557"),
      "Pastel"      = c("#aec6cf","#ffb7b2","#b5ead7","#ffdac1","#c7ceea","#ff9aa2","#e2f0cb","#f8c8d4"),
      "Bold"        = c("#e63946","#2a9d8f","#e9c46a","#264653","#f4a261","#457b9d","#a8dadc","#1d3557"),
      "Publication" = c("#0077b6","#e63946","#2a9d8f","#f4a261","#6a0dad","#e9c46a","#264653","#f77f00"),
      "Colorblind"  = c("#0077bb","#ee7733","#009988","#cc3311","#33bbee","#ee3377","#bbbbbb","#000000"),
      "Viridis"     = c("#440154","#31688e","#35b779","#fde725","#21918c","#5ec962","#90d743","#c8e020"),
      "Grayscale"   = paste0("grey", round(seq(20, 80, length.out = 8)))
    )
    if (!is.null(input$palette_preset_sel) && input$palette_preset_sel %in% names(PALETTE_PRESETS) && has_colourpicker) {
      cols <- rep_len(PALETTE_PRESETS[[input$palette_preset_sel]], n)
      for (i in seq_along(conds)) {
        colourpicker::updateColourInput(session, paste0("color_", i), value = cols[i])
      }
    }
  }, ignoreInit = TRUE)

  cond_colors <- reactive({
    req(plot_data())
    conds <- levels(plot_data()$Condition); n <- length(conds)
    dcols <- if (n == 2) c("#4A90D9", "#D94A4A") else scales::hue_pal()(n)
    cols <- if (has_colourpicker)
      sapply(seq_along(conds), function(i) { v <- input[[paste0("color_", i)]]; if (!is.null(v) && nchar(v) > 0) v else dcols[i] })
    else dcols
    # Names matter: several tabs subset this by condition name (e.g. dropping the
    # reference condition), which silently returns nothing on an unnamed vector.
    setNames(as.character(cols), conds)
  })

  # Replicate summaries
  rep_summ <- reactive({
    req(transformed_data(), input$measure)
    df <- transformed_data()
    m  <- input$measure
    if (!m %in% names(df)) return(NULL)
    df$Time_min <- suppressWarnings(as.numeric(df$Time_min))
    tryCatch(compute_replicate_summaries(df, m), error = function(e) NULL)
  })
  cond_summ <- reactive({
    rs <- rep_summ()
    req(!is.null(rs))
    tryCatch(compute_condition_summaries(rs), error = function(e) NULL)
  })

  # === DATA TAB ===
  output$data_table <- renderDT({ req(raw_data()); datatable(head(raw_data(), 500), options = list(pageLength = 15, scrollX = TRUE), caption = "First 500 rows") })
  output$n_cells_txt  <- renderText({ req(raw_data()); format(nrow(raw_data()), big.mark = ",") })
  output$n_groups_txt <- renderText({ req(raw_data()); nc <- length(unique(raw_data()$Condition)); nt <- length(unique(raw_data()$Time_min)); paste0(nc, " x ", nt, " = ", nc * nt) })
  output$n_reps_txt   <- renderText({ req(raw_data()); paste(length(unique(raw_data()$Replicate)), "per group") })

  output$feature_diagnostics_table <- renderDT({
    req(raw_data(), measure_cols()); df <- raw_data(); mc <- measure_cols()
    diag <- data.frame(Feature = mc,
      N = sapply(mc, function(m) sum(!is.na(df[[m]]))),
      Mean = sapply(mc, function(m) round(mean(df[[m]], na.rm = TRUE), 3)),
      Median = sapply(mc, function(m) round(median(df[[m]], na.rm = TRUE), 3)),
      SD = sapply(mc, function(m) round(sd(df[[m]], na.rm = TRUE), 3)),
      Skewness = sapply(mc, function(m) round(compute_skewness(df[[m]]), 3)), stringsAsFactors = FALSE)
    diag$Suggestion <- ifelse(abs(diag$Skewness) > 1, "Log/Sqrt recommended",
                         ifelse(abs(diag$Skewness) > 0.5, "Slightly skewed", "OK"))
    datatable(diag, options = list(dom = 't', pageLength = 20), caption = "Feature diagnostics")
  })

  # === DISTRIBUTIONS TAB ===
  dist_plot_obj <- reactive({
    req(transformed_data(), input$measure, input$dist_type)
    df  <- transformed_data(); m <- input$measure; rs <- rep_summ(); cols <- cond_colors()
    mlbl <- measure_label()

    # ── Apply sample filter (checkboxes) ────────────────────────────────────
    sel_conds <- input$dist_cond_sel
    if (!is.null(sel_conds) && length(sel_conds) > 0) {
      df <- df[as.character(df$Condition) %in% sel_conds, , drop = FALSE]
      if (!is.null(rs)) rs <- rs[as.character(rs$Condition) %in% sel_conds, , drop = FALSE]
    }

    # ── Apply timepoint filter (per-timepoint checkboxes) ───────────────────
    all_times <- sort(unique(suppressWarnings(as.numeric(transformed_data()$Time_min))))
    all_times <- all_times[!is.na(all_times)]
    sel_times <- all_times[sapply(all_times, function(t) {
      v <- input[[paste0("dist_time_", t)]]
      is.null(v) || isTRUE(v)          # default TRUE if not yet rendered
    })]
    if (length(sel_times) > 0) {
      df <- df[df$Time_min %in% sel_times, , drop = FALSE]
      if (!is.null(rs)) rs <- rs[rs$Time_min %in% sel_times, , drop = FALSE]
    }

    if (nrow(df) == 0) {
      return(ggplot() +
               annotate("text", x=0.5, y=0.5, label="No data — check sample / timepoint filters",
                        size=5, color="grey60") +
               theme_void())
    }

    # Drop unused factor levels after filtering
    if (is.factor(df$Condition)) df$Condition <- droplevels(df$Condition)

    # Recompute cols after possible subsetting
    cond_levels_now <- if (is.factor(df$Condition)) levels(df$Condition) else unique(df$Condition)
    n_now  <- length(cond_levels_now)
    dcols2 <- if (n_now <= 2) c("#4A90D9","#D94A4A") else scales::hue_pal()(n_now)
    cols   <- sapply(seq_along(cond_levels_now), function(i) {
      v <- if (has_colourpicker) input[[paste0("color_", i)]] else NULL
      if (!is.null(v) && nchar(v) > 0) v else dcols2[i]
    })
    cols <- setNames(cols, cond_levels_now)

    df$Time_fac <- factor(fmt_time(df$Time_min),
                          levels = fmt_time(sort(unique(df$Time_min))))
    if (!is.null(rs) && nrow(rs) > 0)
      rs$Time_fac <- factor(fmt_time(rs$Time_min),
                            levels = fmt_time(sort(unique(rs$Time_min))))

    # Grouping mode: by_time = X is time, fill=Condition; by_condition = X is Condition, fill=Time
    grp_mode <- input$dist_grouping_mode %||% "by_time"
    x_is_fill <- FALSE
    if (!has_time()) {
      # No time dimension: plot the samples themselves side by side.
      df$x_var_col <- df$Condition
      df$fill_var  <- df$Condition
      rs$x_var_col <- rs$Condition
      rs$fill_var  <- rs$Condition
      x_levels     <- levels(df$Condition)
      fill_levels  <- levels(df$Condition)
      x_label_base <- "Sample"
      time_fill_cols <- cols
      x_is_fill <- TRUE   # x axis already names each sample; legend would repeat it
    } else if (grp_mode == "by_condition") {
      # Swap X and fill: X=Condition, fill/dodge=Time_fac
      # We temporarily rename for plotting so all downstream code sees Time_fac as the dodge variable
      df$x_var_col  <- df$Condition
      df$fill_var   <- df$Time_fac
      rs$x_var_col  <- rs$Condition
      rs$fill_var   <- rs$Time_fac
      x_levels      <- levels(df$Condition)
      fill_levels   <- levels(df$Time_fac)
      x_label_base  <- "Condition"
      time_fill_cols <- setNames(
        scales::viridis_pal(option="D")(length(fill_levels)),
        fill_levels
      )
    } else {
      df$x_var_col <- df$Time_fac
      df$fill_var  <- df$Condition
      rs$x_var_col <- rs$Time_fac
      rs$fill_var  <- rs$Condition
      x_levels     <- levels(df$Time_fac)
      fill_levels  <- levels(df$Condition)
      x_label_base <- "Time Post-Infection"
      time_fill_cols <- cols
    }

    sample_x_mode <- isTRUE(input$dist_sample_x_mode)
    df$x_group <- as.character(df$x_var_col)
    df$fill_group <- as.character(df$fill_var)
    rs$x_group <- as.character(rs$x_var_col)
    rs$fill_group <- as.character(rs$fill_var)
    combo_df <- unique(df[, c("x_group", "fill_group")])
    combo_df$x_group <- factor(combo_df$x_group, levels = as.character(x_levels))
    combo_df$fill_group <- factor(combo_df$fill_group, levels = as.character(fill_levels))
    combo_df <- combo_df[order(combo_df$x_group, combo_df$fill_group), , drop = FALSE]
    sample_levels <- paste0(as.character(combo_df$x_group), " | ", as.character(combo_df$fill_group))
    df$Sample_Label <- factor(paste0(df$x_group, " | ", df$fill_group), levels = sample_levels)
    rs$Sample_Label <- factor(paste0(rs$x_group, " | ", rs$fill_group), levels = sample_levels)

    # Compute group means/medians for bars
    stat_type <- input$dist_stat_type %||% "mean"
    stat_col <- if (stat_type == "mean") "mean_val" else "median_val"
    group_stats <- rs %>% group_by(Condition, Time_fac) %>%
      summarise(stat_val = mean(.data[[stat_col]]), .groups = "drop")

    # --- SuperPlot options: make biological replicates visually distinguishable ---
    rep_encode   <- input$dist_rep_encode %||% "none"      # summary dots: none | shape
    cells_by_rep <- isTRUE(input$dist_cells_by_rep)        # color raw cells by replicate
    df$Rep_fac <- factor(df$Replicate)
    rs$Rep_fac <- factor(rs$Replicate)
    n_reps_lvl <- max(1, nlevels(rs$Rep_fac))
    rep_shapes <- vapply(seq_len(n_reps_lvl), function(i) {
      v <- suppressWarnings(as.numeric(input[[paste0("dist_rep_shape_", i)]]))
      if (length(v) == 1 && !is.na(v)) v else DEFAULT_REP_SHAPES[((i - 1) %% length(DEFAULT_REP_SHAPES)) + 1]
    }, numeric(1))
    cell_shape <- suppressWarnings(as.numeric(input$dist_cell_shape %||% 19))
    if (is.na(cell_shape)) cell_shape <- 19
    # Okabe-Ito colorblind-safe palette for replicate coloring
    okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
                   "#D55E00", "#CC79A7", "#999999")
    rep_cols <- setNames(rep(okabe_ito, length.out = max(1, nlevels(df$Rep_fac))),
                         levels(df$Rep_fac))

    # --- Build base plot ---
    if (input$dist_type == "violin") {
      v_alpha  <- input$violin_alpha      %||% 0.70
      v_trim   <- isTRUE(input$violin_trim)
      v_linew  <- if (isTRUE(input$violin_outline)) (input$violin_line_w %||% 0.3) else 0
      v_col    <- if (isTRUE(input$violin_outline)) "grey30" else NA
      v_wscale <- input$violin_width_scale %||% 1.0
      rep_sz   <- input$rep_dot_size  %||% 3.5
      rep_al   <- input$rep_dot_alpha %||% 0.9
      rep_col  <- input$rep_dot_color %||% "black"

      if (sample_x_mode) {
        p <- ggplot(df, aes(x = Sample_Label, y = .data[[m]], fill = fill_var, group = Sample_Label)) +
          geom_violin(alpha = v_alpha, scale = "width", linewidth = v_linew, color = v_col,
                      trim = v_trim, width = v_wscale) +
          scale_fill_manual(values = time_fill_cols, guide = "none")
        if (isTRUE(input$show_cell_points)) {
          if (cells_by_rep) {
            p <- p + geom_point(aes(color = Rep_fac), position = position_jitter(width = 0.15, height = 0),
                                size = 0.3, alpha = input$point_alpha, shape = cell_shape) +
              scale_color_manual(values = rep_cols, name = "Replicate",
                                 guide = guide_legend(override.aes = list(size = 3, alpha = 1)))
          } else {
            p <- p + geom_point(position = position_jitter(width = 0.15, height = 0),
                                size = 0.3, alpha = input$point_alpha, color = "grey30", shape = cell_shape)
          }
        }
        if (isTRUE(input$show_group_mean))
          p <- p + stat_summary(aes(group = Sample_Label), fun = if (stat_type == "mean") mean else median,
                                geom = "crossbar", width = 0.35, linewidth = 0.5, color = "black", alpha = 0.6)
        if (isTRUE(input$show_rep_medians)) {
          if (rep_encode == "shape") {
            p <- p + geom_point(data = rs, aes(x = Sample_Label, y = .data[[stat_col]], fill = fill_var, shape = Rep_fac),
                                size = rep_sz, color = rep_col, stroke = 0.8, alpha = rep_al) +
              scale_shape_manual(values = rep_shapes, name = "Replicate")
          } else {
            p <- p + geom_point(data = rs, aes(x = Sample_Label, y = .data[[stat_col]], fill = fill_var),
                                shape = 21, size = rep_sz, color = rep_col, stroke = 0.8, alpha = rep_al, show.legend = FALSE)
          }
        }
        p <- p + labs(x = "Sample", y = mlbl,
                      title = paste("Single-Cell", mlbl, "Distributions"),
                      subtitle = paste0("Violins = cell distributions; dots = replicate ", stat_type, "s; bar = group ", stat_type, ". Sample labels shown on X axis."))
      } else {
        p <- ggplot(df, aes(x = x_var_col, y = .data[[m]], fill = fill_var)) +
          geom_violin(position = position_dodge(0.8), alpha = v_alpha,
                      scale = "width", linewidth = v_linew, color = v_col,
                      trim = v_trim, width = v_wscale) +
          scale_fill_manual(values = time_fill_cols)
        if (isTRUE(input$show_cell_points)) {
          if (cells_by_rep) {
            p <- p + geom_point(aes(color = Rep_fac), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
                                size = 0.3, alpha = input$point_alpha, shape = cell_shape) +
              scale_color_manual(values = rep_cols, name = "Replicate",
                                 guide = guide_legend(override.aes = list(size = 3, alpha = 1)))
          } else {
            p <- p + geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
                                size = 0.3, alpha = input$point_alpha, color = "grey30", shape = cell_shape)
          }
        }
        if (isTRUE(input$show_group_mean))
          p <- p + stat_summary(aes(group = fill_var), fun = if (stat_type == "mean") mean else median,
                                geom = "crossbar", position = position_dodge(0.8), width = 0.35,
                                linewidth = 0.5, color = "black", alpha = 0.6)
        if (isTRUE(input$show_rep_medians)) {
          if (rep_encode == "shape") {
            p <- p + geom_point(data = rs, aes(x = x_var_col, y = .data[[stat_col]], fill = fill_var, group = fill_var, shape = Rep_fac),
                                position = position_dodge(0.8),
                                size = rep_sz, color = rep_col, stroke = 0.8, alpha = rep_al) +
              scale_shape_manual(values = rep_shapes, name = "Replicate")
          } else {
            p <- p + geom_point(data = rs, aes(x = x_var_col, y = .data[[stat_col]], fill = fill_var, group = fill_var),
                                position = position_dodge(0.8), shape = 21,
                                size = rep_sz, color = rep_col, stroke = 0.8, alpha = rep_al)
          }
        }
        p <- p + labs(x = x_label_base, y = mlbl,
                      title = paste("Single-Cell", mlbl, "Distributions"),
                      subtitle = paste0("Violins = cell distributions; dots = replicate ", stat_type, "s; bar = group ", stat_type))
      }

    } else if (input$dist_type == "ridge" && has_ridges) {
      if (has_time()) {
        df$Group <- interaction(df$Condition, df$Time_fac, sep = " | ")
        group_levels <- rev(as.character(interaction(
          rep(levels(df$Condition), length(levels(df$Time_fac))),
          rep(levels(df$Time_fac), each = length(levels(df$Condition))), sep = " | ")))
      } else {
        df$Group <- as.character(df$Condition)
        group_levels <- rev(levels(df$Condition))
      }
      df$Group <- factor(df$Group, levels = group_levels)
      p <- ggplot(df, aes(x = .data[[m]], y = Group, fill = Condition)) +
        ggridges::geom_density_ridges(alpha = 0.7, scale = 1.3, linewidth = 0.3) +
        scale_fill_manual(values = cols) + labs(x = mlbl, y = NULL, title = paste("Ridge Plot --", mlbl))

    } else if (input$dist_type == "strip") {
      if (sample_x_mode) {
        if (cells_by_rep) {
          p <- ggplot(df, aes(x = Sample_Label, y = .data[[m]], color = Rep_fac)) +
            geom_point(position = position_jitter(width = 0.2, height = 0),
                       size = 0.6, alpha = 0.3) +
            scale_color_manual(values = rep_cols, name = "Replicate",
                               guide = guide_legend(override.aes = list(size = 3, alpha = 1)))
        } else {
          p <- ggplot(df, aes(x = Sample_Label, y = .data[[m]], color = fill_var)) +
            geom_point(position = position_jitter(width = 0.2, height = 0),
                       size = 0.6, alpha = 0.15, show.legend = FALSE) +
            scale_color_manual(values = time_fill_cols, guide = "none")
        }
        if (isTRUE(input$show_group_mean))
          p <- p + stat_summary(aes(group = Sample_Label), fun = if (stat_type == "mean") mean else median,
                                geom = "crossbar", width = 0.3, linewidth = 0.5, color = "black", alpha = 0.7)
        if (isTRUE(input$show_rep_medians)) {
          if (rep_encode == "shape") {
            p <- p + geom_point(data = rs, aes(x = Sample_Label, y = .data[[stat_col]], shape = Rep_fac),
                                size = 3.5, fill = "white", color = "black", stroke = 0.8) +
              scale_shape_manual(values = rep_shapes, name = "Replicate")
          } else {
            p <- p + geom_point(data = rs, aes(x = Sample_Label, y = .data[[stat_col]]),
                                shape = 21, size = 3.5, fill = "white", color = "black", stroke = 0.8, show.legend = FALSE)
          }
        }
        p <- p + labs(x = "Sample", y = mlbl,
                      title = paste("Single-Cell", mlbl, "-- Jitter Plot"),
                      subtitle = paste0("Each point = one cell; bar = group ", stat_type, "; circles = replicate ", stat_type, "s."))
      } else {
        if (cells_by_rep) {
          p <- ggplot(df, aes(x = Time_fac, y = .data[[m]])) +
            geom_point(aes(color = Rep_fac, group = Condition),
                       position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7),
                       size = 0.6, alpha = 0.3, shape = cell_shape) +
            scale_color_manual(values = rep_cols, name = "Replicate",
                               guide = guide_legend(override.aes = list(size = 3, alpha = 1)))
        } else {
          p <- ggplot(df, aes(x = Time_fac, y = .data[[m]], color = Condition)) +
            geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7),
                       size = 0.6, alpha = 0.15, shape = cell_shape) + scale_color_manual(values = cols)
        }
        if (isTRUE(input$show_group_mean))
          p <- p + stat_summary(aes(group = Condition), fun = if (stat_type == "mean") mean else median,
                                geom = "crossbar", position = position_dodge(0.7), width = 0.3,
                                linewidth = 0.5, color = "black", alpha = 0.7)
        if (isTRUE(input$show_rep_medians)) {
          if (rep_encode == "shape") {
            p <- p + geom_point(data = rs, aes(x = Time_fac, y = .data[[stat_col]], group = Condition, shape = Rep_fac),
                                position = position_dodge(0.7), size = 3.5, fill = "white", color = "black", stroke = 0.8) +
              scale_shape_manual(values = rep_shapes, name = "Replicate")
          } else {
            p <- p + geom_point(data = rs, aes(x = Time_fac, y = .data[[stat_col]], group = Condition),
                                position = position_dodge(0.7), shape = 21, size = 3.5, fill = "white", color = "black", stroke = 0.8)
          }
        }
        p <- p + labs(x = "Time Post-Infection", y = mlbl,
                      title = paste("Single-Cell", mlbl, "-- Jitter Plot"),
                      subtitle = paste0("Each point = one cell; bar = group ", stat_type, "; circles = replicate ", stat_type, "s"))
      }

    } else if (input$dist_type == "sina") {
      has_ggforce <- requireNamespace("ggforce", quietly = TRUE)
      if (sample_x_mode) {
        cell_aes  <- if (cells_by_rep) aes(color = Rep_fac) else aes(color = fill_var)
        cell_alfa <- if (cells_by_rep) 0.3 else 0.15
        p <- ggplot(df, aes(x = Sample_Label, y = .data[[m]])) + {
          if (has_ggforce) ggforce::geom_sina(cell_aes, size = 0.5, alpha = cell_alfa,
                                              maxwidth = 0.35, show.legend = cells_by_rep,
                                              shape = cell_shape)
          else geom_point(cell_aes, position = position_jitter(width = 0.15, height = 0),
                          size = 0.5, alpha = cell_alfa, show.legend = cells_by_rep,
                          shape = cell_shape)
        }
        if (cells_by_rep) {
          p <- p + scale_color_manual(values = rep_cols, name = "Replicate",
                                      guide = guide_legend(override.aes = list(size = 3, alpha = 1)))
        } else {
          p <- p + scale_color_manual(values = time_fill_cols, guide = "none")
        }
        if (isTRUE(input$show_group_mean))
          p <- p + stat_summary(aes(group = Sample_Label), fun = if (stat_type == "mean") mean else median,
                                geom = "crossbar", width = 0.3, linewidth = 0.5, color = "black", alpha = 0.6)
        if (isTRUE(input$show_rep_medians)) {
          if (rep_encode == "shape") {
            p <- p + geom_point(data = rs, aes(x = Sample_Label, y = .data[[stat_col]], shape = Rep_fac),
                                size = 3.5, fill = "white", color = "black", stroke = 0.8) +
              scale_shape_manual(values = rep_shapes, name = "Replicate")
          } else {
            p <- p + geom_point(data = rs, aes(x = Sample_Label, y = .data[[stat_col]]),
                                shape = 21, size = 3.5, fill = "white", color = "black", stroke = 0.8, show.legend = FALSE)
          }
        }
        p <- p + labs(x = "Sample", y = mlbl,
                      title = paste("Single-Cell", mlbl, "-- Sina Plot"),
                      subtitle = paste0("Density-spread points; dots = replicate ", stat_type, "s; bar = group ", stat_type))
      } else {
        if (cells_by_rep) {
          p <- ggplot(df, aes(x = Time_fac, y = .data[[m]])) + {
            if (has_ggforce) ggforce::geom_sina(aes(color = Rep_fac, group = interaction(Time_fac, Condition)),
                                                position = position_dodge(0.8), size = 0.5, alpha = 0.3,
                                                maxwidth = 0.35, shape = cell_shape)
            else geom_point(aes(color = Rep_fac, group = Condition),
                            position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
                            size = 0.5, alpha = 0.3, shape = cell_shape)
          } + scale_color_manual(values = rep_cols, name = "Replicate",
                                 guide = guide_legend(override.aes = list(size = 3, alpha = 1)))
        } else {
          if (has_ggforce) {
            p <- ggplot(df, aes(x = Time_fac, y = .data[[m]], color = Condition)) +
              ggforce::geom_sina(position = position_dodge(0.8), size = 0.5, alpha = 0.15,
                                 maxwidth = 0.35, shape = cell_shape) +
              scale_color_manual(values = cols)
          } else {
            p <- ggplot(df, aes(x = Time_fac, y = .data[[m]], color = Condition)) +
              geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
                         size = 0.5, alpha = 0.15, shape = cell_shape) +
              scale_color_manual(values = cols)
          }
        }
        if (isTRUE(input$show_group_mean))
          p <- p + stat_summary(aes(group = Condition), fun = if (stat_type == "mean") mean else median,
                                geom = "crossbar", position = position_dodge(0.8), width = 0.3,
                                linewidth = 0.5, color = "black", alpha = 0.6)
        if (isTRUE(input$show_rep_medians)) {
          if (rep_encode == "shape") {
            p <- p + geom_point(data = rs, aes(x = Time_fac, y = .data[[stat_col]], group = Condition, shape = Rep_fac),
                                position = position_dodge(0.8), size = 3.5, fill = "white", color = "black", stroke = 0.8) +
              scale_shape_manual(values = rep_shapes, name = "Replicate")
          } else {
            p <- p + geom_point(data = rs, aes(x = Time_fac, y = .data[[stat_col]], group = Condition),
                                position = position_dodge(0.8), shape = 21, size = 3.5, fill = "white", color = "black", stroke = 0.8)
          }
        }
        p <- p + labs(x = "Time Post-Infection", y = mlbl,
                      title = paste("Single-Cell", mlbl, "-- Sina Plot"),
                      subtitle = paste0("Density-spread points; dots = replicate ", stat_type, "s; bar = group ", stat_type))
      }

    } else if (input$dist_type == "dotplot") {
      if (sample_x_mode) {
        p <- ggplot(rs, aes(x = Sample_Label, y = .data[[stat_col]], color = fill_var)) +
          {if (rep_encode == "shape") geom_point(aes(shape = Rep_fac), size = 4, alpha = 0.85)
           else geom_point(size = 4, alpha = 0.85, show.legend = FALSE)} +
          {if (rep_encode == "shape") scale_shape_manual(values = rep_shapes, name = "Replicate") else NULL} +
          stat_summary(aes(group = Sample_Label), fun = if (stat_type == "mean") mean else median, geom = "crossbar",
                       width = 0.3, linewidth = 0.5, color = "black") +
          scale_color_manual(values = time_fill_cols, guide = "none") +
          labs(x = "Sample", y = paste(tools::toTitleCase(stat_type), mlbl, "(per replicate)"),
               title = paste(mlbl, "-- Replicate Summary"),
               subtitle = paste0("Each point = one replicate ", stat_type, "; crossbar = group ", stat_type))
      } else {
        p <- ggplot(rs, aes(x = Time_fac, y = .data[[stat_col]], color = Condition)) +
          {if (rep_encode == "shape") geom_point(aes(shape = Rep_fac, group = Condition), position = position_dodge(0.5), size = 4, alpha = 0.85)
           else geom_point(position = position_dodge(0.5), size = 4, alpha = 0.85)} +
          {if (rep_encode == "shape") scale_shape_manual(values = rep_shapes, name = "Replicate") else NULL} +
          stat_summary(aes(group = Condition), fun = if (stat_type == "mean") mean else median, geom = "crossbar",
                       position = position_dodge(0.5), width = 0.3, linewidth = 0.5) +
          scale_color_manual(values = cols) +
          labs(x = "Time Post-Infection", y = paste(tools::toTitleCase(stat_type), mlbl, "(per replicate)"),
               title = paste(mlbl, "-- Replicate Summary"),
               subtitle = paste0("Each point = one replicate ", stat_type, "; crossbar = group ", stat_type))
      }

    } else if (input$dist_type == "boxstrip") {
      box_w  <- input$box_notch_w %||% 0.18
      rep_sz <- input$rep_dot_size  %||% 3.5
      rep_al <- input$rep_dot_alpha %||% 0.9
      rep_col <- input$rep_dot_color %||% "black"
      if (sample_x_mode) {
        p <- ggplot(df, aes(x = Sample_Label, y = .data[[m]], fill = fill_var, group = Sample_Label)) +
          geom_boxplot(alpha = 0.6, outlier.shape = NA, linewidth = 0.4, width = box_w * 3, show.legend = FALSE) +
          scale_fill_manual(values = time_fill_cols, guide = "none")
        if (isTRUE(input$show_rep_medians)) {
          if (rep_encode == "shape") {
            p <- p + geom_point(data = rs, aes(x = Sample_Label, y = .data[[stat_col]], shape = Rep_fac),
                                size = rep_sz, fill = "white", color = rep_col, stroke = 0.8, alpha = rep_al) +
              scale_shape_manual(values = rep_shapes, name = "Replicate")
          } else {
            p <- p + geom_point(data = rs, aes(x = Sample_Label, y = .data[[stat_col]]),
                                shape = 18, size = rep_sz, color = rep_col, alpha = rep_al, show.legend = FALSE)
          }
        }
        p <- p + labs(x = "Sample", y = mlbl, title = paste("Box Plot --", mlbl),
                      subtitle = paste0("Replicate ", stat_type, "s overlaid"))
      } else {
        p <- ggplot(df, aes(x = Time_fac, y = .data[[m]], fill = Condition)) +
          geom_boxplot(position = position_dodge(0.8), alpha = 0.6, outlier.shape = NA,
                       linewidth = 0.4, width = box_w * 3) +
          scale_fill_manual(values = cols)
        if (isTRUE(input$show_rep_medians)) {
          if (rep_encode == "shape") {
            p <- p + geom_point(data = rs, aes(x = Time_fac, y = .data[[stat_col]], group = Condition, shape = Rep_fac),
                                position = position_dodge(0.8),
                                size = rep_sz, fill = "white", color = rep_col, stroke = 0.8, alpha = rep_al) +
              scale_shape_manual(values = rep_shapes, name = "Replicate")
          } else {
            p <- p + geom_point(data = rs, aes(x = Time_fac, y = .data[[stat_col]], group = Condition),
                                position = position_dodge(0.8), shape = 18,
                                size = rep_sz, color = rep_col, alpha = rep_al)
          }
        }
        p <- p + labs(x = "Time Post-Infection", y = mlbl, title = paste("Box Plot --", mlbl),
                      subtitle = paste0("Replicate ", stat_type, "s overlaid"))
      }

    } else if (input$dist_type == "ecdf") {
      p <- ggplot(df, aes(x = .data[[m]], color = Condition)) +
        stat_ecdf(linewidth = 0.8) + scale_color_manual(values = cols) +
        {if (has_time()) facet_wrap(~ Time_fac, nrow = 1) else NULL} +
        labs(x = mlbl, y = "Cumulative Proportion", title = paste("ECDF --", mlbl),
             subtitle = "Rightward shift = larger values")
    } else if (input$dist_type == "histogram") {
      bins <- input$dist_hist_bins %||% 40
      if (has_time()) {
        p <- ggplot(df, aes(x = .data[[m]], fill = Time_fac)) +
          geom_histogram(bins = bins, alpha = 0.7, position = "identity", linewidth = 0.2, color = "white") +
          facet_wrap(~ Condition, scales = "free_y") +
          scale_fill_viridis_d(option = "plasma", name = "Time") +
          labs(x = mlbl, y = "Cell Count", title = paste("Histogram --", mlbl),
               subtitle = "Faceted by condition; colored by timepoint")
      } else {
        p <- ggplot(df, aes(x = .data[[m]], fill = Condition)) +
          geom_histogram(bins = bins, alpha = 0.7, position = "identity", linewidth = 0.2, color = "white") +
          facet_wrap(~ Condition, scales = "free_y") +
          scale_fill_manual(values = cols, guide = "none") +
          labs(x = mlbl, y = "Cell Count", title = paste("Histogram --", mlbl),
               subtitle = "One panel per sample")
      }

    } else if (input$dist_type == "stacked_bar") {
      breaks_q <- quantile(df[[m]], probs = seq(0, 1, 0.25), na.rm = TRUE)
      if (length(unique(breaks_q)) < 3) breaks_q <- quantile(df[[m]], probs = c(0, 0.5, 1), na.rm = TRUE)
      df$size_bin <- cut(df[[m]], breaks = unique(breaks_q), include.lowest = TRUE,
                         labels = c("Q1 (Small)", "Q2", "Q3", "Q4 (Large)")[1:(length(unique(breaks_q)) - 1)])
      freq <- df %>% group_by(Condition, Time_fac, size_bin) %>%
        summarise(n = n(), .groups = "drop") %>%
        group_by(Condition, Time_fac) %>% mutate(pct = n / sum(n) * 100) %>% ungroup()
      p <- ggplot(freq, aes(x = Time_fac, y = pct, fill = size_bin)) +
        geom_col(position = "stack", width = 0.75) + facet_wrap(~ Condition) +
        scale_fill_viridis_d(name = paste(mlbl, "Quartile")) +
        labs(x = "Time", y = "% of Cells", title = paste("Population Proportions --", mlbl),
             subtitle = "Cells binned into quartiles of the overall distribution") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

    } else if (input$dist_type == "pie") {
      breaks_q <- quantile(df[[m]], probs = seq(0, 1, 0.25), na.rm = TRUE)
      if (length(unique(breaks_q)) < 3) breaks_q <- quantile(df[[m]], probs = c(0, 0.5, 1), na.rm = TRUE)
      df$size_bin <- cut(df[[m]], breaks = unique(breaks_q), include.lowest = TRUE,
                         labels = c("Q1 (Small)", "Q2", "Q3", "Q4 (Large)")[1:(length(unique(breaks_q)) - 1)])
      freq <- df %>% group_by(Condition, Time_fac, size_bin) %>%
        summarise(n = n(), .groups = "drop") %>%
        group_by(Condition, Time_fac) %>% mutate(pct = n / sum(n) * 100) %>% ungroup()
      p <- ggplot(freq, aes(x = "", y = pct, fill = size_bin)) +
        geom_col(width = 1, color = "white", linewidth = 0.4) +
        coord_polar(theta = "y") +
        facet_grid(Condition ~ Time_fac,
                   labeller = labeller(Condition = function(x) wrap_axis_labels(x, 16))) +
        scale_fill_viridis_d(name = paste(mlbl, "Quartile")) +
        labs(x = NULL, y = NULL, title = paste("Cell Size Composition --", mlbl),
             subtitle = "Proportional pie chart; faceted by condition and time")

    } else if (input$dist_type == "donut") {
      breaks_q <- quantile(df[[m]], probs = seq(0, 1, 0.25), na.rm = TRUE)
      if (length(unique(breaks_q)) < 3) breaks_q <- quantile(df[[m]], probs = c(0, 0.5, 1), na.rm = TRUE)
      df$size_bin <- cut(df[[m]], breaks = unique(breaks_q), include.lowest = TRUE,
                         labels = c("Q1 (Small)", "Q2", "Q3", "Q4 (Large)")[1:(length(unique(breaks_q)) - 1)])
      freq <- df %>% group_by(Condition, Time_fac, size_bin) %>%
        summarise(n = n(), .groups = "drop") %>%
        group_by(Condition, Time_fac) %>% mutate(pct = n / sum(n) * 100) %>% ungroup()
      p <- ggplot(freq, aes(x = 2, y = pct, fill = size_bin)) +
        geom_col(width = 1, color = "white", linewidth = 0.4) +
        coord_polar(theta = "y") + xlim(0.5, 2.5) +
        facet_grid(Condition ~ Time_fac,
                   labeller = labeller(Condition = function(x) wrap_axis_labels(x, 16))) +
        scale_fill_viridis_d(name = paste(mlbl, "Quartile")) +
        labs(x = NULL, y = NULL, title = paste("Cell Composition --", mlbl, "(Donut)"),
             subtitle = "Proportional donut chart; faceted by condition and time")

    } else {
      p <- ggplot(df, aes(x = Time_fac, y = .data[[m]], fill = Condition)) +
        geom_violin(position = position_dodge(0.8), alpha = 0.7) +
        scale_fill_manual(values = cols) + labs(x = "Time", y = mlbl, title = paste(mlbl, "Distributions"))
    }

    # Human-readable legend titles: never leak internal names like "fill_var".
    # (Replicate legends keep their own "Replicate" title from their scales.)
    legend_title <- if (grp_mode == "by_condition") "Time" else "Condition"
    if (input$dist_type %in% c("violin", "boxstrip")) {
      p <- p + labs(fill = legend_title)
    } else if (input$dist_type %in% c("strip", "sina", "dotplot") && !cells_by_rep) {
      p <- p + labs(color = legend_title)
    }

    # --- Significance brackets ---
    # Test + correction follow the Statistics tab settings. The family does too,
    # but only when that tab is left on its default within-time configuration --
    # bracket_family_note() states it in the caption when it is not.
    sig_ymax <- NA_real_
    bracket_types <- c("violin", "strip", "sina", "boxstrip", "dotplot")
    if (isTRUE(input$show_sig) && !sample_x_mode && input$dist_type %in% bracket_types && length(unique(df$Condition)) >= 2) {
      cond_levels <- levels(df$Condition)
      ref_cond <- input$ref_condition %||% cond_levels[1]
      if (!ref_cond %in% cond_levels) ref_cond <- cond_levels[1]
      times_fac <- levels(df$Time_fac)
      sig_yvar <- input$sig_stat_unit %||% "median_val"
      compare_mode <- input$sig_compare_mode %||% "vs_control"
      dodge_w <- switch(input$dist_type, "dotplot" = 0.5, "strip" = 0.7, 0.8)
      y_vals <- if (input$dist_type == "dotplot") rs$median_val else df[[m]]
      y_range <- diff(range(y_vals, na.rm = TRUE)); if (y_range == 0) y_range <- 1

      # Build comparison pairs based on mode
      if (compare_mode == "all_pairwise") {
        pairs_list <- combn(cond_levels, 2, simplify = FALSE)
      } else {
        trt_conds <- setdiff(cond_levels, ref_cond)
        pairs_list <- lapply(trt_conds, function(trt) c(ref_cond, trt))
      }

      sig_res <- compute_sig_table(rs, pairs_list,
        times = sort(unique(rs$Time_min)),
        sig_yvar = sig_yvar,
        test = input$stat_test %||% "welch",
        correction = input$stat_p_adjust %||% "holm",
        label_style = input$sig_label_style %||% "stars",
        show_ns = isTRUE(input$show_ns),
        y_lookup = function(tm, pair, v1, v2) {
          if (input$dist_type == "dotplot") max(c(v1, v2), na.rm = TRUE)
          else max(df[[m]][df$Time_min == tm], na.rm = TRUE)
        })
      sig_annot <- sig_res$annot
      if (nrow(sig_annot) > 0) sig_annot$Time_fac <- fmt_time(sig_annot$Time_min)
      {
        if (nrow(sig_annot) > 0) {
          n_levels <- length(cond_levels)
          bracket_base <- sig_annot$y_pos + y_range * 0.06
          # Taller steps when labels are two lines (stars + p-value)
          bracket_step <- y_range * (if ((input$sig_label_style %||% "stars") == "stars_p") 0.17 else 0.10)
          # Stack brackets at each timepoint so they don't overlap
          for (tf in unique(sig_annot$Time_fac)) {
            idx <- which(sig_annot$Time_fac == tf)
            if (length(idx) > 1) {
              for (j in 2:length(idx)) bracket_base[idx[j]] <- bracket_base[idx[j - 1]] + bracket_step
            }
          }
          for (k in seq_len(nrow(sig_annot))) {
            x_center <- which(times_fac == sig_annot$Time_fac[k])
            cond1_idx <- which(cond_levels == sig_annot$Cond1[k]) - 1
            cond2_idx <- which(cond_levels == sig_annot$Cond2[k]) - 1
            if (isTRUE(x_is_fill)) {
              # Samples occupy whole x slots (no time dimension), so the bracket
              # spans between slots rather than between dodged sub-positions.
              x_left  <- cond1_idx + 1
              x_right <- cond2_idx + 1
            } else if (n_levels == 2) {
              x_left  <- x_center - dodge_w / 4
              x_right <- x_center + dodge_w / 4
            } else {
              x_left  <- x_center + (cond1_idx - (n_levels - 1) / 2) * dodge_w / n_levels
              x_right <- x_center + (cond2_idx - (n_levels - 1) / 2) * dodge_w / n_levels
            }
            if (x_left > x_right) { tmp <- x_left; x_left <- x_right; x_right <- tmp }
            y_bar <- bracket_base[k]; y_tick <- y_bar - y_range * 0.015
            bc <- if (sig_annot$stars[k] == "ns") "grey60" else "black"
            lc <- if (sig_annot$stars[k] == "ns") "grey50" else "black"
            p <- p + annotate("segment", x = x_left, xend = x_left, y = y_tick, yend = y_bar, linewidth = 0.4, color = bc) +
              annotate("segment", x = x_left, xend = x_right, y = y_bar, yend = y_bar, linewidth = 0.4, color = bc) +
              annotate("segment", x = x_right, xend = x_right, y = y_tick, yend = y_bar, linewidth = 0.4, color = bc) +
              annotate("text", x = (x_left + x_right) / 2, y = y_bar + y_range * 0.015,
                       label = sig_annot$label[k], size = (input$dist_base_size %||% 12) / 2.8,
                       vjust = 0, color = lc, fontface = if (sig_annot$stars[k] != "ns") "bold" else "plain")
          }
          label_extra <- if (input$sig_label_style == "stars_p") 0.20 else 0.12
          sig_ymax <- max(bracket_base) + y_range * label_extra
        }
        p <- p + labs(caption = paste0(sig_res$caption, bracket_family_note() %||% ""))
      }
    }

    # Apply axis customizations (single coord_cartesian so user Y limits and
    # bracket headroom cooperate instead of overriding each other)
    if (isTRUE(input$dist_log_y)) {
      if (input$dist_type %in% c("violin", "strip", "sina", "boxstrip", "dotplot")) {
        p <- p + scale_y_log10()
      } else if (input$dist_type %in% c("ridge", "ecdf", "histogram")) {
        p <- p + scale_x_log10()  # measurement sits on the x axis for these
      }
    } else {
      y_min <- if (!is.null(input$dist_y_min) && !is.na(input$dist_y_min)) input$dist_y_min else NA
      y_max <- if (!is.null(input$dist_y_max) && !is.na(input$dist_y_max)) input$dist_y_max else NA
      upper <- suppressWarnings(max(c(y_max, sig_ymax), na.rm = TRUE))
      if (!is.finite(upper)) upper <- NA
      if (!is.na(y_min) || !is.na(upper)) {
        # clip="off" is only needed to let bracket labels sit above the panel;
        # otherwise keep clipping on so data cannot spill past the axis limits.
        p <- p + coord_cartesian(ylim = c(
          if (!is.na(y_min)) y_min else NA,
          if (!is.na(upper)) upper else NA
        ), clip = if (!is.na(sig_ymax)) "off" else "on")
      }
    }

    # Publication theme must be applied BEFORE the per-plot theme tweaks below:
    # it is a complete theme, so adding it later would wipe them.
    p <- p + get_theme(input, "dist")

    # Readable x-axis labels for crowded experiments (discrete-x plots only --
    # ridge/ECDF/histogram put the measurement on a continuous x axis).
    # Wrapping onto multiple lines only helps horizontal labels; angled
    # multi-line labels collide with their neighbors, so wrap only at 0 deg.
    x_angle <- suppressWarnings(as.numeric(input$dist_x_angle %||% 45))
    x_hjust <- ifelse(x_angle == 0, 0.5, 1)
    if (input$dist_type %in% c("violin", "strip", "sina", "boxstrip", "dotplot", "stacked_bar")) {
      # Horizontal labels sit side by side, so the more x groups there are the
      # narrower each label may be before neighbours collide. Angled labels do
      # not have this constraint (and must stay on one line).
      n_x <- {
        xc <- if (sample_x_mode) "Sample_Label" else if (input$dist_type == "violin") "x_var_col" else "Time_fac"
        src <- if (input$dist_type == "dotplot") rs else df
        if (xc %in% names(src)) length(unique(as.character(src[[xc]]))) else 1
      }
      wrap_w <- if (x_angle == 0)
        max(6, min(input$dist_label_wrap %||% 18, floor(130 / max(1, n_x)))) else 500
      p <- p + scale_x_discrete(labels = function(x) {
        lab <- wrap_axis_labels(x, wrap_w)
        if (!is.null(n_lookup)) {
          extra <- unname(n_lookup[as.character(x)])
          lab <- ifelse(is.na(extra), lab, paste0(lab, "\n", extra))
        }
        lab
      })
    }
    left_margin <- if (sample_x_mode && x_angle > 0) 40 else 12
    if (input$dist_type %in% c("pie", "donut")) {
      # Polar charts have no meaningful cartesian axes. get_theme() is a complete
      # theme, so this blanking has to come after it, not in the geom branch.
      p <- p + theme(axis.text = element_blank(), axis.ticks = element_blank(),
                     axis.line = element_blank(), panel.grid = element_blank(),
                     strip.text.y = element_text(angle = 0, hjust = 0, size = rel(0.7),
                                                 lineheight = 0.9),
                     strip.text.x = element_text(size = rel(0.68)),
                     panel.spacing = unit(2, "pt"),
                     plot.margin = margin(12, 12, 12, 12))
    } else if (input$dist_type %in% c("ridge", "ecdf", "histogram")) {
      # x is the continuous measurement here, so leave its numbers horizontal.
      # Ridge y labels are one per group, so shrink them when groups are many.
      n_ygroups <- if (input$dist_type == "ridge") length(unique(df$Condition)) * length(unique(df$Time_min)) else 0
      y_rel <- if (n_ygroups > 24) 0.6 else if (n_ygroups > 14) 0.72 else 1
      p <- p + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
                     axis.text.y = element_text(size = rel(y_rel)),
                     plot.margin = margin(12, 12, 12, 12))
    } else {
      p <- p + theme(axis.text.x = element_text(angle = x_angle, hjust = x_hjust, vjust = 1),
              plot.margin = margin(12, 12, input$dist_bottom_margin %||% 18, left_margin))
    }

    # Optional per-group n annotation. Appended to the x-axis tick labels rather
    # than drawn in the panel, so it can never sit on top of the data.
    n_lookup <- NULL
    if (isTRUE(input$dist_show_n) && input$dist_type %in% c("violin", "strip", "sina", "boxstrip", "dotplot")) {
      x_col <- if (sample_x_mode) "Sample_Label" else if (input$dist_type == "violin") "x_var_col" else "Time_fac"
      n_src <- if (input$dist_type == "dotplot") rs else df
      if (x_col %in% names(n_src)) {
        n_df <- n_src %>% group_by(.data[[x_col]]) %>%
          summarise(n_cells = n(), n_reps = n_distinct(Replicate), .groups = "drop")
        names(n_df)[1] <- "x_grp"
        n_lookup <- setNames(
          if (input$dist_type == "dotplot") paste0("N=", n_df$n_reps)
          else paste0("n=", n_df$n_cells, ", N=", n_df$n_reps),
          as.character(n_df$x_grp))
      }
    }

    # Keep the replicate legend visible in sample-label mode when SuperPlot
    # encoding is active; otherwise hide the (redundant) fill legend as before.
    if (sample_x_mode && rep_encode == "none" && !cells_by_rep) p <- p + theme(legend.position = "none")
    # When the x axis already names each sample, its colour legend is redundant
    # (the Replicate legend, which is a shape scale, is kept).
    if (isTRUE(x_is_fill) && !cells_by_rep) p <- p + guides(fill = "none", color = "none")

    apply_custom_labs(p, input, "dist")
  })

  output$dist_plot <- renderPlot({ dist_plot_obj() })
  output$dl_dist <- downloadHandler(
    filename = function() paste0("distributions_", input$measure, ".", fmt_ext(input$dist_format)),
    content  = function(file) save_figure(file, dist_plot_obj(), input$dist_width, input$dist_height, input$dist_dpi, input$dist_format))

  # --- Keep the figure brackets and the Statistics table in one family --------
  # The brackets always compare conditions WITHIN each timepoint, on replicate
  # summaries, with the Statistics tab's test and correction. The Statistics tab
  # builds its own family from the grouping columns and comparison mode, so the
  # two agree only when it is set to "Pairwise within each time point" over
  # Condition + Time_min -- which is now the shipped default. When it is set to
  # anything else the caption says so, rather than leaving the reader to assume
  # a bracket and a table row with different p-values are the same contrast
  # (audit CM-6; the README used to claim the two "can never disagree").
  bracket_family_note <- function() {
    mode <- input$stat_compare_mode %||% "within_time"
    cols <- input$stat_group_cols %||% c("Condition", "Time_min")
    matches <- identical(mode, "within_time") && all(c("Condition", "Time_min") %in% cols)
    if (matches) return(NULL)
    paste0(" The Statistics tab is currently set to \"", mode, "\" over ",
           paste(cols, collapse = " + "),
           ", which is a different family of comparisons: its adjusted p-values will not match these brackets. ",
           "Set it to \"Pairwise within each time point\" on Condition + Time_min for the two to agree.")
  }

  # === TRAJECTORY TAB ===
  traj_plot_obj <- reactive({
    req(rep_summ(), cond_summ(), input$traj_summary)
    validate(need(has_time(),
      "This dataset has no time column, so there is no trajectory to plot. Use the Distributions tab to compare samples directly."))
    rs <- rep_summ(); cs <- cond_summ(); cols <- cond_colors(); yvar <- input$traj_summary; mlbl <- measure_label()
    ylab <- if (yvar == "median_val") paste("Median", mlbl, "per Replicate") else paste("Mean", mlbl, "per Replicate")
    p <- ggplot()
    if (input$traj_error && input$traj_style == "ribbon") {
      cs$yvar <- cs[[ifelse(yvar == "median_val", "grand_median", "grand_mean")]]
      cs$sem  <- cs[[ifelse(yvar == "median_val", "sem_median", "sem_mean")]]
      p <- p + geom_ribbon(data = cs, aes(x = Time_min, ymin = yvar - sem, ymax = yvar + sem, fill = Condition), alpha = 0.2) +
        scale_fill_manual(values = cols)
    }
    traj_pt_sz <- input$traj_point_size %||% 3
    traj_ln_w  <- input$traj_line_width %||% 1.2
    if (input$traj_connect_reps) p <- p + geom_line(data = rs, aes(x = Time_min, y = .data[[yvar]], group = interaction(Condition, Replicate), color = Condition), alpha = 0.4, linewidth = traj_ln_w * 0.4, linetype = "dashed")
    p <- p + geom_point(data = rs, aes(x = Time_min, y = .data[[yvar]], color = Condition, shape = Replicate), size = traj_pt_sz, alpha = 0.8)
    if (input$traj_smooth) {
      n_times <- length(unique(rs$Time_min))
      if (n_times >= 4) p <- p + geom_smooth(data = rs, aes(x = Time_min, y = .data[[yvar]], color = Condition), method = "loess", se = input$traj_error && input$traj_style == "pointline", linewidth = traj_ln_w, alpha = 0.15, span = 1)
      else p <- p + geom_smooth(data = rs, aes(x = Time_min, y = .data[[yvar]], color = Condition), method = "lm", se = input$traj_error && input$traj_style == "pointline", linewidth = traj_ln_w, alpha = 0.15)
    }
    if (input$traj_error && input$traj_style == "pointline" && !input$traj_smooth) {
      cs$yvar <- cs[[ifelse(yvar == "median_val", "grand_median", "grand_mean")]]
      cs$sem  <- cs[[ifelse(yvar == "median_val", "sem_median", "sem_mean")]]
      p <- p + geom_errorbar(data = cs, aes(x = Time_min, ymin = yvar - sem, ymax = yvar + sem, color = Condition), width = 3, linewidth = 0.6)
    }
    p <- p + scale_color_manual(values = cols) + scale_x_continuous(breaks = sort(unique(rs$Time_min))) +
      labs(x = "Time Post-Infection (min)", y = ylab, title = paste(mlbl, "Trajectory Over Time"),
           subtitle = "Each point = one replicate; dashed lines connect replicates") +
      get_theme(input, "traj")
    apply_custom_labs(p, input, "traj")
  })

  output$traj_plot <- renderPlot({ traj_plot_obj() })
  output$dl_traj <- downloadHandler(
    filename = function() paste0("trajectory_", input$measure, ".", fmt_ext(input$traj_format)),
    content  = function(file) save_figure(file, traj_plot_obj(), input$traj_width, input$traj_height, input$traj_dpi, input$traj_format))

  # === STATISTICS TAB ===
  stat_results <- reactiveVal(NULL)

  observeEvent(input$run_stats, {
    req(flexible_rep_summ(), input$stat_summary)
    rs_flex <- flexible_rep_summ()
    yvar <- input$stat_summary %||% "median_val"
    m <- input$stat_measure %||% input$measure
    results <- list()
    results$method <- paste0("Flexible pairwise statistics on biological replicate-level ",
                             ifelse(yvar == "median_val", "medians", "means"))
    results$measurement <- m
    results$comparison_mode_requested <- input$stat_compare_mode %||% "within_time"
    # Report the mode that was actually built, never the one that was asked for.
    results$comparison_mode <- realised_comparison_mode(comparison_preview(),
                                                       results$comparison_mode_requested)
    results$comparison_mode_message <- attr(comparison_preview(), "mode_message")
    if (!is.null(results$comparison_mode_message)) {
      showNotification(results$comparison_mode_message, type = "warning", duration = 12)
    }
    results$reference_group <- input$stat_reference_group %||% NA_character_
    results$transform_used  <- input$data_transform %||% "none"
    results$transform_label <- describe_transform(transformed_data(), m, input$data_transform)
    results$statistical_unit <- ifelse(yvar == "median_val", "replicate median", "replicate mean")
    results$grouping_variables <- paste(input$stat_group_cols %||% "Condition", collapse = ", ")
    results$p_adjust <- input$stat_p_adjust %||% "holm"
    results$test <- input$stat_test %||% "welch"
    results$test_suggestion <- suggest_statistical_workflow(
      rep_df = rs_flex,
      comparison_tbl = comparison_preview(),
      comparison_mode = results$comparison_mode,
      stat_unit = yvar,
      selected_group_cols = input$stat_group_cols
    )

    desc <- rs_flex %>% group_by(Statistics_Group) %>%
      summarise(n_reps = n_distinct(Replicate),
                n_cells = sum(n_cells, na.rm = TRUE),
                Mean   = round(mean(.data[[yvar]], na.rm = TRUE), 3),
                SD     = round(sd(.data[[yvar]], na.rm = TRUE), 3),
                SEM    = round(sd(.data[[yvar]], na.rm = TRUE) / sqrt(n()), 3),
                Median = round(median(.data[[yvar]], na.rm = TRUE), 3),
                IQR    = round(IQR(.data[[yvar]], na.rm = TRUE), 3),
                Min    = round(min(.data[[yvar]], na.rm = TRUE), 3),
                Max    = round(max(.data[[yvar]], na.rm = TRUE), 3),
                CV_pct = round(sd(.data[[yvar]], na.rm = TRUE) / abs(mean(.data[[yvar]], na.rm = TRUE)) * 100, 1),
                .groups = "drop")
    results$descriptive <- desc
    results$yvar_used <- yvar

    comp_tbl <- comparison_preview()
    pw <- run_pairwise_tests(
      rs_flex, comp_tbl, yvar = "Stat_Value",
      test = input$stat_test %||% "welch",
      correction = input$stat_p_adjust %||% "holm",
      reference_group = input$stat_reference_group %||% NA_character_,
      measurement = m,
      comparison_mode = results$comparison_mode,
      stat_unit = results$statistical_unit
    )
    results$pairwise <- pw
    results$comparison_preview <- comp_tbl
    results$replicate_summary <- rs_flex
    results$pairwise_family_size <- nrow(pw)
    results$pairwise_family <- paste0(
      "Holm-type correction applied across all ", nrow(pw), " comparison",
      if (nrow(pw) == 1) "" else "s", " in this table (one family), using \"",
      results$p_adjust, "\".")
    # Advisory only -- see format_assumption_diagnostics().
    results$diagnostics <- tryCatch(assumption_diagnostics(rs_flex, yvar = "Stat_Value"),
                                    error = function(e) NULL)

    # --- Optional global model -------------------------------------------------
    # Both branches now test each term against the error stratum it belongs to;
    # see fit_rm_anova() and fit_lmm_safe(). Both run on transformed_data(), the
    # same scale as the pairwise tests above, so one report never mixes scales.
    gm <- input$stat_method %||% "anova"
    if (gm == "anova" && all(c("Condition", "Time_min") %in% names(rs_flex))) {
      fit <- fit_rm_anova(rs_flex, yvar = yvar)
      results$global_model <- fit
      results$global_model_kind <- "rm_anova"
      if (isTRUE(fit$ok)) {
        results$model <- fit$model
        results$anova_text <- format_rm_anova(fit)
        results$partial_eta2 <- data.frame(
          Source = fit$table$Term,
          `Tested against` = fit$table$Error_Stratum,
          partial_eta2 = round(fit$table$partial_eta2, 3),
          check.names = FALSE)
        # Residuals of the equivalent fixed-effects fit, for the normality check.
        resid_fit <- tryCatch(residuals(lm(as.formula(paste0(yvar, " ~ Condition * Time_fac + Rep_unique")),
                                           data = fit$data)), error = function(e) NULL)
        results$diagnostics <- tryCatch(assumption_diagnostics(rs_flex, yvar = "Stat_Value",
                                                               residuals = resid_fit),
                                        error = function(e) results$diagnostics)
      } else {
        results$global_model_message <- fit$msg
      }
    } else if (gm == "anova") {
      results$global_model_kind <- "rm_anova"
      results$global_model_message <- paste0(
        "The repeated-measures ANOVA needs Condition and Time_min among the grouping ",
        "columns. Your groups are built from ", results$grouping_variables,
        ", so no global model was fitted. The pairwise results above are unaffected.")
    } else if (gm == "lmm") {
      results$global_model_kind <- "lmm"
      if (!has_lme4) {
        showNotification("lme4 is required for the mixed model. The pairwise results were still generated.",
                         type = "warning")
        results$global_model_message <- "lme4 is not installed, so no mixed model was fitted."
      } else {
        # CM-14: the mixed model reads the SAME frame the pairwise tests read, so
        # one report can no longer hold pairwise tests on log data and a model on
        # raw data. The Data-tab transform is the only transform.
        req(transformed_data(), m)
        lmm_result <- fit_lmm_safe(transformed_data(), m, transform = "none")
        results$lmm <- lmm_result
        results$lmm_warnings <- lmm_result$warnings
        results$method <- paste(results$method, "+ mixed model")
        if (isTRUE(lmm_result$fitted)) {
          results$model <- lmm_result$model
          results$anova_summary <- tryCatch(anova(lmm_result$model), error = function(e) NULL)
          results$lmm_formula <- paste(deparse(lmm_result$formula), collapse = " ")
          results$lmm_df_method <- lmm_result$df_method
          results$ranef_var <- tryCatch(as.data.frame(lme4::VarCorr(lmm_result$model)), error = function(e) NULL)
          if (has_emmeans) {
            # CM-8: `~ Condition | Time_fac` makes Time_fac a `by` variable, so
            # adjust="holm" corrected WITHIN each timepoint -- with two
            # conditions that is one contrast per stratum and no correction at
            # all, while the output was labelled "holm". Cross the factors so the
            # family is every Condition-within-Time contrast, which is the family
            # the documentation has always claimed.
            emm <- tryCatch(emmeans::emmeans(lmm_result$model, ~ Condition * Time_fac),
                            error = function(e) NULL)
            if (!is.null(emm)) {
              results$emmeans <- emm
              # Pull the contrasts unadjusted and apply the correction here, over
              # the whole table, so the family is unambiguous and matches the
              # correction the user chose for the pairwise tests. emmeans' own
              # `adjust=` is applied per `by` group, which is what produced the
              # "holm" label on uncorrected p-values before.
              cons <- tryCatch(as.data.frame(emmeans::contrast(
                        emm, method = "pairwise", simple = "Condition",
                        combine = TRUE, adjust = "none")),
                        error = function(e) NULL)
              if (!is.null(cons) && "p.value" %in% names(cons)) {
                corr <- results$p_adjust
                cons$p_raw <- cons$p.value
                cons$p_adj <- if (is.null(corr) || corr == "none") cons$p_raw
                              else p.adjust(cons$p_raw, method = corr)
                cons$p.value <- NULL
                cons$Significance <- make_significance_label(cons$p_adj)
                results$emm_family <- paste0(
                  "Condition contrasts at every timepoint, corrected as one family of ",
                  nrow(cons), " (\"", corr, "\"), on Kenward-Roger degrees of freedom.")
              }
              results$pairwise_emm <- cons
              results$emm_summary <- tryCatch(as.data.frame(summary(emm)), error = function(e) NULL)
            }
          }
        } else {
          results$global_model_message <- paste(lmm_result$warnings, collapse = " ")
        }
      }
    }

    stat_results(results)
    showNotification("Flexible statistics complete. Pairwise tests used biological replicate summaries, not raw cells.", type = "message")
  })

  # One report body, printed on screen and written to the downloaded .txt, so the
  # two cannot drift apart.
  print_stat_report <- function(r, include_tables = FALSE) {
    cat("=== METHOD ===\n", r$method, "\n\n")
    cat("Measurement:", r$measurement %||% input$stat_measure %||% input$measure, "\n")
    cat("Transform:", r$transform_label %||% "none (raw values)", "\n")
    cat("Grouping variables:", r$grouping_variables %||% "Condition", "\n")
    cat("Comparison mode:", r$comparison_mode %||% NA, "\n")
    if (!is.null(r$comparison_mode_message)) {
      cat("  NOTE:", r$comparison_mode_message, "\n")
    }
    cat("Reference group:", r$reference_group %||% NA, "\n")
    cat("Statistical unit:", r$statistical_unit %||% r$yvar_used, "\n")
    cat("Pairwise test:", r$test %||% NA, "\n")
    cat("P-value correction:", r$p_adjust %||% NA, "\n")
    if (!is.null(r$pairwise_family)) cat("Multiplicity family:", r$pairwise_family, "\n")
    cat("\n")
    if (!is.null(r$test_suggestion)) {
      cat("=== TEST SUGGESTION / RATIONALE ===\n")
      cat(r$test_suggestion$summary, "\n\n")
    }
    cat("Raw cells are visualized, but pairwise hypothesis tests use biological replicate-level summaries.\n\n")
    if (include_tables) {
      if (!is.null(r$descriptive)) { cat("=== DESCRIPTIVE ===\n"); print(as.data.frame(r$descriptive), row.names = FALSE); cat("\n") }
      if (!is.null(r$pairwise))    { cat("=== FLEXIBLE PAIRWISE RESULTS ===\n"); print(as.data.frame(r$pairwise), row.names = FALSE); cat("\n") }
    }
    cat("=== GLOBAL MODEL ===\n")
    if (!is.null(r$anova_text)) {
      cat(r$anova_text, "\n")
    } else if (!is.null(r$anova_summary)) {
      cat(r$lmm$method %||% "Mixed model", "\n")
      if (!is.null(r$lmm_formula)) cat("Formula:", r$lmm_formula, "\n")
      if (!is.null(r$lmm_df_method)) cat("Degrees of freedom:", r$lmm_df_method, "\n")
      print(r$anova_summary)
    } else {
      cat(r$global_model_message %||% "Not run, or not available for the selected grouping.\n")
    }
    if (!is.null(r$lmm_warnings) && length(r$lmm_warnings) > 0) {
      cat("\nModel notes:\n"); for (w in r$lmm_warnings) cat("  -", w, "\n")
    }
    if (!is.null(r$ranef_var)) { cat("\n=== RANDOM EFFECT VARIANCE ===\n"); print(r$ranef_var) }
    if (!is.null(r$partial_eta2)) {
      cat("\n=== EFFECT SIZES (partial eta-squared) ===\n")
      cat("small=0.01, medium=0.06, large=0.14\n")
      print(r$partial_eta2, row.names = FALSE)
    }
    if (!is.null(r$pairwise_emm)) {
      cat("\n=== PAIRWISE (estimated marginal means) ===\n")
      if (!is.null(r$emm_family)) cat(r$emm_family, "\n")
      print(r$pairwise_emm)
    }
    cat("\n=== ASSUMPTION DIAGNOSTICS ===\n")
    cat(format_assumption_diagnostics(r$diagnostics), "\n")
    cat("\nThe Pairwise tab is the primary flexible statistics output.\n")
  }

  output$stat_output <- renderPrint({
    req(stat_results()); print_stat_report(stat_results())
  })

  output$pairwise_table <- renderDT({
    req(stat_results()); r <- stat_results()
    if (!is.null(r$pairwise)) { pw <- r$pairwise; pw[] <- lapply(pw, function(x) if (is.numeric(x)) round(x, 4) else x)
      datatable(pw, options = list(pageLength = 10, dom = 't'))
    } else if (!is.null(r$pairwise_emm)) datatable(as.data.frame(summary(r$pairwise_emm)), options = list(pageLength = 10, dom = 't'))
  })

  output$stat_descriptive_table <- renderDT({
    req(stat_results()); r <- stat_results(); req(!is.null(r$descriptive))
    datatable(r$descriptive,
      options = list(dom = 'tip', pageLength = 20, scrollX = TRUE),
      caption = paste("Descriptive statistics of user-defined groups using per-replicate", r$yvar_used %||% "summaries"))
  })

  emm_plot_obj <- reactive({
    req(stat_results()); r <- stat_results(); req(r$emm_summary)
    emm_df <- as.data.frame(r$emm_summary); cols <- cond_colors()
    # emmeans names its CI columns lower.CL/upper.CL or asymp.LCL/asymp.UCL
    # depending on the fitted model, so resolve them instead of assuming.
    lo <- intersect(c("lower.CL", "asymp.LCL", "lower.HPD"), names(emm_df))[1]
    hi <- intersect(c("upper.CL", "asymp.UCL", "upper.HPD"), names(emm_df))[1]
    validate(need(!is.na(lo) && !is.na(hi) && "emmean" %in% names(emm_df),
                  "Estimated marginal means are not available for this model."))
    p <- ggplot(emm_df, aes(x = Time_fac, y = emmean, color = Condition, group = Condition)) +
      geom_line(position = position_dodge(0.3), linewidth = 0.8) +
      geom_errorbar(aes(ymin = .data[[lo]], ymax = .data[[hi]]),
                    position = position_dodge(0.3), width = 0.2, linewidth = 0.6) +
      geom_point(position = position_dodge(0.3), size = 3) +
      scale_color_manual(values = cols) +
      labs(x = "Time", y = "Estimated Marginal Mean", title = "EMMs with 95% CI",
           subtitle = "From LMM; accounts for replicate random effects") +
      theme_publication(12)
    wrap_legend_labels(p)
  })
  output$emm_plot <- renderPlot({ emm_plot_obj() })
  output$dl_emm_plot <- downloadHandler(
    filename = function() paste0("emm_plot.", fmt_ext(input$emm_format %||% "pdf")),
    content = function(file) save_figure(file, emm_plot_obj(), 7, 4.5, 300, input$emm_format %||% "pdf"))

  # === DOSE-RESPONSE TREND TEST ===
  # Fits lm(replicate summary ~ numeric predictor) within each stratum.
  # Far more powerful than mass pairwise testing for gradient designs,
  # and still anti-pseudoreplication (one value per replicate per group).
  output$trend_predictor_ui <- renderUI({
    req(raw_data())
    df <- raw_data()
    num_cols <- names(df)[vapply(df, is.numeric, TRUE)]
    num_cols <- setdiff(num_cols, c(measure_cols(), "Replicate"))
    selectInput("trend_predictor", "Numeric predictor (e.g., dose)",
      choices = num_cols, selected = grep("dose|conc|nM|uM", num_cols, ignore.case = TRUE, value = TRUE)[1] %||% num_cols[1])
  })
  output$trend_strata_ui <- renderUI({
    req(raw_data())
    df <- raw_data()
    cand <- detect_metadata_columns(df, measure_cols())
    cand <- setdiff(cand, c(input$trend_predictor, "Replicate"))
    selectInput("trend_strata", "Stratify by (optional)",
      choices = cand, selected = intersect(c("Time_min"), cand), multiple = TRUE)
  })

  trend_results <- eventReactive(input$run_trend, {
    req(transformed_data(), input$measure, input$trend_predictor)
    df <- transformed_data()
    m <- input$stat_measure %||% input$measure
    if (!m %in% names(df)) m <- input$measure
    pred <- input$trend_predictor
    strata_cols <- input$trend_strata %||% character(0)
    strata_cols <- strata_cols[strata_cols %in% names(df)]
    validate(need(pred %in% names(df), "Predictor column not found."))
    yvar <- input$stat_summary %||% "median_val"
    # Replicate-level summaries within predictor x strata cells
    grp <- unique(c(pred, strata_cols, "Replicate"))
    rs <- df %>%
      group_by(across(all_of(grp))) %>%
      summarise(median_val = median(.data[[m]], na.rm = TRUE),
                mean_val = mean(.data[[m]], na.rm = TRUE),
                n_cells = sum(!is.na(.data[[m]])), .groups = "drop")
    rs$Stat_Value <- rs[[yvar]]
    rs$Stratum <- if (length(strata_cols) == 0) "All data" else
      do.call(paste, c(lapply(strata_cols, function(cl) paste0(cl, "=", rs[[cl]])), sep = " | "))
    # Order strata by their real values, so numeric levels read 0, 60, 120
    # rather than sorting as text ("0", "120", "60").
    if (length(strata_cols) > 0) {
      key <- unique(rs[, c("Stratum", strata_cols), drop = FALSE])
      ord <- do.call(order, lapply(strata_cols, function(cl) {
        v <- key[[cl]]; nv <- suppressWarnings(as.numeric(as.character(v)))
        if (all(is.na(nv) == is.na(v))) nv else as.character(v)
      }))
      rs$Stratum <- factor(rs$Stratum, levels = key$Stratum[ord])
    }
    out <- rs %>% group_by(Stratum) %>% group_modify(function(d, key) {
      d <- d[is.finite(d$Stat_Value) & is.finite(d[[pred]]), , drop = FALSE]
      n_lvls <- length(unique(d[[pred]]))
      if (nrow(d) < 3 || n_lvls < 3) {
        return(data.frame(Slope = NA_real_, CI95_low = NA_real_, CI95_high = NA_real_,
                          p_value = NA_real_, R2 = NA_real_, Spearman_rho = NA_real_,
                          Spearman_p = NA_real_, N_points = nrow(d), N_levels = n_lvls,
                          Note = "Need >= 3 predictor levels and >= 3 points"))
      }
      fit <- lm(d$Stat_Value ~ d[[pred]])
      sm <- summary(fit)
      ci <- tryCatch(confint(fit)[2, ], error = function(e) c(NA_real_, NA_real_))
      sp <- tryCatch(suppressWarnings(cor.test(d[[pred]], d$Stat_Value, method = "spearman")),
                     error = function(e) NULL)
      data.frame(Slope = coef(fit)[2], CI95_low = ci[1], CI95_high = ci[2],
                 p_value = sm$coefficients[2, 4], R2 = sm$r.squared,
                 Spearman_rho = if (!is.null(sp)) unname(sp$estimate) else NA_real_,
                 Spearman_p = if (!is.null(sp)) sp$p.value else NA_real_,
                 N_points = nrow(d), N_levels = n_lvls, Note = "")
    }) %>% ungroup()
    out$p_adj <- p.adjust(out$p_value, method = "holm")
    out$Significance <- make_significance_label(out$p_adj)
    list(table = out, rs = rs, pred = pred, measure = m, yvar = yvar)
  })

  trend_plot_obj <- reactive({
    req(trend_results())
    tr <- trend_results()
    unit_lbl <- if (tr$yvar == "mean_val") "replicate mean" else "replicate median"
    p <- ggplot(tr$rs, aes(x = .data[[tr$pred]], y = Stat_Value)) +
      geom_point(size = 2.5, alpha = 0.8, color = "#3498DB") +
      geom_smooth(method = "lm", se = TRUE, linewidth = 1, color = "black", alpha = 0.15) +
      # Stratum labels are verbose ("col=value | col=value"); wrap them so the
      # facet strips stay inside their panels instead of being clipped.
      facet_wrap(~ Stratum, scales = "free_y",
                 labeller = labeller(Stratum = function(x) wrap_axis_labels(x, 26))) +
      labs(x = tr$pred, y = paste0(unit_lbl, " of ", tr$measure),
           title = paste("Dose-Response Trend --", tr$measure),
           subtitle = "Each point = one biological replicate; line = linear fit with 95% CI",
           caption = "Trend test on replicate summaries; Holm-corrected across strata.") +
      get_theme(input, "dist")
    p
  })
  output$trend_plot <- renderPlot({ trend_plot_obj() })
  output$trend_table <- renderDT({
    req(trend_results())
    tt <- trend_results()$table
    tt[] <- lapply(tt, function(x) if (is.numeric(x)) round(x, 4) else x)
    datatable(tt, options = list(pageLength = 15, scrollX = TRUE, dom = 'tip'),
      caption = "Linear trend of replicate summaries vs predictor, per stratum (Holm-adjusted)")
  })
  output$dl_trend_plot <- downloadHandler(
    filename = function() paste0("trend_", input$measure, ".pdf"),
    content = function(file) save_figure(file, trend_plot_obj(), 9, 6, 300, "pdf"))
  output$dl_trend_csv <- downloadHandler(
    filename = function() paste0("trend_", input$measure, ".csv"),
    content = function(file) write.csv(trend_results()$table, file, row.names = FALSE))

  output$dl_stats <- downloadHandler(
    filename = function() paste0("stats_", input$stat_measure %||% input$measure, "_", Sys.Date(), ".txt"),
    content  = function(file) {
      r <- stat_results(); sink(file)
      cat("CellMorphR ", APP_VERSION, " -- Flexible Statistical Results\n", sep = "")
      cat("Date:", as.character(Sys.time()), "\n\n")
      print_stat_report(r, include_tables = TRUE)
      sink()
    })

  output$dl_stats_csv <- downloadHandler(
    filename = function() paste0("pairwise_stats_", input$stat_measure %||% input$measure, "_", Sys.Date(), ".csv"),
    content  = function(file) {
      r <- stat_results()
      out <- if (!is.null(r$pairwise)) r$pairwise else data.frame()
      write.csv(out, file, row.names = FALSE)
    })

  # === EFFECT SIZES TAB ===
  es_data <- reactive({
    req(rep_summ(), input$es_summary, input$ref_condition)
    rs <- rep_summ(); yvar <- input$es_summary; ref <- input$ref_condition
    conds <- setdiff(levels(rs$Condition), ref)
    if (length(conds) == 0) conds <- setdiff(unique(as.character(rs$Condition)), ref)
    times <- sort(unique(rs$Time_min))
    trt_rows <- lapply(times, function(tm) {
      cv <- rs[[yvar]][rs$Condition == ref & rs$Time_min == tm]
      lapply(conds, function(cd) {
        tv <- rs[[yvar]][rs$Condition == cd & rs$Time_min == tm]
        cm <- mean(cv); tm2 <- mean(tv)
        # Hedges' g with a noncentral-t interval, and NA (never an error) when a
        # condition is absent or has a single replicate at this timepoint.
        g <- hedges_g(cv, tv)
        data.frame(Time_min = tm, Condition = cd,
          n_reps_ref = length(cv[is.finite(cv)]), n_reps_cond = length(tv[is.finite(tv)]),
          ctrl_mean = cm, trt_mean = tm2,
          abs_diff = tm2 - cm,
          pct_change = if (isTRUE(abs(cm) > 1e-10)) (tm2 - cm) / cm * 100 else NA_real_,
          hedges_g = g$g, hedges_g_low = g$ci_low, hedges_g_high = g$ci_high,
          note = g$note %||% NA_character_, is_control = FALSE,
          stringsAsFactors = FALSE)
      }) %>% bind_rows()
    }) %>% bind_rows()
    # Add control baseline rows (effect size = 0)
    ctrl_rows <- data.frame(Time_min = times, Condition = ref,
      n_reps_ref = NA_integer_, n_reps_cond = NA_integer_, ctrl_mean = NA, trt_mean = NA,
      abs_diff = 0, pct_change = 0, hedges_g = 0, hedges_g_low = NA_real_,
      hedges_g_high = NA_real_, note = NA_character_, is_control = TRUE,
      stringsAsFactors = FALSE)
    rbind(ctrl_rows, trt_rows)
  })

  es_plot_obj <- reactive({
    req(es_data(), input$es_type)
    esd <- es_data(); cols <- cond_colors(); yvar <- input$es_type; mlbl <- measure_label()
    ylab <- switch(yvar, abs_diff = paste("Delta", mlbl), pct_change = paste("% Change in", mlbl),
                   hedges_g = "Hedges' g (95% CI)")
    # Filter control if not requested
    esd_plot <- if (isTRUE(input$es_show_control)) esd else esd[!esd$is_control, , drop = FALSE]
    validate(need(nrow(esd_plot) > 0 && any(!esd$is_control),
      "Effect sizes need at least one condition besides the reference condition."))
    cols_es <- if (isTRUE(input$es_show_control)) cols else {
      ref <- input$ref_condition %||% names(cols)[1]
      cols[!names(cols) %in% ref]
    }
    # Keep the legend in the app's condition order (i.e. dose order) rather than
    # letting a character column sort alphabetically (100nM before 12.5nM).
    esd_plot$Condition <- factor(as.character(esd_plot$Condition),
                                 levels = intersect(names(cols), unique(as.character(esd_plot$Condition))))
    style <- input$es_plot_style %||% "linepoint"
    p <- ggplot(esd_plot, aes(x = Time_min, y = .data[[yvar]], color = Condition)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")
    if (style == "bubble") {
      esd_plot$abs_mag <- abs(esd_plot[[yvar]])
      p <- ggplot(esd_plot, aes(x = Time_min, y = .data[[yvar]], color = Condition, size = abs_mag)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
        geom_point(alpha = 0.75) + scale_size_area(max_size = 18, name = "|Effect|")
    } else if (style == "points") {
      p <- p + geom_point(size = 4)
    } else {
      p <- p + geom_line(linewidth = 1.2) + geom_point(size = 4)
    }
    p <- p + scale_color_manual(values = cols_es) + scale_x_continuous(breaks = sort(unique(esd$Time_min))) +
      labs(x = "Time (min)", y = ylab, title = paste("Effect on", mlbl, "Over Time"),
           subtitle = if (isTRUE(input$es_show_control)) "Control shown as baseline (0)" else "Treatment vs control") +
      get_theme(input, "es")
    if (yvar == "hedges_g") {
      # The guides are Cohen's conventional benchmarks, not thresholds for this
      # experiment; they are only honest next to an interval, so the interval is
      # drawn with them (audit CM-3).
      ci <- esd_plot[!esd_plot$is_control & is.finite(esd_plot$hedges_g_low), , drop = FALSE]
      if (nrow(ci) > 0) p <- p + geom_errorbar(
        data = ci, aes(x = Time_min, ymin = hedges_g_low, ymax = hedges_g_high, color = Condition),
        width = 2, linewidth = 0.6, alpha = 0.8, inherit.aes = FALSE)
      p <- p + geom_hline(yintercept = c(-0.8, 0.8), linetype = "dotted", color = "grey70") +
        labs(caption = paste0("Hedges' g (small-sample corrected) with 95% noncentral-t intervals, ",
                              "computed on biological replicate summaries. Dotted lines are Cohen's ",
                              "conventional |0.8| benchmark, not a significance threshold."))
    }
    dropped <- esd[!esd$is_control & is.na(esd$hedges_g) & !is.na(esd$note), , drop = FALSE]
    if (yvar == "hedges_g" && nrow(dropped) > 0) {
      p <- p + labs(subtitle = paste0(
        (if (isTRUE(input$es_show_control)) "Control shown as baseline (0)" else "Treatment vs control"),
        ". Not estimable at ", nrow(dropped), " condition-timepoint combination",
        if (nrow(dropped) == 1) "" else "s", " (fewer than 2 replicates); see the Table tab."))
    }
    apply_custom_labs(p, input, "es")
  })

  output$es_plot <- renderPlot({ es_plot_obj() })
  output$es_table <- renderDT({ req(es_data()); esd <- es_data(); esd[] <- lapply(esd, function(x) if (is.numeric(x)) round(x, 3) else x); datatable(esd, options = list(dom = 't')) })
  output$dl_es <- downloadHandler(
    filename = function() paste0("effect_sizes_", input$measure, ".", fmt_ext(input$es_format)),
    content  = function(file) save_figure(file, es_plot_obj(), input$es_width, input$es_height, input$es_dpi, input$es_format))

  # === MORPHOSPACE (PCA) TAB ===
  output$pca_exclude_features_ui <- renderUI({
    req(measure_cols()); mc <- measure_cols()
    selectInput("pca_exclude_features", "Exclude features from PCA",
      choices = mc, selected = NULL, multiple = TRUE,
      selectize = TRUE)
  })

  pca_features <- reactive({
    mc <- measure_cols(); req(mc)
    excl <- input$pca_exclude_features
    if (!is.null(excl) && length(excl) > 0) mc <- setdiff(mc, excl)
    req(length(mc) >= 2)
    mc
  })

  pca_result <- reactive({ req(transformed_data(), pca_features()); run_pca(transformed_data(), pca_features(), scale = isTRUE(input$pca_scale)) })

  pca_scatter_obj <- reactive({
    req(pca_result()); res <- pca_result(); req(res$ok)
    scores <- as.data.frame(res$pca$x[, 1:min(ncol(res$pca$x), 2)])
    scores$Condition <- res$df$Condition; scores$Time_min <- factor(res$df$Time_min)
    scores$cond_time <- interaction(res$df$Condition, res$df$Time_min, sep = " | ")
    var_pct <- round(res$var_explained[2, 1:min(ncol(res$pca$x), 2)] * 100, 1)
    color_var <- input$pca_color_by; cols <- cond_colors()
    p <- ggplot(scores, aes(x = PC1, y = PC2, color = .data[[color_var]])) +
      geom_point(alpha = input$pca_point_alpha, size = input$pca_point_size %||% 1) +
      labs(x = paste0("PC1 (", var_pct[1], "%)"), y = paste0("PC2 (", var_pct[2], "%)"),
           title = "Single-Cell Morphospace", subtitle = "Each point = one cell") +
      get_theme(input, "pca")
    # Apply appropriate color scale based on color_by selection
    if (color_var == "Condition") {
      p <- p + scale_color_manual(values = cols)
    } else if (color_var == "Time_min") {
      p <- p + scale_color_viridis_d(option = "plasma", name = "Time (min)")
    } else if (color_var == "cond_time") {
      p <- p + scale_color_viridis_d(option = "turbo", name = "Condition | Time")
    }
    # Centroids use same color variable to avoid scale conflicts
    if (isTRUE(input$pca_show_centroids)) {
      cen <- scores %>% group_by(Condition, Time_min, cond_time) %>%
        summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")
      # Black outline ring for visibility
      p <- p + geom_point(data = cen, aes(x = PC1, y = PC2), size = 6, shape = 5,
                          color = "black", stroke = 1.2, inherit.aes = FALSE)
      # Colored centroid point — use same color variable as main scatter
      p <- p + geom_point(data = cen, aes(x = PC1, y = PC2, color = .data[[color_var]]),
                          size = 4.5, shape = 18, inherit.aes = FALSE)
    }
    # Density contours
    if (isTRUE(input$pca_show_density)) {
      if (color_var == "Condition") {
        p <- p + geom_density_2d(aes(color = Condition), linewidth = 0.3, alpha = 0.5)
      } else {
        p <- p + geom_density_2d(color = "grey50", linewidth = 0.3, alpha = 0.5)
      }
    }
    if (isTRUE(input$pca_show_loadings)) {
      ld <- as.data.frame(res$pca$rotation[, 1:2]); ld$Feature <- rownames(ld)
      # Scale the loading arrows to the spread of the scores. A fixed multiplier
      # leaves every arrow (and its label) bunched at the origin whenever the
      # scores span a wide range, which made the feature names unreadable.
      score_rng <- max(abs(c(scores$PC1, scores$PC2)), na.rm = TRUE)
      ld_max <- max(sqrt(ld$PC1^2 + ld$PC2^2), na.rm = TRUE)
      auto_sc <- if (is.finite(ld_max) && ld_max > 0 && is.finite(score_rng)) 0.6 * score_rng / ld_max else 1
      sc <- auto_sc * ((input$pca_loading_scale %||% 3) / 3)
      p <- p + geom_segment(data = ld, aes(x = 0, y = 0, xend = PC1 * sc, yend = PC2 * sc),
                            arrow = arrow(length = unit(0.2, "cm")), color = "grey30", inherit.aes = FALSE) +
        geom_text(data = ld, aes(x = PC1 * sc * 1.08, y = PC2 * sc * 1.08, label = Feature),
                  size = 3, color = "grey20", inherit.aes = FALSE,
                  check_overlap = TRUE, fontface = "bold")
    }
    apply_custom_labs(p, input, "pca")
  })

  output$pca_scatter <- renderPlot({
    res <- pca_result()
    if (!is.null(res) && res$ok) pca_scatter_obj()
    else { plot.new(); text(0.5, 0.5, if (!is.null(res)) res$msg else "Need >= 2 features.", cex = 1.2) }
  })
  pca_scree_obj <- reactive({
    req(pca_result()); res <- pca_result(); req(res$ok)
    n <- min(ncol(res$var_explained), 10)
    vdf <- data.frame(PC = factor(paste0("PC", 1:n), levels = paste0("PC", 1:n)), Var = res$var_explained[2, 1:n] * 100)
    ggplot(vdf, aes(x = PC, y = Var)) + geom_col(fill = "#3498DB", alpha = 0.8) +
      geom_text(aes(label = paste0(round(Var, 1), "%")), vjust = -0.5, size = 3) +
      labs(x = NULL, y = "% Variance", title = "Scree Plot") + theme_publication(10)
  })
  output$pca_scree <- renderPlot({ pca_scree_obj() })
  output$dl_pca_scree <- downloadHandler(
    filename = function() paste0("pca_scree.", fmt_ext(input$pca_format)),
    content = function(file) save_figure(file, pca_scree_obj(), input$pca_width, input$pca_height, input$pca_dpi, input$pca_format))

  output$pca_loadings_table <- renderDT({
    req(pca_result()); res <- pca_result(); req(res$ok)
    n <- min(3, ncol(res$pca$rotation)); ld <- as.data.frame(res$pca$rotation[, 1:n])
    ld$Feature <- rownames(ld); ld <- ld[, c("Feature", setdiff(names(ld), "Feature"))]
    ld[, -1] <- round(ld[, -1], 3); datatable(ld, options = list(dom = 't', pageLength = 20))
  })

  pca_traj_obj <- reactive({
    req(pca_result(), input$pca_show_trajectories); res <- pca_result(); req(res$ok)
    scores <- as.data.frame(res$pca$x[, 1:2]); scores$Condition <- res$df$Condition; scores$Time_min <- res$df$Time_min
    cen <- scores %>% group_by(Condition, Time_min) %>% summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop") %>% arrange(Condition, Time_min)
    vpct <- round(res$var_explained[2, 1:2] * 100, 1); cols <- cond_colors()
    ggplot(cen, aes(x = PC1, y = PC2, color = Condition)) +
      geom_path(aes(group = Condition), linewidth = 1.2, arrow = arrow(length = unit(0.25, "cm"), type = "closed")) +
      geom_point(size = 4) + geom_text(aes(label = fmt_time(Time_min)), vjust = -1, size = 3.5) +
      scale_color_manual(values = cols) +
      labs(x = paste0("PC1 (", vpct[1], "%)"), y = paste0("PC2 (", vpct[2], "%)"),
           title = "Morphological Trajectory", subtitle = "Centroids connected over time") +
      get_theme(input, "pca")
  })

  output$pca_trajectory <- renderPlot({
    if (isTRUE(input$pca_show_trajectories)) pca_traj_obj()
    else { plot.new(); text(0.5, 0.5, "Enable 'Show centroid trajectories' in sidebar.", cex = 1.2) }
  })

  # --- PCA Loadings Bar Plot ---
  pca_loadings_bar_obj <- reactive({
    req(pca_result()); res <- pca_result(); req(res$ok)
    rot <- as.data.frame(res$pca$rotation[, 1:min(2, ncol(res$pca$rotation))])
    rot$Feature <- rownames(rot)
    vpct <- round(res$var_explained[2, 1:min(2, ncol(res$var_explained))] * 100, 1)
    # PC1 loadings
    rot1 <- rot[order(abs(rot$PC1), decreasing = TRUE), ]
    rot1$Feature <- factor(rot1$Feature, levels = rev(rot1$Feature))
    rot1$Direction <- ifelse(rot1$PC1 >= 0, "Positive", "Negative")
    p1 <- ggplot(rot1, aes(x = Feature, y = PC1, fill = Direction)) +
      geom_col(alpha = 0.85, width = 0.7) +
      coord_flip() + geom_hline(yintercept = 0, linewidth = 0.3) +
      scale_fill_manual(values = c("Positive" = "#3498DB", "Negative" = "#E74C3C")) +
      labs(x = NULL, y = paste0("Loading on PC1 (", vpct[1], "%)"),
           title = "Feature Contributions to PC1",
           subtitle = "Which morphological features drive the main axis of variation") +
      get_theme(input, "pca") + theme(legend.position = "none")
    if (ncol(rot) >= 3) {
      rot2 <- rot[order(abs(rot$PC2), decreasing = TRUE), ]
      rot2$Feature <- factor(rot2$Feature, levels = rev(rot2$Feature))
      rot2$Direction <- ifelse(rot2$PC2 >= 0, "Positive", "Negative")
      p2 <- ggplot(rot2, aes(x = Feature, y = PC2, fill = Direction)) +
        geom_col(alpha = 0.85, width = 0.7) +
        coord_flip() + geom_hline(yintercept = 0, linewidth = 0.3) +
        scale_fill_manual(values = c("Positive" = "#3498DB", "Negative" = "#E74C3C")) +
        labs(x = NULL, y = paste0("Loading on PC2 (", vpct[2], "%)"),
             title = "Feature Contributions to PC2") +
        get_theme(input, "pca") + theme(legend.position = "none")
      if (has_patchwork) p1 / p2 else p1
    } else p1
  })
  output$pca_loadings_bar <- renderPlot({ pca_loadings_bar_obj() })
  output$dl_pca_loadings <- downloadHandler(
    filename = function() paste0("pca_loadings.", fmt_ext(input$pca_format)),
    content = function(file) save_figure(file, pca_loadings_bar_obj(), input$pca_width, input$pca_height, input$pca_dpi, input$pca_format))

  # --- Morphological Diversity Metrics ---
  pca_diversity_data <- reactive({
    req(pca_result()); res <- pca_result(); req(res$ok)
    scores <- as.data.frame(res$pca$x[, 1:min(2, ncol(res$pca$x))])
    scores$Condition <- res$df$Condition
    scores$Time_min <- res$df$Time_min

    div_list <- list()
    groups <- unique(interaction(scores$Condition, scores$Time_min, sep = " | "))
    for (g in groups) {
      parts <- strsplit(as.character(g), " \\| ")[[1]]
      cond_val <- parts[1]; time_val <- as.numeric(parts[2])
      idx <- scores$Condition == cond_val & scores$Time_min == time_val
      sub <- scores[idx, , drop = FALSE]
      if (nrow(sub) < 3) next
      pc1_vals <- sub$PC1; pc2_vals <- sub$PC2
      cen_pc1 <- mean(pc1_vals); cen_pc2 <- mean(pc2_vals)
      # Convex hull area using shoelace formula
      hull_area <- tryCatch({
        hi <- chull(pc1_vals, pc2_vals)
        hx <- pc1_vals[hi]; hy <- pc2_vals[hi]
        n <- length(hx)
        0.5 * abs(sum(hx[1:(n-1)] * hy[2:n]) + hx[n] * hy[1] -
                  sum(hx[2:n] * hy[1:(n-1)]) - hx[1] * hy[n])
      }, error = function(e) NA_real_)
      # Mean distance to centroid
      dists <- sqrt((pc1_vals - cen_pc1)^2 + (pc2_vals - cen_pc2)^2)
      div_list[[length(div_list) + 1]] <- data.frame(
        Condition = cond_val, Time_min = time_val, N_cells = nrow(sub),
        Hull_area = round(hull_area, 3),
        PC1_var = round(var(pc1_vals), 3), PC2_var = round(var(pc2_vals), 3),
        Mean_dist_centroid = round(mean(dists), 3),
        Centroid_PC1 = round(cen_pc1, 3), Centroid_PC2 = round(cen_pc2, 3),
        stringsAsFactors = FALSE)
    }
    do.call(rbind, div_list)
  })

  output$pca_diversity_table <- renderDT({
    req(pca_diversity_data()); dd <- pca_diversity_data()
    dd$Time_label <- fmt_time(dd$Time_min)
    show_cols <- c("Condition", "Time_label", "N_cells", "Hull_area", "PC1_var",
                   "PC2_var", "Mean_dist_centroid", "Centroid_PC1", "Centroid_PC2")
    datatable(dd[, show_cols], options = list(pageLength = 20, scrollX = TRUE, dom = 'tip'),
      caption = "Morphological diversity: convex hull area, PC variance, mean distance to centroid",
      colnames = c("Condition", "Time", "N Cells", "Hull Area", "Var(PC1)",
                   "Var(PC2)", "Mean Dist.", "Centroid PC1", "Centroid PC2"))
  })

  pca_diversity_plot_obj <- reactive({
    req(pca_diversity_data()); dd <- pca_diversity_data(); cols <- cond_colors()
    dd$Time_fac <- factor(fmt_time(dd$Time_min),
                          levels = fmt_time(sort(unique(dd$Time_min))))
    p1 <- ggplot(dd, aes(x = Time_fac, y = Hull_area, fill = Condition)) +
      geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.85) +
      scale_fill_manual(values = cols) +
      labs(x = "Time", y = "Convex Hull Area (PC1-PC2)",
           title = "Morphospace Expansion",
           subtitle = "Larger hull = more morphological diversity") +
      get_theme(input, "pca")
    p2 <- ggplot(dd, aes(x = Time_fac, y = Mean_dist_centroid, fill = Condition)) +
      geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.85) +
      scale_fill_manual(values = cols) +
      labs(x = "Time", y = "Mean Distance to Centroid",
           title = "Morphological Heterogeneity",
           subtitle = "Higher = more phenotypic spread from group center") +
      get_theme(input, "pca")
    if (has_patchwork) p1 + p2 else p1
  })
  output$pca_diversity_plot <- renderPlot({ pca_diversity_plot_obj() })
  output$dl_pca_diversity <- downloadHandler(
    filename = function() paste0("pca_diversity.", fmt_ext(input$pca_format)),
    content = function(file) save_figure(file, pca_diversity_plot_obj(), input$pca_width, input$pca_height, input$pca_dpi, input$pca_format))

  output$dl_pca <- downloadHandler(filename = function() paste0("morphospace.", fmt_ext(input$pca_format)),
    content = function(file) save_figure(file, pca_scatter_obj(), input$pca_width, input$pca_height, input$pca_dpi, input$pca_format))
  output$dl_pca_traj <- downloadHandler(filename = function() paste0("trajectory.", fmt_ext(input$pca_format)),
    content = function(file) { req(input$pca_show_trajectories); save_figure(file, pca_traj_obj(), input$pca_width, input$pca_height, input$pca_dpi, input$pca_format) })

  # === PHENOTYPES (CLUSTERING) TAB ===
  cluster_result <- reactiveVal(NULL)
  observeEvent(input$run_clustering, {
    req(pca_result()); res <- pca_result(); req(res$ok)
    if (isTRUE(input$clust_on_pca)) { n_pcs <- min(input$clust_n_pcs, ncol(res$pca$x)); cdata <- res$pca$x[, 1:n_pcs, drop = FALSE]
    } else cdata <- scale(res$df[res$features_used])
    if (input$clust_k_method == "manual") { k <- input$clust_k; ki <- NULL
    } else { withProgress(message = "Finding optimal k...", value = 0.3, { ki <- find_optimal_k(cdata, max_k = 10) }); k <- ki$optimal_k }
    withProgress(message = paste("Clustering k =", k, "..."), value = 0.7, {
      if (input$clust_method == "kmeans") { km <- tryCatch(kmeans(cdata, centers = k, nstart = 25, iter.max = 100), error = function(e) NULL)
        if (!is.null(km)) clusters <- km$cluster else { showNotification("K-means failed.", type = "error"); return() }
      } else { hc <- tryCatch({ d <- dist(cdata); hclust(d, method = "ward.D2") }, error = function(e) NULL)
        if (!is.null(hc)) clusters <- cutree(hc, k = k) else { showNotification("hclust failed.", type = "error"); return() } }
    })
    df_out <- res$df; df_out$Phenotype <- factor(paste0("P", clusters), levels = paste0("P", 1:k))
    ps <- as.data.frame(res$pca$x[, 1:min(2, ncol(res$pca$x))]); df_out$PC1 <- ps$PC1; df_out$PC2 <- ps$PC2
    cluster_result(list(df = df_out, k = k, k_info = ki, features_used = res$features_used, var_explained = res$var_explained))
    showNotification(paste(k, "phenotypes identified."), type = "message")
  })

  clust_k_plot_obj <- reactive({
    req(cluster_result()); cr <- cluster_result(); req(cr$k_info)
    ki <- cr$k_info; pdf <- data.frame(k = ki$k_range, Sil = ki$sil)
    pdf <- pdf[!is.na(pdf$Sil), ]
    ggplot(pdf, aes(x = k, y = Sil)) + geom_line(linewidth = 1, color = "#3498DB") + geom_point(size = 3, color = "#3498DB") +
      geom_vline(xintercept = cr$k, linetype = "dashed", color = "#D94A4A") +
      annotate("text", x = cr$k + 0.3, y = max(pdf$Sil, na.rm = TRUE), label = paste("k =", cr$k), color = "#D94A4A") +
      scale_x_continuous(breaks = ki$k_range) + labs(x = "k", y = "Silhouette", title = "Optimal k") + theme_publication(10)
  })
  output$clust_k_plot <- renderPlot({ clust_k_plot_obj() })
  output$dl_clust_k <- downloadHandler(
    filename = function() paste0("cluster_optimal_k.", fmt_ext(input$clust_format)),
    content = function(file) save_figure(file, clust_k_plot_obj(), input$clust_width, input$clust_height, input$clust_dpi, input$clust_format))

  clust_pca_plot_obj <- reactive({
    req(cluster_result()); cr <- cluster_result()
    vpct <- round(cr$var_explained[2, 1:2] * 100, 1)
    ggplot(cr$df, aes(x = PC1, y = PC2, color = Phenotype)) + geom_point(alpha = 0.3, size = 1) +
      labs(x = paste0("PC1 (", vpct[1], "%)"), y = paste0("PC2 (", vpct[2], "%)"), title = "By Phenotype") + theme_publication(10)
  })
  output$clust_pca_plot <- renderPlot({ clust_pca_plot_obj() })
  output$dl_clust_pca <- downloadHandler(
    filename = function() paste0("cluster_pca.", fmt_ext(input$clust_format)),
    content = function(file) save_figure(file, clust_pca_plot_obj(), input$clust_width, input$clust_height, input$clust_dpi, input$clust_format))

  clust_comp_obj <- reactive({
    req(cluster_result()); cr <- cluster_result()
    comp <- cr$df %>% group_by(Condition, Time_min, Phenotype) %>% summarise(n = n(), .groups = "drop") %>%
      group_by(Condition, Time_min) %>% mutate(pct = n / sum(n) * 100) %>% ungroup()
    comp$Time_fac <- factor(fmt_time(comp$Time_min), levels = fmt_time(sort(unique(comp$Time_min))))
    p <- ggplot(comp, aes(x = Time_fac, y = pct, fill = Phenotype)) + geom_col(position = "stack", width = 0.7) +
      facet_wrap(~ Condition) + labs(x = "Time", y = "% Cells", title = "Phenotype Composition Over Time") +
      get_theme(input, "clust") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    apply_custom_labs(p, input, "clust")
  })

  output$clust_composition <- renderPlot({ clust_comp_obj() })

  output$clust_comp_table <- renderDT({
    req(cluster_result()); cr <- cluster_result()
    comp <- cr$df %>% group_by(Condition, Time_min, Phenotype) %>%
      summarise(Count = n(), .groups = "drop") %>%
      group_by(Condition, Time_min) %>%
      mutate(Total = sum(Count), Percent = round(Count / Total * 100, 1)) %>% ungroup()
    comp$Time_label <- fmt_time(comp$Time_min)
    datatable(comp[, c("Condition", "Time_label", "Phenotype", "Count", "Total", "Percent")],
      options = list(pageLength = 30, scrollX = TRUE, dom = 'tip'),
      caption = "Phenotype composition: exact counts and percentages per Condition × Time",
      colnames = c("Condition", "Time", "Phenotype", "Count", "Total Cells", "% of Group"))
  })

  output$clust_char_table <- renderDT({
    req(cluster_result()); cr <- cluster_result()
    feats <- cr$features_used
    # Compute mean, median, SD, IQR, CV for each phenotype × feature
    char_list <- lapply(feats, function(f) {
      cr$df %>% group_by(Phenotype) %>%
        summarise(
          Feature = f,
          N      = n(),
          Mean   = round(mean(.data[[f]], na.rm = TRUE), 3),
          Median = round(median(.data[[f]], na.rm = TRUE), 3),
          SD     = round(sd(.data[[f]], na.rm = TRUE), 3),
          IQR    = round(IQR(.data[[f]], na.rm = TRUE), 3),
          CV_pct = round(if (abs(mean(.data[[f]], na.rm = TRUE)) > 1e-10)
            sd(.data[[f]], na.rm = TRUE) / abs(mean(.data[[f]], na.rm = TRUE)) * 100 else NA_real_, 1),
          .groups = "drop"
        )
    })
    ch <- do.call(rbind, char_list)
    ch <- ch[, c("Phenotype", "Feature", "N", "Mean", "Median", "SD", "IQR", "CV_pct")]
    datatable(ch, options = list(pageLength = 30, scrollX = TRUE, dom = 'tip'),
      caption = "Detailed characterization of each phenotype by feature",
      colnames = c("Phenotype", "Feature", "N", "Mean", "Median", "SD", "IQR", "CV (%)"))
  })
  output$dl_clust <- downloadHandler(filename = function() paste0("phenotypes.", fmt_ext(input$clust_format)),
    content = function(file) save_figure(file, clust_comp_obj(), input$clust_width, input$clust_height, input$clust_dpi, input$clust_format))

  # === MULTI-FEATURE SUMMARY TAB ===
  mf_results <- reactiveVal(NULL)

  observeEvent(input$run_multifeature, {
    req(transformed_data(), measure_cols(), input$ref_condition)
    df <- transformed_data()
    mcols <- if (!is.null(input$mf_features)) input$mf_features else measure_cols()
    n_features <- length(mcols)
    alpha <- input$mf_alpha %||% 0.05
    results_list <- list(); div_list <- list()
    withProgress(message = "Analyzing...", value = 0, {
      for (i in seq_along(mcols)) {
        incProgress(1 / n_features, detail = mcols[i])
        r <- tryCatch(run_anova_for_feature(df, mcols[i], input$mf_summary), error = function(e) NULL)
        if (!is.null(r)) results_list[[i]] <- r
        d <- tryCatch(detect_earliest_divergence(df, mcols[i], input$ref_condition, alpha), error = function(e) NULL)
        if (!is.null(d)) div_list[[i]] <- d
      }
    })
    mf_results(list(results = do.call(rbind, results_list), divergence = do.call(rbind, div_list)))
    showNotification("Multi-feature analysis complete.", type = "message")
  })

  mf_heatmap_obj <- reactive({
    req(mf_results()); res <- mf_results()$results; req(nrow(res) > 0)
    heat <- res %>% select(Feature, Condition_p, Time_p, Interaction_p) %>%
      pivot_longer(-Feature, names_to = "Effect", values_to = "p_value") %>%
      mutate(nlp = -log10(pmax(p_value, 1e-20)),
             sig = ifelse(p_value < 0.001, "***", ifelse(p_value < 0.01, "**", ifelse(p_value < 0.05, "*", "ns"))),
             pval_txt = ifelse(p_value < 0.001, "<0.001", formatC(p_value, format = "f", digits = 3)),
             Effect = recode(Effect, Condition_p = "Condition", Time_p = "Time", Interaction_p = "Interaction"))
    heat$Effect <- factor(heat$Effect, levels = c("Condition", "Time", "Interaction"))

    # Choose label
    cell_label <- if (isTRUE(input$mf_show_pvals)) heat$pval_txt
                  else if (isTRUE(input$mf_show_stars)) heat$sig
                  else ""

    # Choose palette
    pal <- input$mf_palette %||% "wyr"
    base_sz <- input$mf_base_size %||% 12

    p <- ggplot(heat, aes(x = Effect, y = Feature, fill = nlp)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = cell_label), size = base_sz / 3, fontface = "bold") +
      labs(x = NULL, y = NULL, fill = "-log10(p)",
           title = "Multi-Feature Statistical Summary",
           subtitle = paste0("Repeated-measures ANOVA (dish as error stratum) on replicate ",
                             if ((input$mf_summary %||% "median_val") == "mean_val") "means" else "medians",
                             ". * p<0.05, ** p<0.01, *** p<0.001"))

    p <- p + switch(pal,
      "wyr" = scale_fill_gradient2(low = "white", mid = "#FFDD57", high = "#D94A4A", midpoint = 1.3),
      "wb"  = scale_fill_gradient(low = "white", high = "#2C3E50"),
      "viridis" = scale_fill_viridis_c(option = "viridis"),
      "wgb" = scale_fill_gradient2(low = "white", mid = "#18BC9C", high = "#1a1a1a", midpoint = 1.3),
      scale_fill_gradient2(low = "white", mid = "#FFDD57", high = "#D94A4A", midpoint = 1.3))

    p <- p + theme_publication(base_sz) + theme(axis.text.x = element_text(angle = 0))
    apply_custom_labs(p, input, "mf")
  })

  output$mf_heatmap <- renderPlot({ mf_heatmap_obj() })

  mf_importance_obj <- reactive({
    req(mf_results()); res <- mf_results()$results; req(nrow(res) > 0)
    imp <- res %>% select(Feature, Condition_partial_eta2, Time_partial_eta2, Interaction_partial_eta2) %>%
      pivot_longer(-Feature, names_to = "Effect", values_to = "eta2") %>%
      mutate(Effect = recode(Effect, Condition_partial_eta2 = "Condition",
                             Time_partial_eta2 = "Time", Interaction_partial_eta2 = "Interaction"))
    imp$Effect <- factor(imp$Effect, levels = c("Condition", "Time", "Interaction"))
    ggplot(imp, aes(x = reorder(Feature, eta2), y = eta2, fill = Effect)) +
      geom_col(position = "dodge", alpha = 0.85) + coord_flip() +
      labs(x = NULL, y = "Partial Eta-Squared", title = "Feature Importance",
           subtitle = "small=0.01, medium=0.06, large=0.14") + theme_publication(10)
  })
  output$mf_importance <- renderPlot({ mf_importance_obj() })
  output$dl_mf_importance <- downloadHandler(
    filename = function() paste0("feature_importance.", fmt_ext(input$mf_format)),
    content = function(file) save_figure(file, mf_importance_obj(), input$mf_width, input$mf_height, input$mf_dpi, input$mf_format))

  output$mf_divergence_table <- renderDT({
    req(mf_results()); div <- mf_results()$divergence
    if (!is.null(div) && nrow(div) > 0) datatable(div, options = list(dom = 't', pageLength = 20),
      caption = "First significant timepoint per feature") })

  output$mf_results_table <- renderDT({
    req(mf_results()); res <- mf_results()$results
    if (!is.null(res) && nrow(res) > 0) { res[] <- lapply(res, function(x) if (is.numeric(x)) round(x, 4) else x)
      datatable(res, options = list(pageLength = 20, scrollX = TRUE, dom = 'tip')) }
  })

  output$dl_mf_heatmap <- downloadHandler(filename = function() paste0("heatmap.", fmt_ext(input$mf_format)),
    content = function(file) save_figure(file, mf_heatmap_obj(), input$mf_width, input$mf_height, input$mf_dpi, input$mf_format))
  output$dl_mf_table <- downloadHandler(filename = function() paste0("multifeature_", Sys.Date(), ".csv"),
    content = function(file) { req(mf_results()); write.csv(mf_results()$results, file, row.names = FALSE) })

  # --- RAW VALUES HEATMAP ---
  mf_raw_heatmap_obj <- reactive({
    req(raw_data()); df <- transformed_data()
    mcols <- if (!is.null(input$mf_features)) input$mf_features else measure_cols()
    req(length(mcols) >= 1)
    stat_fn <- input$mf_raw_stat %||% "mean"
    heat <- compute_raw_heatmap_data(df, mcols, stat_fn, time_fmt = fmt_time)
    req(nrow(heat) > 0)

    # Z-score scaling per feature if requested
    if (isTRUE(input$mf_raw_zscore)) {
      heat <- heat %>% group_by(Feature) %>%
        mutate(value = (value - mean(value, na.rm = TRUE)) /
                         (sd(value, na.rm = TRUE) + 1e-10)) %>% ungroup()
    }

    # Order features by variability (most variable on top)
    feat_order <- heat %>% group_by(Feature) %>%
      summarise(var = var(value, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(var)) %>% pull(Feature)
    heat$Feature <- factor(heat$Feature, levels = rev(feat_order))

    # Ordered group labels
    heat$Group <- factor(heat$Group,
      levels = unique(heat[order(heat$Condition, heat$Time_min), "Group", drop = TRUE]))

    base_sz <- input$mf_base_size %||% 12
    pal_choice <- input$mf_raw_palette %||% "bwr"
    scale_label <- if (isTRUE(input$mf_raw_zscore)) "Z-score" else paste0(stat_fn, " value")

    pal_colors <- switch(pal_choice,
      "bwr" = c("#2166AC", "white", "#B2182B"),
      "gwr" = c("#1B7837", "white", "#E31A1C"),
      "pwo" = c("#7B3294", "white", "#E66101"),
      c("#440154", "#21908C", "#FDE725")  # viridis fallback
    )

    p <- ggplot(heat, aes(x = Group, y = Feature, fill = value)) +
      geom_tile(color = "grey90", linewidth = 0.4)

    if (pal_choice == "viridis") {
      p <- p + scale_fill_viridis_c(name = scale_label, option = "viridis")
    } else {
      p <- p + scale_fill_gradient2(low = pal_colors[1], mid = pal_colors[2],
        high = pal_colors[3], midpoint = if (isTRUE(input$mf_raw_zscore)) 0
          else median(heat$value, na.rm = TRUE), name = scale_label)
    }

    if (isTRUE(input$mf_raw_show_values)) {
      p <- p + geom_text(aes(label = round(value, 2)),
        size = base_sz / 4, color = "black")
    }
    p + labs(title = "Raw Values Heatmap",
             subtitle = paste0("Feature × (Condition × Time) — ", stat_fn,
               if (isTRUE(input$mf_raw_zscore)) " (z-scored)" else ""),
             x = NULL, y = NULL) +
      theme_publication(base_sz, show_grid = FALSE) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = rel(0.8)),
            panel.border = element_blank())
  })

  output$mf_raw_heatmap <- renderPlot({ mf_raw_heatmap_obj() })
  output$dl_mf_raw_heatmap <- downloadHandler(
    filename = function() paste0("raw_heatmap.", fmt_ext(input$mf_format)),
    content = function(file) save_figure(file, mf_raw_heatmap_obj(), input$mf_width, input$mf_height, input$mf_dpi, input$mf_format))

  # --- SPIDER / RADAR PLOT ---
  output$mf_radar_time_ui <- renderUI({
    req(raw_data()); times <- sort(unique(raw_data()$Time_min))
    selectInput("mf_radar_time", "Timepoint",
      choices = c("All (aggregated)" = "all", setNames(times, paste0(times, " min"))),
      selected = "all")
  })

  mf_radar_obj <- reactive({
    req(raw_data()); df <- transformed_data()
    mcols <- if (!is.null(input$mf_features)) input$mf_features else measure_cols()
    req(length(mcols) >= 3)
    stat_fn <- input$mf_radar_stat %||% "mean"
    sel_time <- input$mf_radar_time %||% "all"
    if (sel_time != "all") {
      df <- df[df$Time_min == as.numeric(sel_time), , drop = FALSE]
      time_lbl <- paste0(sel_time, " min")
    } else {
      time_lbl <- "All timepoints (aggregated)"
    }
    radar_data <- prep_radar_data(df, mcols, stat_fn)
    cols <- cond_colors()
    base_sz <- input$mf_base_size %||% 12
    show_pts <- isTRUE(input$mf_radar_points)
    show_grd <- isTRUE(input$mf_radar_grid)
    alpha_f  <- input$mf_radar_alpha %||% 0.25
    p <- plot_radar(radar_data, cols, base_sz, show_pts, show_grd, alpha_f)
    p + labs(title = "Spider / Radar Plot",
             subtitle = paste0("Features scaled [0-1], ", stat_fn, " — ", time_lbl))
  })

  output$mf_radar_plot <- renderPlot({ mf_radar_obj() })
  output$dl_mf_radar <- downloadHandler(
    filename = function() paste0("radar_plot.", fmt_ext(input$mf_format)),
    content = function(file) save_figure(file, mf_radar_obj(), input$mf_width, input$mf_height, input$mf_dpi, input$mf_format))

  # === CELL METRICS TAB ===
  cm_data <- reactive({
    req(raw_data()); df <- transformed_data(); m <- input$measure
    compute_cell_metrics(df, m)
  })

  # Dynamic UI: scatter feature selectors
  # Dynamic UI: measurement selector for Cell Metrics tab
  output$cm_measure_ui <- renderUI({
    req(measure_cols()); mc <- measure_cols()
    global_m <- input$measure %||% mc[1]
    selectInput("cm_measure", "Measurement feature",
      choices = mc, selected = if (global_m %in% mc) global_m else mc[1])
  })

  output$cm_scatter_features <- renderUI({
    req(measure_cols()); mc <- measure_cols()
    tagList(
      selectInput("cm_scatter_x", "X-axis measurement", choices = mc, selected = mc[1]),
      selectInput("cm_scatter_y", "Y-axis measurement", choices = mc,
                  selected = if (length(mc) >= 2) mc[2] else mc[1])
    )
  })

  # Dynamic UI: radar timepoint selector
  output$cm_radar_time_ui <- renderUI({
    req(raw_data()); times <- sort(unique(raw_data()$Time_min))
    selectInput("cm_radar_time", "Timepoint",
      choices = c("All (aggregated)" = "all", setNames(times, paste0(times, " min"))),
      selected = "all")
  })

  cm_plot_obj <- reactive({
    req(raw_data(), input$cm_measure)
    df <- transformed_data(); cols <- cond_colors()
    # Use Cell Metrics-specific measure (not global input$measure)
    m <- input$cm_measure; req(m %in% names(df))
    mlbl <- m
    dat <- compute_cell_metrics(df, m)
    rs <- compute_replicate_summaries(df, m)
    ctype <- input$cm_chart_type %||% "violin"
    base_sz <- input$cm_base_size %||% 12
    legend_pos <- input$cm_legend_pos %||% "top"
    grouping <- input$cm_grouping %||% "by_time"
    metric <- input$cm_metric %||% "mean_val"
    metric_label <- switch(metric,
      mean_val = "Mean", median_val = "Median", sd_val = "SD",
      iqr_val = "IQR", cv_val = "CV (%)", skewness = "Skewness",
      n_cells = "N Cells", q25 = "Q25", q75 = "Q75", metric)

    # Prepare factor columns (consistent " min" suffix format)
    all_times <- sort(unique(df$Time_min))
    time_levels <- fmt_time(all_times)
    dat$Time_fac <- factor(fmt_time(dat$Time_min), levels = time_levels)
    df$Time_fac <- factor(fmt_time(df$Time_min), levels = time_levels)
    rs$Time_fac <- factor(fmt_time(rs$Time_min), levels = time_levels)

    # --- Grouping axis logic ---
    if (grouping == "by_cond") {
      x_var <- "Condition"; dodge_var <- "Time_fac"; x_lab <- "Condition"
    } else {
      x_var <- "Time_fac"; dodge_var <- "Condition"; x_lab <- "Time"
    }

    # ==== SUMMARY-BASED CHARTS (bar, line, heatmap, bubble) ====
    if (ctype == "bar") {
      if (grouping == "by_cond") {
        p <- ggplot(dat, aes(x = Condition, y = .data[[metric]], fill = Time_fac)) +
          geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.85) +
          labs(x = "Condition", y = metric_label, fill = "Time",
               title = paste(metric_label, "by Condition"))
      } else {
        p <- ggplot(dat, aes(x = Time_fac, y = .data[[metric]], fill = Condition)) +
          geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.85) +
          scale_fill_manual(values = cols) +
          labs(x = "Time (min)", y = metric_label, title = paste(metric_label, "across Conditions"))
      }
      if (isTRUE(input$cm_show_values))
        p <- p + geom_text(aes(label = round(.data[[metric]], 2)),
          position = position_dodge(0.8), vjust = -0.3, size = base_sz / 4)

    } else if (ctype == "line") {
      if (grouping == "by_cond") {
        p <- ggplot(dat, aes(x = Condition, y = .data[[metric]],
                             color = Time_fac, group = Time_fac)) +
          geom_line(linewidth = 1) + geom_point(size = 3) +
          labs(x = "Condition", y = metric_label, color = "Time",
               title = paste(metric_label, "by Condition"))
      } else {
        p <- ggplot(dat, aes(x = Time_min, y = .data[[metric]], color = Condition, group = Condition)) +
          geom_line(linewidth = 1) + geom_point(size = 3) +
          scale_color_manual(values = cols) +
          scale_x_continuous(breaks = sort(unique(dat$Time_min))) +
          labs(x = "Time (min)", y = metric_label, title = paste(metric_label, "Trajectory"))
      }
      if (isTRUE(input$cm_show_values))
        p <- p + geom_text(aes(label = round(.data[[metric]], 2)), vjust = -1, size = base_sz / 4)

    } else if (ctype == "heatmap") {
      # A true 2-D grid (conditions x timepoints) rather than every group
      # squeezed onto one axis, which made the labels unreadable.
      if (grouping == "by_cond") {
        p <- ggplot(dat, aes(x = Condition, y = Time_fac, fill = .data[[metric]]))
      } else {
        dat$Condition <- factor(dat$Condition, levels = rev(levels(factor(dat$Condition))))
        p <- ggplot(dat, aes(x = Time_fac, y = Condition, fill = .data[[metric]]))
      }
      p <- p + geom_tile(color = "grey90", linewidth = 0.5) +
        scale_fill_viridis_c(option = "plasma", name = metric_label) +
        labs(x = NULL, y = NULL, title = paste(metric_label, "Heatmap"),
             subtitle = "Each tile = one condition at one timepoint")
      if (isTRUE(input$cm_show_values))
        p <- p + geom_text(aes(label = round(.data[[metric]], 2)), size = base_sz / 3.5,
                           color = "white", fontface = "bold")

    } else if (ctype == "bubble") {
      size_metric <- input$cm_bubble_size %||% "n_cells"
      size_label <- switch(size_metric, n_cells = "N Cells", sd_val = "SD", cv_val = "CV (%)", size_metric)
      if (grouping == "by_cond") {
        p <- ggplot(dat, aes(x = Condition, y = .data[[metric]], color = Time_fac,
                             size = .data[[size_metric]])) +
          geom_point(alpha = 0.75) +
          scale_size_area(max_size = 15, name = size_label) +
          labs(x = "Condition", y = metric_label, color = "Time",
               title = paste(metric_label, "Bubble Plot"),
               subtitle = paste("Bubble size =", size_label))
      } else {
        p <- ggplot(dat, aes(x = Time_fac, y = .data[[metric]], color = Condition,
                             size = .data[[size_metric]])) +
          geom_point(alpha = 0.75) + scale_color_manual(values = cols) +
          scale_size_area(max_size = 15, name = size_label) +
          labs(x = "Time (min)", y = metric_label,
               title = paste(metric_label, "Bubble Plot"),
               subtitle = paste("Bubble size =", size_label))
      }

    # ==== CELL-LEVEL CHARTS (violin, boxplot, ridgeline, ecdf) ====
    } else if (ctype == "violin") {
      if (grouping == "by_cond") {
        p <- ggplot(df, aes(x = Condition, y = .data[[m]], fill = Time_fac)) +
          geom_violin(position = position_dodge(0.8), alpha = 0.7, scale = "width",
                      linewidth = 0.3, color = "grey30") +
          stat_summary(aes(group = Time_fac), fun = median, geom = "crossbar",
                       position = position_dodge(0.8), width = 0.35, linewidth = 0.4, color = "black") +
          labs(x = "Condition", y = mlbl, fill = "Time", title = paste(mlbl, "Violin Plot"))
      } else {
        p <- ggplot(df, aes(x = Time_fac, y = .data[[m]], fill = Condition)) +
          geom_violin(position = position_dodge(0.8), alpha = 0.7, scale = "width",
                      linewidth = 0.3, color = "grey30") +
          stat_summary(aes(group = Condition), fun = median, geom = "crossbar",
                       position = position_dodge(0.8), width = 0.35, linewidth = 0.4, color = "black") +
          scale_fill_manual(values = cols) +
          labs(x = "Time", y = mlbl, title = paste(mlbl, "Violin Plot"))
      }

    } else if (ctype == "boxplot") {
      if (grouping == "by_cond") {
        p <- ggplot(df, aes(x = Condition, y = .data[[m]], fill = Time_fac)) +
          geom_boxplot(position = position_dodge(0.8), alpha = 0.65, outlier.alpha = 0.15,
                       outlier.size = 0.5, linewidth = 0.4) +
          labs(x = "Condition", y = mlbl, fill = "Time", title = paste(mlbl, "Boxplot"))
      } else {
        p <- ggplot(df, aes(x = Time_fac, y = .data[[m]], fill = Condition)) +
          geom_boxplot(position = position_dodge(0.8), alpha = 0.65, outlier.alpha = 0.15,
                       outlier.size = 0.5, linewidth = 0.4) +
          scale_fill_manual(values = cols) +
          labs(x = "Time", y = mlbl, title = paste(mlbl, "Boxplot"))
      }

    } else if (ctype == "ridgeline") {
      if (has_ridges) {
        if (grouping == "by_cond") {
          p <- ggplot(df, aes(x = .data[[m]], y = Condition, fill = Time_fac)) +
            ggridges::geom_density_ridges(alpha = 0.55, scale = 1.2, linewidth = 0.3) +
            facet_wrap(~ Time_fac) +
            labs(x = mlbl, y = "Condition", fill = "Time",
                 title = paste(mlbl, "Ridgeline Density"),
                 subtitle = "Distribution per Condition, faceted by Timepoint")
        } else {
          p <- ggplot(df, aes(x = .data[[m]], y = Time_fac, fill = Condition)) +
            ggridges::geom_density_ridges(alpha = 0.55, scale = 1.2, linewidth = 0.3) +
            scale_fill_manual(values = cols) + facet_wrap(~ Condition) +
            labs(x = mlbl, y = "Time", title = paste(mlbl, "Ridgeline Density"),
                 subtitle = "Distribution per Condition per Timepoint")
        }
      } else {
        p <- ggplot() + annotate("text", x = 0.5, y = 0.5,
          label = "Install ggridges: install.packages('ggridges')", size = 5) + theme_void()
      }

    } else if (ctype == "ecdf") {
      if (grouping == "by_cond") {
        p <- ggplot(df, aes(x = .data[[m]], color = Time_fac)) +
          stat_ecdf(linewidth = 0.8) +
          facet_wrap(~ Condition, nrow = 1) +
          labs(x = mlbl, y = "Cumulative Proportion", color = "Time",
               title = paste(mlbl, "ECDF"),
               subtitle = "Faceted by Condition, colored by Timepoint")
      } else {
        p <- ggplot(df, aes(x = .data[[m]], color = Condition)) +
          stat_ecdf(linewidth = 0.8) + scale_color_manual(values = cols) +
          facet_wrap(~ Time_fac, nrow = 1) +
          labs(x = mlbl, y = "Cumulative Proportion", title = paste(mlbl, "ECDF"),
               subtitle = "Rightward shift = larger values")
      }

    # ==== TWO-FEATURE CHARTS (scatter, hexbin) ====
    } else if (ctype == "scatter") {
      req(input$cm_scatter_x, input$cm_scatter_y)
      fx <- input$cm_scatter_x; fy <- input$cm_scatter_y
      p <- ggplot(df, aes(x = .data[[fx]], y = .data[[fy]], color = Condition)) +
        geom_point(alpha = 0.15, size = 0.8) + scale_color_manual(values = cols) +
        facet_wrap(~ Time_fac, nrow = 1) +
        labs(x = fx, y = fy, title = paste(fy, "vs", fx),
             subtitle = "Each point = one cell, faceted by timepoint")

    } else if (ctype == "hexbin") {
      req(input$cm_scatter_x, input$cm_scatter_y)
      fx <- input$cm_scatter_x; fy <- input$cm_scatter_y
      hex_geom <- if (has_hexbin) geom_hex(bins = 40) else geom_bin2d(bins = 40)
      p <- ggplot(df, aes(x = .data[[fx]], y = .data[[fy]])) +
        hex_geom + scale_fill_viridis_c(option = "inferno", name = "Count") +
        facet_grid(Condition ~ Time_fac) +
        labs(x = fx, y = fy, title = paste(fy, "vs", fx, "- Density"),
             subtitle = "Color = cell density per bin")

    # ==== SPIDER / RADAR (in Cell Metrics context) ====
    } else if (ctype == "radar") {
      mcols <- measure_cols(); req(length(mcols) >= 3)
      sel_time <- input$cm_radar_time %||% "all"
      if (sel_time != "all") {
        df_filt <- df[df$Time_min == as.numeric(sel_time), , drop = FALSE]
        radar_data <- prep_radar_data(df_filt, mcols, "mean")
        time_lbl <- paste0(sel_time, " min")
      } else {
        radar_data <- prep_radar_data(df, mcols, "mean")
        time_lbl <- "All timepoints (aggregated)"
      }
      p <- plot_radar(radar_data, cols, base_sz, show_points = TRUE,
                      show_grid = TRUE, alpha_fill = 0.25) +
        labs(title = "Morphological Profile (Spider / Radar)",
             subtitle = paste0("Features scaled [0-1], ", time_lbl))
    }

    # --- Apply significance brackets for bracket-compatible charts ---
    bracket_types <- c("bar", "violin", "boxplot")
    if (ctype %in% bracket_types && isTRUE(input$cm_show_sig) && grouping == "by_time" &&
        !is.null(rs) && nrow(rs) > 0) {
      compare_mode <- input$cm_sig_mode %||% "vs_control"
      sig_yvar <- input$cm_sig_stat %||% "median_val"
      cond_levels <- levels(df$Condition)
      if (is.null(cond_levels)) cond_levels <- unique(df$Condition)
      if (length(cond_levels) >= 2) {
        # Determine pairs
        if (compare_mode == "all_pairwise") {
          pairs_list <- combn(cond_levels, 2, simplify = FALSE)
        } else {
          ref <- cond_levels[1]
          pairs_list <- lapply(cond_levels[-1], function(x) c(ref, x))
        }
        # Compute significance per timepoint (test/correction follow Statistics tab)
        times_fac_levels <- levels(df$Time_fac)
        if (is.null(times_fac_levels)) times_fac_levels <- unique(as.character(df$Time_fac))
        sig_res <- compute_sig_table(rs, pairs_list,
          times = sort(unique(rs$Time_min)),
          sig_yvar = sig_yvar,
          test = input$stat_test %||% "welch",
          correction = input$stat_p_adjust %||% "holm",
          label_style = input$cm_sig_label %||% "stars",
          show_ns = isTRUE(input$cm_show_ns),
          y_lookup = function(tm, pair, v1, v2) {
            if (ctype %in% c("violin", "boxplot")) max(df[[m]][df$Time_min == tm], na.rm = TRUE)
            else max(dat[[metric]][dat$Time_min == tm], na.rm = TRUE)
          })
        sig_annot <- sig_res$annot
        if (nrow(sig_annot) > 0) sig_annot$Time_fac <- fmt_time(sig_annot$Time_min)
        {
          if (nrow(sig_annot) > 0) {
            if (ctype %in% c("violin", "boxplot")) {
              y_range <- diff(range(df[[m]], na.rm = TRUE))
            } else {
              y_range <- diff(range(dat[[metric]], na.rm = TRUE))
            }
            bracket_base <- sig_annot$y_pos + y_range * 0.06
            bracket_step <- y_range * (if ((input$cm_sig_label %||% "stars") == "stars_p") 0.17 else 0.10)
            for (tf in unique(sig_annot$Time_fac)) {
              idx <- which(sig_annot$Time_fac == tf)
              if (length(idx) > 1) {
                for (j in 2:length(idx))
                  bracket_base[idx[j]] <- bracket_base[idx[j - 1]] + bracket_step
              }
            }
            # Render brackets
            n_levels <- length(cond_levels)
            for (i in seq_len(nrow(sig_annot))) {
              c1i <- match(sig_annot$Cond1[i], cond_levels)
              c2i <- match(sig_annot$Cond2[i], cond_levels)
              ti <- match(sig_annot$Time_fac[i], times_fac_levels)
              if (n_levels == 2) {
                x1 <- ti - 0.2; x2 <- ti + 0.2
              } else {
                dw <- 0.8; step <- dw / n_levels
                x1 <- ti - dw / 2 + (c1i - 0.5) * step
                x2 <- ti - dw / 2 + (c2i - 0.5) * step
              }
              y_bar <- bracket_base[i]; y_tick <- y_bar - y_range * 0.015
              p <- p +
                annotate("segment", x = x1, xend = x1, y = y_tick, yend = y_bar, linewidth = 0.35) +
                annotate("segment", x = x1, xend = x2, y = y_bar, yend = y_bar, linewidth = 0.35) +
                annotate("segment", x = x2, xend = x2, y = y_tick, yend = y_bar, linewidth = 0.35) +
                annotate("text", x = (x1 + x2) / 2, y = y_bar + y_range * 0.012,
                         label = sig_annot$label[i], size = base_sz / 4, vjust = 0)
            }
            label_extra <- if ((input$cm_sig_label %||% "stars") == "stars_p") 0.20 else 0.12
            max_bracket <- max(bracket_base) + y_range * label_extra
            p <- p + coord_cartesian(ylim = c(NA, max_bracket), clip = "off")
            p <- p + labs(caption = paste0(sig_res$caption, bracket_family_note() %||% ""))
          }
        }
      }
    }

    # Apply theme
    p <- p + get_theme(input, "cm") +
      if (ctype %in% c("scatter", "hexbin", "radar", "ecdf", "ridgeline")) labs() else
        labs(subtitle = paste0(mlbl, if (ctype %in% bracket_types) "" else paste0(" — ", metric_label)))
    apply_custom_labs(p, input, "cm")
  })

  output$cm_chart <- renderPlot({ cm_plot_obj() })

  output$cm_table <- renderDT({
    req(raw_data(), input$cm_measure)
    m_cm <- input$cm_measure; req(m_cm %in% names(transformed_data()))
    dat <- compute_cell_metrics(transformed_data(), m_cm)
    dat[] <- lapply(dat, function(x) if (is.numeric(x)) round(x, 4) else x)
    datatable(dat, options = list(pageLength = 20, scrollX = TRUE, dom = 'tip'),
      caption = paste("Cell Metrics for", m_cm))
  })

  output$cm_about <- renderUI({
    tags$div(class = "p-3",
      tags$h4("Cell Metrics Explorer"),
      tags$p("Population-level morphometric analysis for cell biology. Compute summary statistics,
              visualize distributions, compare conditions, and test for significance.
              Charts follow best practices from high-impact cell morphology publications."),
      tags$h5("Summary Charts (population statistics)"),
      tags$ul(
        tags$li(tags$strong("Bar Chart:"), " Grouped bars. Clear comparison of summary metrics at each timepoint."),
        tags$li(tags$strong("Line Chart:"), " Trajectories showing metric evolution over time."),
        tags$li(tags$strong("Heatmap:"), " Color intensity grid for quick multi-group overview."),
        tags$li(tags$strong("Bubble Plot:"), " Size encodes a second metric (N Cells, SD, CV).")
      ),
      tags$h5("Distribution Charts (cell-level data)"),
      tags$ul(
        tags$li(tags$strong("Violin Plot:"), " Best overall summary. Shows distribution shape, density,
                 median, and variability. Reveals multimodal populations and skewness."),
        tags$li(tags$strong("Boxplot:"), " Classic statistical summary: median, Q1/Q3, whiskers, outliers.
                 Good for clean publication-ready comparisons."),
        tags$li(tags$strong("Ridgeline:"), " Stacked density plots. Great for visualizing population
                 dynamics across many timepoints or strains."),
        tags$li(tags$strong("ECDF:"), " Empirical cumulative distribution. Very powerful for detecting
                 distributional shifts across the full range. Underused but highly informative.")
      ),
      tags$h5("Feature Relationship Charts"),
      tags$ul(
        tags$li(tags$strong("Scatter (2 features):"), " Shows morphological relationships (e.g., length vs width).
                 Reveals constraints, scaling, and subpopulations."),
        tags$li(tags$strong("Hexbin (2 features):"), " 2D density for large populations. Color = cell count
                 per hex. Reveals morphological regimes without overplotting.")
      ),
      tags$h5("Multi-Feature"),
      tags$ul(
        tags$li(tags$strong("Spider / Radar:"), " Morphological profile fingerprint. Each axis is a feature
                 (scaled 0-1). Each condition creates a shape signature for quick phenotype comparison.
                 Can filter to specific timepoints.")
      ),
      tags$h5("Significance Testing"),
      tags$p("Brackets show Welch's t-test on biological replicate summaries (not individual cells)
              to avoid pseudoreplication. Holm-Bonferroni corrected for multiple comparisons.
              Available on violin, boxplot, and bar chart types in 'By timepoint' grouping."),
      tags$h5("Grouping Control"),
      tags$ul(
        tags$li(tags$strong("By timepoint:"), " X-axis = Time, compare conditions side-by-side at each timepoint."),
        tags$li(tags$strong("By condition:"), " X-axis = Condition, compare timepoints within each condition.
                 Controls on one side, treatments on the other.")
      ),
      tags$p(class = "small", style = "color:#8888a8; margin-top:1rem;",
        tags$em("Tip: Always show distribution plots, not just means. Means hide heterogeneity,
                 bimodal populations, and rare phenotypes."))
    )
  })

  output$dl_cm_plot <- downloadHandler(
    filename = function() paste0("cell_metrics.", fmt_ext(input$cm_format)),
    content = function(file) save_figure(file, cm_plot_obj(), input$cm_width, input$cm_height, input$cm_dpi, input$cm_format))
  output$dl_cm_table <- downloadHandler(
    filename = function() paste0("cell_metrics_", input$cm_measure, "_", Sys.Date(), ".csv"),
    content = function(file) {
      req(input$cm_measure, transformed_data())
      dat <- compute_cell_metrics(transformed_data(), input$cm_measure)
      write.csv(dat, file, row.names = FALSE)
    })

  # === PUBLICATION COMPOSITE FIGURE ===
  pub_panels_list <- reactive({
    req(input$pub_panels); panels <- list()
    base <- input$pub_base_size %||% 11; df <- transformed_data(); rs <- rep_summ()
    cols <- cond_colors(); m <- input$measure; mlbl <- measure_label()
    if ("dist" %in% input$pub_panels) {
      df$Time_fac <- factor(fmt_time(df$Time_min), levels = fmt_time(sort(unique(df$Time_min))))
      rs$Time_fac <- factor(fmt_time(rs$Time_min), levels = fmt_time(sort(unique(rs$Time_min))))
      panels[["dist"]] <- ggplot(df, aes(x = Time_fac, y = .data[[m]], fill = Condition)) +
        geom_violin(position = position_dodge(0.8), alpha = 0.7, scale = "width", linewidth = 0.3, color = "grey30") +
        geom_point(data = rs, aes(x = Time_fac, y = median_val, fill = Condition, group = Condition),
                   position = position_dodge(0.8), shape = 21, size = 2.5, color = "black", stroke = 0.7) +
        scale_fill_manual(values = cols) + labs(x = NULL, y = mlbl, fill = "Condition") +
        theme_publication(base) + theme(legend.position = "bottom")
    }
    if ("traj" %in% input$pub_panels) {
      ptj <- ggplot(rs, aes(x = Time_min, y = median_val, color = Condition)) +
        geom_line(aes(group = interaction(Condition, Replicate)), alpha = 0.35, linewidth = 0.4, linetype = "dashed") +
        geom_point(aes(shape = Replicate), size = 2.5, alpha = 0.8) + scale_color_manual(values = cols) +
        scale_x_continuous(breaks = sort(unique(rs$Time_min))) + labs(x = "Time (min)", y = paste("Median", mlbl)) + theme_publication(base)
      if (length(unique(rs$Time_min)) >= 4) ptj <- ptj + geom_smooth(method = "loess", se = TRUE, alpha = 0.12, linewidth = 1.1, span = 1)
      else ptj <- ptj + stat_summary(fun = mean, geom = "line", linewidth = 1.1)
      # No in-panel legend: the composite carries one shared legend at the bottom
      panels[["traj"]] <- ptj + theme(legend.position = "none")
    }
    if ("effect" %in% input$pub_panels) {
      esd <- es_data()
      esd_trt <- esd[!esd$is_control, , drop = FALSE]
      ref <- input$ref_condition %||% names(cols)[1]
      cols_trt <- cols[!names(cols) %in% ref]
      panels[["effect"]] <- ggplot(esd_trt, aes(x = Time_min, y = pct_change, color = Condition)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
        geom_line(linewidth = 1.1) + geom_point(size = 3.5) + scale_color_manual(values = cols_trt) +
        scale_x_continuous(breaks = sort(unique(esd$Time_min))) + labs(x = "Time (min)", y = paste("% Delta", mlbl)) +
        theme_publication(base) + theme(legend.position = "none")
    }
    if ("morphospace" %in% input$pub_panels) {
      res <- pca_result()
      if (!is.null(res) && res$ok) {
        sc <- as.data.frame(res$pca$x[, 1:min(ncol(res$pca$x), 2)]); sc$Condition <- res$df$Condition
        vpct <- round(res$var_explained[2, 1:min(ncol(res$pca$x), 2)] * 100, 1)
        panels[["morphospace"]] <- ggplot(sc, aes(x = PC1, y = PC2, color = Condition)) +
          geom_point(alpha = 0.15, size = 0.5) + scale_color_manual(values = cols) +
          labs(x = paste0("PC1 (", vpct[1], "%)"), y = paste0("PC2 (", vpct[2], "%)")) +
          theme_publication(base) + theme(legend.position = "none")
      }
    }
    if ("phenotype" %in% input$pub_panels) {
      cr <- cluster_result()
      if (!is.null(cr)) {
        comp <- cr$df %>% group_by(Condition, Time_min, Phenotype) %>% summarise(n = n(), .groups = "drop") %>%
          group_by(Condition, Time_min) %>% mutate(pct = n / sum(n) * 100) %>% ungroup()
        comp$Time_fac <- factor(fmt_time(comp$Time_min), levels = fmt_time(sort(unique(comp$Time_min))))
        # Facet strips shrink in a composite panel, so scale their text down to
        # stop neighbouring strip labels from colliding.
        panels[["phenotype"]] <- ggplot(comp, aes(x = Time_fac, y = pct, fill = Phenotype)) +
          geom_col(position = "stack", width = 0.7) +
          facet_wrap(~ Condition, labeller = labeller(Condition = function(x) wrap_axis_labels(x, 14))) +
          labs(x = NULL, y = "% Cells") + theme_publication(base) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none",
                strip.text = element_text(size = rel(0.62), lineheight = 0.9),
                panel.spacing = unit(4, "pt"))
      }
    }
    if ("cellmetrics" %in% input$pub_panels) {
      cm <- compute_cell_metrics(df, m)
      cm$Time_fac <- factor(cm$Time_min)
      panels[["cellmetrics"]] <- ggplot(cm, aes(x = Time_min, y = mean_val, color = Condition, group = Condition)) +
        geom_line(linewidth = 1) + geom_point(size = 3) + scale_color_manual(values = cols) +
        labs(x = "Time (min)", y = paste("Mean", mlbl)) + theme_publication(base) + theme(legend.position = "none")
    }
    if ("radar" %in% input$pub_panels) {
      mcols <- measure_cols()
      if (length(mcols) >= 3) {
        radar_data <- prep_radar_data(df, mcols, "mean")
        panels[["radar"]] <- plot_radar(radar_data, cols, base_sz = base, show_points = TRUE,
          show_grid = TRUE, alpha_fill = 0.2) + theme(legend.position = "none")
      }
    }
    if ("morph_traj" %in% input$pub_panels) {
      res <- pca_result()
      if (!is.null(res) && res$ok) {
        sc2 <- as.data.frame(res$pca$x[, 1:2]); sc2$Condition <- res$df$Condition; sc2$Time_min <- res$df$Time_min
        cen2 <- sc2 %>% group_by(Condition, Time_min) %>% summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop") %>% arrange(Condition, Time_min)
        vpct2 <- round(res$var_explained[2, 1:2] * 100, 1)
        panels[["morph_traj"]] <- ggplot(cen2, aes(x = PC1, y = PC2, color = Condition)) +
          geom_path(aes(group = Condition), linewidth = 1.2, arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
          geom_point(size = 3) + geom_text(aes(label = paste0(Time_min, "m")), vjust = -1, size = 2.5) +
          scale_color_manual(values = cols) + labs(x = paste0("PC1 (", vpct2[1], "%)"), y = paste0("PC2 (", vpct2[2], "%)")) +
          theme_publication(base) + theme(legend.position = "none")
      }
    }
    if ("morph_div" %in% input$pub_panels) {
      dd <- tryCatch(pca_diversity_data(), error = function(e) NULL)
      if (!is.null(dd) && nrow(dd) > 0) {
        dd$Time_fac <- factor(fmt_time(dd$Time_min), levels = fmt_time(sort(unique(dd$Time_min))))
        panels[["morph_div"]] <- ggplot(dd, aes(x = Time_fac, y = Hull_area, fill = Condition)) +
          geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.85) +
          scale_fill_manual(values = cols) + labs(x = "Time", y = "Hull Area") +
          theme_publication(base) + theme(legend.position = "none")
      }
    }
    panels
  })

  pub_plot_obj <- reactive({
    req(has_patchwork)
    panels <- pub_panels_list()
    if (length(panels) == 0) return(NULL)
    base <- input$pub_base_size %||% 11; mlbl <- measure_label()
    # Exactly one panel keeps its legend (the first); patchwork then collects it
    # into a single shared legend below the figure. Whichever panels the user
    # picks, no panel ends up with a legend spilling over its neighbours.
    labeled <- lapply(seq_along(panels), function(i) {
      panels[[i]] + ggtitle(LETTERS[i]) +
        theme(legend.position = if (i == 1) "bottom" else "none")
    })
    n_panels <- length(labeled)
    pub_title <- {
      ct <- input$pub_custom_title %||% ""
      if (nzchar(trimws(ct))) ct else paste(mlbl, "-- Time Course")
    }
    # guides = "collect" lifts the single kept legend out of its panel and into
    # one shared legend under the whole figure, so it cannot overflow a panel.
    patchwork::wrap_plots(labeled, ncol = if (n_panels == 1) 1 else 2, guides = "collect") +
      patchwork::plot_annotation(
        title = pub_title,
        theme = theme_publication(base) +
          theme(legend.position = "bottom", legend.box = "horizontal",
                legend.justification = "center"))
  })

  output$pub_plot <- renderPlot({ if (has_patchwork) pub_plot_obj() })
  output$dl_pub <- downloadHandler(filename = function() paste0("figure_", input$measure, ".", fmt_ext(input$pub_format)),
    content = function(file) save_figure(file, pub_plot_obj(), input$pub_width, input$pub_height, input$pub_dpi, input$pub_format))

  # === GUIDE / SCRIPT BUILDER =================================================
  output$guide_status <- renderUI({
    if (is.null(raw_data()))
      return(tags$div(class = "alert alert-info py-2 px-3",
        icon("circle-info"), " No data loaded yet. Start at step 1 below, or click ",
        tags$strong("Load Demo Data"), " on the Data tab."))
    tags$div(class = "alert alert-success py-2 px-3",
      icon("circle-check"), sprintf(" %s cells loaded across %d sample(s)%s.",
        format(nrow(raw_data()), big.mark = ","),
        length(unique(raw_data()$Condition)),
        if (has_time()) sprintf(" and %d timepoint(s)", length(time_values())) else ""))
  })

  guide_report_data <- eventReactive(input$guide_check, {
    req(raw_data())
    df <- raw_data()
    reps <- df %>% dplyr::group_by(Condition, Time_min) %>%
      dplyr::summarise(n_reps = dplyr::n_distinct(Replicate),
                       n_cells = dplyr::n(), .groups = "drop")
    list(df = df, reps = reps,
         min_reps = min(reps$n_reps), max_reps = max(reps$n_reps),
         n_cond = length(unique(df$Condition)),
         measures = measure_cols(), has_time = has_time(),
         n_times = length(time_values()))
  })

  output$guide_report <- renderUI({
    r <- guide_report_data()
    item <- function(ok, title, detail) {
      cls <- if (ok == "good") "alert-success" else if (ok == "warn") "alert-warning" else "alert-danger"
      ic  <- if (ok == "good") "circle-check" else if (ok == "warn") "triangle-exclamation" else "circle-xmark"
      tags$div(class = paste("alert py-2 px-3 mb-2", cls),
               icon(ic), tags$strong(paste0(" ", title)), tags$br(), tags$span(class = "small", detail))
    }
    tagList(
      item("good", sprintf("%d samples, %d measurement columns", r$n_cond, length(r$measures)),
           paste("Measurements found:", paste(utils::head(r$measures, 12), collapse = ", "),
                 if (length(r$measures) > 12) "..." else "")),
      if (r$has_time)
        item("good", sprintf("Time course with %d timepoints", r$n_times),
             "Trajectories, time-resolved statistics and the trend test are all available.")
      else
        item("good", "No time column - single snapshot design",
             "Samples are compared directly. The Trajectories tab does not apply."),
      if (r$min_reps >= 3)
        item("good", sprintf("%d-%d biological replicates per group", r$min_reps, r$max_reps),
             "Enough for pairwise tests. Welch's t-test with Holm correction is a sound default.")
      else if (r$min_reps == 2)
        item("warn", "Only 2 biological replicates in at least one group",
             paste("Pairwise p-values are weak here (a t-test has about 1 degree of freedom, and a",
                   "Wilcoxon test cannot reach p < 0.05 at all). Prefer the linear mixed model, a",
                   "dose-response trend test, or present a SuperPlot without p-values."))
      else
        item("bad", "Only 1 biological replicate in at least one group",
             paste("No hypothesis test is valid for that group - there is nothing to estimate",
                   "sample-to-sample variability from. You can still plot the distributions.")),
      if (r$max_reps >= 2 * max(1, r$min_reps))
        item("warn", "Replicate counts are unbalanced across groups",
             "Interpret comparisons cautiously; Welch's t-test is the right choice here.") else NULL,
      tags$p(tags$strong("Cells and replicates per group"), class = "mt-3 mb-1"),
      renderTable(head(as.data.frame(r$reps), 40))
    )
  })

  output$script_measure_ui <- renderUI({
    req(measure_cols())
    selectInput("script_measure", "Measurement", choices = measure_cols(),
                selected = input$measure %||% measure_cols()[1])
  })
  output$script_group_ui <- renderUI({
    req(raw_data())
    meta <- detect_metadata_columns(raw_data(), measure_cols())
    sel <- intersect(c("Condition", if (has_time()) "Time_min"), meta)
    checkboxGroupInput("script_groups", "Group comparisons by",
                       choices = meta, selected = if (length(sel)) sel else meta[1])
  })
  output$script_trend_pred_ui <- renderUI({
    req(raw_data())
    num <- names(raw_data())[vapply(raw_data(), is.numeric, logical(1))]
    num <- setdiff(num, c(measure_cols(), "Replicate"))
    if (!length(num)) return(tags$p(class = "small text-muted",
      "No numeric predictor column yet - use 'Split Condition into factors' on the Data tab to create one."))
    selectInput("script_trend_pred", "Numeric predictor",
                choices = num,
                selected = grep("dose|conc", num, ignore.case = TRUE, value = TRUE)[1] %||% num[1])
  })

  build_analysis_script <- function() {
    m      <- input$script_measure %||% (measure_cols() %||% "Length")[1]
    groups <- input$script_groups %||% "Condition"
    unit   <- input$script_stat_unit %||% "median_val"
    test   <- input$script_test %||% "welch"
    corr   <- input$script_corr %||% "holm"
    fmt    <- input$script_fmt %||% "pdf"
    plot_t <- input$script_plot %||% "violin"
    ref    <- input$ref_condition %||% (unique(as.character(raw_data()$Condition))[1])
    all_m  <- isTRUE(input$script_all_measures)
    q  <- function(x) paste0('"', x, '"')
    qv <- function(x) paste0("c(", paste(vapply(x, q, character(1)), collapse = ", "), ")")
    unit_lab <- if (unit == "mean_val") "mean" else "median"
    geom <- switch(plot_t,
      violin = '    geom_violin(position = position_dodge(0.8), alpha = 0.7, scale = "width",\n                linewidth = 0.3, color = "grey30") +',
      sina   = '    geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),\n                size = 0.4, alpha = 0.12) +',
      box    = '    geom_boxplot(position = position_dodge(0.8), alpha = 0.6, outlier.shape = NA,\n                 linewidth = 0.4) +',
      dot    = "")

    L <- c(
      "# =============================================================================",
      paste0("# CellMorphR analysis script - generated ", format(Sys.Date())),
      "#",
      "# Run it with:   Rscript my_analysis.R",
      "# Keep this file next to CellMorphR_v3.R (it reuses the app's own functions,",
      "# so the numbers and the figure match what you saw in the app exactly).",
      "# =============================================================================",
      "",
      'DATA_FILE <- "REPLACE_WITH_YOUR_FILE.csv"   # <- your data',
      paste0('OUT_DIR   <- "cellmorphr_output"'),
      "",
      "suppressPackageStartupMessages(invisible(source(\"CellMorphR_v3.R\")))",
      "dir.create(OUT_DIR, showWarnings = FALSE)",
      "",
      "# ---- 1. Load and validate ---------------------------------------------------",
      "raw <- read.csv(DATA_FILE, check.names = FALSE)",
      "v   <- validate_data(raw)",
      'if (!isTRUE(v$ok)) stop("Data problem: ", v$msg)',
      "d <- v$df",
      'd$Replicate <- as.factor(d$Replicate)',
      'cat("Loaded", nrow(d), "cells,", length(v$measures), "measurements\\n")',
      "",
      paste0("MEASURES <- ", if (all_m) "v$measures" else qv(m)),
      paste0("GROUP_COLS <- ", qv(groups)),
      paste0("REFERENCE  <- ", q(ref), "   # control / baseline condition"),
      "",
      "for (MEASURE in MEASURES) {",
      '  cat("\\n=== ", MEASURE, " ===\\n", sep = "")',
      "",
      "  # ---- 2. Replicate-level summaries (this is your n) ------------------------",
      "  rs <- compute_flexible_replicate_summary(d, MEASURE, GROUP_COLS,",
      paste0("                                          yvar = ", q(unit), ")"),
      "",
      "  # ---- 3. Pairwise statistics on replicate summaries ------------------------",
      '  groups_present <- sort(unique(as.character(rs$Statistics_Group)))',
      "  ref_group <- grep(REFERENCE, groups_present, fixed = TRUE, value = TRUE)[1]",
      "  if (is.na(ref_group)) ref_group <- groups_present[1]",
      '  cmp <- build_comparison_table(rs, mode = "vs_reference", reference_group = ref_group)',
      "  res <- run_pairwise_tests(rs, cmp,",
      paste0("                            test = ", q(test), ", correction = ", q(corr), ","),
      "                            reference_group = ref_group, measurement = MEASURE)",
      "  print(res[, c(\"Group_1\", \"Group_2\", \"N_BioRep_Group_1\", \"N_BioRep_Group_2\",",
      "                \"Difference_Group2_minus_Group1\", \"CI95_low\", \"CI95_high\",",
      "                \"p_adj\", \"Significance\")])",
      "  write.csv(res, file.path(OUT_DIR, paste0(\"stats_\", MEASURE, \".csv\")), row.names = FALSE)",
      ""
    )

    if (isTRUE(input$script_trend) && !is.null(input$script_trend_pred)) {
      L <- c(L,
        "  # ---- 3b. Dose-response trend on replicate summaries -----------------------",
        paste0("  PRED <- ", q(input$script_trend_pred)),
        "  tr <- d %>% dplyr::group_by(dplyr::across(dplyr::all_of(c(PRED, \"Replicate\")))) %>%",
        paste0("    dplyr::summarise(val = ", unit_lab, "(.data[[MEASURE]], na.rm = TRUE), .groups = \"drop\")"),
        "  fit <- lm(val ~ tr[[PRED]], data = tr)",
        '  cat("trend slope =", signif(coef(fit)[2], 4),',
        '      " p =", signif(summary(fit)$coefficients[2, 4], 3),',
        '      " R2 =", signif(summary(fit)$r.squared, 3), "\\n")',
        "")
    }

    L <- c(L,
      "  # ---- 4. SuperPlot figure --------------------------------------------------",
      '  # reference condition first, so it sits at the left of each group',
      '  d$Condition <- factor(d$Condition,',
      '                        levels = c(REFERENCE, setdiff(sort(unique(as.character(d$Condition))), REFERENCE)))',
      '  d$Time_fac <- factor(d$Time_min)',
      '  rs_plot <- compute_replicate_summaries(d, MEASURE)',
      '  rs_plot$Time_fac <- factor(rs_plot$Time_min)',
      "  p <- ggplot(d, aes(x = Time_fac, y = .data[[MEASURE]], fill = Condition)) +")
    if (nzchar(geom)) L <- c(L, geom)
    L <- c(L,
      "    geom_point(data = rs_plot,",
      paste0("               aes(x = Time_fac, y = ", unit, ", fill = Condition,"),
      "                   group = Condition, shape = factor(Replicate)),",
      "               position = position_dodge(0.8), size = 2.6, colour = \"black\", stroke = 0.7) +",
      '    scale_shape_manual(values = c(21, 22, 24, 23, 25), name = "Replicate") +',
      "    scale_fill_manual(values = scales::hue_pal()(length(unique(d$Condition)))) +",
      paste0("    labs(x = ", q(if (isolate(has_time())) "Time" else "Sample"),
             ", y = ", q(m), ", fill = \"Condition\","),
      paste0("         title = paste(MEASURE, \"by condition\"),"),
      paste0("         subtitle = \"Points = biological replicate ", unit_lab, "s\","),
      paste0("         caption = \"", if (test == "wilcox") "Wilcoxon rank-sum" else "Welch's t-test",
             " on replicate ", unit_lab, "s, ", corr, "-corrected.\") +"),
      "    theme_publication(12)",
      "",
      paste0("  ggsave(file.path(OUT_DIR, paste0(\"figure_\", MEASURE, \".", fmt, "\")),"),
      "         p, width = 9, height = 5.5, dpi = 300)",
      "}",
      "",
      'cat("\\nDone. Figures and statistics are in ", normalizePath(OUT_DIR), "\\n", sep = "")')
    paste(L, collapse = "\n")
  }

  script_text <- reactiveVal("")
  observeEvent(input$script_build, {
    req(raw_data())
    txt <- tryCatch(build_analysis_script(),
                    error = function(e) paste("# Could not build script:", conditionMessage(e)))
    script_text(txt)
    showNotification("Script generated. Copy it, or use Download .R file.", type = "message")
  })
  output$script_preview <- renderText({
    if (!nzchar(script_text()))
      "# Choose your options on the left, then click 'Generate script'."
    else script_text()
  })
  output$script_download <- downloadHandler(
    filename = function() paste0("cellmorphr_analysis_", Sys.Date(), ".R"),
    content = function(file) {
      txt <- script_text()
      if (!nzchar(txt)) txt <- build_analysis_script()
      writeLines(txt, file)
    })

  # --- Export each composite panel as an editable PPTX slide ---
  output$dl_pub_deck <- downloadHandler(
    filename = function() paste0("figure_panels_", input$measure, ".pptx"),
    content = function(file) {
      if (!(has_officer && has_rvg)) {
        showNotification("PPTX deck export needs the 'officer' and 'rvg' packages (run run_cellmorphr.bat to install).", type = "warning", duration = 10)
        return(save_figure(file, pub_plot_obj(), input$pub_width, input$pub_height, input$pub_dpi, "pdf"))
      }
      panels <- pub_panels_list()
      validate(need(length(panels) > 0, "No panels selected."))
      doc <- officer::read_pptx()
      for (i in seq_along(panels)) {
        doc <- officer::add_slide(doc, layout = "Title and Content", master = "Office Theme")
        doc <- tryCatch({
          d <- officer::ph_with(doc, paste0("Panel ", LETTERS[i], " -- ", names(panels)[i]),
                                location = officer::ph_location_type(type = "title"))
          officer::ph_with(d, rvg::dml(ggobj = panels[[i]]),
                           location = officer::ph_location(left = 0.5, top = 1.2, width = 9, height = 5.8))
        }, error = function(e) doc)
      }
      print(doc, target = file)
    })
}

# =============================================================================
# RUN
# =============================================================================
shinyApp(ui, server)
