# =============================================================================
# CellMorphR launcher
# Started by run_cellmorphr.bat. Installs any missing packages into the user
# library on first run, then launches the Shiny app in the default browser.
# =============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))

# --- Resolve the folder this script lives in and work from there -------------
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) > 0) {
  self_path <- normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")
  setwd(dirname(self_path))
}
message("CellMorphR launcher -- working directory: ", getwd())

# --- Make sure a writable user library exists and is first on the path -------
ulib <- Sys.getenv("R_LIBS_USER")
if (nzchar(ulib)) {
  if (!dir.exists(ulib)) dir.create(ulib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(ulib, .libPaths()))
}

# --- Install anything missing (binary installs; no Rtools needed) ------------
required <- c("shiny", "bslib", "ggplot2", "dplyr", "tidyr", "readxl", "DT",
              "scales", "writexl", "cluster")
optional <- c("lme4", "lmerTest", "emmeans", "ggridges", "ggbeeswarm", "ggforce",
              "patchwork", "colourpicker", "hexbin", "svglite",
              "officer", "rvg")
pkgs <- c(required, optional)
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing) > 0) {
  message("First run: installing ", length(missing), " package(s): ",
          paste(missing, collapse = ", "))
  message("This can take several minutes. Please leave this window open...")
  tryCatch(
    install.packages(missing, lib = if (nzchar(ulib)) ulib else NULL, type = "binary"),
    error = function(e) message("Package installation problem: ", conditionMessage(e))
  )
  still_missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  hard_missing <- intersect(still_missing, required)
  soft_missing <- intersect(still_missing, optional)
  if (length(soft_missing) > 0) {
    message("Optional packages not installed (app still works, some features disabled): ",
            paste(soft_missing, collapse = ", "))
  }
  if (length(hard_missing) > 0) {
    stop("Required packages could not be installed: ", paste(hard_missing, collapse = ", "),
         "\nCheck your internet connection and try again.")
  }
}

# --- Find the app file (survives future renames like CellMorphR_v3_1...) -----
app_candidates <- sort(Sys.glob("CellMorphR*.R"), decreasing = TRUE)
if (length(app_candidates) == 0) stop("No CellMorphR*.R app file found in ", getwd())
app_file <- app_candidates[1]
message("Launching: ", app_file)

# --- Pick a free port (4321-4330) --------------------------------------------
port <- 4321
for (try_port in 4321:4330) {
  ok <- tryCatch({
    srv <- httpuv::startServer("127.0.0.1", try_port, list())
    httpuv::stopServer(srv)
    TRUE
  }, error = function(e) FALSE)
  if (ok) { port <- try_port; break }
}
message("Opening http://127.0.0.1:", port, " in your default browser...")
message("Keep this window open while using CellMorphR. Close it (or press Ctrl+C) to stop.")

tryCatch(
  shiny::runApp(shiny::shinyAppFile(app_file), port = port, host = "127.0.0.1",
                launch.browser = TRUE),
  error = function(e) {
    message("\nERROR while running CellMorphR: ", conditionMessage(e))
    message("If this mentions a missing package, re-run this launcher to retry installation.")
  }
)
