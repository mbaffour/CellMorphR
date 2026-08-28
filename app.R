# =============================================================================
# CellMorphR - RStudio entry point
#
# Open this file in RStudio and click "Run App" (top right of the editor),
# or from any R console in this folder:  shiny::runApp()
#
# On Windows you can instead just double-click run_cellmorphr.bat - no RStudio
# needed. Both routes run exactly the same app (CellMorphR_v3.R).
# =============================================================================

local({
  # RStudio's "Run App" and shiny::runApp() both set the working directory to
  # this folder. Fail with a clear message if the app file is not beside us.
  if (!file.exists("CellMorphR_v3.R")) {
    stop("CellMorphR_v3.R was not found in ", getwd(),
         ".
Open app.R from inside the CellMorphR folder (or set that folder as the working directory).")
  }

  # Install anything missing into your personal library (first run only)
  required <- c("shiny", "bslib", "ggplot2", "dplyr", "tidyr", "readxl", "DT",
                "scales", "writexl", "cluster")
  optional <- c("lme4", "lmerTest", "emmeans", "ggridges", "ggbeeswarm", "ggforce",
                "patchwork", "colourpicker", "hexbin", "svglite", "officer", "rvg")
  pkgs <- c(required, optional)
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    ul <- Sys.getenv("R_LIBS_USER")
    if (nzchar(ul)) {
      ul <- strsplit(ul, .Platform$path.sep, fixed = TRUE)[[1]][1]
      if (!dir.exists(ul)) dir.create(ul, recursive = TRUE, showWarnings = FALSE)
      .libPaths(c(ul, .libPaths()))
    }
    message("First run: installing ", length(missing), " package(s): ",
            paste(missing, collapse = ", "), "\nThis can take several minutes...")
    try(install.packages(missing, lib = if (nzchar(ul)) ul else NULL))
    still <- intersect(required, pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)])
    if (length(still))
      stop("Required packages could not be installed: ", paste(still, collapse = ", "))
  }
})

# The app itself lives in CellMorphR_v3.R; this returns it for RStudio to run.
shiny::shinyAppFile("CellMorphR_v3.R")
