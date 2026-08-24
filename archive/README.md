# Archive

Superseded releases, kept only for reference and reproducibility. These are
**not maintained** and do not receive fixes.

| File | Version | Notes |
|------|---------|-------|
| `CellMorphR_v1.0.R` | 1.0.0 | The original single-file app (6 analysis tabs). Shipped in this repo as `R-code` until v3.0.0 replaced it. |

To run an archived version:

```bash
Rscript -e "shiny::runApp(shiny::shinyAppFile('archive/CellMorphR_v1.0.R'), launch.browser = TRUE)"
```

Note: v2.8 was a local-only working copy that was edited directly into v3.0.0,
so it was never committed and cannot be recovered.

For current work use `CellMorphR_v3.R` in the repository root.
