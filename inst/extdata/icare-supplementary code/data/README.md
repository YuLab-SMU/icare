# data/

Optional local inputs. The repository runs without anything here — this folder only
exists so a local copy can override the shipped default.

## `PMID37633276_DIA_plasmaproomic.xlsx`

The UTUC plasma-proteome workbook (n = 601) used by `01_leakage/` and
`02_head2head/`. Both scripts look here first:

```r
xlsx <- file.path(ROOT, "data", "PMID37633276_DIA_plasmaproomic.xlsx")
if (!file.exists(xlsx))   # fall back to the copy shipped inside the icare package
  xlsx <- system.file("extdata", "PMID37633276_DIA_plasmaproomic.xlsx", package = "icare")
```

So a normal install needs nothing in this folder. Drop the workbook here only if
you want to pin a specific version, in which case it takes precedence over the
package copy (and is git-ignored, see `.gitignore`).

Sheets used: sheet 3 for the 23 routine laboratory features (cases `KC*` dropped, as
in the manuscript), `Table S1B` / `Table S1E` for the plasma proteome.

## NHANES

Not here. Each scenario keeps its own cache in
`03_clinical_scenarios/<scenario>/data/*.rds` (≈ 3 MB per scenario), written by
`01_download.R`. Delete a scenario's `data/` folder to force a fresh download from
the CDC with `nhanesA`.
