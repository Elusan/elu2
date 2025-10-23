#' Export color-coded percent-change tables to LaTeX (SDM & GLM)
#'
#' Takes the output of [`compute_param_pct_change()`]—specifically `pct$wide$SDM`
#' and/or `pct$wide$GLM`—and writes **two LaTeX tables** (one per `type`) with
#' compact formatting, alternating row shading within families, and colored cells
#' for positive/negative/near-zero changes. Files are saved as
#' `file.path(output_dir, paste0(file_stub, "_SDM.tex"))` and
#' `file.path(output_dir, paste0(file_stub, "_GLM.tex"))`.
#'
#' @param pct A list returned by [`compute_param_pct_change()`], containing a
#'   `wide` element with data frames `wide$SDM` and/or `wide$GLM`. Each wide
#'   table must have columns: `Family`, `Param`, and one column per scenario
#'   (excluding the base), whose values are formatted strings like `"+3.2%"`.
#' @param output_dir Directory where the `.tex` files will be written.
#'   Default `"FIG"`. Created if it does not exist.
#' @param file_stub File name stub (prefix) used for both outputs. Default
#'   `"pct_change_vs_S4"`.
#' @param caption_base Base caption text used for both tables; the function
#'   appends the type (SDM/GLM). Default
#'   `"Percent change of estimates relative to base scenario S4"`.
#' @param zero_tol Numeric tolerance in percent units; cells with
#'   `abs(% change) <= zero_tol` are colored gray to indicate “no meaningful
#'   change.” Default `0.05` (i.e., 0.05\\%).
#' @param sdm_wrap_footnotesize Logical; if `TRUE`, wraps the **SDM** table in
#'   a `\\begin{footnotesize} ... \\end{footnotesize}` block for extra compactness.
#'   Default `TRUE`. The GLM table is not wrapped by default.
#'
#' @return (Invisibly) `TRUE` on success. Side effect: writes up to two `.tex`
#'   files to `output_dir`.
#'
#' @details
#' **Expected input structure:** `pct$wide$SDM` and/or `pct$wide$GLM` as produced
#' by [`compute_param_pct_change()`], with:
#' - `Family` (model family label), `Param` (parameter label),
#' - one column per scenario (e.g., `S1`, `S2`, `S5`, ...), whose cells are
#'   strings like `"+1.2%"`, `"-0.4%"`, or `""` for missing values.
#'
#' **Color rules:**
#' - `|value| <= zero_tol` → `\{\\textcolor{gray}{...}\}` (neutral).
#' - `value > 0` → `\{\\textcolor{teal}{...}\}` (positive).
#' - `value < 0` → `\{\\textcolor{red!70!black}{...}\}` (negative).
#'
#' **Row shading:** Within each `Family` block, rows alternate `\\cellcolor{gray!10}`.
#'
#' **LaTeX requirements:** Ensure your preamble includes
#' ```latex
#' \usepackage[table]{xcolor}  % for \cellcolor and \textcolor
#' \usepackage{booktabs}       % optional, for \addlinespace if you use it elsewhere
#' ```
#'
#' @seealso [compute_param_pct_change()]
#'
#' @examples
#' \dontrun{
#' pct <- compute_param_pct_change(param_df, base_scenario = "S4")
#' export_pct_change_tables_awesome22(
#'   pct,
#'   output_dir   = "FIG",
#'   file_stub    = "pct_change_vs_S4",
#'   caption_base = "Percent change of estimates relative to base scenario S4",
#'   zero_tol     = 0.05
#' )
#' # Files written:
#' #   FIG/pct_change_vs_S4_SDM.tex
#' #   FIG/pct_change_vs_S4_GLM.tex
#' }
#'
#' @export
export_pct_change_tables_awesome22 <- function(
    pct,
    output_dir   = "FIG",                          # where .tex files go
    file_stub    = "pct_change_vs_S4",             # filename prefix
    caption_base = "Percent change of estimates relative to base scenario S4",
    zero_tol     = 0.05,                           # |%| ≤ tol → gray (≈ no change)
    sdm_wrap_footnotesize = TRUE                   # wrap SDM table in \begin{footnotesize}
) {
  # Validate input structure: must be the list returned by compute_param_pct_change()
  if (!is.list(pct) || is.null(pct$wide)) stop("pct must be the output of compute_param_pct_change().")
  # Ensure target directory exists (quietly)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # --- helpers -------------------------------------------------------------

  # Escape LaTeX special chars in text cells (%, &, #, _).
  # Numbers/s already formatted like "+1.2%" should be left as-is except for escaping.
  esc <- function(x) {
    if (is.factor(x)) x <- as.character(x)
    if (!is.character(x)) return(x)
    x <- gsub("([%&#_])", "\\\\\\1", x, perl = TRUE)  # → \% \& \# \_
    x
  }

  # Color a percentage string based on magnitude:
  #   |value| ≤ tol  → gray (neutral)
  #   value  > 0     → teal (positive change)
  #   value  < 0     → red!70!black (negative change)
  # Non-numeric / empty strings are returned escaped without coloring.
  colorize_pct <- function(s, tol = zero_tol) {
    if (!nzchar(s)) return("")                                # empty cell stays empty
    num <- suppressWarnings(as.numeric(gsub("%", "", s)))     # parse "+1.2%" → 1.2
    if (!is.finite(num)) return(esc(s))                       # not a number → plain escaped
    s_tex <- esc(s)
    if (abs(num) <= tol) {
      sprintf("{\\textcolor{gray}{%s}}", s_tex)
    } else if (num > 0) {
      sprintf("{\\textcolor{teal}{%s}}", s_tex)
    } else {
      sprintf("{\\textcolor{red!70!black}{%s}}", s_tex)
    }
  }

  # Build a single LaTeX table for a given type ("SDM" or "GLM").
  # Expects `tab` columns: Family, Param, then one column per scenario (S1,S2,... excluding base).
  build_table_tex <- function(tab, type_label, wrap_footnotesize = FALSE) {
    if (is.null(tab) || !nrow(tab)) return("")               # nothing to write

    fam_col <- "Family"
    par_col <- "Param"
    sc_cols <- setdiff(names(tab), c(fam_col, par_col))      # dynamic scenario columns

    # Preserve grouping per Family; we also alternate light row shading within each family
    families <- unique(tab[[fam_col]])
    # For each family, start shading on the first row for odd families (1-based index)
    fam_starts_shaded <- setNames(seq_along(families) %% 2 == 1, families)

    lines <- character(0)
    add <- function(...) { lines <<- c(lines, paste0(...)) } # small appender

    if (wrap_footnotesize) add("\\begin{footnotesize}")

    add("\\begin{table}[!h]")
    add("\\centering")
    add(sprintf("\\caption{%s (Type: %s)}", esc(caption_base), type_label))
    add("\\vspace{-0.23cm}")                                 # compact vertical spacing
    add("\\centering")
    add("\\fontsize{7}{9}\\selectfont")                      # small, tight table
    add(sprintf("\\begin{tabular}[t]{ll%s}",                 # 2 left cols + right-aligned scenarios
                paste(rep("r", length(sc_cols)), collapse = "")))
    add("\\hline")

    # Column headers: Model, Param., then scenario IDs
    hdr <- c("\\textbf{Model}", "\\textbf{Param.}", sprintf("\\textbf{%s}", esc(sc_cols)))
    add(paste(hdr, collapse = " & "), "\\\\")
    add("\\hline")
    add("\\addlinespace[0.3em]")

    # ---- table body (grouped by Family) -----------------------------------
    for (fi in seq_along(families)) {
      fam <- families[fi]
      # Family group label spanning the full width
      add(sprintf("\\multicolumn{%d}{l}{\\textbf{%s}}\\\\", 2 + length(sc_cols), esc(fam)))

      # Subset rows for this family in their current order
      sub <- tab[tab[[fam_col]] == fam, , drop = FALSE]
      # Whether the first row in this family starts shaded
      shaded <- isTRUE(fam_starts_shaded[[fam]])

      for (ri in seq_len(nrow(sub))) {
        row <- sub[ri, , drop = FALSE]

        # Build left-hand cells (Model/Family and Param) with optional shading
        shade_prefix <- if (shaded) "\\cellcolor{gray!10}" else ""
        lhs_fam <- paste0("\\hspace{1em}",
                          if (shaded) paste0(shade_prefix, "{", esc(row[[fam_col]]), "}")
                          else esc(row[[fam_col]]))
        lhs_par <- if (shaded) paste0(shade_prefix, "{", esc(row[[par_col]]), "}")
        else esc(row[[par_col]])

        # Scenario cells: colorized, with same row shading if toggled
        cells <- character(length(sc_cols))
        for (j in seq_along(sc_cols)) {
          val <- as.character(row[[sc_cols[j]]])             # e.g., "+1.2%" or ""
          colstr <- colorize_pct(val, zero_tol)
          cells[j] <- if (shaded) sprintf("%s{%s}", shade_prefix, colstr) else colstr
        }

        # Emit the row
        add(paste(c(lhs_fam, lhs_par, cells), collapse = " & "), "\\\\")
        shaded <- !shaded                                     # alternate shading per row
      }

      # Light spacer between families (except after the last one)
      if (fi < length(families)) add("\\addlinespace[0.3em]")
    }

    add("\\hline")
    add("\\end{tabular}")
    add("\\end{table}")
    if (wrap_footnotesize) add("\\end{footnotesize}")

    paste(lines, collapse = "\n")
  }

  # --- writer --------------------------------------------------------------
  # Write .tex to disk and message the path (useBytes avoids locale quirks)
  write_tex <- function(tex, path) { writeLines(tex, path, useBytes = TRUE); message("Saved: ", path) }

  # SDM table (optional footnotesize wrapper)
  if (!is.null(pct$wide$SDM)) {
    tex_sdm <- build_table_tex(pct$wide$SDM, "SDM", wrap_footnotesize = sdm_wrap_footnotesize)
    write_tex(tex_sdm, file.path(output_dir, sprintf("%s_SDM.tex", file_stub)))
  }
  # GLM table
  if (!is.null(pct$wide$GLM)) {
    tex_glm <- build_table_tex(pct$wide$GLM, "GLM", wrap_footnotesize = FALSE)
    write_tex(tex_glm, file.path(output_dir, sprintf("%s_GLM.tex", file_stub)))
  }

  invisible(TRUE)  # quiet success for piping
}
