#' Export a management table to a LaTeX **tabular** (no float)
#'
#' @description
#' Writes a LaTeX `tabular` (not a floating table) summarizing SPiCT management
#' scenarios to a `.tex` file. The table is produced with **pure base R**
#' (robust escaping/formatting, no external formatting dependencies) and assumes
#' that `rep$man` has been created (i.e., you have run `manage(rep)`).
#'
#' The function internally calls \code{sumspict.manage()} (from **SPiCT**) to
#' build the estimates, then formats and escapes them for LaTeX. Column types
#' are handled robustly: numeric columns are rounded with per-column digits and
#' missing/non-finite values are printed as `--`. Scenario names are UTF-8
#' encoded and LaTeX-escaped in a safe order (backslash first, then specials).
#'
#' @details
#' - The output is a **bare tabular** (no `table` environment, no caption, no label),
#'   intended for inclusion via `\\input{...}` in your LaTeX document.
#' - The column spec uses `\\RaggedRight` from **ragged2e** and `p{<width>}` from
#'   **array** for the scenario column. Ensure these packages are loaded in your
#'   LaTeX preamble, e.g.:
#'
#'   \\preformatted{
#'   \\usepackage{array}
#'   \\usepackage{ragged2e}
#'   }
#'
#' - The function writes bytes with UTF-8 encoding; XeLaTeX/LuaLaTeX are recommended.
#' - Required columns from \code{sumspict.manage()} are: \code{"C"}, \code{"B/Bmsy"},
#'   and \code{"F/Fmsy"}; an informative error is thrown if any are missing.
#'
#' @param rep A fitted SPiCT report object (class \code{spictcls}) **with**
#'   management scenarios present in \code{rep$man} (i.e., after \code{manage(rep)}).
#' @param file Character scalar. Output path to the `.tex` file to write. Parent
#'   directories are created as needed.
#' @param digits Named vector (or list) of per-column rounding digits for numeric
#'   columns. Defaults to \code{c("C" = 1, "B/Bmsy" = 2, "F/Fmsy" = 2)}.
#'   Columns not listed fall back to 2 decimal places.
#' @param select_rows Optional integer vector of row indices to subset the
#'   scenarios from the internal estimates table (after the function builds it).
#'   Indices are validated and must be within \code{1..nrow(df)}.
#' @param scenario_colwidth Character scalar giving the LaTeX width for the
#'   first column (scenario names), e.g. \code{"3cm"} (default).
#'
#' @return
#' (Invisibly) returns the data frame that was formatted and written, with
#' columns: \code{"Management scenario"}, \code{"C"}, \code{"B/Bmsy"},
#' \code{"F/Fmsy"}. On error, throws with an informative message.
#'
#' @section LaTeX requirements:
#' The generated tabular assumes the following packages in your preamble:
#' \itemize{
#'   \item \strong{array} – for \code{p\{\}} column types and \code{\\arraybackslash}
#'   \item \strong{ragged2e} – for \code{\\RaggedRight}
#' }
#'
#' @seealso \code{\link[spict]{manage}}, \code{\link[spict]{sumspict.manage}}
#'
#' @examples
#' \dontrun{
#' # After fitting and running management:
#' # rep <- fit.spict(inp)
#' # rep <- manage(rep)
#' elu2_export_manage_tabular_only_tex3(
#'   rep,
#'   file = "FIG/manage_table_S1.tex",
#'   digits = c("C" = 0, "B/Bmsy" = 2, "F/Fmsy" = 2),
#'   select_rows = c(1, 3, 4),
#'   scenario_colwidth = "3.2cm"
#' )
#' # Later in LaTeX:
#' # \\input{FIG/manage_table_S1.tex}
#' }
#'
#' @export
#' @encoding UTF-8
#' @importFrom utils capture.output
elu2_export_manage_tabular_only_tex3 <- function(
    rep,
    file,
    digits = c("C" = 1, "B/Bmsy" = 2, "F/Fmsy" = 2),
    select_rows = NULL,
    scenario_colwidth = "3cm"
) {
  ## ---- Basic checks ----
  if (missing(file) || !is.character(file) || !nzchar(file[1L])) {
    stop("Provide an output `.tex` path in `file`.")
  }
  if (is.null(rep) || is.null(rep$man)) {
    stop("Run manage(rep) first: `rep$man` not found in `rep`.")
  }
  # Create parent dir if needed (silently)
  fdir <- dirname(file)
  if (is.na(fdir) || !nzchar(fdir)) fdir <- "."
  if (!file.exists(fdir)) {
    dir.create(fdir, recursive = TRUE, showWarnings = FALSE)
    if (!file.exists(fdir)) stop("Cannot create output directory: ", fdir)
  }

  ## ---- Quietly build the management estimates table (requires SPiCT) ----
  .elu2_get_manage_est_table_quiet <- function(rep, include_EBinf) {
    # capture without importing utils at top-level
    junk <- utils::capture.output({
      res <- sumspict.manage(
        rep           = rep,
        include.EBinf = isTRUE(include_EBinf),
        include.unc   = FALSE,
        include.abs   = FALSE,
        timeline      = TRUE,
        verbose       = FALSE
      )
    })
    if (is.null(res) || is.null(res$est)) {
      stop("sumspict.manage() returned no $est; check that `rep$man` exists and is valid.")
    }
    df <- data.frame("Management scenario" = rownames(res$est), res$est,
                     row.names = NULL, check.names = FALSE)
    df
  }

  # You can keep include_EBinf = TRUE safely; we strip unused cols later
  df <- .elu2_get_manage_est_table_quiet(rep, include_EBinf = TRUE)

  ## ---- Validate required columns and optionally subset rows ----
  needed <- c("Management scenario", "C", "B/Bmsy", "F/Fmsy")
  miss <- setdiff(needed, names(df))
  if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))
  df <- df[, needed, drop = FALSE]

  if (!is.null(select_rows)) {
    # tolerate double/numeric; coerce safely to integer indices
    select_rows <- as.integer(select_rows)
    if (any(!is.finite(select_rows))) stop("`select_rows` must be finite integers.")
    if (length(select_rows) == 0L) stop("`select_rows` is empty after coercion.")
    if (any(select_rows < 1L | select_rows > nrow(df))) {
      stop("`select_rows` out of range 1..", nrow(df), ".")
    }
    df <- df[select_rows, , drop = FALSE]
  }

  ## ---- Encode text as UTF-8 (XeLaTeX-friendly) ----
  df[["Management scenario"]] <- enc2utf8(df[["Management scenario"]])

  ## ---- Robust LaTeX escaping (no double-escape) ----
  # Escape order matters: backslash first, then other specials.
  escape_latex <- function(x) {
    x <- as.character(x)
    # Keep NA as-is for now; we map later
    nas <- !nzchar(ifelse(is.na(x), "", x))
    # Backslash -> \textbackslash{}
    x <- gsub("\\\\", "\\\\textbackslash{}", x, perl = TRUE)
    # Then specials not already escaped:
    x <- gsub("(?<!\\\\)%",  "\\\\%",                 x, perl = TRUE)
    x <- gsub("(?<!\\\\)_",  "\\\\_",                 x, perl = TRUE)
    x <- gsub("(?<!\\\\)&",  "\\\\&",                 x, perl = TRUE)
    x <- gsub("(?<!\\\\)#",  "\\\\#",                 x, perl = TRUE)
    x <- gsub("(?<!\\\\)\\$", "\\\\$",                x, perl = TRUE)
    x <- gsub("(?<!\\\\)\\{", "\\\\{",                x, perl = TRUE)
    x <- gsub("(?<!\\\\)\\}", "\\\\}",                x, perl = TRUE)
    # ~ and ^ (use text-mode ascii symbols)
    x <- gsub("(?<!\\\\)~",  "\\\\textasciitilde{}",  x, perl = TRUE)
    x <- gsub("(?<!\\\\)\\^","\\\\textasciicircum{}", x, perl = TRUE)
    # Keep original NA empty state
    x[nas & is.na(x)] <- NA_character_
    x
  }

  df[["Management scenario"]] <- escape_latex(df[["Management scenario"]])

  ## ---- Numeric rounding & stable formatting ----
  # Map NA/NaN/Inf to NA; then print with format() (non-scientific, trimmed).
  # Per-column digits from named vector/list `digits`.
  format_numeric_col <- function(v, nm) {
    v <- suppressWarnings(as.numeric(v))
    v[!is.finite(v)] <- NA_real_
    d <- if (!is.null(digits) && !is.null(digits[[nm]])) {
      as.integer(digits[[nm]])
    } else 2L
    # round first for reproducibility, then format
    vv <- round(v, d)
    out <- format(vv, trim = TRUE, scientific = FALSE, nsmall = d)
    out[is.na(v)] <- "--"
    out
  }

  # Apply to the three numeric columns if numeric; if not numeric, just escape.
  for (nm in c("C", "B/Bmsy", "F/Fmsy")) {
    if (is.numeric(df[[nm]]) || is.integer(df[[nm]]) || is.double(df[[nm]])) {
      df[[nm]] <- format_numeric_col(df[[nm]], nm)
    } else {
      # For character factors etc., encode and escape; keep as-is if they already contain symbols.
      df[[nm]] <- escape_latex(enc2utf8(as.character(df[[nm]])))
      # Replace empty strings with NA marker for table consistency
      df[[nm]][!nzchar(df[[nm]])] <- "--"
    }
  }

  ## ---- Build LaTeX rows (type-stable) ----
  make_row <- function(i) {
    paste0(
      df[i, "Management scenario"], " & ",
      df[i, "C"],                  " & ",
      df[i, "B/Bmsy"],             " & ",
      df[i, "F/Fmsy"],             "\\\\"
    )
  }
  rows_tex <- paste(vapply(seq_len(nrow(df)), make_row, character(1L)), collapse = "\n")

  ## ---- Header / Footer (note: requires ragged2e + array in preamble) ----
  header <- paste0(
    "\\setlength{\\tabcolsep}{2.5pt}%\n",
    "\\renewcommand{\\arraystretch}{1.05}%\n",
    "\\begin{tabular}{@{} >{\\RaggedRight\\arraybackslash}p{", scenario_colwidth, "} lll @{} }\n",
    "\\hline\n",
    "\\textbf{Management scenario} & \\textbf{C} & \\textbf{B/Bmsy} & \\textbf{F/Fmsy}\\\\\n",
    "\\hline\n"
  )
  footer <- "\n\\hline\n\\end{tabular}%\n"  # keep the % to kill trailing whitespace

  tex <- paste0(header, rows_tex, footer)

  ## ---- Write as bytes, preserving UTF-8 ----
  con <- file(file, open = "wb", encoding = "UTF-8")
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  # writeChar preserves bytes; set eos = NULL for raw write
  writeChar(tex, con, eos = NULL, useBytes = TRUE)

  invisible(df)
}
