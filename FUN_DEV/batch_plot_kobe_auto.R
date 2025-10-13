#' Batch-export Kobe 2×3 grids for an arbitrary set of models
#'
#' @description
#' A flexible batch runner around [plot_kobe_grid_2x3_last()] that **does not
#' assume any fixed set of scenario names**. It iterates over the supplied
#' `models` list (any length, any names), creates output directories as needed,
#' and renders/saves each model’s 2×3 Kobe-style grid.
#'
#' The arguments `out_dir` and `file_basename` are deliberately flexible so you
#' can control per-model destinations and filenames:
#'
#' - **`out_dir`** can be:
#'   - a single string (all files saved in that folder),
#'   - a **named** character vector keyed by model names,
#'   - a function `function(name)` returning the directory for that model.
#'
#' - **`file_basename`** can be:
#'   - `NULL` (defaults to `"<name>_custom_grid"`),
#'   - a **named** character vector keyed by model names,
#'   - a function `function(name)` returning the basename for that model,
#'   - a single un-named string (use the same basename for every model).
#'
#' @param models A list of fitted model objects. The list can be unnamed; if so,
#'   names are generated as `"scenario_1"`, `"scenario_2"`, … by position.
#' @param out_dir Destination directory (see flexibility notes above). Default
#'   `"FIG/KobePhasesNEW"`.
#' @param file_basename Basename for output files (see flexibility notes above).
#'   Default `NULL` meaning `"<name>_custom_grid"`.
#' @param width,height Plot size in inches passed to [plot_kobe_grid_2x3_last()].
#'   Defaults `width = 20`, `height = 9`.
#' @param dpi Resolution passed to [plot_kobe_grid_2x3_last()]. Default `400`.
#' @param man.legend Logical, forwarded to [plot_kobe_grid_2x3_last()] to show/hide
#'   the management legend inside the Kobe grid. Default `FALSE`.
#' @param save Logical, forwarded to [plot_kobe_grid_2x3_last()] indicating whether
#'   to write image files to disk. Default `TRUE`.
#' @param verbose Logical; if `TRUE`, prints progress messages. Default `TRUE`.
#' @param ... Additional arguments forwarded to [plot_kobe_grid_2x3_last()].
#'
#' @return
#' A `data.frame` with one row per model containing:
#' \itemize{
#'   \item `model` — model name used,
#'   \item `out_dir` — final output directory,
#'   \item `file_basename` — final basename used,
#'   \item `ok` — `TRUE` if the call succeeded, `FALSE` otherwise,
#'   \item `error` — error message when `ok = FALSE`, otherwise `NA`,
#'   \item `path` — return value from [plot_kobe_grid_2x3_last()] when available.
#' }
#'
#' @details
#' - If `models` contains duplicated names, they are made unique via
#'   [base::make.unique()] and a warning is issued.
#' - Any directories required by `out_dir` are created (recursively) if missing.
#'
#' @seealso [plot_kobe_grid_2x3_last()]
#'
#' @examples
#' \dontrun{
#' # 1) All to a single folder, default basenames "<name>_custom_grid"
#' batch_plot_kobe_auto(
#'   models   = all_models,
#'   out_dir  = "FIG/KobePhasesNEW",
#'   width = 20, height = 9, dpi = 400,
#'   man.legend = FALSE
#' )
#'
#' # 2) Per-model folders via a function:
#' batch_plot_kobe_auto(
#'   models   = all_models,
#'   out_dir  = function(nm) file.path("FIG", "KobePhasesNEW", nm)
#' )
#'
#' # 3) Mixed folders via a named vector:
#' batch_plot_kobe_auto(
#'   models  = all_models,
#'   out_dir = c(S1 = "FIG/KobePhasesNew", S2 = "FIG/KobePhasesNEW")
#' )
#'
#' # 4) Force a single basename for all:
#' batch_plot_kobe_auto(
#'   models        = all_models,
#'   file_basename = "mygrid"
#' )
#' }
#'
#' @export
batch_plot_kobe_auto <- function(models,
                                 out_dir = "FIG/KobePhasesNEW",
                                 file_basename = NULL,
                                 width = 20, height = 9, dpi = 400,
                                 man.legend = FALSE,
                                 save = TRUE,
                                 verbose = TRUE,
                                 ...) {
  stopifnot(is.list(models))

  # Derive model names (fallback to generic if missing)
  nm <- names(models)
  if (is.null(nm)) nm <- paste0("scenario_", seq_along(models))
  # ensure unique
  if (anyDuplicated(nm)) {
    nm <- make.unique(nm)
    warning("Duplicate model names detected; made unique with 'make.unique()'.")
  }

  # Helpers to resolve per-model out_dir and basename
  resolve_out_dir <- function(name) {
    if (is.character(out_dir) && length(out_dir) == 1) return(out_dir)
    if (is.character(out_dir) && !is.null(names(out_dir)) && name %in% names(out_dir)) return(out_dir[[name]])
    if (is.function(out_dir)) return(out_dir(name))
    if (is.character(out_dir) && length(out_dir) >= 1) return(out_dir[[1]])
    "FIG/KobePhasesNEW"
  }

  resolve_basename <- function(name) {
    if (is.null(file_basename)) return(paste0(name, "_custom_grid"))
    if (is.character(file_basename) && length(file_basename) == 1 && is.null(names(file_basename)))
      return(file_basename)
    if (is.character(file_basename) && !is.null(names(file_basename)) && name %in% names(file_basename))
      return(file_basename[[name]])
    if (is.function(file_basename)) return(file_basename(name))
    paste0(name, "_custom_grid")
  }

  # Ensure all required dirs exist
  dirs <- unique(vapply(nm, resolve_out_dir, character(1)))
  invisible(lapply(dirs, function(d) dir.create(d, recursive = TRUE, showWarnings = FALSE)))

  # Run
  res <- lapply(seq_along(models), function(i) {
    name   <- nm[[i]]
    model  <- models[[i]]
    odir   <- resolve_out_dir(name)
    base   <- resolve_basename(name)

    if (verbose) message(sprintf("[%-s] -> %s/%s.*", name, odir, base))

    ok <- TRUE; err <- NA_character_; out_path <- NA_character_

    tryCatch({
      out_path <- plot_kobe_grid_2x3_last(
        model,
        save          = save,
        out_dir       = odir,
        width         = width,
        height        = height,
        dpi           = dpi,
        file_basename = base,
        man.legend    = man.legend,
        ...
      )
    }, error = function(e) {
      ok  <<- FALSE
      err <<- conditionMessage(e)
      if (verbose) message(sprintf("  ✗ %s failed: %s", name, err))
    })

    data.frame(
      model        = name,
      out_dir      = odir,
      file_basename= base,
      ok           = ok,
      error        = if (ok) NA_character_ else err,
      path         = if (!is.null(out_path)) as.character(out_path) else NA_character_,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, res)
  rownames(out) <- NULL
  out
}
