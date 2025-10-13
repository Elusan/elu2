# ------------------------------------------------------------------------------
# Helpers + Dynamic Scenario Runner (keeps your existing functions untouched)
# ------------------------------------------------------------------------------

#' @title (Internal) Safe adapter for `priors_as_named_list()`
#' @description Returns a named list of ggplots or an empty list on error.
#' @keywords internal
priors_as_named_list_safe <- function(model, ...) {
  out <- try(priors_as_named_list(model, ...), silent = TRUE)
  if (inherits(out, "try-error") || is.null(out)) return(list())
  if (!is.list(out)) return(list())
  out
}

#' @title (Internal) Infer/repair `inp$priorsuseflags` from plotted prior names
#' @description
#' Given the names produced by `priors_as_named_list()` for a model (e.g.,
#' "logsdi1", "logq2", "logsdc"), ensure the model has usable
#' `inp$priorsuseflags` by turning on the corresponding base names
#' ("logsdi", "logq", "logsdc") if missing or zeroed.
#' @param model A fitted ELU/SPiCT model object with `$inp$priors`.
#' @param plot_names Character vector of prior plot names (from safe adapter).
#' @return The same model object with `$inp$priorsuseflags` ensured/updated.
#' @keywords internal
infer_and_set_priorsuseflags <- function(model, plot_names) {
  if (!length(plot_names)) return(model)

  base_of <- function(x) sub("\\d+$", "", x)  # "logsdi1" -> "logsdi"
  active_bases <- unique(base_of(plot_names))

  if (is.null(model$inp) || is.null(model$inp$priors)) return(model)

  priors_names <- names(model$inp$priors)
  if (is.null(priors_names) || !length(priors_names)) return(model)

  need_new_flags <- is.null(model$inp$priorsuseflags) ||
    length(model$inp$priorsuseflags) != length(priors_names) ||
    is.null(names(model$inp$priorsuseflags)) ||
    !all(names(model$inp$priorsuseflags) %in% priors_names)

  if (need_new_flags) {
    flags <- integer(length(priors_names))
    names(flags) <- priors_names
    on <- intersect(priors_names, active_bases)
    if (length(on)) flags[on] <- 1L
    model$inp$priorsuseflags <- flags
    return(model)
  }

  # Existing flags: ensure active bases are turned on
  flags <- model$inp$priorsuseflags
  names(flags) <- priors_names
  on <- intersect(priors_names, active_bases)
  if (length(on)) flags[on] <- 1L
  model$inp$priorsuseflags <- flags
  model
}

#' @title Plot Prior–Posterior Grids for Any Number of Scenarios
#' @description
#' Iterates over \code{all_models} (a named list of scenarios), and for each
#' scenario calls your existing \code{elu_prior_posterior_grid_999()} to build
#' the grid. It automatically:
#' \itemize{
#'   \item Works with \strong{any} number of scenarios and models (no fixed counts).
#'   \item Skips scenarios where no model has active priors.
#'   \item Repairs missing/zeroed \code{inp$priorsuseflags} using the names
#'         returned by \code{priors_as_named_list()}.
#'   \item Optionally autosizes the output figure based on rows/cols; you can
#'         keep fixed size by setting \code{autosize = FALSE}.
#' }
#' Your functions \code{theme_minimal_compact_this()} and
#' \code{elu_prior_posterior_grid_999()} are \emph{not modified}.
#'
#' @param all_models Named list of scenarios; each element is a named list
#'   of fitted models for that scenario.
#' @param outdir Output directory for PNG files (passed to your grid function).
#' @param file_prefix_prefix Prefix used before the scenario name to compose the
#'   filename stem.
#' @param dpi Output DPI passed through to your grid function.
#' @param autosize Logical; if \code{TRUE}, width/height are inferred from the
#'   number of columns/rows. If \code{FALSE}, use \code{fixed_width}/\code{fixed_height}.
#' @param base_width_per_col Inches per prior-column when \code{autosize=TRUE}.
#' @param base_height_per_row Inches per model-row when \code{autosize=TRUE}.
#' @param min_width Minimum width (inches) when \code{autosize=TRUE}.
#' @param min_height Minimum height (inches) when \code{autosize=TRUE}.
#' @param fixed_width Fixed width (inches) when \code{autosize=FALSE}.
#' @param fixed_height Fixed height (inches) when \code{autosize=FALSE}.
#'
#' @return \code{invisible(NULL)} (files are saved by your existing function).
#'
#' @examples
#' \dontrun{
#' plot_all_scenarios_dynamic(
#'   all_models = all_models,
#'   outdir = "FIG",
#'   file_prefix_prefix = "elu_prior_grid_",
#'   dpi = 300
#' )
#' }
#' @importFrom patchwork wrap_plots plot_spacer
#' @importFrom ggplot2 theme_void
#' @export
plot_all_scenarios_dynamic <- function(all_models,
                                       outdir = "FIG",
                                       file_prefix_prefix = "elu_prior_grid_",
                                       dpi = 300,
                                       autosize = TRUE,
                                       base_width_per_col = 2.8,
                                       base_height_per_row = 2.6,
                                       min_width = 10,
                                       min_height = 6,
                                       fixed_width = 12,
                                       fixed_height = 11) {
  if (!length(all_models)) {
    message("No scenarios in all_models.")
    return(invisible(NULL))
  }

  for (s in names(all_models)) {
    message("Running elu_prior_posterior_grid_999() for ", s)

    models_s <- all_models[[s]]
    if (is.null(models_s) || !length(models_s)) {
      message("  -> No models in scenario ", s, ". Skipping.")
      next
    }

    # 1) Probe active priors per model WITHOUT plotting
    np_lists <- lapply(seq_along(models_s), function(i) {
      priors_as_named_list_safe(
        models_s[[i]],
        model_id = names(models_s)[i],
        do.plot  = NULL,
        stamp    = NULL
      )
    })
    names(np_lists) <- names(models_s)

    # 2) Drop models that truly have no active priors
    keep <- vapply(np_lists, function(x) is.list(x) && length(x) > 0, logical(1))
    if (!any(keep)) {
      message("No active priors across the supplied models in scenario ", s, ". Nothing to plot.")
      next
    }
    models_s <- models_s[keep]
    np_lists  <- np_lists[keep]

    # 3) Repair priorsuseflags based on what priors are actually plotted
    models_s <- mapply(function(m, nplist) {
      infer_and_set_priorsuseflags(m, names(nplist))
    }, models_s, np_lists, SIMPLIFY = FALSE)
    names(models_s) <- names(np_lists)

    # 4) Determine size
    if (isTRUE(autosize)) {
      plot_name_union <- unique(unlist(lapply(np_lists, names)))
      n_cols <- max(1L, length(plot_name_union))
      n_rows <- max(1L, length(models_s))
      width  <- max(min_width,  base_width_per_col  * n_cols)
      height <- max(min_height, base_height_per_row * n_rows)
    } else {
      width  <- fixed_width
      height <- fixed_height
    }

    # 5) Call YOUR existing function (UNCHANGED)
    elu_prior_posterior_grid_999(
      models           = models_s,
      outdir           = outdir,
      file_prefix      = paste0(file_prefix_prefix, s),
      width            = width,
      height           = height,
      dpi              = dpi,
      return_patchwork = FALSE
    )
  }

  invisible(NULL)
}
