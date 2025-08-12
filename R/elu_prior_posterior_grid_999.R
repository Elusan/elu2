#' Plot and Export a Grid of Prior vs. Posterior Distributions for Multiple Models
#'
#' Generates and saves a grid of prior vs. posterior plots for all active prior parameters across multiple SPiCT or ELU models. Each row represents a model, and each column a prior parameter; only priors actively used (estimated) in each model are shown. Empty panels are included where a prior is not used in a particular model.
#'
#' @param models A named list of model fit objects (e.g., SPiCT/ELU fits), each with \code{$inp$} and active priors.
#' @param priors_fun A function (default: \code{priors.elu6}) that generates a list of ggplot prior vs posterior plots for a single model.
#' @param width Numeric. Width (in inches) of the saved grid plot. Default is 13.
#' @param height Numeric. Height (in inches) of the saved grid plot. Default is 10.
#' @param dpi Numeric. Resolution (in dots per inch) for the saved image. Default is 300.
#' @param outdir Character. Output directory where the grid plot image will be saved. Default is \code{"FIG"}.
#' @param file_prefix Character. Optional file name prefix for the saved plot. If \code{NULL} (default), will use the scenario name extracted from model names.
#' @param return_patchwork Logical. If \code{TRUE} (default), returns the patchwork grid object. If \code{FALSE}, returns \code{NULL} invisibly.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Automatically determines the set of prior parameters that are active (used/estimated) in any of the models provided.
#'   \item Calls \code{priors_fun} for each model to generate the relevant prior vs posterior plots.
#'   \item Arranges all plots in a grid (rows = models, columns = priors), filling empty panels if a given prior is not estimated for a particular model.
#'   \item Saves the grid as a high-resolution PNG image to \code{outdir}, with a file name based on \code{file_prefix} or scenario.
#'   \item Optionally returns the \code{patchwork} object for further customization or immediate display.
#' }
#' This is especially useful for scenario-based sensitivity analysis and manuscript reporting.
#'
#' @return
#' A patchwork grid object if \code{return_patchwork = TRUE}; otherwise, \code{NULL} invisibly.
#'
#' @seealso \code{\link{priors.elu6}}
#'
#' @examples
#' \dontrun{
#' grid_plot <- elu_prior_posterior_grid_999(
#'   models = list(S1P = fit1, S1S = fit2, S1F = fit3),
#'   outdir = "FIG",
#'   file_prefix = "prior_posterior_grid_S1"
#' )
#' }
#'
#' @export
elu_prior_posterior_grid_999 <- function(models,
                                         priors_fun = priors.elu6,
                                         width = 13, height = 10, dpi = 300,
                                         outdir = "FIG/Priors_posterios_grids", file_prefix = NULL,
                                         return_patchwork = TRUE) {

  # -- Automatically set scenario prefix for filename --
  scenario_prefix <- ""
  model_names <- names(models)
  if (length(model_names) > 0) {
    scenario_prefix <- sub("([A-Za-z]+)$", "", model_names[1]) # remove trailing letters (P/S/F)
    scenario_prefix <- gsub("_+$", "", scenario_prefix) # clean trailing underscores
  }
  if (is.null(file_prefix) || file_prefix == "") {
    file_prefix <- paste0("prior_posterior_grid_", scenario_prefix)
  }
  # ----------------------------------------------------

  # Determine all priors that are active (used = 1) in any model
  all_priors <- unique(unlist(lapply(models, function(rep) {
    inp <- rep$inp
    useflags <- inp$priorsuseflags
    inds <- which(useflags == 1)
    names(inp$priors)[inds]
  })))
  n_models <- length(models)
  n_priors <- length(all_priors)
  model_ids <- names(models)

  # Helper: Get ggplot for one prior, or return a blank plot if not available
  get_prior_plot_or_empty <- function(rep, prior, model_id) {
    inp <- rep$inp
    useflags <- inp$priorsuseflags
    priors_avail <- names(inp$priors)[which(useflags == 1)]
    if (prior %in% priors_avail) {
      priors_row <- priors_fun(rep, model_id = model_id, do.plot = NULL, stamp = NULL)
      which_col <- which(priors_avail == prior)
      if (length(which_col) == 1 && length(priors_row) >= which_col) {
        return(priors_row[[which_col]])
      }
    }
    patchwork::plot_spacer() + ggplot2::theme_void()
  }

  # Generate all individual plots in a grid (rows = models, columns = priors)
  plots_grid <- lapply(seq_along(models), function(i) {
    model_id <- model_ids[i]
    rep <- models[[i]]
    lapply(all_priors, function(prior) {
      get_prior_plot_or_empty(rep, prior, model_id)
    })
  })

  plots_vec <- unlist(plots_grid, recursive = FALSE)
  final_patchwork <- patchwork::wrap_plots(plots_vec, nrow = n_models, ncol = n_priors)

  # Save to disk
  if (!is.null(outdir)) {
    if (!dir.exists(outdir)) dir.create(outdir)
    out_file <- file.path(outdir, paste0(file_prefix, ".png"))
    ggplot2::ggsave(out_file, final_patchwork, width = width, height = height, dpi = dpi)
  }

  if (return_patchwork) {
    return(final_patchwork)
  } else {
    invisible(NULL)
  }
}
