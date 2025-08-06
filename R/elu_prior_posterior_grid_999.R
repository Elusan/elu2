#' Create and Save a Grid of Prior vs Posterior Plots for Multiple Models
#'
#' This function generates a grid of prior vs posterior plots for a list of SPiCT models,
#' organizing them with models as rows and priors as columns. Missing priors are replaced
#' with empty spacer plots to maintain consistent layout.
#'
#' @param models A named list of fitted SPiCT model objects (e.g., from `fit.elu2()`).
#'               The names should follow a pattern like "S1P", "S1F", etc., where the
#'               scenario prefix (e.g., "S1") is automatically extracted.
#' @param priors_fun A function used to generate the prior vs posterior plots for a single model.
#'                   Must return a list of ggplot objects (one per active prior).
#'                   Default is \code{\link{priors.elu6}}.
#' @param width Width of the final PNG figure in inches. Default: \code{15}.
#' @param height Height of the final PNG figure in inches. Default: \code{10}.
#' @param dpi Resolution (dots per inch) for the saved figure. Default: \code{300}.
#' @param outdir Directory to save the figure in. Created if it does not exist. Default: \code{"FIG"}.
#'               If \code{NULL}, no file is saved.
#' @param file_prefix Optional prefix for the filename. If \code{NULL}, the prefix is extracted
#'                    from the first model name by removing trailing model identifiers (e.g., "P").
#' @param return_patchwork Logical. If \code{TRUE}, returns the patchwork object containing the full grid.
#'                         If \code{FALSE}, the function runs for its side effects (i.e., saving the figure).
#'                         Default: \code{TRUE}.
#'
#' @details
#' The function determines which priors are active across all models and creates a grid of plots where
#' each row corresponds to a model and each column to a prior. If a prior is not used in a model,
#' an empty placeholder is inserted in its position to preserve alignment. The final grid can be saved
#' as a PNG and/or returned as a patchwork object for further use.
#'
#' @return A patchwork object if \code{return_patchwork = TRUE}, otherwise invisibly returns \code{NULL}.
#'
#' @importFrom patchwork wrap_plots plot_spacer
#' @importFrom ggplot2 ggsave theme_void
#'
#' @examples
#' \dontrun{
#' models <- list(S1P = fit1, S1F = fit2, S1S = fit3)
#' elu_prior_posterior_grid_999(models, priors_fun = priors.elu6)
#' }
#'
#' @export
elu_prior_posterior_grid_999 <- function(models,
                                         priors_fun = priors.elu6,
                                         width = 15, height = 10, dpi = 300,
                                         outdir = "FIG", file_prefix = NULL,
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
