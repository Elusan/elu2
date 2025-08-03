#' Grid of Prior and Posterior Distributions Across Multiple Models
#'
#' Generates a grid of plots comparing prior and posterior distributions for all active priors across multiple fitted ELU models. Each row in the grid corresponds to a model and each column to a prior parameter, enabling easy comparison of prior/posterior overlap and informativeness by model and parameter.
#'
#' @param models A named list of fitted model objects (e.g., from \code{fit.elu2()}), where each element represents a scenario or model variant.
#' @param priors_fun Function to generate prior/posterior plots for a single model. Should take arguments (\code{rep}, \code{model_id}, etc.) and return a list of ggplot objects, one for each prior (default: \code{priors.elu6}).
#' @param width Numeric; width (in inches) of the saved plot grid image.
#' @param height Numeric; height (in inches) of the saved plot grid image.
#' @param dpi Numeric; dots per inch for saved image (default: 300).
#' @param outdir Directory to save the image (default: \code{"FIG"}). If NULL, no file is written.
#' @param file_prefix Character prefix for the saved file name. Default is to generate one from the first model name.
#' @param return_patchwork Logical; if TRUE, returns the combined patchwork plot object. If FALSE, returns nothing (useful for batch saving).
#'
#' @details
#' The function collects all unique active priors across the provided models. For each model (row) and prior (column), it generates a prior/posterior plot using \code{priors_fun}. If a prior is not available for a model, the corresponding cell is left blank. The resulting grid is saved as a PNG if \code{outdir} is not NULL and is also returned as a patchwork object if \code{return_patchwork = TRUE}.
#'
#' @return A patchwork plot object representing the grid of prior/posterior plots, or invisibly NULL if \code{return_patchwork = FALSE}.
#'
#' @importFrom patchwork wrap_plots plot_spacer
#' @importFrom ggplot2 theme_void ggsave
#' @examples
#' \dontrun{
#' grid_plot <- elu_prior_posterior_grid(list(Fox = fit1, Schaefer = fit2, Pella = fit3))
#' }
#' @export
elu_prior_posterior_grid <- function(models,
                                     priors_fun = priors.elu6,
                                     width = 15, height = 10, dpi = 300,
                                     outdir = "FIG", file_prefix = NULL,  # <- NULL by default
                                     return_patchwork = TRUE) {
  require(patchwork)
  require(ggplot2)

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

  all_priors <- unique(unlist(lapply(models, function(rep) {
    inp <- rep$inp
    useflags <- inp$priorsuseflags
    inds <- which(useflags == 1)
    names(inp$priors)[inds]
  })))
  n_models <- length(models)
  n_priors <- length(all_priors)
  model_ids <- names(models)

  # TO KEEP HERE ---------
  # get_prior_plot_or_empty <- function(rep, prior, model_id) {
  #   inp <- rep$inp
  #   useflags <- inp$priorsuseflags
  #   priors_avail <- names(inp$priors)[which(useflags == 1)]
  #   if (prior %in% priors_avail) {
  #     priors_row <- priors_fun(rep, model_id = model_id, do.plot = NULL, stamp = NULL)
  #     which_col <- which(priors_avail == prior)
  #     if (length(which_col) == 1 && length(priors_row) >= which_col) {
  #       return(priors_row[[which_col]])
  #     }
  #   }
  #   patchwork::plot_spacer() + theme_void()
  # }

  get_prior_plot_or_empty <- function(rep, prior, model_id) {
    inp <- rep$inp
    useflags <- inp$priorsuseflags
    priors_avail <- names(inp$priors)[which(useflags == 1)]

    # Call priors.elu6() to get all plots with names
    plot_list <- priors_fun(rep, model_id = model_id, do.plot = NULL, stamp = NULL, return_list = TRUE)

    # Safely return the matching plot if available
    if (!is.null(plot_list) && prior %in% names(plot_list)) {
      return(plot_list[[prior]])
    }

    # Return blank panel if not available
    patchwork::plot_spacer() + ggplot2::theme_void()
  }

  plots_grid <- lapply(seq_along(models), function(i) {
    model_id <- model_ids[i]
    rep <- models[[i]]
    lapply(all_priors, function(prior) {
      get_prior_plot_or_empty(rep, prior, model_id)
    })
  })

  plots_vec <- unlist(plots_grid, recursive = FALSE)

  final_patchwork <- wrap_plots(plots_vec, nrow = n_models, ncol = n_priors)


  #final_patchwork <- wrap_plots(plots_vec, nrow = n_models, ncol = n_priors) +
    #plot_layout(guides = "collect")

  if (!is.null(outdir)) {
    if (!dir.exists(outdir)) dir.create(outdir)
    out_file <- file.path(outdir, paste0(file_prefix, ".png"))
    ggsave(out_file, final_patchwork, width = width, height = height, dpi = dpi)
  }

  if (return_patchwork) {
    return(final_patchwork)
  } else {
    invisible(NULL)
  }
}
