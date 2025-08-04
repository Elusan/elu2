#' Grid of Prior and Posterior Distributions Across Multiple Models
#'
#' Generates a grid of plots comparing prior and posterior distributions for all active priors across multiple fitted ELU models. Each row in the grid corresponds to a model and each column to a prior parameter, enabling easy comparison of prior/posterior overlap and informativeness by model and parameter.
#'
#' @param models A named list of fitted model objects (e.g., from \code{fit.elu2()}), where each element represents a scenario or model variant.
#' @param priors_fun Function to generate prior/posterior plots for a single model. Should return a **named list** of ggplot objects. Default: \code{priors.elu6} or \code{priors.elu66}.
#' @param width Numeric; width (in inches) of the saved plot grid image.
#' @param height Numeric; height (in inches) of the saved plot grid image.
#' @param dpi Numeric; dots per inch for saved image (default: 300).
#' @param outdir Directory to save the image (default: \code{"FIG"}). If NULL, no file is written.
#' @param file_prefix Character prefix for the saved file name. Default is to generate one from the first model name.
#' @param return_patchwork Logical; if TRUE, returns the combined patchwork plot object. If FALSE, returns nothing (useful for batch saving).
#'
#' @return A patchwork plot object representing the grid of prior/posterior plots, or invisibly NULL if \code{return_patchwork = FALSE}.
#'
#' @importFrom patchwork wrap_plots plot_spacer
#' @importFrom ggplot2 theme_void ggsave
#' @export
elu_prior_posterior_grid2 <- function(models,
                                      priors_fun = priors.elu66,
                                      width = 15, height = 10, dpi = 300,
                                      outdir = "FIG", file_prefix = NULL,
                                      return_patchwork = TRUE) {
  # Dependencies
  requireNamespace("patchwork")
  requireNamespace("ggplot2")

  # Auto file prefix from first model
  model_names <- names(models)
  scenario_prefix <- if (length(model_names) > 0) {
    tmp <- sub("([A-Za-z]+)$", "", model_names[1])
    gsub("_+$", "", tmp)
  } else ""
  if (is.null(file_prefix) || file_prefix == "") {
    file_prefix <- paste0("prior_posterior_grid_", scenario_prefix)
  }

  # 1. Determine unique prior names across all models
  all_priors <- unique(unlist(lapply(models, function(r) {
    names(r$inp$priors)[which(r$inp$priorsuseflags == 1)]
  })))

  n_models  <- length(models)
  n_priors  <- length(all_priors)
  model_ids <- names(models)

  # 2. Helper: get each cell by name
  get_cell <- function(rep, model_id, prior) {
    plots_list <- priors_fun(rep, model_id = model_id,
                             do.plot = NULL, stamp = NULL)
    if (is.list(plots_list) && prior %in% names(plots_list)) {
      return(plots_list[[prior]])
    } else {
      return(patchwork::plot_spacer() + ggplot2::theme_void())
    }
  }

  # 3. Build grid: row = models, col = priors
  plots_grid <- unlist(
    lapply(seq_along(models), function(i) {
      rep <- models[[i]]
      id  <- model_ids[i]
      lapply(all_priors, function(p) get_cell(rep, id, p))
    }), recursive = FALSE
  )

  # 4. Combine using patchwork
  final_pw <- patchwork::wrap_plots(
    plots_grid, nrow = n_models, ncol = n_priors
  )

  # 5. Save if requested
  if (!is.null(outdir)) {
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
    out_file <- file.path(outdir, paste0(file_prefix, ".png"))
    ggplot2::ggsave(out_file, final_pw,
                    width = width, height = height, dpi = dpi)
  }

  # 6. Return or invisible
  if (return_patchwork) {
    return(final_pw)
  } else {
    invisible(NULL)
  }
}
