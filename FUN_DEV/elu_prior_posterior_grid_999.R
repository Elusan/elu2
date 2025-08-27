elu_prior_posterior_grid_999 <- function(models,
                                     priors_fun = priors.elu10,
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
    patchwork::plot_spacer() + theme_void()
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
