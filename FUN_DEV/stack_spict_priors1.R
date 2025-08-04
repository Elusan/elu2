#' @title Stack SPiCT prior–posterior plots for multiple submodels
#'
#' @param reps        Named list of SPiCT fits (e.g. list(Pella=S8_Pella, Fox=S8_Fox, Schaefer=S8_Schaefer))
#' @param do.plot     Integer: max number of priors to plot (passed to plotspict.priors.ggplot98)
#' @param stamp       Character to stamp on each plot (passed to plotspict.priors.ggplot98)
#' @param CI          Confidence level for posterior CIs (passed to plotspict.priors.ggplot98)
#' @param save_path   If non‐NULL, path to write out a PNG of the composite (via ggplot2::ggsave)
#' @param width       Width in inches for the saved PNG (default 12)
#' @param height      Height in inches for the saved PNG (default 8)
#' @param dpi         Dots per inch for the saved PNG (default 300)
#' @return A patchwork object (invisible if \code{save_path} is non‐NULL)
#' @importFrom patchwork wrap_plots
#' @importFrom ggplot2 ggsave
#' @export
stack_spict_priors1 <- function(reps,
                                do.plot   = NULL,
                                stamp      = get.version(),
                                CI         = 0.95,
                                save_path  = NULL,
                                width      = 12,
                                height     = 10,
                                dpi        = 300) {
  # require patchwork for arranging
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Please install the patchwork package to use stack_spict_priors().")
  }
  # build one full plot per submodel
  plots <- lapply(names(reps), function(model_id) {
    # call your existing function (it will print internally, but we capture the returned gg object)
    p <- plotspict.priors.ggplot98(
      reps[[model_id]],
      do.plot = do.plot,
      stamp   = stamp,
      CI      = CI
    )
    # add a row‐header so we know which is which
    p + ggplot2::labs(subtitle = model_id)
  })
  # stack them
  combo <- patchwork::wrap_plots(plots, ncol = 1)

  # save if requested
  if (!is.null(save_path)) {
    dir <- dirname(save_path)
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
    ggplot2::ggsave(save_path, combo,
                    width  = width,
                    height = height,
                    dpi    = dpi)
    invisible(combo)
  } else {
    return(combo)
  }
}
