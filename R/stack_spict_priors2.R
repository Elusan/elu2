#' @title Stack SPiCT prior–posterior plots for multiple submodels
#'
#' @param reps        Named list of SPiCT fits, e.g. list(S8P = S8_Pella, S8F = S8_Fox, S8S = S8_Schaefer)
#' @param do.plot     Integer: max number of priors to plot (passed to plotspict.priors.ggplot98)
#' @param stamp       Character to stamp (passed to plotspict.priors.ggplot98)
#' @param CI          Confidence level for posterior CIs (passed to plotspict.priors.ggplot98)
#' @param save_path   If non-NULL, full path to write a PNG of the composite
#' @param width       Width in inches for the saved PNG (default 12)
#' @param height      Height in inches for the saved PNG (default 8)
#' @param dpi         Dots per inch for the saved PNG (default 300)
#' @return A patchwork object (invisible if save_path is non-NULL)
#' @importFrom patchwork wrap_plots
#' @importFrom ggplot2 ggsave annotation_custom
#' @importFrom grid textGrob unit gpar
#' @export
stack_spict_priors2 <- function(reps,
                                do.plot   = NULL,
                                stamp      = get.version(),
                                CI         = 0.95,
                                save_path  = NULL,
                                width      = 12,
                                height     = 8,
                                dpi        = 300) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Please install the {patchwork} package to use stack_spict_priors().")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install ggplot2 to use stack_spict_priors().")
  }
  # For grid::textGrob etc.
  requireNamespace("grid", quietly = TRUE)

  # 1. build one full priors plot per submodel, then inset the label
  plots <- lapply(names(reps), function(model_id) {
    p <- plotspict.priors.ggplot98(
      reps[[model_id]],
      do.plot = do.plot,
      stamp   = stamp,
      CI      = CI
    )
    # insert a bold model_id inside top-left of the panel area
    p + ggplot2::annotation_custom(
      grob = grid::textGrob(
        label = model_id,
        x     = grid::unit(0.02, "npc"),
        y     = grid::unit(0.98, "npc"),
        just  = c("left", "top"),
        gp    = grid::gpar(fontface = "bold", col = "grey30")
      ),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    )
  })

  # 2. stack them vertically
  combo <- patchwork::wrap_plots(plots, ncol = 1)

  # 3. save if requested
  if (!is.null(save_path)) {
    dir <- dirname(save_path)
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
    ggplot2::ggsave(
      filename  = save_path,
      plot      = combo,
      width     = width,
      height    = height,
      dpi       = dpi,
      limitsize = FALSE
    )
    return(invisible(combo))
  }

  # 4. otherwise return for interactive use
  combo
}
