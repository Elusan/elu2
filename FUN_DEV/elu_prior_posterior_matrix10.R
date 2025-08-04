#' Save Prior–Posterior Matrix per Scenario
#'
#' For each scenario (e.g. "S1"), calls `priors_fun()` on its three submodels
#' ("S1P","S1F","S1S"), each of which must return a **named** list of ggplots
#' (one plot per prior parameter).  It then arranges those into a matrix:
#'   - Rows: P, F, S
#'   - Columns: every unique prior name across the three
#' missing cells get a blank spacer.
#'
#' @param models Named list of fitted‐model objects.  Names must follow
#'   `<Scenario><P|F|S>` (e.g. `"S1P"`, `"S1F"`, `"S1S"`, …).
#' @param priors_fun Function that, given a model and its `model_id`, returns
#'   a **named** list of ggplot objects, one per prior parameter.
#'   Default: `priors.elu6`.
#' @param outdir Directory to save PNGs.  Created if needed.  Default: `"FIG"`.
#' @param panel_width  Width (in inches) per column.  Default: 4.
#' @param panel_height Height (in inches) per row.  Default: 4.
#' @param dpi           PNG resolution. Default: 300.
#' @return Invisibly `NULL`.  Side‐effect: saves one file per scenario named
#'   `priors_<Scenario>_matrix.png`.
#' @importFrom patchwork wrap_plots
#' @export
elu_prior_posterior_matrix10 <- function(models,
                                       priors_fun    = priors.elu6,
                                       outdir        = "FIG",
                                       panel_width   = 4,
                                       panel_height  = 4,
                                       dpi           = 300) {
  # ensure output dir
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  # figure out scenario tags
  scenarios <- unique(sub("(P|F|S)$", "", names(models)))

  for (scn in scenarios) {
    # expected submodel keys
    sub_keys <- paste0(scn, c("P", "F", "S"))
    # call priors_fun for each submodel (may be NULL or error)
    priors_lists <- lapply(sub_keys, function(key) {
      if (!key %in% names(models)) return(NULL)
      obj <- models[[key]]
      # guard: skip if no result
      if (is.null(obj) || inherits(obj, "error") || is.null(obj$result)) return(NULL)
      priors_fun(obj$result, model_id = key)
    })
    names(priors_lists) <- sub_keys

    # collect all prior names
    all_priors <- sort(unique(unlist(lapply(priors_lists, names))))
    if (length(all_priors) == 0) {
      warning("No priors found for scenario ", scn)
      next
    }

    # build a list of plots in row-major order
    plot_list <- vector("list", length(sub_keys) * length(all_priors))
    idx <- 1
    for (key in sub_keys) {
      pri_list <- priors_lists[[key]]
      for (pr in all_priors) {
        if (!is.null(pri_list) && pr %in% names(pri_list)) {
          plot_list[[idx]] <- pri_list[[pr]] +
            ggplot2::labs(title = paste0(key, " — ", pr))
        } else {
          plot_list[[idx]] <- ggplot2::ggplot() + ggplot2::theme_void()
        }
        idx <- idx + 1
      }
    }

    # arrange into a 3×K grid
    grid_plot <- patchwork::wrap_plots(plot_list, ncol = length(all_priors))

    # compute overall dimensions
    w <- panel_width  * length(all_priors)
    h <- panel_height * length(sub_keys)

    # save
    outfile <- file.path(outdir, paste0("priors_", scn, "_matrix.png"))
    message("Writing: ", outfile)
    ggplot2::ggsave(outfile,
                    plot      = grid_plot,
                    width     = w,
                    height    = h,
                    dpi       = dpi,
                    limitsize = FALSE)
  }

  invisible(NULL)
}
