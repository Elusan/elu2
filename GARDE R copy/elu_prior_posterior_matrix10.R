#' Save Prior–Posterior Matrix per Scenario (v11)
#'
#' Like elu_prior_posterior_matrix10(), but will also
#' extract individual plots from a single patchwork return.
#'
#' @param models Named list of models. Names must end in "P","F","S" (e.g. "S8P","S8F","S8S").
#' @param priors_fun Function(model_obj, model_id) → either
#'        * a **named** list of ggplots, or
#'        * a **patchwork** object whose `$patches` each contain `$plot` and `$tag`.
#'        Default: `priors.elu6`.
#' @param outdir Directory for output PNGs. Default `"FIG"`.
#' @param panel_width Width per column (in). Default `4`.
#' @param panel_height Height per row (in). Default `4`.
#' @param dpi PNG resolution. Default `300`.
#' @return Invisibly `NULL`; side‐effect: one `<outdir>/priors_<scn>_matrix.png` per scenario.
#' @importFrom patchwork wrap_plots
#' @export
elu_prior_posterior_matrix10 <- function(models,
                                         priors_fun    = priors.elu6,
                                         outdir        = "FIG",
                                         panel_width   = 4,
                                         panel_height  = 4,
                                         dpi           = 300) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  # helper: turn patchwork → named list of ggplots
  extract_from_patchwork <- function(pw) {
    patches <- pw$patches
    if (is.null(patches) || length(patches)==0) {
      stop("Cannot extract from empty patchwork")
    }
    plots <- lapply(patches, "[[", "plot")
    tags  <- vapply(patches, function(x) {
      if (!is.null(x$tag) && nzchar(x$tag)) x$tag else NA_character_
    }, character(1), USE.NAMES=FALSE)
    # fallback names
    names(plots) <- ifelse(is.na(tags), paste0("prior", seq_along(plots)), tags)
    plots
  }

  # get unique scenario IDs
  scenarios <- unique(sub("(P|F|S)$", "", names(models)))

  for (scn in scenarios) {
    sub_keys <- paste0(scn, c("P","F","S"))
    priors_lists <- vector("list", length(sub_keys))
    names(priors_lists) <- sub_keys

    # 1) call priors_fun for each submodel
    for (key in sub_keys) {
      if (!key %in% names(models)) next
      obj <- models[[key]]
      fit_obj <- if (is.list(obj) && !is.null(obj$result)) obj$result else obj

      pr_raw <- try(priors_fun(fit_obj, model_id = key), silent=TRUE)
      if (inherits(pr_raw, "try-error") || is.null(pr_raw)) {
        warning("priors_fun() failed for ", key); next
      }

      # 2) normalize into a named list
      if (inherits(pr_raw, "patchwork")) {
        priors_lists[[key]] <- extract_from_patchwork(pr_raw)
      } else if (is.list(pr_raw) && all(vapply(pr_raw, inherits, logical(1), "ggplot"))) {
        # assume already named
        priors_lists[[key]] <- pr_raw
      } else {
        warning("priors_fun() returned unexpected type for ", key)
      }
    }

    # 3) gather all prior names
    all_priors <- sort(unique(unlist(lapply(priors_lists, names))))
    if (length(all_priors)==0) {
      warning("No priors found for scenario ", scn); next
    }

    # 4) build a 3×K list of plots
    plot_list <- vector("list", length(sub_keys) * length(all_priors))
    idx <- 1
    for (key in sub_keys) {
      pl <- priors_lists[[key]]
      for (pr in all_priors) {
        if (!is.null(pl) && pr %in% names(pl)) {
          plot_list[[idx]] <- pl[[pr]] +
            ggplot2::labs(title = paste0(key, " — ", pr))
        } else {
          plot_list[[idx]] <- ggplot2::ggplot() + ggplot2::theme_void()
        }
        idx <- idx + 1
      }
    }

    # 5) stitch into 3×K grid
    grid_plot <- patchwork::wrap_plots(plot_list, ncol = length(all_priors))

    # 6) save with correct dimensions
    w <- panel_width  * length(all_priors)
    h <- panel_height * length(sub_keys)
    outfile <- file.path(outdir, paste0("priors_", scn, "_matrix.png"))
    message("Writing matrix for ", scn, " → ", outfile)
    ggplot2::ggsave(
      filename   = outfile,
      plot       = grid_plot,
      width      = w,
      height     = h,
      dpi        = dpi,
      limitsize  = FALSE
    )
  }

  invisible(NULL)
}
