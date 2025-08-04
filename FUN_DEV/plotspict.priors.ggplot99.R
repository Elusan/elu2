#' @title Plot priors and posterior distributions using ggplot2 + theme_sleek
#' @param rep     A result from fit.spict(), or a list with `$result`
#' @param do.plot Maximum number of priors to display (NULL = all)
#' @param stamp   Optional version label
#' @param CI      Posterior confidence level
#' @return ggplot object (invisibly)
#' @export
plotspict.priors.ggplot99 <- function(rep,
                                      do.plot = NULL,
                                      stamp   = if (exists("get.version")) get.version() else "v0.0",
                                      CI      = 0.95) {
  if (!is.list(rep)) stop("`rep` must be a list.")
  if (!is.null(rep$result)) rep <- rep$result
  if (is.null(rep$inp)) stop("Missing $inp in SPiCT result.")

  inp   <- rep$inp
  flags <- inp$priorsuseflags
  inds  <- which(flags == 1)
  if (length(inds) == 0) {
    message("No priors used.")
    return(invisible(NULL))
  }

  df_list <- list()
  for (j in inds) {
    orig_nm  <- names(inp$priors)[j]
    label_nm <- sub("gamma", "", sub("^log", "", orig_nm))
    pv       <- inp$priors[[j]]
    parmat   <- get.par(orig_nm, rep, exp = FALSE, CI = CI)

    if (orig_nm %in% c("logB","logF","logBBmsy","logFFmsy")) {
      sel      <- pv[5]
      parmat   <- parmat[sel, , drop = FALSE]
      label_nm <- paste0(sub("^log", "", orig_nm), fd(pv[4]))
    }

    for (rr in seq_len(nrow(parmat))) {
      this_label <- if (nrow(parmat) > 1) paste0(label_nm, rr) else label_nm
      prvec <- if (is.list(pv)) pv[[rr]] else pv
      mu <- ifelse(is.na(parmat[rr, 4]), prvec[1], parmat[rr, 2])
      sd <- ifelse(is.na(parmat[rr, 4]), prvec[2], parmat[rr, 4])
      isG <- grepl("gamma", orig_nm)

      if (isG && is.na(parmat[rr, 4])) {
        xmin <- 1e-12
        xmax <- stats::qgamma(0.99, prvec[1], prvec[2])
      } else {
        xmin <- mu - 3 * sd
        xmax <- mu + 3 * sd
      }

      xpr <- seq(xmin, xmax, length.out = 200)
      xpo <- if (isG && !is.na(parmat[rr, 4])) seq(mu - 3*sd, mu + 3*sd, length.out = 200) else xpr

      prior_vals <- if (isG) stats::dgamma(xpr, prvec[1], prvec[2]) else stats::dnorm(xpr, prvec[1], prvec[2])
      post_vals <- if (is.na(parmat[rr, 4]) || parmat[rr, 4] <= 0) rep(NA_real_, length(xpr)) else stats::dnorm(xpo, parmat[rr, 2], parmat[rr, 4])

      df_list[[length(df_list) + 1]] <- data.frame(
        param     = this_label,
        x         = exp(xpr),
        prior     = prior_vals,
        posterior = if (length(post_vals) == length(xpr)) post_vals else rep(NA_real_, length(xpr)),
        stringsAsFactors = FALSE
      )
    }
  }

  df <- do.call(rbind, df_list)
  if (!is.null(do.plot)) {
    keep <- head(unique(df$param), do.plot)
    df <- df[df$param %in% keep, ]
  }

  dfl <- reshape2::melt(df,
                        id.vars = c("param", "x"),
                        measure.vars = c("prior", "posterior"),
                        variable.name = "type",
                        value.name = "density"
  )
  dfl$param <- factor(dfl$param, levels = unique(dfl$param))

  p <- ggplot2::ggplot(dfl, ggplot2::aes(x = x, y = density, color = type, fill = type)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = density), alpha = 0.14, show.legend = FALSE) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = 0, colour = "white", linewidth = 0.5) +
    ggplot2::facet_wrap(~param, scales = "free") +
    ggplot2::scale_x_log10() +
    ggplot2::scale_color_manual(values = c("prior" = "black", "posterior" = "red"),
                                labels = c("Prior", "Post.")) +
    ggplot2::scale_fill_manual(values = c("prior" = "black", "posterior" = "red"), guide = FALSE) +
    ggplot2::labs(x = NULL, y = "Density", color = NULL) +
    theme_sleek() +
    ggplot2::theme(
      legend.position = c(0.98, 0.98),
      legend.justification = c("right", "top")
    )

  if (!is.null(stamp)) {
    p <- p +
      ggplot2::labs(tag = stamp) +
      ggplot2::theme(
        plot.tag.position = c(0, 0),
        plot.tag = ggplot2::element_text(hjust = 0, vjust = 0, size = 3, color = "grey30")
      )
  }

  print(p)
  invisible(p)
}
