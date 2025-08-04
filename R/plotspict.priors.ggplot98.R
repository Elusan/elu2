#' @name plotspict.priors.ggplot
#' @title Plot priors and posterior distributions using ggplot2 + theme_sleek
#' @param rep     A result from fit.spict, or a list with `$result`
#' @param do.plot Integer: max number of priors to plot (NULL = all)
#' @param stamp   Character to stamp on the plot (default get.version())
#' @param CI      Confidence level for posterior CIs (default 0.95)
#' @return Invisibly, the ggplot object
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot geom_ribbon geom_line facet_wrap scale_x_log10
#' @importFrom ggplot2 scale_color_manual scale_fill_manual labs theme element_rect element_text
#' @importFrom ggplot2 element_line margin unit
#' @export
plotspict.priors.ggplot98 <- function(rep,
                                      do.plot = NULL,
                                      stamp   = get.version(),
                                      CI      = 0.95) {
  if (!is.null(rep$result)) rep <- rep$result
  inp   <- rep$inp
  flags <- inp$priorsuseflags
  inds  <- which(flags == 1)
  if (length(inds) == 0) {
    message("No priors are used")
    return(invisible(NULL))
  }

  df_list <- list()
  for (j in inds) {
    orig_nm  <- names(inp$priors)[j]
    label_nm <- sub("gamma","", sub("^log","", orig_nm))
    pv       <- inp$priors[[j]]
    parmat   <- get.par(orig_nm, rep, exp = FALSE, CI = CI)

    repriors <- c("logB","logF","logBBmsy","logFFmsy")
    if (orig_nm %in% repriors) {
      sel      <- pv[5]
      parmat   <- parmat[sel, , drop = FALSE]
      label_nm <- paste0(sub("^log","", orig_nm), fd(pv[4]))
    }

    for (rr in seq_len(nrow(parmat))) {
      this_label <- if (nrow(parmat)>1) paste0(label_nm, rr) else label_nm
      prvec <- if (is.list(pv)) pv[[rr]] else pv

      mu <- ifelse(is.na(parmat[rr,4]), prvec[1], parmat[rr,2])
      sd <- ifelse(is.na(parmat[rr,4]), prvec[2], parmat[rr,4])
      isG <- grepl("gamma", orig_nm)

      if (isG && is.na(parmat[rr,4])) {
        xmin <- 1e-12
        xmax <- stats::qgamma(0.99, prvec[1], prvec[2])
      } else {
        xmin <- mu - 3*sd
        xmax <- mu + 3*sd
      }

      xpr <- seq(xmin, xmax, length.out=200)
      xpo <- if (isG && !is.na(parmat[rr,4])) {
        seq(mu - 3*sd, mu + 3*sd, length.out=200)
      } else xpr

      prior_vals <- if (isG)
        stats::dgamma(xpr, prvec[1], prvec[2]) else
          stats::dnorm(xpr, prvec[1], prvec[2])

      ## — avoid sd == 0 or negative —
      post_vals <- if (is.na(parmat[rr,4]) || parmat[rr,4] <= 0) {
        rep(NA_real_, length(xpo))
      } else {
        stats::dnorm(xpo, parmat[rr,2], parmat[rr,4])
      }

      ## — build the data frame as before —
      df_list[[length(df_list)+1]] <-
        data.frame(
          param     = this_label,
          x         = exp(xpr),
          prior     = prior_vals,
          posterior = c(post_vals, rep(NA_real_, length(xpr) - length(post_vals))),
          stringsAsFactors = FALSE
        )

    }
  }

  df <- do.call(rbind, df_list)
  if (!is.null(do.plot)) {
    keep <- unique(df$param)[1:do.plot]
    df   <- df[df$param %in% keep, ]
  }

  library(reshape2)
  dfl <- melt(df,
              id.vars       = c("param","x"),
              measure.vars  = c("prior","posterior"),
              variable.name = "type",
              value.name    = "density")

  library(ggplot2)
  p <- ggplot(dfl, aes(x=x, y=density, color=type, fill=type)) +
    geom_ribbon(aes(ymin=0, ymax=density), alpha=0.14, show.legend=FALSE) +
    geom_line(linewidth=0.8, show.legend=TRUE) +
    geom_hline(yintercept = 0, colour = "white", linewidth = 0.5, show.legend = FALSE) +
    facet_wrap(~ param, scales="free") +
    scale_x_log10() +
    scale_color_manual(
      values = c("prior"="black","posterior"="red"),
      labels = c("Prior","Post.")
    ) +
    scale_fill_manual(values = c("prior"="black","posterior"="red"), guide=FALSE) +
    labs(x=NULL, y="Density", color=NULL) +
    theme_sleek() +
    theme(
      legend.position      = c(0.98, 0.98),
      legend.justification = c("right", "top")
    )

  #–– Remove the old annotate block entirely and replace with:

  if (!is.null(stamp)) {
    p <- p +
      labs(
        tag = stamp
      ) +
      theme(
        plot.tag.position = c(0, 0),
        plot.tag           = element_text(
          hjust = 0, vjust = 0,
          size  = 3,
          color = "grey30"
        )
      )
  }

  print(p)
  invisible(p)
}
