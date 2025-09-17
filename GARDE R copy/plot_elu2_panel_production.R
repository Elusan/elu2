#' Production curve (ggplot2) — theoretical P(B) + time-ordered trajectory
#'
#' @description
#' Draws the theoretical production curve \eqn{P(B)} from a fitted SPiCT model,
#' optionally overlaying the time-ordered production/biomass trajectory
#' (if `inp$reportall` is `TRUE`). When time-varying growth or log-m-covariance
#' is enabled, the y-axis is normalised against the maximum of the theoretical
#' curves (as in the base plotting logic).
#'
#' @param rep A fitted SPiCT report object (`spictcls`).
#' @param n.plotyears Integer; if the number of trajectory points is less than
#'   `n.plotyears`, year labels are drawn (first year is **bold**).
#'   Default `40`.
#' @param main Plot title. Default `"Production curve"`.
#' @param stamp Small stamp (e.g. version string) placed bottom-right.
#'   Default `spict::get.version()`.
#' @param CI Confidence level (0,1) for parameter extraction. Default `0.95`.
#'
#' @return Invisibly returns the `ggplot` object after drawing it.
#'
#' @details
#' Minimal equivalents of external helpers are used via \pkg{spict}:
#' `get.par()` and `get.version()`. A compact theme helper
#' is inlined (file-local) and applied at the end.
#'
#' @examples
#' \dontrun{
#' rep <- fit.spict(inp)
#' plot_elu2_panel_plotspict.production_gg(rep)
#' }
#'
#' @export
plot_elu2_panel_production <- function(rep, n.plotyears = 40,
                                                    main = "Production curve",
                                                    stamp = get.version(),
                                                    CI = 0.95) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  # ---- file-local minimal validator (avoids non-exported spict::check.rep) ----
  .elu2_check_rep <- function(x) {
    if (!inherits(x, "spictcls") || is.null(x$inp) || is.null(x$opt))
      stop("The argument 'rep' must be a fitted spict object (class 'spictcls').")
    invisible(TRUE)
  }
  .elu2_check_rep(rep)

  if ("sderr" %in% names(rep)) return(invisible(NULL))

  inp     <- rep$inp
  tvgflag <- isTRUE(inp$timevaryinggrowth) || isTRUE(inp$logmcovflag)

  Kest  <- get.par("logK",  rep, exp = TRUE, CI = CI)
  mest  <- get.par("logm",  rep, exp = TRUE, CI = CI)
  nr    <- nrow(mest)
  gamma <- get.par("gamma", rep,           CI = CI)
  npar  <- get.par("logn",  rep, exp = TRUE, CI = CI)
  Pest  <- get.par("P",     rep,           CI = CI)

  # indices that tie P to B(t)
  binds <- inp$ic[1:nrow(Pest)]

  # y scaling under time-varying growth
  if (tvgflag) {
    yscal <- get.par("logMSYvec", rep, exp = TRUE, CI = CI)[binds, 2]
  } else {
    yscal <- rep(1, length(binds))
  }

  # ---- file-local helpers (plain comments; not roxygen) ----
  # Compact minimal theme applied at the end
  .spict_theme_minimal_compact22 <- function(base_size = 10, base_family = "") {
    ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12),
        axis.title = ggplot2::element_text(face = "bold", size = 12),
        axis.text  = ggplot2::element_text(size = 10, face = "bold"),
        legend.position = "none",
        panel.grid = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(fill = NA, colour = "grey35", linewidth = 1),
        axis.ticks = ggplot2::element_line(linewidth = 0.5, color = "grey35"),
        axis.ticks.length = grid::unit(3, "pt"),
        strip.background = ggplot2::element_rect(fill = "grey35", color = "grey35", linewidth = 0.5),
        strip.text = ggplot2::element_text(face = "bold", size = ggplot2::rel(1)),
        text = ggplot2::element_text(face = "bold", size = 10),
        plot.margin = ggplot2::margin(3, 3, 3, 3)
      )
  }

  # Grid for theoretical curves
  nBplot <- 200
  Bplot  <- seq(0.5 * 1e-8, Kest[2], length.out = nBplot)

  pfun <- function(gamma, m, K, n, B) gamma * m / K * B * (1 - (B / K)^(n - 1))

  Pst_list <- lapply(seq_len(nr), function(i) {
    pfun(gamma[2], mest[i, 2], Kest[2], npar[2], Bplot)
  })
  Pstscal <- if (tvgflag) max(unlist(Pst_list)) else 1

  have_report <- isTRUE(inp$reportall)
  if (have_report) {
    Best  <- get.par("logB", rep, exp = TRUE, CI = CI)
    Bplot <- seq(0.5 * min(c(1e-8, Best[, 2])),
                 1.0 * max(c(Kest[2], Best[, 2])),
                 length.out = nBplot)
    Pst_list <- lapply(seq_len(nr), function(i) {
      pfun(gamma[2], mest[i, 2], Kest[2], npar[2], Bplot)
    })
    Bvec <- Best[binds, 2]
  }

  # Theoretical curves data (x = B/K; y normalized if needed)
  df_curves <- do.call(
    rbind,
    lapply(seq_len(nr), function(i) {
      data.frame(
        curve_id = i,
        x = Bplot / Kest[2],
        y = Pst_list[[i]] / Pstscal
      )
    })
  )

  # Time-ordered path (guarantees continuity across years)
  if (have_report) {
    t_years <- inp$time[binds]
    df_path <- data.frame(
      t = t_years,
      x = Bvec / Kest[2],
      y = Pest[, 2] / yscal
    )
    df_path <- df_path[order(df_path$t), , drop = FALSE]
  } else {
    df_path <- NULL
  }

  # Axis limits (mirror base logic)
  if (have_report) {
    xlim <- range(c(df_path$x, 0, 1), na.rm = TRUE)
    ylim <- c(
      min(0, df_path$y, na.rm = TRUE),
      max(df_path$y, df_curves$y, na.rm = TRUE)
    )
  } else {
    xlim <- range(df_curves$x, na.rm = TRUE)
    ylim <- c(0, max(df_curves$y, na.rm = TRUE))
  }

  # Tick labels (only if we have report)
  if (have_report) {
    x_breaks <- seq(0, 1, by = 0.2)
    x_labels <- formatC(x_breaks, format = "f", digits = 1)
  } else {
    x_breaks <- ggplot2::waiver()
    x_labels <- ggplot2::waiver()
  }

  # Y label
  if (tvgflag) {
    ylab_txt <- "Production (normalised)"
  } else {
    cu <- as.character(inp$catchunit)
    ylab_txt <- if (nzchar(cu)) paste0("Production, ", cu) else "Production"
  }

  # x at maximum production
  mx <- (1 / npar[2])^(1 / (npar[2] - 1))

  p <- ggplot2::ggplot() +
    { if (nr > 1) ggplot2::geom_line(
      data = subset(df_curves, curve_id < nr),
      ggplot2::aes(x = x, y = y), color = "gray"
    ) } +
    ggplot2::geom_line(
      data = subset(df_curves, curve_id == nr),
      ggplot2::aes(x = x, y = y), color = "black"
    ) +
    {
      if (!is.null(df_path)) {
        list(
          ggplot2::geom_path(data = df_path, ggplot2::aes(x = x, y = y),
                             color = "blue", linewidth = 0.3),
          ggplot2::geom_point(data = df_path, ggplot2::aes(x = x, y = y),
                              color = "blue", size = 0.5)
        )
      }
    } +
    {
      if (!is.null(df_path) && length(inp$ic) < n.plotyears) {
        inds <- unique(c(1, nrow(df_path), seq(1, nrow(df_path), by = 2)))
        ggplot2::geom_text(
          data = transform(df_path[inds, , drop = FALSE], lab = round(t, 2)),
          ggplot2::aes(x = x, y = y, label = lab),
          hjust = 0, nudge_x = 0.01, size = 1.8
        )
      }
    } +
    {
      if (!is.null(df_path) && length(inp$ic) < n.plotyears) {
        ggplot2::geom_text(
          data = transform(df_path[1, , drop = FALSE], lab = round(t, 2)),
          ggplot2::aes(x = x, y = y, label = lab),
          fontface = "bold", hjust = 0, nudge_x = 0.01, size = 1.8
        )
      }
    } +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
    ggplot2::geom_vline(xintercept = mx, linetype = "dotted") +
    ggplot2::labs(title = main, x = "B/K", y = ylab_txt) +
    ggplot2::scale_x_continuous(limits = xlim,
                                breaks = x_breaks, labels = x_labels,
                                expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    ggplot2::scale_y_continuous(limits = ylim,
                                expand = ggplot2::expansion(mult = c(0.04, 0.04))
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(hjust = 0.5),
      panel.border = ggplot2::element_rect(color = "grey35", fill = NA, linewidth = 1)
    ) +
    ggplot2::coord_cartesian(clip = "off")

  if (!is.null(stamp) && !is.na(stamp) && nzchar(stamp)) {
    p <- p + ggplot2::annotate("text", x = Inf, y = -Inf, label = stamp,
                               hjust = 1.1, vjust = -0.8, size = 1.5)
  }
  if (!is.null(rep$opt$convergence) && rep$opt$convergence != 0) {
    p <- p +
      ggplot2::annotate("point", x = xlim[1], y = ylim[2], shape = 24, size = 3,
                        fill = "yellow", color = "black") +
      ggplot2::annotate("text",  x = xlim[1], y = ylim[2], label = "!",
                        vjust = -0.3, size = 3)
  }

  p <- p +
    .spict_theme_minimal_compact22() +
    ggplot2::theme(
      legend.position      = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.title         = ggplot2::element_blank(),
      legend.background    = ggplot2::element_rect(
        fill   = grDevices::adjustcolor("white", alpha.f = 0.90),
        colour = "grey35", linewidth = 0.5
      ),
      legend.key           = ggplot2::element_rect(fill = NA, colour = NA),
      legend.key.size      = grid::unit(10, "pt"),
      legend.text          = ggplot2::element_text(size = 9, face = "bold"),
      panel.border         = ggplot2::element_rect(color = "grey35", fill = NA, linewidth = 1),
      plot.title           = ggplot2::element_text(hjust = 0.5)
    )

  print(p)
  invisible(p)
}
