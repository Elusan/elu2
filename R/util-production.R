# ===========================
# Production curve
# ===========================

#' Production curve (ggplot2, theoretical + trajectory)
#'
#' Draws the theoretical production curve \eqn{P(B)} from a fitted SPiCT model
#' (optionally time-varying growth) and overlays the observed production–biomass
#' trajectory. Places dotted guides at \eqn{B/K} maximizing production and at
#' zero production.
#'
#' @param rep A fitted SPiCT result (\code{spictcls}).
#' @param n.plotyears If the number of labeled points is small, annotate years.
#' @param main Plot title.
#' @param stamp Stamp string placed in the plot corner.
#' @param CI Confidence level for parameter extraction (\code{0.95} by default).
#'
#' @return A printed \code{ggplot} object (invisibly).
#' @family ELU2-utilities
#' @export
#' @examples
#' \dontrun{
#' plotspict.production_gg(all_models$S1$S1P.SDM)
#' }
plotspict.production_gg <- function(rep, n.plotyears = 40,
                                    main = "Production curve",
                                    stamp = get.version(),
                                    CI = 0.95) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
  check.rep(rep)
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
    # widen Bplot to cover both theoretical & observed support
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
    # sort by time to connect chronologically (e.g., 192 → 2022)
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

  # ----- Only add base-style tick VALUES; do not change axis settings -----
  if (have_report) {
    x_breaks <- seq(0, 1, by = 0.2)
    x_labels <- formatC(x_breaks, format = "f", digits = 1)  # 0.0, 0.2, ..., 1.0
  } else {
    x_breaks <- ggplot2::waiver()
    x_labels <- ggplot2::waiver()
  }
  # -----------------------------------------------------------------------

  # Y label
  if (tvgflag) {
    ylab_txt <- "Production (normalised)"
  } else {
    cu <- as.character(inp$catchunit)
    ylab_txt <- if (nzchar(cu)) paste0("Production, ", cu) else "Production"
  }

  # x at maximum production
  mx <- (1 / npar[2])^(1 / (npar[2] - 1))

  library(ggplot2)
  p <- ggplot() +
    # Previous theoretical curves (gray)
    {
      if (nr > 1) {
        geom_line(
          data = subset(df_curves, curve_id < nr),
          aes(x = x, y = y),
          color = "gray"
        )
      }
    } +
    # Current theoretical curve (black, thin)
    geom_line(
      data = subset(df_curves, curve_id == nr),
      aes(x = x, y = y),
      color = "black"
    ) +
    # Production trajectory (time-ordered): continuous blue path + points
    {
      if (!is.null(df_path)) {
        list(
          geom_path(data = df_path, aes(x = x, y = y), color = "blue", linewidth = 0.3),
          geom_point(data = df_path, aes(x = x, y = y), color = "blue", size = 0.5)
        )
      }
    } +
    # Year labels if few points
    {
      if (!is.null(df_path) && length(inp$ic) < n.plotyears) {
        inds <- unique(c(1, nrow(df_path), seq(1, nrow(df_path), by = 2)))
        geom_text(
          data = transform(df_path[inds, , drop = FALSE],
                           lab = round(t, 2)),
          aes(x = x, y = y, label = lab),
          hjust = 0, nudge_x = 0.01, size = 1.8
        )
      }
    } +
    # --- Bold the starting year label (added layer; everything else unchanged) ---
    {
      if (!is.null(df_path) && length(inp$ic) < n.plotyears) {
        geom_text(
          data = transform(df_path[1, , drop = FALSE],
                           lab = round(t, 2)),
          aes(x = x, y = y, label = lab),
          fontface = "bold",
          hjust = 0, nudge_x = 0.01, size = 1.8
        )
      }
    } +
    # Guides
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = mx, linetype = "dotted") +
    # Labels & scales (unchanged; only breaks/labels added)
    labs(title = main, x = "B/K", y = ylab_txt) +
    scale_x_continuous(limits = xlim,
                       breaks = x_breaks, labels = x_labels,
                       expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(limits = ylim, expand = expansion(mult = c(0.04, 0.04))) +
    theme_classic(base_size = 11) +
    theme(
      plot.title   = element_text(hjust = 0.5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    coord_cartesian(clip = "off")

  # Stamp and convergence warning
  if (!is.null(stamp) && !is.na(stamp) && nzchar(stamp)) {
    p <- p + annotate("text", x = Inf, y = -Inf, label = stamp,
                      hjust = 1.1, vjust = -0.8, size = 1.5)
  }
  if (!is.null(rep$opt$convergence) && rep$opt$convergence != 0) {
    p <- p +
      annotate("point", x = xlim[1], y = ylim[2], shape = 24, size = 3,
               fill = "yellow", color = "black") +
      annotate("text",  x = xlim[1], y = ylim[2], label = "!", vjust = -0.3, size = 3)
  }

  # --- NEW: apply your compact theme and in-panel legend styling (no legend will show unless mapped) ---
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
      # restore this plot’s original black border & centered title
      panel.border         = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title           = ggplot2::element_text(hjust = 0.5)
    )

  print(p)
  invisible(p)

}  # (body unchanged)
