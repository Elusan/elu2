#' Time-to-Bmsy plot (ggplot2; with visual pads around y = 1 and x = 0)
#'
#' @description
#' Simulates biomass trajectories under several multiples of \eqn{F_{MSY}}
#' and plots the time it would take for biomass to reach (or return to)
#' B\eqn{_{MSY}} from the current biomass level. Uses SPiCT helpers
#' `check.rep()`, `get.par()`, and `calc.gamma()` for inputs/parameters.
#'
#' @param rep A fitted SPiCT report object (class `spictcls`).
#' @param main Plot title. Default `"Time to Bmsy"`.
#' @param stamp Small stamp (e.g., version string) placed bottom-right.
#'   Default `get.version()` from SPiCT.
#' @param CI Confidence level (0,1) for parameter extraction. Default `0.95`.
#'
#' @return Invisibly returns the `ggplot` object after drawing it.
#'
#' @details
#' The function draws the simulated trajectories for a set of
#' \eqn{F = \{0, 0.75, 0.95, 1\}\times F_{MSY}} if B\eqn{_0} < B\eqn{_{MSY}}
#' (or \eqn{F = \{2, 1.25, 1.05, 1\}\times F_{MSY}} if B\eqn{_0} > B\eqn{_{MSY}}),
#' and marks the crossing time (triangle marker on the x-axis). A compact theme
#' is applied and the legend is locked to the **top-right inside** the panel.
#'
#' @section Dependencies:
#' Requires \pkg{ggplot2}. Uses \pkg{spict} helpers:
#' `check.rep()`, `get.par()`, `calc.gamma()`, and `get.version()`.
#'
#' @examples
#' \dontrun{
#'   rep <- fit.spict(inp)
#'   plot_elu2_panel_plotspict.tc_gg(rep)
#' }
#'
#' @import ggplot2
#' @importFrom grid unit
#' @importFrom grDevices adjustcolor
#' @importFrom spict check.rep get.par calc.gamma get.version
#' @export
plot_elu2_panel_time_to_bmsy <- function(rep,
                                            main  = "Time to Bmsy",
                                            stamp = get.version(),
                                            CI    = 0.95) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotspict.tc_gg(). Please install it.")
  }

  check.rep(rep)
  if ('sderr' %in% names(rep) || rep$opt$convergence != 0) return(invisible(NULL))

  inp   <- rep$inp
  B0cur <- get.par('logBl',  rep, exp = TRUE, CI = CI)[2]
  Kest  <- get.par('logK',   rep, exp = TRUE, CI = CI)
  m     <- get.par('logm',   rep, exp = TRUE, CI = CI)
  mmean <- tail(m, 1)[2]
  n     <- get.par('logn',   rep, exp = TRUE, CI = CI)
  gamma <- calc.gamma(n[2])
  sdb   <- get.par('logsdb', rep, exp = TRUE, CI = CI)
  Fmsy  <- tail(get.par('logFmsy', rep, exp = TRUE, CI = CI), 1)
  Bmsy  <- tail(get.par('logBmsy', rep, exp = TRUE, CI = CI), 1)

  if (is.na(Bmsy[2])) return(invisible(NULL))

  do.flag <- TRUE
  if (B0cur < Bmsy[2])  do.flag <- (B0cur / Bmsy[2]) <= 0.95
  if (B0cur > Bmsy[2])  do.flag <- (Bmsy[2] / B0cur) <= 0.95
  if (!do.flag) return(invisible(NULL))

  if (B0cur < Bmsy[2]) facvec <- c(0, 0.75, 0.95, 1) else facvec <- c(2, 1.25, 1.05, 1)
  cols   <- c('green3', 'blue', 'red', 'orange', 5:8)
  Fvec   <- round(facvec * Fmsy[2], 4)
  nFvec  <- length(Fvec)

  gstep <- function(F, K, m, n, sdb, B0, dt) {
    exp(log(B0) + (gamma * m / K - gamma * m / K * (B0 / K)^(n - 1) - F - 0.5 * sdb^2) * dt)
  }
  simdt <- 0.01
  nt    <- 10000
  Bsim  <- matrix(0, nFvec, nt)
  time  <- matrix(0, nFvec, nt)
  for (i in 1:nFvec) {
    time[i, ]  <- seq(0, simdt * (nt - 1), by = simdt)
    Bsim[i, 1] <- B0cur
    for (j in 2:nt) {
      Bsim[i, j] <- gstep(Fvec[i], Kest[2], mmean, n[2], sdb[2], Bsim[i, j - 1], simdt)
    }
  }
  Bsim <- Bsim / Bmsy[2]

  frac <- 0.95
  if (B0cur < Bmsy[2]) inds <- which(Bsim[nFvec, ] < 0.99) else inds <- which(Bsim[nFvec, ] > (1 / 0.99))
  if (!length(inds)) return(invisible(NULL))
  ylim <- range(Bsim[nFvec, ], na.rm = TRUE)
  xlim <- range(time[nFvec, inds]); xlim[2] <- min(xlim[2], 15)

  yr <- diff(ylim); if (!is.finite(yr) || yr <= 0) yr <- 1
  xr <- diff(xlim); if (!is.finite(xr) || xr <= 0) xr <- 1
  ypad <- max(0.02 * yr, 0.02)
  xpad <- max(0.02 * xr, 0.10)
  if (ylim[1] >= 1) {
    ylim[1] <- 1 - ypad
  } else if (abs(ylim[1] - 1) < 1e-8) {
    ylim[1] <- ylim[1] - ypad
  }
  if (xlim[1] >= 0) xlim[1] <- 0 - xpad

  vt <- numeric(nFvec)
  if (B0cur < Bmsy[2]) {
    for (i in 1:nFvec) vt[i] <- time[i, max(which(Bsim[i, ] < frac))]
  } else {
    for (i in 1:nFvec) vt[i] <- time[i, max(which(1 / Bsim[i, ] < frac))]
  }

  lab_vec <- paste('F =', facvec, 'x Fmsy')
  line_df <- do.call(rbind, lapply(seq_len(nFvec), function(i) {
    data.frame(scenario = lab_vec[i], time = time[i, ], prop = Bsim[i, ], stringsAsFactors = FALSE)
  }))
  vt_df <- data.frame(scenario = lab_vec, vt = vt, y0 = ylim[1])

  # original dynamic positioning (kept, but overridden by compact theme below)
  lg_pos <- if (B0cur > Bmsy[2]) c(0.98, 0.98) else c(0.98, 0.02)
  lg_jst <- if (B0cur > Bmsy[2]) c(1, 1) else c(1, 0)

  col_map <- setNames(cols[seq_len(nFvec)], lab_vec)
  hvals <- c(frac, 1/frac); hvals <- hvals[hvals >= ylim[1] & hvals <= ylim[2]]

  p <- ggplot2::ggplot() +
    (if (length(hvals))
      ggplot2::geom_hline(yintercept = hvals, colour = "lightgray", linewidth = 0.6)
     else NULL) +
    ggplot2::geom_hline(yintercept = 1, linetype = 3, linewidth = 0.7, colour = "black") +
    ggplot2::geom_line(
      data = line_df,
      ggplot2::aes(x = time, y = prop, colour = scenario, group = scenario),
      linewidth = 0.7, lineend = "butt", na.rm = TRUE
    ) +
    ggplot2::geom_vline(
      data = vt_df,
      ggplot2::aes(xintercept = vt, colour = scenario),
      linetype = 2, linewidth = 0.5, show.legend = FALSE, na.rm = TRUE
    ) +
    ggplot2::geom_point(
      data = vt_df,
      ggplot2::aes(x = vt, y = y0, colour = scenario),
      shape = 4, size = 1.2, stroke = 0.8, show.legend = FALSE, na.rm = TRUE
    ) +
    ggplot2::labs(title = main, x = "Years to Bmsy", y = "Proportion of Bmsy") +
    ggplot2::scale_colour_manual(values = col_map, breaks = names(col_map), labels = names(col_map)) +
    ggplot2::theme_light(base_size = 11) +
    ggplot2::theme(
      legend.title         = ggplot2::element_blank(),
      legend.position      = lg_pos,
      legend.justification = lg_jst,
      panel.border         = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1),
      plot.title           = ggplot2::element_text(hjust = 0.5)
    ) +
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE, clip = "on") +
    ggplot2::annotate("text", x = Inf, y = -Inf, label = stamp, hjust = 1.02, vjust = -0.8, size = 3) +

    # ---- Compact in-panel legend and styling override ----
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
      legend.text          = ggplot2::element_text(size = 9, face = "bold")
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        override.aes = list(linewidth = 1.1),
        ncol = 1, byrow = TRUE, keywidth = 1.2, keyheight = 0.6
      )
    )

  print(p)
  invisible(p)
}

# ---- Internal compact theme (file-local helper) ----

#' A compact minimal theme for SPiCT faceting/legends used in ELU2 panels
#'
#' @keywords internal
#' @noRd
.spict_theme_minimal_compact22 <- function(base_size = 10, base_family = "") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12),
      axis.title = ggplot2::element_text(face = "bold", size = 12),
      axis.text  = ggplot2::element_text(size = 10, face = "bold"),
      legend.position = "none",
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(fill = NA, colour = "grey35", linewidth = 2),
      axis.ticks = ggplot2::element_line(linewidth = 0.5, color = "grey35"),
      axis.ticks.length = grid::unit(3, "pt"),
      strip.background = ggplot2::element_rect(fill = "grey35", color = "grey35", linewidth = 0.5),
      strip.text = ggplot2::element_text(face = "bold", size = ggplot2::rel(1)),
      text = ggplot2::element_text(face = "bold", size = 10),
      plot.margin = ggplot2::margin(3, 3, 3, 3)
    )
}
