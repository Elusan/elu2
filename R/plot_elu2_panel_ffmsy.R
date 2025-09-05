#' F/Fmsy panel (ribbon + outlines; shows sub-annual Ft/Fmsy if seasonal)
#'
#' @param model A fitted SPiCT object (`spictcls`).
#' @param line_color Line color for estimates.
#' @param show_CIs Logical; draw ribbons & dotted CI edges.
#' @param CI Confidence level for intervals (default 0.95).
#' @return A ggplot object.
#' @export
#' @import ggplot2
plot_elu2_panel_ffmsy <- function(model,
                                   line_color = "blue",
                                   show_CIs   = TRUE,
                                   CI         = 0.95) {
  stopifnot(inherits(model, "spictcls"))
  vline_col  <- "grey50"; vline_size <- 0.2
  manflag <- ("man" %in% names(model))
  inp <- model$inp

  FF <- get.par("logFFmsy", model, exp = TRUE, CI = CI)
  df <- data.frame(time = .spict_time_from_par(model, FF), lwr = FF[,1], est = FF[,2], upr = FF[,3])

  ind_in <- inp$indest
  if (length(ind_in) > 0) ind_in <- ind_in[seq_len(length(ind_in) - 1)]
  ind_pr <- if (!manflag) inp$indpred else integer(0)

  df_in <- if (length(ind_in)) df[ind_in, , drop = FALSE] else df[0, ]
  df_pr <- if (length(ind_pr)) df[ind_pr, , drop = FALSE] else df[0, ]

  rib <- .spict_make_ribbon_df(df$time, df$lwr, df$upr)

  p <- ggplot()
  if (isTRUE(show_CIs) && nrow(rib)) {
    p <- p + geom_ribbon(data = rib, aes(x = time, ymin = ymin, ymax = ymax), fill = grDevices::rgb(0,0,1,0.10)) +
      geom_line(data = transform(rib, y = ymin), aes(x = time, y = y), color = grDevices::rgb(0,0,1,0.20), linewidth = 0.6) +
      geom_line(data = transform(rib, y = ymax), aes(x = time, y = y), color = grDevices::rgb(0,0,1,0.20), linewidth = 0.6)
  }
  if (nrow(df_in)) p <- p + geom_line(data = df_in, aes(x = time, y = est), color = line_color, linewidth = 0.8)
  if (nrow(df_pr)) p <- p + geom_line(data = df_pr, aes(x = time, y = est), linetype = "dashed", color = line_color, linewidth = 0.8)

  obs_end <- .spict_obs_end_overall(model)
  if (is.finite(obs_end)) p <- p + geom_vline(xintercept = obs_end, color = vline_col, linewidth = vline_size)

  p <- p +
    geom_hline(yintercept = 1, linetype = "solid", color = "black", linewidth = 0.8) +
    labs(title = "Relative fishing mortality", x = "Year", y = expression(bold(F/F[MSY]))) +
    .spict_theme_minimal_compact2()

  # Seasonal overlay (unchanged)
  if (min(inp$dtc, na.rm = TRUE) < 1 && nrow(df_in)) {
    FFs <- get.par("logFFmsy", model, exp = TRUE, CI = CI)
    p <- p + geom_line(
      data = data.frame(time = df_in$time, est = FFs[seq_len(nrow(df_in)), 2]),
      aes(x = time, y = est),
      color = grDevices::rgb(0,0,1,0.4), linewidth = 0.6
    )
  }

  p
}

# ---- Inlined helpers (only those used by this panel) ----

#' Minimal compact ggplot2 theme (ELU2)
#' @keywords internal
#' @noRd
.spict_theme_minimal_compact2 <- function(base_size = 10, base_family = "") {
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

#' End of observation time across components
#' @keywords internal
#' @noRd
.spict_obs_end_overall <- function(model) {
  inp <- model$inp
  tr <- inp$timerangeObs
  if (!is.null(tr) && length(tr) >= 2 && is.finite(tr[2])) return(tr[2])
  if (!is.null(inp$timeC) && length(inp$timeC))       return(tail(inp$timeC, 1))
  if (!is.null(inp$timeI) && length(inp$timeI) && length(inp$timeI[[1]]) > 0) return(tail(inp$timeI[[1]], 1))
  if (!is.null(inp$time)  && length(inp$time))        return(max(inp$time, na.rm = TRUE))
  NA_real_
}

#' Build a ribbon dataframe from lwr/upr vectors
#' @keywords internal
#' @noRd
.spict_make_ribbon_df <- function(time, lwr, upr) {
  ok <- is.finite(time) & is.finite(lwr) & is.finite(upr)
  data.frame(time = time[ok], ymin = lwr[ok], ymax = upr[ok])
}

#' Prefer rownames(time) else fallback to inp$time
#' @keywords internal
#' @noRd
.spict_time_from_par <- function(model, par_matrix) {
  rn <- suppressWarnings(as.numeric(rownames(par_matrix)))
  if (length(rn) == nrow(par_matrix) && !all(is.na(rn))) return(rn)
  t_full <- as.numeric(model$inp$time)
  if (length(t_full) >= nrow(par_matrix)) return(t_full[seq_len(nrow(par_matrix))])
  return(rep_len(t_full, nrow(par_matrix)))
}
