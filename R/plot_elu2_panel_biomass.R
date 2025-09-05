#' Biomass panel (with Bmsy band and B/Bmsy-on-B ribbon)
#'
#' @description
#' Plot absolute biomass \eqn{B_t} with:
#' \itemize{
#'   \item a light gray band for the uncertainty on \eqn{B_{MSY}} over time,
#'   \item a semi-transparent blue ribbon for \eqn{B/B_{MSY}} mapped onto the biomass scale
#'         with solid blue boundary lines,
#'   \item an overlaid solid line for \eqn{B_{MSY}(t)} (or constant if not time-varying),
#'   \item optional q-scaled index points,
#'   \item a thin grey vertical line at the end of observations.
#' }
#'
#' @param model A fitted SPiCT object (`spictcls`).
#' @param line_color Line color for estimates. Default `"blue"`.
#' @param show_CIs Logical; draw ribbons & CI edges. Default `TRUE`.
#' @param CI Confidence level for intervals (default `0.95`).
#'
#' @details
#' Assumes the following SPiCT helper functions are available in the search path:
#' `get.par()`, `get.manmax()`, and `get.msyvec()`.
#'
#' @return A ggplot object.
#' @export
#' @import ggplot2
plot_elu2_panel_biomass <- function(model,
                                    line_color = "blue",
                                    show_CIs   = TRUE,
                                    CI         = 0.95) {
  stopifnot(inherits(model, "spictcls"))
  vline_col  <- "grey50"; vline_size <- 0.2
  manflag <- ("man" %in% names(model))

  Best <- get.par("logB",     model, exp = TRUE, CI = CI)
  BB   <- get.par("logBBmsy", model, exp = TRUE, CI = CI)

  repmax   <- if (manflag) get.manmax(model) else model
  Bmsy_all <- get.par("logBmsy", repmax, exp = TRUE, CI = CI)
  Bmsyvec  <- get.msyvec(repmax$inp, Bmsy_all)
  Bmsy     <- if (!is.null(nrow(Bmsy_all))) Bmsy_all[1, ] else Bmsy_all

  df_B <- data.frame(
    time = .spict_time_from_par(model, Best),
    lwr  = Best[, 1], est = Best[, 2], upr = Best[, 3]
  )

  inp <- model$inp
  ind_in <- inp$indest
  if (length(ind_in) > 0) ind_in <- ind_in[seq_len(length(ind_in) - 1)]
  ind_pr <- if (!manflag) inp$indpred else integer(0)

  df_B_in <- if (length(ind_in)) df_B[ind_in, , drop = FALSE] else df_B[0, ]
  df_B_pr <- if (length(ind_pr)) df_B[ind_pr, , drop = FALSE] else df_B[0, ]

  Bmsy_est <- Bmsy[2]
  df_BB_rib <- .spict_make_ribbon_df(
    time = .spict_time_from_par(model, BB),
    lwr  = BB[, 1] * Bmsy_est,
    upr  = BB[, 3] * Bmsy_est
  )

  df_Bmsy_band <- data.frame(time = repmax$inp$time, ymin = Bmsyvec$ll, ymax = Bmsyvec$ul)
  df_Bmsy_line <- data.frame(time = repmax$inp$time, y = Bmsyvec$msy)

  p <- ggplot()

  # --- Bands / ribbons ---
  if (isTRUE(show_CIs)) {
    # BMSY band: lightgray (matches base R polygon 'lightgray')
    p <- p + geom_ribbon(
      data = df_Bmsy_band,
      aes(x = time, ymin = ymin, ymax = ymax),
      fill = "lightgray", colour = NA
    )

    # B/BMSY-on-B ribbon: fill rgb(0,0,1,0.10), solid boundary lines rgb(0,0,1,0.20)
    if (nrow(df_BB_rib)) {
      p <- p +
        geom_ribbon(
          data = df_BB_rib,
          aes(x = time, ymin = ymin, ymax = ymax),
          fill = grDevices::rgb(0, 0, 1, 0.10), colour = NA
        ) +
        geom_line(
          data = transform(df_BB_rib, y = ymin),
          aes(x = time, y = y),
          color = grDevices::rgb(0, 0, 1, 0.20), linewidth = 0.6, linetype = "solid"
        ) +
        geom_line(
          data = transform(df_BB_rib, y = ymax),
          aes(x = time, y = y),
          color = grDevices::rgb(0, 0, 1, 0.20), linewidth = 0.6, linetype = "solid"
        )
    }
  }

  # BMSY mean line
  p <- p + geom_line(data = df_Bmsy_line, aes(x = time, y = y), color = "black", linewidth = 0.7)

  # --- Bt (absolute) ---
  # Estimation window: mean solid; CI dashed (matches lty=2)
  if (nrow(df_B_in)) {
    if (isTRUE(show_CIs)) {
      p <- p +
        geom_line(data = df_B_in, aes(x = time, y = lwr),
                  linetype = "dashed", linewidth = 0.6, color = line_color) +
        geom_line(data = df_B_in, aes(x = time, y = upr),
                  linetype = "dashed", linewidth = 0.6, color = line_color)
    }
    p <- p + geom_line(data = df_B_in, aes(x = time, y = est),
                       linewidth = 0.8, color = line_color) # solid
  }

  # Prediction window: mean dotted (matches lty=3); CI dashed (same as base)
  if (!manflag && nrow(df_B_pr)) {
    p <- p +
      geom_line(data = df_B_pr, aes(x = time, y = est),
                linetype = "dotted", linewidth = 0.8, color = line_color)
    if (isTRUE(show_CIs)) {
      p <- p +
        geom_line(data = df_B_pr, aes(x = time, y = lwr),
                  linetype = "dashed", linewidth = 0.6, color = line_color) +
        geom_line(data = df_B_pr, aes(x = time, y = upr),
                  linetype = "dashed", linewidth = 0.6, color = line_color)
    }
  }

  # Thin grey vertical line at end of observations
  obs_end <- .spict_obs_end_overall(model)
  if (is.finite(obs_end)) {
    p <- p + geom_vline(xintercept = obs_end, color = vline_col, linewidth = vline_size)
  }

  # q-scaled index points (optional visuals kept)
  idx <- .spict_index_points(model)
  if (length(idx) >= 1) p <- p + geom_point(
    data = idx[[1]], aes(x = time, y = obs),
    color = "blue", shape = 16, size = 2, inherit.aes = FALSE
  )
  if (length(idx) >= 2) p <- p + geom_point(
    data = idx[[2]], aes(x = time, y = obs),
    shape = 22, color = "black", fill = "green", size = 2, stroke = 0.5, inherit.aes = FALSE
  )

  p <- p +
    labs(title = "Absolute biomass", x = "Year", y = expression(B[t])) +
    .spict_theme_minimal_compact2() +
    scale_y_continuous(sec.axis = sec_axis(~ . / Bmsy_est, name = expression(B[t]/B[MSY])))

  p
}

# ---- Inlined helpers (only those used by this panel) ----

#' Minimal compact ggplot2 theme (ELU2)
#' @param base_size Base font size
#' @param base_family Base font family
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
#' @param model Fitted SPiCT object
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
#' @param time Time vector
#' @param lwr Lower bound
#' @param upr Upper bound
#' @keywords internal
#' @noRd
.spict_make_ribbon_df <- function(time, lwr, upr) {
  ok <- is.finite(time) & is.finite(lwr) & is.finite(upr)
  data.frame(time = time[ok], ymin = lwr[ok], ymax = upr[ok])
}

#' q-scaled index points for first two indices (if present)
#' @param model Fitted SPiCT object
#' @keywords internal
#' @noRd
.spict_index_points <- function(model) {
  out <- list()
  inp <- model$inp
  if (!is.null(inp$timeI) && length(inp$timeI) >= 1) {
    qest <- try(get.par("logq", model, exp = TRUE), silent = TRUE)
    if (!inherits(qest, "try-error")) {
      for (i in seq_along(inp$timeI)) {
        qrow <- if (!is.null(inp$mapq) && length(inp$mapq) >= i) inp$mapq[i] else i
        qfac <- if (is.finite(qest[qrow, 2])) qest[qrow, 2] else 1
        out[[i]] <- data.frame(
          time = inp$timeI[[i]],
          obs  = inp$obsI[[i]] / qfac,
          idx  = i
        )
      }
    }
  }
  out
}

#' Prefer rownames(time) else fallback to inp$time
#' @param model Fitted SPiCT object
#' @param par_matrix Matrix returned by get.par(...)
#' @keywords internal
#' @noRd
.spict_time_from_par <- function(model, par_matrix) {
  rn <- suppressWarnings(as.numeric(rownames(par_matrix)))
  if (length(rn) == nrow(par_matrix) && !all(is.na(rn))) return(rn)
  t_full <- as.numeric(model$inp$time)
  if (length(t_full) >= nrow(par_matrix)) return(t_full[seq_len(nrow(par_matrix))])
  return(rep_len(t_full, nrow(par_matrix)))
}
