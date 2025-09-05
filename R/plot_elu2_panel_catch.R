#' Catch panel (MSY band; dashed CI; dotted pred mean; no true lines)
#'
#' @param model A fitted SPiCT object (`spictcls`).
#' @param line_color Line color for estimates. Default `"blue"`.
#' @param show_CIs Logical; draw CI edges. Default `TRUE`.
#' @param CI Confidence level (default 0.95).
#' @param qlegend Logical; keep for compatibility (unused in ggplot).
#' @return A ggplot object.
#' @export
#' @import ggplot2
plot_elu2_panel_catch <- function(model,
                                  line_color = "blue",
                                  show_CIs   = TRUE,
                                  CI         = 0.95,
                                  qlegend    = TRUE) {
  stopifnot(inherits(model, "spictcls"))
  inp <- model$inp
  manflag <- ("man" %in% names(model))
  vline_col  <- "grey50"; vline_size <- 0.2

  # --- MSY time series (time-varying or constant) ---
  repmax <- if (manflag) get.manmax(model) else model
  tvgflag <- isTRUE(repmax$inp$timevaryinggrowth) || isTRUE(repmax$inp$logmcovflag)

  if (tvgflag) {
    MSY <- get.par("logMSYvec", repmax, exp = TRUE, CI = CI)
    MSYvec <- as.data.frame(MSY)
    MSYvec$msy <- MSYvec$est
    df_msy_band <- data.frame(
      time = .spict_time_from_par(repmax, MSY),
      ymin = MSY[, 1],
      ymax = MSY[, 3]
    )
    df_msy_line <- data.frame(
      time = .spict_time_from_par(repmax, MSY),
      y    = MSYvec$msy
    )
  } else {
    MSY  <- get.par("logMSY", repmax, exp = TRUE, CI = CI)
    mvec <- get.msyvec(repmax$inp, MSY)
    df_msy_band <- data.frame(time = repmax$inp$time, ymin = mvec$ll, ymax = mvec$ul)
    df_msy_line <- data.frame(time = repmax$inp$time, y = mvec$msy)
  }

  # --- Predicted catches (logCpred, exp=TRUE) ---
  Cpred <- get.par("logCpred", model, exp = TRUE, CI = CI)
  Cpred[Cpred < 0] <- 0
  if ("Cp" %in% names(model)) model$Cp[model$Cp < 0] <- 0

  ind_est <- which(inp$timeCpred <= tail(inp$timeC, 1))
  ind_prd <- which(inp$timeCpred >= tail(inp$timeC, 1))

  # Annualize if needed (matches base plotspict.catch)
  if (min(inp$dtc) < 1) {
    # Observed
    alo  <- annual(inp$timeC, inp$obsC / inp$dtc)
    timeo <- alo$anntime
    obs   <- alo$annvec

    # Estimation window (<= last obs)
    al1 <- annual(inp$timeCpred[ind_est], Cpred[ind_est, 1] / inp$dtcp[ind_est])
    al2 <- annual(inp$timeCpred[ind_est], Cpred[ind_est, 2] / inp$dtcp[ind_est])
    al3 <- annual(inp$timeCpred[ind_est], Cpred[ind_est, 3] / inp$dtcp[ind_est])
    inds <- which(!is.na(al2$annvec))
    time_est <- al2$anntime[inds]
    c_est    <- al2$annvec[inds]
    cl_est   <- al1$annvec[inds]
    cu_est   <- al3$annvec[inds]

    # Prediction window (>= last obs)
    al1p <- annual(inp$timeCpred[ind_prd], Cpred[ind_prd, 1] / inp$dtcp[ind_prd])
    al2p <- annual(inp$timeCpred[ind_prd], Cpred[ind_prd, 2] / inp$dtcp[ind_prd])
    al3p <- annual(inp$timeCpred[ind_prd], Cpred[ind_prd, 3] / inp$dtcp[ind_prd])
    inds <- which(!is.na(al2p$annvec))
    time_prd <- al2p$anntime[inds]
    c_prd    <- al2p$annvec[inds]
    cl_prd   <- al1p$annvec[inds]
    cu_prd   <- al3p$annvec[inds]

    # Merge in any dtc==1 exact points (same as base)
    if (any(inp$dtc == 1)) {
      inds1 <- which(inp$dtc == 1)
      timeo <- c(timeo, inp$timeC[inds1])
      obs   <- c(obs,   inp$obsC[inds1])

      time_uns <- c(time_est, inp$timeCpred[inds1])
      o <- sort(time_uns, index = TRUE)
      time_est <- o$x
      c_est    <- c(c_est, Cpred[inds1, 2])[o$ix]
      cl_est   <- c(cl_est, Cpred[inds1, 1])[o$ix]
      cu_est   <- c(cu_est, Cpred[inds1, 3])[o$ix]
    }
  } else {
    # No annualization needed
    timeo  <- inp$timeC
    obs    <- inp$obsC / inp$dtc

    time_est <- inp$timeCpred[ind_est]
    c_est    <- Cpred[ind_est, 2] / inp$dtcp[ind_est]
    cl_est   <- Cpred[ind_est, 1]
    cu_est   <- Cpred[ind_est, 3]

    time_prd <- inp$timeCpred[ind_prd]
    c_prd    <- Cpred[ind_prd, 2] / inp$dtcp[ind_prd]
    cl_prd   <- Cpred[ind_prd, 1]
    cu_prd   <- Cpred[ind_prd, 3]
  }

  # Data frames
  df_obs <- data.frame(time = timeo, obs = obs)
  df_est <- data.frame(time = time_est, est = c_est, lwr = cl_est, upr = cu_est)
  df_prd <- data.frame(time = time_prd, est = c_prd, lwr = cl_prd, upr = cu_prd)

  # Build plot
  p <- ggplot()

  # MSY band (light gray) + mean line (black, solid)
  p <- p +
    geom_ribbon(data = df_msy_band, aes(x = time, ymin = ymin, ymax = ymax),
                fill = "lightgray", colour = NA) +
    geom_line(data = df_msy_line, aes(x = time, y = y),
              color = "black", linewidth = 0.7)

  # Estimation window: mean solid; CI dashed
  if (nrow(df_est)) {
    if (isTRUE(show_CIs)) {
      p <- p +
        geom_line(data = df_est, aes(x = time, y = lwr),
                  linetype = "dashed", linewidth = 0.6, color = line_color) +
        geom_line(data = df_est, aes(x = time, y = upr),
                  linetype = "dashed", linewidth = 0.6, color = line_color)
    }
    p <- p +
      geom_line(data = df_est, aes(x = time, y = est),
                linewidth = 0.8, color = line_color)
  }

  # Prediction window: mean dotted; CI dashed (only if not in man mode and predictions exist)
  if (!manflag && isTRUE(inp$dtpredc > 0) && nrow(df_prd)) {
    p <- p +
      geom_line(data = df_prd, aes(x = time, y = est),
                linetype = "dotted", linewidth = 0.8, color = line_color)
    if (isTRUE(show_CIs)) {
      p <- p +
        geom_line(data = df_prd, aes(x = time, y = lwr),
                  linetype = "dashed", linewidth = 0.6, color = line_color) +
        geom_line(data = df_prd, aes(x = time, y = upr),
                  linetype = "dashed", linewidth = 0.6, color = line_color)
    }
  }

  # Observations (points)
  if (nrow(df_obs)) {
    p <- p + geom_point(data = df_obs, aes(x = time, y = obs),
                        size = 1.8, shape = 16, color = "black")
  }

  # Thin grey vertical line at the end of observations (aligned to last integer year)
  vline_time <- .spict_obs_end_for_catch(inp)
  if (is.finite(vline_time) && !manflag) {
    p <- p + geom_vline(xintercept = vline_time, color = vline_col, linewidth = vline_size)
  }

  # Labels and theme
  y_lab <- .spict_add_catchunit_label("Catch", inp$catchunit)
  p <- p +
    labs(title = "Catch", x = "Year", y = y_lab) +
    .spict_theme_minimal_compact2()

  p
}

# ----------------------- helpers (inlined) -----------------------

# minimal compact theme
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
      text = ggplot2::element_text(face = "bold", size = 10),
      plot.margin = ggplot2::margin(3, 3, 3, 3)
    )
}

# prefer rownames(time) else fallback to inp$time
#' @keywords internal
#' @noRd
.spict_time_from_par <- function(model, par_matrix) {
  rn <- suppressWarnings(as.numeric(rownames(par_matrix)))
  if (length(rn) == nrow(par_matrix) && !all(is.na(rn))) return(rn)
  t_full <- as.numeric(model$inp$time)
  if (length(t_full) >= nrow(par_matrix)) return(t_full[seq_len(nrow(par_matrix))])
  return(rep_len(t_full, nrow(par_matrix)))
}

# last integer-year boundary for obs end (matches base abline logic)
#' @keywords internal
#' @noRd
.spict_obs_end_for_catch <- function(inp) {
  if (is.null(inp$timeC) || !length(inp$timeC)) return(NA_real_)
  ix <- which(inp$timeC %% 1 == 0)
  if (!length(ix)) return(NA_real_)
  tail(inp$timeC[ix], 1)
}

# add catch unit to label (simple, ggplot-friendly)
#' @keywords internal
#' @noRd
.spict_add_catchunit_label <- function(lab, cu) {
  cu <- as.character(cu %||% "")
  if (nzchar(cu)) paste0(lab, ", ", cu) else lab
}

# safe pipe for NULL coalesce
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (!is.null(a)) a else b
