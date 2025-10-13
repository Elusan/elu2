# ---- Local helpers (tiny, internal) ------------------------------------------

# Truth overlay color (semi-transparent orange) — same as in plotspict.*
true.col <- function() grDevices::rgb(1, 165/255, 0, alpha = 0.7)

# Prefer rownames(time) else fallback to inp$time
.spict_time_from_par <- function(model, par_matrix) {
  rn <- suppressWarnings(as.numeric(rownames(par_matrix)))
  if (length(rn) == nrow(par_matrix) && any(is.finite(rn))) return(rn)
  t_full <- as.numeric(if (!is.null(model$inp$time)) model$inp$time else numeric(0))
  if (length(t_full) >= nrow(par_matrix)) return(t_full[seq_len(nrow(par_matrix))])
  if (length(t_full)) return(rep_len(t_full, nrow(par_matrix)))
  seq_len(nrow(par_matrix)) # last resort
}

# Overall end-of-observation time
.spict_obs_end_overall <- function(model) {
  inp <- model$inp
  tr  <- if (!is.null(inp$timerangeObs)) inp$timerangeObs else c(NA_real_, NA_real_)
  if (length(tr) >= 2 && is.finite(tr[2])) return(tr[2])
  if (!is.null(inp$timeC) && length(inp$timeC)) return(utils::tail(inp$timeC, 1))
  if (!is.null(inp$timeI) && length(inp$timeI) && length(inp$timeI[[1]]) > 0)
    return(utils::tail(inp$timeI[[1]], 1))
  if (!is.null(inp$time) && length(inp$time)) return(max(inp$time, na.rm = TRUE))
  NA_real_
}

# Build a ribbon df from lwr/upr; clamp to >= 0
.spict_make_ribbon_df <- function(time, lwr, upr) {
  ok <- is.finite(time) & is.finite(lwr) & is.finite(upr)
  data.frame(time = time[ok], ymin = pmax(lwr[ok], 0), ymax = pmax(upr[ok], 0))
}

# ---- Main panel ---------------------------------------------------------------

#' Absolute F[t] panel with Fmsy band, scaled ribbons, and scenarios
#'
#' @description
#' Plots absolute fishing mortality \eqn{F_t} with an Fmsy reference band
#' (constant or time-varying), plus ribbons derived from \eqn{(F/F_{MSY})\times F_{MSY}}.
#' If management scenarios exist in \code{rep_man$man}, their trajectories are
#' overlaid with a stable, repeatable colour mapping. For raw fits (no scenarios),
#' in-sample \eqn{F_t} is solid and one-step-ahead is dotted; a bridge segment is drawn.
#' When simulation truth is available in \code{rep_man$inp$true}, the true
#' \eqn{F_t} path and horizontal Fmsy are overlaid in semi-transparent orange,
#' matching \strong{plotspict} styling.
#'
#' @param rep_man A fitted SPiCT report object (\code{spictcls}), optionally with \code{$man}.
#' @param scenario_color Optional named character vector for scenario colours; if \code{NULL},
#'   a repeatable internal palette is used.
#' @param show_CIs Logical; draw CI bands/edges where applicable. Default \code{TRUE}.
#' @param CI Confidence level in (0, 1). Default \code{0.95}.
#' @param show_legend Logical; show scenario legend when scenarios exist. Default \code{FALSE}.
#'
#' @return A \code{ggplot} object.
#'
#' @details
#' - Fmsy handling supports both constant and time-varying growth configurations.
#' - Sub-annual \code{Fs} is drawn (semi-transparent) when \code{dtc < 1}.
#' - Vertical management lines are added using \code{inp$maninterval} when present,
#'   otherwise a line at the last observation time is drawn.
#'
#' @examples
#' \dontrun{
#' p <- my_plot_manage_f_panel_sim(fit)
#' print(p)
#' }
#'
#' @import ggplot2
#' @importFrom grDevices rgb col2rgb
#' @importFrom stats approx median
#' @importFrom utils tail
#' @export
my_plot_manage_f_panel_sim <- function(rep_man,
                                       scenario_color = NULL,
                                       show_CIs = TRUE,
                                       CI = 0.95,
                                       show_legend = FALSE) {
  stopifnot(inherits(rep_man, "spictcls"))

  theme_minimal_compact2_good_local <- function(base_size = 10, base_family = "") {
    ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12),
        axis.title = ggplot2::element_text(face = "bold", size = 12),
        axis.text  = ggplot2::element_text(size = 10, face = "bold"),
        legend.position = c(0.4, 0.98),
        legend.justification = c("left", "top"),
        legend.background = ggplot2::element_blank(),
        legend.box.background = ggplot2::element_blank(),
        legend.title = ggplot2::element_blank(),
        legend.text  = ggplot2::element_text(size = 10),
        legend.key.size = grid::unit(1, "lines"),
        panel.grid = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(fill = NA, colour = "grey45", linewidth = 1.2),
        axis.ticks = ggplot2::element_line(linewidth = 0.5, color = "grey45"),
        axis.ticks.length = grid::unit(3, "pt"),
        strip.background = ggplot2::element_rect(fill = "grey45", color = "grey45", linewidth = 0.5),
        strip.text = ggplot2::element_text(face = "bold", size = ggplot2::rel(1)),
        text = ggplot2::element_text(face = "bold", size = 10),
        plot.margin = ggplot2::margin(3, 3, 3, 3)
      )
  }

  vline_col  <- "grey75"; vline_size <- 0.4
  inp <- rep_man$inp
  manflag <- isTRUE(("man" %in% names(rep_man)) && length(rep_man$man) > 0)

  # --- F[t] absolute (notS)
  Fest <- elu2::get.par("logFnotS", rep_man, exp = TRUE, CI = CI)
  dfF  <- data.frame(time = .spict_time_from_par(rep_man, Fest),
                     lwr = Fest[,1], est = Fest[,2], upr = Fest[,3])

  # In-sample vs one-ahead split (raw fit only)
  ind_in <- if (!is.null(inp$indest)) inp$indest else integer(0)
  if (length(ind_in) > 0) ind_in <- ind_in[seq_len(length(ind_in) - 1)]
  ind_pr <- if (!manflag && !is.null(inp$indpred)) inp$indpred else integer(0)
  dfF_in <- if (length(ind_in)) dfF[ind_in, , drop = FALSE] else dfF[0, ]
  dfF_pr <- if (length(ind_pr)) dfF[ind_pr, , drop = FALSE] else dfF[0, ]

  # --- Fmsy band & mean
  repmax <- if (exists("get.manmax", where = asNamespace("spict"), inherits = FALSE))
    try(elu2::get.manmax(rep_man), silent = TRUE) else rep_man
  if (inherits(repmax, "try-error") || is.null(repmax$inp)) repmax <- rep_man

  tvgflag <- isTRUE(repmax$inp$timevaryinggrowth) || isTRUE(repmax$inp$logmcovflag)

  get_Fmsy_constant <- function() {
    Fmsy_all <- try(elu2::get.par("logFmsy", repmax, exp = TRUE, CI = CI), silent = TRUE)
    if (inherits(Fmsy_all, "try-error") || any(!is.finite(as.numeric(Fmsy_all)))) {
      Fmsy_all <- try(elu2::get.par("logFmsyd", repmax, exp = TRUE, CI = CI), silent = TRUE)
    }
    if (inherits(Fmsy_all, "try-error")) return(NULL)
    as.numeric(Fmsy_all[1, 2]) # est
  }

  if (tvgflag) {
    Fmsyvec <- try(elu2::get.par("logFmsyvec", repmax, exp = TRUE, CI = CI), silent = TRUE)
    if (inherits(Fmsyvec, "try-error")) {
      # fallback: pretend constant if vec unavailable
      Fmsy_est <- get_Fmsy_constant()
      Fmsy_band <- data.frame(time = repmax$inp$time,
                              ymin = rep(Fmsy_est, length(repmax$inp$time)),
                              ymax = rep(Fmsy_est, length(repmax$inp$time)))
      Fmsy_line <- data.frame(time = repmax$inp$time, y = rep(Fmsy_est, length(repmax$inp$time)))
      sec_axis_obj <- if (is.finite(Fmsy_est)) ggplot2::sec_axis(~ . / Fmsy_est, name = expression(F[t]/F[MSY])) else ggplot2::waiver()
      Fmsy_series <- function(x) rep(Fmsy_est, length(x))
    } else {
      tt_fmsy <- .spict_time_from_par(repmax, Fmsyvec)
      Fmsy_band <- data.frame(time = tt_fmsy, ymin = Fmsyvec[,1], ymax = Fmsyvec[,3])
      Fmsy_line <- data.frame(time = tt_fmsy, y = Fmsyvec[,2])
      sec_axis_obj <- ggplot2::waiver() # time-varying: no static sec axis
      # Build a time-matched series function for (F/Fmsy) scaling
      Fmsy_series <- function(x) {
        stats::approx(x = tt_fmsy, y = Fmsyvec[,2], xout = x, rule = 2)$y
      }
    }
  } else {
    Fmsy_est <- get_Fmsy_constant()
    # Try to use msyvec for CI band; if missing, use a flat band
    Fmsy_all <- try(elu2::get.par("logFmsy", repmax, exp = TRUE, CI = CI), silent = TRUE)
    if (inherits(Fmsy_all, "try-error")) Fmsy_all <- try(spict::get.par("logFmsyd", repmax, exp = TRUE, CI = CI), silent = TRUE)
    if (!inherits(Fmsy_all, "try-error")) {
      mv <- try(elu2::get.msyvec(repmax$inp, Fmsy_all), silent = TRUE)
    } else mv <- NULL
    if (!inherits(mv, "try-error") && !is.null(mv)) {
      Fmsy_band <- data.frame(time = repmax$inp$time, ymin = mv$ll, ymax = mv$ul)
      Fmsy_line <- data.frame(time = repmax$inp$time, y = mv$msy)
    } else {
      Fmsy_band <- data.frame(time = repmax$inp$time,
                              ymin = rep(Fmsy_est, length(repmax$inp$time)),
                              ymax = rep(Fmsy_est, length(repmax$inp$time)))
      Fmsy_line <- data.frame(time = repmax$inp$time, y = rep(Fmsy_est, length(repmax$inp$time)))
    }
    sec_axis_obj <- if (is.finite(Fmsy_est)) ggplot2::sec_axis(~ . / Fmsy_est, name = expression(F[t]/F[MSY])) else ggplot2::waiver()
    Fmsy_series <- function(x) rep(Fmsy_est, length(x))
  }

  # --- (F/Fmsy) × Fmsy ribbons (pre & post)
  FFrel_all <- try(spict::get.par("logFFmsynotS", rep_man, exp = TRUE, CI = CI), silent = TRUE)
  if (inherits(FFrel_all, "try-error")) {
    FFrel_all <- spict::get.par("logFFmsy", rep_man, exp = TRUE, CI = CI)
  }
  FFrel_all <- FFrel_all[, 1:3, drop = FALSE]
  tt_rel    <- .spict_time_from_par(rep_man, FFrel_all)

  last_obs_time <- .spict_obs_end_overall(rep_man)
  man_left <- if (!is.null(inp$maninterval) && length(inp$maninterval) >= 1L) inp$maninterval[1] else NA_real_
  t_pre_end <- if (manflag && is.finite(man_left)) man_left else last_obs_time

  # PRE mask: vectorized
  pre_mask <- is.finite(tt_rel) & (tt_rel <= t_pre_end)

  # POST mask: guard with a scalar, then make a vector mask (no && on vectors)
  has_post <- (!manflag) && (length(ind_pr) > 0)
  post_mask <- if (has_post) {
    is.finite(tt_rel) & (tt_rel %in% tt_rel[ind_pr])
  } else {
    rep(FALSE, length(tt_rel))
  }

  # scale FFrel by matched Fmsy at the same times
  FFabs_all <- FFrel_all
  FFabs_all[,1] <- FFrel_all[,1] * Fmsy_series(tt_rel)
  FFabs_all[,2] <- FFrel_all[,2] * Fmsy_series(tt_rel)
  FFabs_all[,3] <- FFrel_all[,3] * Fmsy_series(tt_rel)

  ribFF_pre  <- if (any(pre_mask))  .spict_make_ribbon_df(tt_rel[pre_mask],  FFabs_all[pre_mask,1],  FFabs_all[pre_mask,3]) else .spict_make_ribbon_df(numeric(0),numeric(0),numeric(0))
  ribFF_post <- if (any(post_mask)) .spict_make_ribbon_df(tt_rel[post_mask], FFabs_all[post_mask,1], FFabs_all[post_mask,3]) else .spict_make_ribbon_df(numeric(0),numeric(0),numeric(0))


  # --- Scenario Ft paths
  sc_order <- c("currentCatch","currentF","Fmsy","noF","reduceF25","increaseF25","msyHockeyStick","ices")
  sc_present <- if (manflag) names(rep_man$man) else character(0)
  sc_core   <- intersect(sc_order, sc_present)
  sc_other  <- setdiff(sc_present, sc_core)
  scs_final <- c(sc_core, sort(sc_other))

  df_scen_F <- NULL
  if (length(scs_final)) {
    lst <- vector("list", length(scs_final))
    for (i in seq_along(scs_final)) {
      sc <- scs_final[i]
      rp <- rep_man$man[[sc]]
      Fi <- elu2::get.par("logFnotS", rp, exp = TRUE, CI = CI)
      ti <- .spict_time_from_par(rp, Fi)
      lst[[i]] <- data.frame(time = ti, est = Fi[,2], scenario = sc)
    }
    df_scen_F <- do.call(rbind, lst)
    df_scen_F$scenario <- factor(df_scen_F$scenario, levels = scs_final)
  }

  # --- Colors
  man_cols <- function(n) {
    base <- c('darkmagenta','cyan3','darkgreen','coral1','black',
              'magenta','gold','green','cadetblue3','chocolate3',
              'darkolivegreen3','cyan','darkred')
    rep(base, length.out = n)
  }
  if (is.null(scenario_color)) {
    cols <- if (length(scs_final)) { x <- man_cols(length(scs_final)); names(x) <- scs_final; x } else NULL
  } else {
    cols <- scenario_color
    if (length(scs_final)) {
      if (is.null(names(cols))) names(cols) <- scs_final
      cols <- cols[match(scs_final, names(cols))]; names(cols) <- scs_final
    }
  }
  spict_blue_mean    <- "#0000FF"
  spict_blue_ci_line <- grDevices::rgb(0, 0, 1, 0.20)
  spict_blue_ci_fill <- grDevices::rgb(0, 0, 1, 0.10)

  # --- Build plot
  p <- ggplot2::ggplot()

  # Fmsy band & mean
  if (isTRUE(show_CIs) && nrow(Fmsy_band)) {
    p <- p + ggplot2::geom_ribbon(
      data = Fmsy_band, ggplot2::aes(x = time, ymin = pmax(ymin,0), ymax = pmax(ymax,0)),
      fill = "lightgray", colour = NA
    )
  }
  if (nrow(Fmsy_line)) {
    p <- p + ggplot2::geom_line(
      data = Fmsy_line, ggplot2::aes(x = time, y = y),
      color = "black", linewidth = 0.7
    )
  }

  # (F/Fmsy)×Fmsy ribbons (pre & post)
  if (isTRUE(show_CIs) && nrow(ribFF_pre)) {
    p <- p +
      ggplot2::geom_ribbon(data = ribFF_pre,  ggplot2::aes(x = time, ymin = ymin, ymax = ymax),
                           fill = spict_blue_ci_fill, colour = NA) +
      ggplot2::geom_line(  data = transform(ribFF_pre,  y = ymin), ggplot2::aes(x = time, y = y),
                           color = spict_blue_ci_line, linewidth = 0.6) +
      ggplot2::geom_line(  data = transform(ribFF_pre,  y = ymax), ggplot2::aes(x = time, y = y),
                           color = spict_blue_ci_line, linewidth = 0.6)
  }
  if (!manflag && isTRUE(show_CIs) && nrow(ribFF_post)) {
    p <- p +
      ggplot2::geom_ribbon(data = ribFF_post, ggplot2::aes(x = time, ymin = ymin, ymax = ymax),
                           fill = spict_blue_ci_fill, colour = NA) +
      ggplot2::geom_line(  data = transform(ribFF_post, y = ymin), ggplot2::aes(x = time, y = y),
                           color = spict_blue_ci_line, linewidth = 0.6) +
      ggplot2::geom_line(  data = transform(ribFF_post, y = ymax), ggplot2::aes(x = time, y = y),
                           color = spict_blue_ci_line, linewidth = 0.6)
  }

  # V-lines
  add_management_vlines_BF_good <- function(rep,
                                            color = "grey30",
                                            linetype = "dashed",
                                            linewidth = 0.6,
                                            lineend = "butt") {
    mi <- if (!is.null(rep$inp$maninterval)) rep$inp$maninterval else numeric(0)
    if (!length(mi) || !is.finite(mi[1])) return(list())
    left  <- mi[1]
    right <- if (length(mi) >= 2L && is.finite(mi[2])) mi[2] else NA_real_
    out <- list(ggplot2::geom_vline(
      xintercept = left, color = color, linetype = linetype,
      linewidth = linewidth, lineend = lineend
    ))
    if (is.finite(right)) {
      out[[length(out) + 1L]] <- ggplot2::geom_vline(
        xintercept = right, color = color, linetype = linetype,
        linewidth = linewidth, lineend = lineend
      )
    }
    out
  }
  if (manflag) {
    p <- p + add_management_vlines_BF_good(rep_man, color = vline_col,
                                           linetype = "solid", linewidth = vline_size, lineend = "butt")
  } else {
    obs_end <- .spict_obs_end_overall(rep_man)
    if (is.finite(obs_end)) {
      p <- p + ggplot2::geom_vline(xintercept = obs_end, color = vline_col,
                                   linetype = "solid", linewidth = vline_size)
    }
  }

  # Scenario F paths
  if (!is.null(df_scen_F) && nrow(df_scen_F)) {
    p <- p + ggplot2::geom_line(
      data = df_scen_F, ggplot2::aes(x = time, y = est, color = scenario),
      linewidth = 0.8, na.rm = TRUE
    )
  }

  # Sub-annual Ft (semi-transparent) if available and meaningful
  has_dtc <- !is.null(inp$dtc) && length(inp$dtc)
  if (has_dtc && suppressWarnings(min(inp$dtc, na.rm = TRUE)) < 1) {
    sF <- try(elu2::get.par("logFs", rep_man, exp = TRUE, CI = CI), silent = TRUE)
    if (!inherits(sF, "try-error")) {
      tt_sF <- .spict_time_from_par(rep_man, sF)
      mask  <- is.finite(tt_sF) & tt_sF <= t_pre_end
      if (any(mask)) {
        p <- p + ggplot2::geom_line(
          data = data.frame(time = tt_sF[mask], est = sF[mask, 2]),
          ggplot2::aes(x = time, y = est),
          color = grDevices::rgb(0,0,1,0.4), linewidth = 0.6
        )
      }
    }
  }

  # In-sample Ft mean/CI
  if (nrow(dfF_in)) {
    if (isTRUE(show_CIs)) {
      p <- p +
        ggplot2::geom_line(data = dfF_in, ggplot2::aes(x = time, y = lwr),
                           linetype = "dashed", linewidth = 0.6, colour = spict_blue_mean) +
        ggplot2::geom_line(data = dfF_in, ggplot2::aes(x = time, y = upr),
                           linetype = "dashed", linewidth = 0.6, colour = spict_blue_mean)
    }
    p <- p + ggplot2::geom_line(
      data = dfF_in, ggplot2::aes(x = time, y = est),
      linewidth = 0.9, colour = spict_blue_mean
    )
  }

  # Raw fit only: one-ahead dotted mean + bridge
  if (!manflag && nrow(dfF_pr)) {
    p <- p + ggplot2::geom_line(
      data = dfF_pr, ggplot2::aes(x = time, y = est),
      linetype = "dotted", linewidth = 0.8, colour = spict_blue_mean
    )
    obs_end <- .spict_obs_end_overall(rep_man)
    ttF <- dfF$time
    i0 <- if (any(ttF <= obs_end)) max(which(ttF <= obs_end)) else NA_integer_
    i1 <- if (any(ttF  > obs_end)) min(which(ttF  > obs_end)) else NA_integer_
    if (is.finite(i0) && is.finite(i1)) {
      p <- p + ggplot2::geom_segment(
        ggplot2::aes(x = ttF[i0], xend = ttF[i1], y = dfF$est[i0], yend = dfF$est[i1]),
        linetype = "dotted", color = spict_blue_mean, linewidth = 0.8, lineend = "round"
      )
    }
  }

  # ---------- Simulation overlays (exactly like plotspict.f) ----------
  if ("true" %in% names(inp) && !is.null(inp$true)) {
    tr <- inp$true
    # True F(t) trajectory
    if (!is.null(tr$time) && !is.null(tr$Fs) && length(tr$time) == length(tr$Fs)) {
      p <- p + ggplot2::geom_line(
        data = data.frame(time = tr$time, y = tr$Fs),
        ggplot2::aes(x = time, y = y),
        inherit.aes = FALSE, color = true.col(), linewidth = 1.0, show.legend = FALSE
      )
    }
    # True Fmsy horizontal (solid orange + black dotted echo)
    if (!is.null(tr$Fmsy) && is.finite(tr$Fmsy)) {
      p <- p +
        ggplot2::geom_hline(yintercept = tr$Fmsy, color = true.col(), linewidth = 0.9, show.legend = FALSE) +
        ggplot2::geom_hline(yintercept = tr$Fmsy, color = "black", linetype = "dotted", linewidth = 0.6, show.legend = FALSE)
    }
  }
  # --------------------------------------------------------------------

  # Legend, labels, axis, theme
  p +
    ggplot2::labs(title = "Absolute fishing mortality", x = "Year", y = expression(F[t])) +
    { if (length(scs_final)) ggplot2::scale_color_manual(
      values = cols, guide = if (show_legend) "legend" else "none"
    ) else ggplot2::scale_color_discrete(guide = "none") } +
    ggplot2::scale_fill_manual(values = if (is.null(cols)) NA else cols, guide = "none") +
    ggplot2::scale_y_continuous(sec.axis = sec_axis_obj,
                                expand = ggplot2::expansion(mult = c(0.02, 0.02))) +
    theme_minimal_compact2_good_local()
}
