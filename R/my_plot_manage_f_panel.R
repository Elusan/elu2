# ---- Local helpers (tiny, internal) ------------------------------------------

# Accept a fitted spictcls OR a rep$man-style list of spictcls
.as_spict_like <- function(x) {
  if (inherits(x, "spictcls")) return(x)
  if (is.list(x) && length(x) && all(vapply(x, function(z) inherits(z, "spictcls"), logical(1)))) {
    first <- x[[1]]
    inp_fallback <- if (!is.null(first$inp)) first$inp else list()
    out <- list(inp = inp_fallback, man = x)
    class(out) <- "spictcls"
    return(out)
  }
  if (is.list(x) && !length(x)) {
    out <- list(inp = list(), man = list()); class(out) <- "spictcls"; return(out)
  }
  stop("Input must be a fitted spictcls, or a list like rep$man containing spictcls objects.")
}

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
  inp <- model$inp %||% list()
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

# Small `%||%` helper
`%||%` <- function(x, y) if (is.null(x)) y else x

# ---- Palette (unchanged) -----------------------------------------------------
man_cols <- function(n) {
  base <- c('darkmagenta','cyan3','darkgreen','coral1','black',
            'magenta','gold','green','cadetblue3','chocolate3',
            'darkolivegreen3','cyan','darkred')
  rep(base, length.out = n)
}

# ---- Main panel ---------------------------------------------------------------

#' Absolute F[t] panel (with Fmsy band, (F/Fmsy)×Fmsy ribbon, scenarios)
#' @export
my_plot_manage_f_panel <- function(rep_man,
                                   scenario_color = NULL,
                                   show_CIs = TRUE,
                                   CI = 0.95,
                                   show_legend = FALSE) {
  # accept fitted or rep$man
  rep_man <- .as_spict_like(rep_man)
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
  inp <- rep_man$inp %||% list()
  manflag <- isTRUE(("man" %in% names(rep_man)) && length(rep_man$man) > 0)

  # --- F[t] absolute (notS)
  Fest <- try(get.par("logFnotS", rep_man, exp = TRUE, CI = CI), silent = TRUE)
  if (inherits(Fest, "try-error")) Fest <- try(get.par("logF", rep_man, exp = TRUE, CI = CI), silent = TRUE)
  if (inherits(Fest, "try-error")) stop("Could not extract F time series from object.")

  dfF  <- data.frame(time = .spict_time_from_par(rep_man, Fest),
                     lwr = Fest[,1], est = Fest[,2], upr = Fest[,3])

  # In-sample vs one-ahead split (raw fit only)
  ind_in <- inp$indest %||% integer(0)
  if (length(ind_in) > 0) ind_in <- ind_in[seq_len(max(1, length(ind_in) - 1))]
  ind_pr <- if (!manflag) (inp$indpred %||% integer(0)) else integer(0)
  dfF_in <- if (length(ind_in)) dfF[ind_in, , drop = FALSE] else dfF[0, ]
  dfF_pr <- if (length(ind_pr)) dfF[ind_pr, , drop = FALSE] else dfF[0, ]

  # --- Fmsy band & mean (robust to time-varying / alt names)
  repmax <- try(if (exists("get.manmax", where = asNamespace("spict"), inherits = FALSE))
    get.manmax(rep_man) else rep_man, silent = TRUE)
  if (inherits(repmax, "try-error") || is.null(repmax$inp)) repmax <- rep_man

  tvgflag <- isTRUE(repmax$inp$timevaryinggrowth) || isTRUE(repmax$inp$logmcovflag)

  get_Fmsy_constant <- function() {
    Fmsy_all <- try(get.par("logFmsy", repmax, exp = TRUE, CI = CI), silent = TRUE)
    if (inherits(Fmsy_all, "try-error") || any(!is.finite(as.numeric(Fmsy_all)))) {
      Fmsy_all <- try(get.par("logFmsyd", repmax, exp = TRUE, CI = CI), silent = TRUE)
    }
    if (inherits(Fmsy_all, "try-error")) return(NA_real_)
    as.numeric(Fmsy_all[1, 2])
  }

  if (tvgflag) {
    Fmsyvec <- try(get.par("logFmsyvec", repmax, exp = TRUE, CI = CI), silent = TRUE)
    if (!inherits(Fmsyvec, "try-error")) {
      tt_fmsy <- .spict_time_from_par(repmax, Fmsyvec)
      Fmsy_band <- data.frame(time = tt_fmsy, ymin = Fmsyvec[,1], ymax = Fmsyvec[,3])
      Fmsy_line <- data.frame(time = tt_fmsy, y = Fmsyvec[,2])
      sec_axis_obj <- ggplot2::waiver()  # time-varying: no static sec axis
      Fmsy_series <- function(x) stats::approx(x = tt_fmsy, y = Fmsyvec[,2], xout = x, rule = 2)$y
    } else {
      Fmsy_est <- get_Fmsy_constant()
      tvec <- repmax$inp$time %||% rep_man$inp$time %||% numeric(0)
      Fmsy_band <- data.frame(time = tvec,
                              ymin = rep(Fmsy_est, length(tvec)),
                              ymax = rep(Fmsy_est, length(tvec)))
      Fmsy_line <- data.frame(time = tvec, y = rep(Fmsy_est, length(tvec)))
      sec_axis_obj <- if (is.finite(Fmsy_est)) ggplot2::sec_axis(~ . / Fmsy_est, name = expression(F[t]/F[MSY])) else ggplot2::waiver()
      Fmsy_series <- function(x) rep(Fmsy_est, length(x))
    }
  } else {
    Fmsy_est <- get_Fmsy_constant()
    Fmsy_all <- try(get.par("logFmsy", repmax, exp = TRUE, CI = CI), silent = TRUE)
    if (inherits(Fmsy_all, "try-error")) Fmsy_all <- try(get.par("logFmsyd", repmax, exp = TRUE, CI = CI), silent = TRUE)
    mv <- try(if (!inherits(Fmsy_all, "try-error")) get.msyvec(repmax$inp, Fmsy_all) else NULL, silent = TRUE)

    if (!inherits(mv, "try-error") && !is.null(mv)) {
      tvec <- repmax$inp$time %||% rep_man$inp$time %||% numeric(0)
      Fmsy_band <- data.frame(time = tvec, ymin = mv$ll, ymax = mv$ul)
      Fmsy_line <- data.frame(time = tvec, y = mv$msy)
    } else {
      tvec <- repmax$inp$time %||% rep_man$inp$time %||% numeric(0)
      Fmsy_band <- data.frame(time = tvec,
                              ymin = rep(Fmsy_est, length(tvec)),
                              ymax = rep(Fmsy_est, length(tvec)))
      Fmsy_line <- data.frame(time = tvec, y = rep(Fmsy_est, length(tvec)))
    }
    sec_axis_obj <- if (is.finite(Fmsy_est)) ggplot2::sec_axis(~ . / Fmsy_est, name = expression(F[t]/F[MSY])) else ggplot2::waiver()
    Fmsy_series <- function(x) rep(Fmsy_est, length(x))
  }

  # --- (F/Fmsy) × Fmsy ribbons (pre & post)
  FFrel_all <- try(get.par("logFFmsynotS", rep_man, exp = TRUE, CI = CI), silent = TRUE)
  if (inherits(FFrel_all, "try-error")) {
    FFrel_all <- try(get.par("logFFmsy", rep_man, exp = TRUE, CI = CI), silent = TRUE)
  }
  if (inherits(FFrel_all, "try-error")) {
    FFrel_all <- matrix(numeric(0), ncol = 3)
  } else {
    FFrel_all <- FFrel_all[, 1:3, drop = FALSE]
  }
  tt_rel <- if (nrow(FFrel_all)) .spict_time_from_par(rep_man, FFrel_all) else numeric(0)

  last_obs_time <- .spict_obs_end_overall(rep_man)
  man_left <- if (!is.null(inp$maninterval) && length(inp$maninterval) >= 1L) inp$maninterval[1] else NA_real_
  t_pre_end <- if (manflag && is.finite(man_left)) man_left else last_obs_time

  pre_mask  <- if (length(tt_rel)) is.finite(tt_rel) & (tt_rel <= t_pre_end) else logical(0)
  has_post  <- (!manflag) && length(ind_pr) > 0 && length(tt_rel) > 0
  post_mask <- if (has_post) is.finite(tt_rel) & (tt_rel %in% tt_rel[ind_pr]) else rep(FALSE, length(tt_rel))

  if (length(tt_rel)) {
    fmsy_at <- Fmsy_series(tt_rel)
    FFabs_all <- cbind(FFrel_all[,1] * fmsy_at,
                       FFrel_all[,2] * fmsy_at,
                       FFrel_all[,3] * fmsy_at)
  } else {
    FFabs_all <- matrix(numeric(0), ncol = 3)
  }

  ribFF_pre  <- if (any(pre_mask))  .spict_make_ribbon_df(tt_rel[pre_mask],  FFabs_all[pre_mask,1],  FFabs_all[pre_mask,3]) else .spict_make_ribbon_df(numeric(0),numeric(0),numeric(0))
  ribFF_post <- if (any(post_mask)) .spict_make_ribbon_df(tt_rel[post_mask], FFabs_all[post_mask,1], FFabs_all[post_mask,3]) else .spict_make_ribbon_df(numeric(0),numeric(0),numeric(0))

  # --- Scenario Ft paths (keep canonical order, then extras alpha)
  sc_order   <- c("currentCatch","currentF","Fmsy","noF","reduceF25","increaseF25","msyHockeyStick","ices")
  sc_present <- if (manflag) names(rep_man$man) else character(0)
  sc_core    <- intersect(sc_order, sc_present)
  sc_other   <- setdiff(sc_present, sc_core)
  scs_final  <- c(sc_core, sort(sc_other))

  df_scen_F <- NULL
  if (length(scs_final)) {
    lst <- vector("list", length(scs_final))
    for (i in seq_along(scs_final)) {
      sc <- scs_final[i]
      rp <- rep_man$man[[sc]]
      Fi <- try(get.par("logFnotS", rp, exp = TRUE, CI = CI), silent = TRUE)
      if (inherits(Fi, "try-error")) Fi <- try(get.par("logF", rp, exp = TRUE, CI = CI), silent = TRUE)
      if (inherits(Fi, "try-error")) next
      ti <- .spict_time_from_par(rp, Fi)
      lst[[i]] <- data.frame(time = ti, est = Fi[,2], scenario = sc)
    }
    lst <- Filter(NROW, lst)
    if (length(lst)) {
      df_scen_F <- do.call(rbind, lst)
      df_scen_F$scenario <- factor(df_scen_F$scenario, levels = scs_final, ordered = TRUE)
    } else {
      df_scen_F <- NULL
    }
  }

  # --- Colors (unchanged palette, but locked to order)
  if (is.null(scenario_color)) {
    cols <- if (length(scs_final)) { x <- man_cols(length(scs_final)); names(x) <- scs_final; x } else NULL
  } else {
    cols <- scenario_color
    if (length(scs_final)) {
      if (is.null(names(cols))) names(cols) <- scs_final
      cols <- cols[match(scs_final, names(cols))]
      # backfill any missing with your palette to keep look stable
      if (anyNA(cols)) {
        fill <- man_cols(length(scs_final))
        cols[is.na(cols)] <- fill[is.na(cols)]
      }
      names(cols) <- scs_final
    }
  }

  spict_blue_mean    <- "#0000FF"
  spict_blue_ci_line <- grDevices::rgb(0, 0, 1, 0.20)
  spict_blue_ci_fill <- grDevices::rgb(0, 0, 1, 0.10)

  # --- Build plot
  p <- ggplot2::ggplot()

  # Fmsy band & mean
  if (isTRUE(show_CIs) && NROW(Fmsy_band)) {
    p <- p + ggplot2::geom_ribbon(
      data = Fmsy_band, ggplot2::aes(x = time, ymin = pmax(ymin,0), ymax = pmax(ymax,0)),
      fill = "lightgray", colour = NA
    )
  }
  if (NROW(Fmsy_line)) {
    p <- p + ggplot2::geom_line(
      data = Fmsy_line, ggplot2::aes(x = time, y = y),
      color = "black", linewidth = 0.7
    )
  }

  # (F/Fmsy)×Fmsy ribbons (pre & post)
  if (isTRUE(show_CIs) && NROW(ribFF_pre)) {
    p <- p +
      ggplot2::geom_ribbon(data = ribFF_pre,  ggplot2::aes(x = time, ymin = ymin, ymax = ymax),
                           fill = spict_blue_ci_fill, colour = NA) +
      ggplot2::geom_line(  data = transform(ribFF_pre,  y = ymin), ggplot2::aes(x = time, y = y),
                           color = spict_blue_ci_line, linewidth = 0.6) +
      ggplot2::geom_line(  data = transform(ribFF_pre,  y = ymax), ggplot2::aes(x = time, y = y),
                           color = spict_blue_ci_line, linewidth = 0.6)
  }
  if (!manflag && isTRUE(show_CIs) && NROW(ribFF_post)) {
    p <- p +
      ggplot2::geom_ribbon(data = ribFF_post, ggplot2::aes(x = time, ymin = ymin, ymax = ymax),
                           fill = spict_blue_ci_fill, colour = NA) +
      ggplot2::geom_line(  data = transform(ribFF_post, y = ymin), ggplot2::aes(x = time, y = y),
                           color = spict_blue_ci_line, linewidth = 0.6) +
      ggplot2::geom_line(  data = transform(ribFF_post, y = ymax), ggplot2::aes(x = time, y = y),
                           color = spict_blue_ci_line, linewidth = 0.6)
  }

  # V-lines
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
  if (!is.null(df_scen_F) && NROW(df_scen_F)) {
    p <- p + ggplot2::geom_line(
      data = df_scen_F, ggplot2::aes(x = time, y = est, color = scenario),
      linewidth = 0.8, na.rm = TRUE
    )
  }

  # Sub-annual Ft (semi-transparent) if available and meaningful
  has_dtc <- !is.null(inp$dtc) && length(inp$dtc)
  if (has_dtc && suppressWarnings(min(inp$dtc, na.rm = TRUE)) < 1) {
    sF <- try(get.par("logFs", rep_man, exp = TRUE, CI = CI), silent = TRUE)
    if (!inherits(sF, "try-error")) {
      tt_sF <- .spict_time_from_par(rep_man, sF)
      mask  <- is.finite(tt_sF) & tt_sF <= (if (isTRUE(manflag) && !is.null(inp$maninterval)) inp$maninterval[1] else .spict_obs_end_overall(rep_man))
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
  if (NROW(dfF_in)) {
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
  if (!manflag && NROW(dfF_pr)) {
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

  # Canonical scenario order & locked colours (prevents fallback)
  scs_final <- if (exists("scs_final")) scs_final else character(0)
  if (length(scs_final)) {
    p <- p + ggplot2::scale_color_manual(
      values  = cols,
      limits  = scs_final,
      drop    = FALSE,
      guide   = if (show_legend) "legend" else "none",
      na.translate = FALSE
    )
  } else {
    p <- p + ggplot2::scale_color_discrete(guide = "none")
  }

  p +
    ggplot2::labs(title = "Absolute fishing mortality", x = "Year", y = expression(F[t])) +
    ggplot2::scale_fill_manual(values = if (is.null(cols)) NA else cols, guide = "none") +
    ggplot2::scale_y_continuous(sec.axis = sec_axis_obj,
                                expand = ggplot2::expansion(mult = c(0.02, 0.02))) +
    theme_minimal_compact2_good_local()
}
