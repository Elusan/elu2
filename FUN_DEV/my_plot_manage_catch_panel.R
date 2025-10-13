#' Reusable colour palette for management scenarios
#'
#' @description Returns a repeatable vector of colours to map across an
#' arbitrary number of management scenarios. Colours recycle to the requested
#' length.
#' @param n Integer. Number of colours required.
#' @return Character vector of hex colour strings.
#' @examples
#' man_cols(5)
#' @keywords internal
#' @noRd
man_cols <- function(n) {
  base <- c('darkmagenta','cyan3','darkgreen','coral1','black',
            'magenta','gold','green','cadetblue3','chocolate3',
            'darkolivegreen3','cyan','darkred')
  rep(base, length.out = n)
}

# ---- tiny safe-field helpers (internal) --------------------------------------
.sf  <- function(x, nm, default = NULL) if (!is.null(x) && nm %in% names(x)) x[[nm]] else default
.sfb <- function(x, nm) isTRUE(.sf(x, nm, FALSE))
.sfi <- function(x, nm, default = integer(0)) .sf(x, nm, default)
.sfn <- function(x, nm, default = numeric(0)) .sf(x, nm, default)

#' Observed indices on the relative B/Bmsy scale
#' @keywords internal
#' @noRd
obs_indices_rel_df <- function(rep_main) {
  stopifnot(inherits(rep_main, "spictcls"))

  qest <- try(elu2::get.par("logq", rep_main, exp = TRUE), silent = TRUE)
  Bmsy <- try(elu2::get.par("logBmsy", rep_main, exp = TRUE), silent = TRUE)
  if (inherits(qest, "try-error") || inherits(Bmsy, "try-error")) return(NULL)

  Bmsy2 <- if (is.null(dim(Bmsy))) {
    as.numeric(Bmsy[2])
  } else {
    as.numeric(utils::tail(Bmsy[, 2], 1))
  }
  if (!is.finite(Bmsy2) || Bmsy2 <= 0) return(NULL)

  inp  <- rep_main$inp
  mapq <- .sfi(inp, "mapq", integer(0))

  .one_index <- function(i) {
    ti <- .sf(inp, "timeI", list())[[i]]
    oi <- .sf(inp, "obsI",  list())[[i]]
    if (is.null(ti) || is.null(oi) || length(oi) != length(ti)) return(NULL)

    qi_idx <- if (length(mapq) >= i) mapq[i] else i
    if (is.null(dim(qest))) {
      qhat <- as.numeric(qest[2])
    } else {
      if (!is.finite(qi_idx) || qi_idx < 1 || qi_idx > nrow(qest)) return(NULL)
      qhat <- as.numeric(qest[qi_idx, 2])
    }
    if (!is.finite(qhat) || qhat <= 0) return(NULL)

    data.frame(index  = paste0("Index", i),
               time   = ti,
               obsrel = (oi / qhat) / Bmsy2,
               row.names = NULL)
  }

  out <- Filter(NROW, list(.one_index(1L), .one_index(2L)))
  if (!length(out)) NULL else do.call(rbind, out)
}

#' Historical B/Bmsy and F/Fmsy up to the last observed time
#' @keywords internal
#' @noRd
get_base_BB_FF_pre <- function(rep_man, CI = 0.95) {
  stopifnot(inherits(rep_man, "spictcls"))

  t_last_obs <- .sfn(rep_man$inp, "timerangeObs", c(NA_real_, NA_real_))[2]
  tt         <- .sfn(rep_man$inp, "time", numeric(0))
  if (!length(tt)) {
    return(list(t_last_obs = NA_real_, BB = data.frame(), FF = data.frame()))
  }

  keep <- which(tt <= t_last_obs)

  BB <- elu2::get.par("logBBmsy", rep_man, exp = TRUE)
  FF <- elu2::get.par("logFFmsy", rep_man, exp = TRUE)

  if (NROW(BB) >= length(tt)) BB <- BB[seq_along(tt), , drop = FALSE]
  if (NROW(FF) >= length(tt)) FF <- FF[seq_along(tt), , drop = FALSE]

  BB <- BB[keep, , drop = FALSE]
  FF <- FF[keep, , drop = FALSE]

  list(
    t_last_obs = t_last_obs,
    BB = data.frame(time = tt[keep], lwr = BB[,1], est = BB[,2], upr = BB[,3], row.names = NULL),
    FF = data.frame(time = tt[keep], lwr = FF[,1], est = FF[,2], upr = FF[,3], row.names = NULL)
  )
}

#' Prepare data frames for management panels (fit or manage)
#' @keywords internal
#' @noRd
prepare_manage_panel_data <- function(rep_man, CI = 0.95) {
  stopifnot(inherits(rep_man, "spictcls"))

  has_man <- ("man" %in% names(rep_man)) && length(rep_man$man)

  co_df <- if (!is.null(rep_man$inp$timeC) && !is.null(rep_man$inp$obsC)) {
    data.frame(time = rep_man$inp$timeC,
               catch = rep_man$inp$obsC,
               catch_type = "Observed",
               row.names = NULL)
  } else NULL

  if (has_man) {
    scenarios <- names(rep_man$man)
    ns <- length(scenarios)

    bb_list <- vector("list", ns)
    ff_list <- vector("list", ns)
    cp_list <- vector("list", ns)
    eval_list <- vector("list", ns)

    for (i in seq_len(ns)) {
      sc <- scenarios[i]
      rp <- rep_man$man[[sc]]

      BB <- elu2::get.par("logBBmsy", rp, exp = TRUE)
      FF <- elu2::get.par("logFFmsy", rp, exp = TRUE)

      rnt <- suppressWarnings(as.numeric(rownames(BB)))
      if (!length(rnt) || any(!is.finite(rnt))) {
        ti <- .sfn(rp$inp, "time", numeric(0))
        rnt <- if (length(ti) >= NROW(BB)) ti[seq_len(NROW(BB))] else seq_len(NROW(BB))
      }

      bb_list[[i]] <- data.frame(time = rnt, lwr = BB[,1], est = BB[,2], upr = BB[,3],
                                 scenario = sc, row.names = NULL)
      ff_list[[i]] <- data.frame(time = rnt, lwr = FF[,1], est = FF[,2], upr = FF[,3],
                                 scenario = sc, row.names = NULL)

      CP  <- elu2::get.par("logCpred", rp, exp = TRUE)
      tcp <- .sfn(rp$inp, "timeCpred", numeric(0))
      if (!length(tcp)) tcp <- seq_len(NROW(CP))
      cp_list[[i]] <- data.frame(time = tcp, lwr = CP[,1], catch = CP[,2], upr = CP[,3],
                                 scenario = sc, catch_type = "Predicted", row.names = NULL)

      ti  <- .sfn(rp$inp, "time", numeric(0))
      tE  <- .sfn(rp$inp, "maneval", NA_real_)
      if (length(ti) && is.finite(tE)) {
        idx <- which.min(abs(ti - tE))
        if (idx >= 1 && idx <= NROW(FF)) {
          eval_list[[i]] <- data.frame(scenario = sc, maneval = tE,
                                       BB_est = BB[idx,2], FF_est = FF[idx,2],
                                       row.names = NULL)
        }
      }
    }

    t0 <- if (!is.null(rep_man$inp$maninterval) && length(rep_man$inp$maninterval) >= 1L)
      rep_man$inp$maninterval[1] else NA_real_

    list(
      bbmsy      = do.call(rbind, bb_list),
      ffmsy      = do.call(rbind, ff_list),
      catch_pred = do.call(rbind, cp_list),
      catch_obs  = co_df,
      eval_pts   = if (length(Filter(NROW, eval_list))) do.call(rbind, eval_list) else NULL,
      obsI_rel   = obs_indices_rel_df(rep_man),
      t0         = t0
    )

  } else {
    list(
      bbmsy      = data.frame(time = numeric(0), lwr = numeric(0), est = numeric(0),
                              upr = numeric(0), scenario = character(0)),
      ffmsy      = data.frame(time = numeric(0), lwr = numeric(0), est = numeric(0),
                              upr = numeric(0), scenario = character(0)),
      catch_pred = data.frame(time = numeric(0), lwr = numeric(0), catch = numeric(0),
                              upr = numeric(0), scenario = character(0), catch_type = character(0)),
      catch_obs  = co_df,
      eval_pts   = NULL,
      obsI_rel   = obs_indices_rel_df(rep_man),
      t0         = NA_real_
    )
  }
}

#' Robust time step for catch series
#' @keywords internal
#' @noRd
.safe_dt_step <- function(rep) {
  if (!is.null(rep$inp$timeC) && length(rep$inp$timeC) > 1L) {
    dt <- stats::median(diff(rep$inp$timeC))
    if (is.finite(dt) && dt > 0) return(dt)
  }
  if (!is.null(rep$inp$dtc)) {
    dt <- suppressWarnings(min(rep$inp$dtc, na.rm = TRUE))
    if (is.finite(dt) && dt > 0) return(dt)
  }
  1
}

#' Management interval vertical lines for B/Bmsy and F/Fmsy panels
#' @keywords internal
#' @noRd
add_management_vlines_BF_good <- function(rep,
                                          color = "grey30",
                                          linetype = "dashed",
                                          linewidth = 0.6,
                                          lineend = "butt") {
  mi <- .sfn(rep$inp, "maninterval", numeric(0))
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

#' Management vertical lines for Catch panel
#' @keywords internal
#' @noRd
add_management_vlines_catch_good <- function(rep,
                                             color = "grey30",
                                             linetype = "dashed",
                                             linewidth = 0.6,
                                             lineend = "butt") {
  mi <- .sfn(rep$inp, "maninterval", numeric(0))
  last_obs_start <- if (!is.null(rep$inp$timeC) && length(rep$inp$timeC)) {
    utils::tail(rep$inp$timeC, 1)
  } else NA_real_

  if (!length(mi) || !is.finite(mi[1])) {
    if (is.finite(last_obs_start)) {
      return(list(ggplot2::geom_vline(
        xintercept = last_obs_start,
        color = color, linetype = linetype,
        linewidth = linewidth, lineend = lineend
      )))
    } else {
      return(list())
    }
  }

  old_left <- mi[1]
  dt_step  <- .safe_dt_step(rep)
  new_left  <- if (is.finite(last_obs_start)) last_obs_start else (old_left - dt_step)
  new_right <- old_left

  list(
    ggplot2::geom_vline(xintercept = new_left,
                        color = color, linetype = linetype,
                        linewidth = linewidth, lineend = lineend),
    ggplot2::geom_vline(xintercept = new_right,
                        color = color, linetype = linetype,
                        linewidth = linewidth, lineend = lineend)
  )
}

#' Catch panel (ggplot2) with MSY band and optional scenario overlays
#'
#' @description
#' Plots observed and predicted catch time series with an MSY reference band.
#' If management scenarios are present in `rep_man$man`, their predicted catch
#' trajectories are overlaid with a stable colour mapping and canonical ordering.
#' For raw fits (no scenarios), base predicted catch is solid up to the last
#' observed catch start and dotted thereafter; CI edges can be shown before the
#' management start (`t0`) when available.
#'
#' @param rep_man A fitted SPiCT report object (`spictcls`), optionally with `$man`.
#' @param scenario_color Optional named character vector of colours for scenarios.
#' @param show_CIs Logical; if `TRUE`, draws dashed CI edges for base Cpred up to `t0`.
#' @param CI Confidence level in (0, 1). Default `0.95`.
#' @param show_legend Logical; show scenario legend when scenarios exist. Default `FALSE`.
#' @return A `ggplot` object.
#' @export
my_plot_manage_catch_panel <- function(rep_man,
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

  dat <- prepare_manage_panel_data(rep_man, CI = CI)

  sc_order <- c("currentCatch","currentF","Fmsy","noF",
                "reduceF25","increaseF25","msyHockeyStick","ices")
  sc_present <- if (("man" %in% names(rep_man)) && length(rep_man$man)) names(rep_man$man) else character(0)
  sc_core   <- intersect(sc_order, sc_present)
  sc_other  <- setdiff(sc_present, sc_core)
  scs_final <- c(sc_core, sort(sc_other))
  has_man   <- length(scs_final) > 0

  if (!is.null(dat$catch_pred) && nrow(dat$catch_pred)) {
    dat$catch_pred$scenario <- factor(dat$catch_pred$scenario, levels = scs_final)
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

  spict_blue_mean    <- "blue"
  spict_blue_ci_line <- grDevices::rgb(0, 0, 1, 0.20)

  # Management start (mask for base Cpred)
  t0 <- dat$t0

  # Base Cpred up to t0 (if available)
  CP0 <- try(spict::get.par("logCpred", rep_man, exp = TRUE), silent = TRUE)
  base_C_t0 <- NULL
  if (!inherits(CP0, "try-error") && !is.null(rep_man$inp$timeCpred)) {
    tc <- rep_man$inp$timeCpred
    maskC <- is.finite(tc)
    if (is.finite(t0)) maskC <- maskC & (tc < t0)
    if (any(maskC)) {
      base_C_t0 <- data.frame(
        time = tc[maskC],
        lwr  = CP0[maskC, 1],
        est  = CP0[maskC, 2],
        upr  = CP0[maskC, 3]
      )
    }
  }

  # MSY series (time-varying vs constant)
  tvgflag <- .sfb(rep_man$inp, "timevaryinggrowth") || .sfb(rep_man$inp, "logmcovflag")
  msy_df <- NULL
  if (isTRUE(tvgflag)) {
    MSYmat <- try(elu2::get.par("logMSYvec", rep_man, exp = TRUE, CI = CI), silent = TRUE)
    if (!inherits(MSYmat, "try-error")) {
      t_msy  <- suppressWarnings(as.numeric(rownames(MSYmat)))
      if (!length(t_msy) || any(!is.finite(t_msy))) t_msy <- .sfn(rep_man$inp, "time", numeric(0))
      msy_df <- data.frame(time = t_msy, lwr = MSYmat[,1], est = MSYmat[,2], upr = MSYmat[,3])
    }
  }
  if (is.null(msy_df)) {
    MSYone <- try(elu2::get.par("logMSY", rep_man, exp = TRUE, CI = CI), silent = TRUE) # c(ll, est, ul)
    t_msy  <- .sfn(rep_man$inp, "time", numeric(0))
    if (!inherits(MSYone, "try-error") && length(t_msy)) {
      msy_df <- data.frame(
        time = t_msy,
        lwr = rep(MSYone[1], length(t_msy)),
        est = rep(MSYone[2], length(t_msy)),
        upr = rep(MSYone[3], length(t_msy))
      )
    } else {
      msy_df <- data.frame(time = numeric(0), lwr = numeric(0), est = numeric(0), upr = numeric(0))
    }
  }
  # Clamp negatives
  msy_df$lwr[msy_df$lwr < 0] <- 0
  msy_df$est[msy_df$est < 0] <- 0
  msy_df$upr[msy_df$upr < 0] <- 0

  darker_color <- function(col, f = 0.6) {
    m <- grDevices::col2rgb(col) / 255
    grDevices::rgb(m[1]*f, m[2]*f, m[3]*f)
  }
  dot_blue <- darker_color(spict_blue_mean, 0.6)

  # Last observed catch start (for styling & vline handling in raw fits)
  obs_end <- if (!is.null(rep_man$inp$timeC) && length(rep_man$inp$timeC))
    utils::tail(rep_man$inp$timeC, 1) else NA_real_

  # Y label with units if available
  cu <- .sf(rep_man$inp, "catchunit", "")
  ylab <- if (nzchar(as.character(cu))) paste0("Catch (", cu, ")") else "Catch"

  pC <- ggplot2::ggplot() +
    # MSY ribbon + mean at the back
    { if (nrow(msy_df)) ggplot2::geom_ribbon(
      data = msy_df,
      ggplot2::aes(x = time, ymin = lwr, ymax = upr),
      inherit.aes = FALSE, fill = "lightgray"
    ) } +
    { if (nrow(msy_df)) ggplot2::geom_line(
      data = msy_df,
      ggplot2::aes(x = time, y = est),
      inherit.aes = FALSE, color = "black", linewidth = 0.8
    ) } +
    # Vlines: with scenarios (two), without scenarios (maybe one)
    { if (has_man)
      add_management_vlines_catch_good(rep_man, color = "grey75", linetype = "solid",
                                       linewidth = 0.4, lineend = "butt")
      else
        add_management_vlines_catch_good(rep_man, color = "grey75", linetype = "solid",
                                         linewidth = 0.4, lineend = "butt")[1]
    } +
    # Scenario predicted catch (no-op if none)
    ggplot2::geom_line(
      data = dat$catch_pred,
      ggplot2::aes(x = time, y = catch, color = scenario),
      linewidth = 0.6, na.rm = TRUE
    ) +
    # Observed catch
    { if (!is.null(dat$catch_obs) && nrow(dat$catch_obs)) ggplot2::geom_point(
      data = dat$catch_obs,
      ggplot2::aes(x = time, y = catch),
      color = dot_blue, shape = 16, size = 3, na.rm = TRUE
    ) } +
    # Base predicted catch (≤ t0 if present): CI edges (dashed)
    { if (show_CIs && !is.null(base_C_t0) && nrow(base_C_t0)) list(
      ggplot2::geom_line(
        data = base_C_t0, ggplot2::aes(x = time, y = lwr),
        color = spict_blue_ci_line, linetype = "dashed", linewidth = 0.7, alpha = 0.95, inherit.aes = FALSE
      ),
      ggplot2::geom_line(
        data = base_C_t0, ggplot2::aes(x = time, y = upr),
        color = spict_blue_ci_line, linetype = "dashed", linewidth = 0.7, alpha = 0.95, inherit.aes = FALSE
      )
    ) } +
    # Mean line: managed = solid; raw fit = solid <= obs_end, dotted >= obs_end
    { if (!is.null(base_C_t0) && nrow(base_C_t0)) {
      if (has_man || !is.finite(obs_end)) {
        ggplot2::geom_line(
          data = base_C_t0, ggplot2::aes(x = time, y = est),
          color = spict_blue_mean, linewidth = 0.6, inherit.aes = FALSE
        )
      } else {
        list(
          ggplot2::geom_line(
            data = subset(base_C_t0, time <= obs_end),
            ggplot2::aes(x = time, y = est),
            color = spict_blue_mean, linewidth = 0.6, inherit.aes = FALSE
          ),
          ggplot2::geom_line(
            data = subset(base_C_t0, time >= obs_end),
            ggplot2::aes(x = time, y = est),
            color = spict_blue_mean, linewidth = 0.6, inherit.aes = FALSE,
            linetype = "dotted"
          )
        )
      }
    } } +
    ggplot2::labs(title = "Catch", x = "Year", y = ylab) +
    { if (length(scs_final)) ggplot2::scale_color_manual(values = cols, guide = if (show_legend) "legend" else "none")
      else ggplot2::scale_color_discrete(guide = "none") } +
    ggplot2::scale_fill_manual(values = if (is.null(cols)) NA else cols, guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.04))) +
    theme_minimal_compact2_good_local()

  pC
}
