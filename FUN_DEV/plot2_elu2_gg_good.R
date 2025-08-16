#' Multi-panel SPiCT management plot with Kobe panel (ggplot2)
#'
#' @description
#' Produces a 2×2 layout of management diagnostics for a fitted **SPiCT** object
#' with scenarios generated via `spict::manage()`:
#' **B/Bmsy**, **Catch**, **F/Fmsy**, and a **Kobe** plot.
#'
#' - Time-series panels show a thin **solid grey** vertical line at the start
#'   and end of the management interval.
#' - In the **Catch** panel, the **left** vertical line is aligned with the
#'   **last observed catch start** (`timeC` tail), while the **right** vertical
#'   line is at the management start.
#' - Historical base (≤ last observation) B/Bmsy and F/Fmsy appear in SPiCT-blue
#'   with CIs; scenario lines are overlaid in user-supplied or default colors.
#' - The Kobe panel uses **absolute** primary axes (Bt, Ft) with **relative**
#'   secondary axes (Bt/Bmsy, Ft/Fmsy). When growth is time-varying, axes switch
#'   to relative.
#'
#' @param rep_man A fitted `spictcls` object for which `manage(rep)` has been run.
#' @param scenario_color Optional named character vector of colors for scenarios.
#'   If `NULL`, a repeatable default palette is used. Names must match
#'   `names(rep_man$man)`. Unnamed vectors will be matched in order.
#' @param show_CIs Logical; show SPiCT-style CIs for the historical base series
#'   and for base predicted catch up to the management start. Default `TRUE`.
#' @param CI Numeric in (0,1); confidence level for intervals. Default `0.95`.
#' @param return_patchwork Logical; if `TRUE` (default), returns a 2×2 patchwork
#'   layout. If `FALSE`, returns a named list of ggplot objects:
#'   `list(bbmsy, catch, ffmsy, kobe)`.
#'
#' @return Either a patchwork object (default) or a named list of ggplot objects.
#'
#' @section Notes:
#' - Observed catches are plotted at `timeC` (start of interval); predicted catches
#'   at `timeCpred` (SPiCT convention).
#' - Vertical lines in time-series panels are **solid grey** (thin).
#' - In the Catch panel, the left line is at the **last observed** catch start
#'   (tail of `timeC`), and the right line is at the **management start**.
#'
#' @examples
#' \dontrun{
#'   library(spict)
#'   # fit <- fit.spict(inp)             # your input list
#'   # rep <- fit                         # spictcls
#'   # rep <- manage(rep)                 # create scenarios
#'   # p   <- plot2_elu2_gg_good(rep)
#'   # p
#' }
#'
#' @seealso [spict::manage()]
#'
#' @import ggplot2
#' @import patchwork
#' @import grid
#' @importFrom stats median qnorm
#' @importFrom utils tail
#' @importFrom grDevices rgb adjustcolor
#' @export
plot2_elu2_gg_good <- function(rep_man,
                                  scenario_color = NULL,
                                  show_CIs = TRUE,
                                  CI = 0.95,
                                  return_patchwork = TRUE) {
  stopifnot(inherits(rep_man, "spictcls"))
  if (!"man" %in% names(rep_man)) stop("Run manage(rep) before plot2_elu2_gg_good().")

  # Prepared data (scenarios + observed/predicted catch + indices)
  dat <- prepare_manage_panel_data(rep_man, CI = CI)
  scs <- unique(dat$bbmsy$scenario)

  # ---- Scenario color mapping (shared across panels & Kobe) ----
  if (is.null(scenario_color)) {
    cols <- man_cols(length(scs)); names(cols) <- scs
  } else {
    cols <- scenario_color
    if (is.null(names(cols))) names(cols) <- scs
    cols <- cols[match(scs, names(cols))]; names(cols) <- scs
  }

  # ---- SPiCT-style fixed blues (base historical series) ----
  spict_blue_mean    <- "blue"
  spict_blue_ci_line <- grDevices::rgb(0, 0, 1, 0.20)
  spict_blue_ci_fill <- grDevices::rgb(0, 0, 1, 0.10)

  darker_color <- function(col, f = 0.8) {
    m <- grDevices::col2rgb(col) / 255
    grDevices::rgb(m[1]*f, m[2]*f, m[3]*f)
  }
  dot_blue <- darker_color(spict_blue_mean, 0.6)

  # Historical base (≤ last observation) BB & FF (SPiCT-style)
  base_pre <- get_base_BB_FF_pre(rep_man, CI = CI)

  # Management start (used to truncate base catch predictions for CI display)
  t0 <- dat$t0

  # Base (rep) Cpred up to t0 (for blue overlay in Catch panel)
  CP0 <- try(get.par("logCpred", rep_man, exp = TRUE), silent = TRUE)
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

  # MSY series for Catch panel (time-varying vs constant)
  tvgflag <- isTRUE(rep_man$inp$timevaryinggrowth) || isTRUE(rep_man$inp$logmcovflag)
  if (tvgflag) {
    MSYmat <- get.par("logMSYvec", rep_man, exp = TRUE, CI = CI)
    t_msy  <- suppressWarnings(as.numeric(rownames(MSYmat)))
    if (!length(t_msy) || any(!is.finite(t_msy))) t_msy <- rep_man$inp$time
    msy_df <- data.frame(time = t_msy, lwr = MSYmat[, 1], est = MSYmat[, 2], upr = MSYmat[, 3])
  } else {
    MSYone <- get.par("logMSY", rep_man, exp = TRUE, CI = CI) # c(ll, est, ul)
    t_msy  <- rep_man$inp$time
    msy_df <- data.frame(time = t_msy,
                         lwr = rep(MSYone[1], length(t_msy)),
                         est = rep(MSYone[2], length(t_msy)),
                         upr = rep(MSYone[3], length(t_msy)))
  }
  msy_df$lwr[msy_df$lwr < 0] <- 0
  msy_df$est[msy_df$est < 0] <- 0
  msy_df$upr[msy_df$upr < 0] <- 0

  # ---- Compact no-frills theme (panel frames bold, legend inside) ----
  theme_minimal_compact2 <- function(base_size = 10, base_family = "") {
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
        legend.text  = ggplot2::element_text(size = 14),
        legend.key.size = grid::unit(0.8, "lines"),
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

  ## --------------------
  ## Panel: B / Bmsy
  ## --------------------
  pB <- ggplot2::ggplot() +
    { if (show_CIs && nrow(base_pre$BB)) ggplot2::geom_ribbon(
      data = base_pre$BB,
      ggplot2::aes(x = time, ymin = lwr, ymax = upr),
      inherit.aes = FALSE, fill = spict_blue_ci_fill, color = NA
    ) } +
    { if (show_CIs && nrow(base_pre$BB)) list(
      ggplot2::geom_line(data = base_pre$BB, ggplot2::aes(x = time, y = lwr),
                         inherit.aes = FALSE, color = spict_blue_ci_line, linewidth = 0.6),
      ggplot2::geom_line(data = base_pre$BB, ggplot2::aes(x = time, y = upr),
                         inherit.aes = FALSE, color = spict_blue_ci_line, linewidth = 0.6)
    ) } +

    { if (nrow(base_pre$BB)) ggplot2::geom_line(
      data = base_pre$BB,
      ggplot2::aes(x = time, y = est),
      inherit.aes = FALSE,
      color = spict_blue_mean,   # <- fixed blue
      linewidth = 0.8
    ) } +
    ggplot2::geom_line(
      data = dat$bbmsy,
      ggplot2::aes(x = time, y = est, color = scenario),
      linewidth = 0.8, na.rm = TRUE
    ) +
    { if (!is.null(dat$obsI_rel)) {
      list(
        ggplot2::geom_point(
          data = dat$obsI_rel[dat$obsI_rel$index == "Index1", , drop = FALSE],
          ggplot2::aes(x = time, y = obsrel), inherit.aes = FALSE,
          shape = 16, size = 3, color = dot_blue
        ),
        ggplot2::geom_point(
          data = dat$obsI_rel[dat$obsI_rel$index == "Index2", , drop = FALSE],
          ggplot2::aes(x = time, y = obsrel), inherit.aes = FALSE,
          shape = 22, size = 3, fill = "green", color = "black", stroke = 0.8
        )
      )
    } } +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
    add_management_vlines_BF(rep_man, color = "grey75", linetype = "solid", linewidth = 0.4) +
    ggplot2::labs(x = "Year", y = expression(bold(B/B[MSY]))) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::scale_fill_manual(values = cols, guide = "none") +
    ggplot2::geom_line(
      data = base_pre$BB,
      mapping = ggplot2::aes(x = time, y = est, group = 1),
      inherit.aes = FALSE,
      colour = I(spict_blue_mean),   # force fixed blue (immune to scales)
      linewidth = 0.9,
      show.legend = FALSE
    ) +

    theme_minimal_compact2()

  ## --------------------
  ## Panel: F / Fmsy
  ## --------------------
  pF <- ggplot2::ggplot() +
    { if (show_CIs && nrow(base_pre$FF)) ggplot2::geom_ribbon(
      data = base_pre$FF,
      ggplot2::aes(x = time, ymin = lwr, ymax = upr),
      inherit.aes = FALSE, fill = spict_blue_ci_fill, color = NA
    ) } +
    { if (show_CIs && nrow(base_pre$FF)) list(
      ggplot2::geom_line(data = base_pre$FF, ggplot2::aes(x = time, y = lwr),
                         inherit.aes = FALSE, color = spict_blue_ci_line, linewidth = 0.7),
      ggplot2::geom_line(data = base_pre$FF, ggplot2::aes(x = time, y = upr),
                         inherit.aes = FALSE, color = spict_blue_ci_line, linewidth = 0.7)
    ) } +
    { if (nrow(base_pre$FF)) ggplot2::geom_line(
      data = base_pre$FF, ggplot2::aes(x = time, y = est),
      inherit.aes = FALSE,
      color = spict_blue_mean,
      linewidth = 0.8
    ) } +
    ggplot2::geom_line(
      data = dat$ffmsy,
      ggplot2::aes(x = time, y = est, color = scenario),
      linewidth = 0.8, na.rm = TRUE
    ) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
    add_management_vlines_BF(rep_man, color = "grey75", linetype = "solid", linewidth = 0.4) +
    ggplot2::labs(x = "Year", y = expression(bold(F/F[MSY]))) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::scale_fill_manual(values = cols, guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.02))) +
    ggplot2::guides(color = "none") +
    ggplot2::geom_line(
      data = base_pre$FF,
      mapping = ggplot2::aes(x = time, y = est, group = 1),
      inherit.aes = FALSE,
      colour = I(spict_blue_mean),   # force fixed blue (immune to scales)
      linewidth = 0.9,
      show.legend = FALSE
    ) +

    theme_minimal_compact2()

  ## --------------------
  ## Panel: Catch
  ## --------------------
  pC <- ggplot2::ggplot() +
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
    ggplot2::geom_line(
      data = dat$catch_pred,
      ggplot2::aes(x = time, y = catch, color = scenario),
      linewidth = 0.8, na.rm = TRUE
    ) +
    ggplot2::geom_point(
      data = dat$catch_obs,
      ggplot2::aes(x = time, y = catch),
      color = dot_blue, shape = 16, size = 3, na.rm = TRUE
    ) +
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
    { if (!is.null(base_C_t0) && nrow(base_C_t0)) ggplot2::geom_line(
      data = base_C_t0, ggplot2::aes(x = time, y = est),
      color = spict_blue_mean, linewidth = 0.8, inherit.aes = FALSE
    ) } +
    add_management_vlines_catch(rep_man, color = "grey75", linetype = "solid", linewidth = 0.4) +
    ggplot2::labs(x = "Year", y = "Catch (tons)") +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::scale_fill_manual(values = cols, guide = "none") +
    ggplot2::guides(color = "none") +
    theme_minimal_compact2()

  ## --------------------
  ## Panel: Kobe (absolute primary axes + relative secondary axes)
  ## --------------------
  pK <- kobe_all_in_one_gg(
    rep                 = rep_man,
    logax               = FALSE,
    plot.legend         = TRUE,
    man.legend          = FALSE,   # hide scenario legend here
    ext                 = TRUE,
    rel.axes            = FALSE,
    CI                  = CI,
    scenario_colors     = cols,
    origin_offset_frac  = 0.04
  )

  if (isTRUE(return_patchwork)) {
    ((pB | pC) / (pF | pK)) + patchwork::plot_layout(heights = c(1, 1))
  } else {
    list(bbmsy = pB, catch = pC, ffmsy = pF, kobe = pK)
  }
}

# ----- helpers (internal) -----------------------------------------------------

#' Repeatable scenario palette
#' @keywords internal
#' @noRd
man_cols <- function(n) {
  base <- c('darkmagenta','cyan3','darkgreen','coral1','black',
            'magenta','gold','green','cadetblue3','chocolate3',
            'darkolivegreen3','cyan','darkred')
  rep(base, length.out = n)
}

#' Observed indices on B/Bmsy scale (Index1/Index2), SPiCT-style
#' @keywords internal
#' @noRd
obs_indices_rel_df <- function(rep_main) {
  qest  <- get.par("logq", rep_main, exp = TRUE)
  Bmsy  <- get.par("logBmsy", rep_main, exp = TRUE)
  Bmsy2 <- if (is.null(nrow(Bmsy))) Bmsy[2] else Bmsy[1, 2]
  inp   <- rep_main$inp

  out_list <- list()
  k <- 0L

  if (!is.null(inp$timeI) && length(inp$timeI) >= 1L &&
      !is.null(inp$obsI)  && length(inp$obsI)  >= 1L) {
    k <- k + 1L
    out_list[[k]] <- data.frame(
      index  = "Index1",
      time   = inp$timeI[[1]],
      obsrel = (inp$obsI[[1]] / qest[inp$mapq[1], 2]) / Bmsy2
    )
  }
  if (!is.null(inp$timeI) && length(inp$timeI) >= 2L &&
      !is.null(inp$obsI)  && length(inp$obsI)  >= 2L) {
    k <- k + 1L
    out_list[[k]] <- data.frame(
      index  = "Index2",
      time   = inp$timeI[[2]],
      obsrel = (inp$obsI[[2]] / qest[inp$mapq[2], 2]) / Bmsy2
    )
  }

  if (length(out_list) == 0L) NULL else do.call(rbind, out_list)
}

#' Historical (≤ last obs) B/Bmsy and F/Fmsy in SPiCT style
#' @keywords internal
#' @noRd
get_base_BB_FF_pre <- function(rep_man, CI = 0.95) {
  stopifnot(inherits(rep_man, "spictcls"))
  t_last_obs <- rep_man$inp$timerangeObs[2]
  tt <- rep_man$inp$time
  keep <- which(tt <= t_last_obs)

  BB <- get.par("logBBmsy", rep_man, exp = TRUE)
  FF <- get.par("logFFmsy", rep_man, exp = TRUE)

  BB <- BB[keep, , drop = FALSE]
  FF <- FF[keep, , drop = FALSE]

  list(
    t_last_obs = t_last_obs,
    BB = data.frame(
      time = tt[keep], lwr = BB[, 1], est = BB[, 2], upr = BB[, 3],
      row.names = NULL
    ),
    FF = data.frame(
      time = tt[keep], lwr = FF[, 1], est = FF[, 2], upr = FF[, 3],
      row.names = NULL
    )
  )
}

#' Construct tidy scenario data for B/Bmsy, F/Fmsy, Catch panels
#' @keywords internal
#' @noRd
prepare_manage_panel_data <- function(rep_man, CI = 0.95) {
  stopifnot(inherits(rep_man, "spictcls"))
  if (!"man" %in% names(rep_man)) stop("Run manage(rep) before calling.")

  scenarios <- names(rep_man$man)
  ns        <- length(scenarios)

  bb_list <- vector("list", ns)
  ff_list <- vector("list", ns)
  cp_list <- vector("list", ns)
  eval_list <- vector("list", ns)

  # Observed catch (shared) — plotted at timeC (start-of-interval)
  co_df <- data.frame(
    time       = rep_man$inp$timeC,
    catch      = rep_man$inp$obsC,
    catch_type = "Observed"
  )

  for (i in seq_len(ns)) {
    sc <- scenarios[i]
    rp <- rep_man$man[[sc]]

    BB <- get.par("logBBmsy", rp, exp = TRUE)
    FF <- get.par("logFFmsy", rp, exp = TRUE)

    bb_list[[i]] <- data.frame(
      time = as.numeric(rownames(BB)),
      lwr  = BB[, 1],
      est  = BB[, 2],
      upr  = BB[, 3],
      scenario = sc
    )
    ff_list[[i]] <- data.frame(
      time = as.numeric(rownames(FF)),
      lwr  = FF[, 1],
      est  = FF[, 2],
      upr  = FF[, 3],
      scenario = sc
    )

    CP <- get.par("logCpred", rp, exp = TRUE)
    cp_list[[i]] <- data.frame(
      time       = rp$inp$timeCpred,
      lwr        = CP[, 1],
      catch      = CP[, 2],
      upr        = CP[, 3],
      scenario   = sc,
      catch_type = "Predicted"
    )

    ti  <- rp$inp$time
    tE  <- rp$inp$maneval
    idx <- which(ti == tE)
    if (length(idx) == 1L) {
      eval_list[[i]] <- data.frame(
        scenario = sc, maneval = tE,
        BB_est   = BB[idx, 2], FF_est = FF[idx, 2]
      )
    }
  }

  # management start (for CI masking / alignment)
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
}

#' Robust step size (annual or sub-annual)
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

#' Add management interval vertical lines to B/Bmsy and F/Fmsy panels
#' @keywords internal
#' @noRd
add_management_vlines_BF <- function(rep,
                                     color = "grey30",
                                     linetype = "dashed",
                                     linewidth = 0.6) {
  mi <- rep$inp$maninterval
  if (is.null(mi) || length(mi) < 1L || !is.finite(mi[1])) return(list())

  left  <- mi[1]
  right <- if (length(mi) >= 2L && is.finite(mi[2])) mi[2] else NA_real_

  out <- list(ggplot2::geom_vline(xintercept = left,
                                  color = color, linetype = linetype, linewidth = linewidth))
  if (is.finite(right)) {
    out[[length(out) + 1L]] <- ggplot2::geom_vline(
      xintercept = right, color = color, linetype = linetype, linewidth = linewidth
    )
  }
  out
}

#' Add management interval vertical lines to Catch panel
#'
#' Left line = last observed catch **start** (tail of `timeC`);
#' Right line = original management start.
#' @keywords internal
#' @noRd
add_management_vlines_catch <- function(rep,
                                        color = "grey30",
                                        linetype = "dashed",
                                        linewidth = 0.6) {
  mi <- rep$inp$maninterval
  if (is.null(mi) || length(mi) < 1L || !is.finite(mi[1])) return(list())

  old_left <- mi[1]
  dt_step  <- .safe_dt_step(rep)

  last_obs_start <- if (!is.null(rep$inp$timeC) && length(rep$inp$timeC)) {
    utils::tail(rep$inp$timeC, 1)
  } else NA_real_

  new_left  <- if (is.finite(last_obs_start)) last_obs_start else (old_left - dt_step)
  new_right <- old_left

  list(
    ggplot2::geom_vline(xintercept = new_left,
                        color = color, linetype = linetype, linewidth = linewidth),
    ggplot2::geom_vline(xintercept = new_right,
                        color = color, linetype = linetype, linewidth = linewidth)
  )
}

#' Extract catch series for a single scenario (SPiCT convention)
#' @keywords internal
#' @noRd
elu_extract_catch_data_spict <- function(model, scenario_name = "Scenario") {
  stopifnot(inherits(model, "spictcls"))
  CP <- get.par("logCpred", model, exp = TRUE)

  pred <- data.frame(
    time       = model$inp$timeCpred,
    lwr        = CP[, 1],
    catch      = CP[, 2],
    upr        = CP[, 3],
    catch_type = "Predicted",
    scenario   = scenario_name,
    stringsAsFactors = FALSE
  )

  obs <- data.frame(
    time       = model$inp$timeC,
    catch      = model$inp$obsC,
    lwr        = NA_real_,
    upr        = NA_real_,
    catch_type = "Observed",
    scenario   = scenario_name,
    stringsAsFactors = FALSE
  )

  rbind(pred, obs)
}

#' Kobe plot (all-in-one, absolute primary axes + relative secondary axes)
#' @keywords internal
#' @noRd
kobe_all_in_one_gg <- function(rep,
                               logax = FALSE,
                               plot.legend = TRUE,
                               man.legend = TRUE,
                               ext = TRUE,
                               rel.axes = FALSE,
                               xlim = NULL,
                               ylim = NULL,
                               labpos = c(1, 1),
                               xlabel = NULL,
                               stamp = NULL,
                               verbose = TRUE,
                               CI = 0.95,
                               scenario_colors = NULL,
                               origin_offset_frac = 0.04) {

  # ----- small helpers (local) -----
  check_rep <- function(rep, reportmode0 = TRUE) {
    if (!inherits(rep, "spictcls") || !"opt" %in% names(rep)) stop("The argument 'rep' must be a fitted spict object (fit.spict()).")
    if (reportmode0 && rep$inp$reportmode != 0) stop("All states must be reported! Set 'inp$reportmode <- 0' and refit.")
    TRUE
  }
  get_par_loc <- function(parname, rep, exp = FALSE, CI = 0.95) {
    z <- stats::qnorm(CI + (1 - CI) / 2)
    ind_ran <- which(names(rep$par.random) == parname)
    ind_fix <- which(names(rep$par.fixed)  == parname)
    ind_sdr <- which(names(rep$value)      == parname)
    ind_opt <- which(names(rep$opt$par)    == parname)
    est <- ll <- ul <- sdv <- NULL
    if (length(ind_ran)) { est <- rep$par.random[ind_ran]; sdv <- sqrt(rep$diag.cov.random[ind_ran]); ll <- est - z*sdv; ul <- est + z*sdv }
    if (length(ind_fix)) { est <- rep$par.fixed[ind_fix];  sdv <- sqrt(diag(rep$cov.fixed))[ind_fix]; ll <- est - z*sdv; ul <- est + z*sdv }
    if (length(ind_sdr)) { est <- rep$value[ind_sdr];       sdv <- rep$sd[ind_sdr];                    ll <- est - z*sdv; ul <- est + z*sdv }
    if (length(est) == 0) {
      if (length(ind_opt)) { est <- rep$opt$par[ind_opt]; sdv <- rep(NA_real_, length(est)); ll <- ul <- rep(NA_real_, length(est))
      } else { est <- NA_real_; sdv <- ll <- ul <- NA_real_ }
    }
    n <- length(est); if (length(sdv)==0) sdv <- rep(NA_real_,n); if (length(ll)==0) ll <- rep(NA_real_,n); if (length(ul)==0) ul <- rep(NA_real_,n)
    if (isTRUE(exp)) { ll <- exp(ll); ul <- exp(ul); est <- exp(est); ul[is.infinite(ul)] <- exp(705) }
    cbind(ll, est, ul, sdv, if (isTRUE(exp)) NA_real_ else sdv/est)
  }
  add_catchunit_loc <- function(lab, cu) { cu <- as.character(cu); if (nzchar(cu)) eval(bquote(.(lab[[1]]) * ',' ~ .(cu))) else lab }
  annual_avg <- function(intime, vec, type = "mean") {
    fun <- match.fun(type)
    anntime <- unique(floor(intime))
    floortime <- floor(intime)
    nstepvec <- vapply(anntime, function(a) sum(a == floortime), numeric(1))
    anntime <- anntime[which(nstepvec == max(nstepvec))]
    annvec  <- vapply(anntime, function(a) fun(vec[which(a == floortime)]), numeric(1))
    list(anntime = anntime, annvec = annvec)
  }
  calc_EBinf <- function(K, n, Fl, Fmsy, sdb2) {
    base <- 1 - (n - 1) / n * (Fl / Fmsy); base <- max(0, base)
    corr <- 1 - n / 2 / (1 - (1 - n * Fmsy + (n - 1) * Fl))
    max(0, K * (base)^(1/(n-1)) * (1 - corr * sdb2))
  }
  get_EBinf_loc <- function(rep) {
    K <- get_par_loc("logK", rep, exp=TRUE)[2]
    n <- get_par_loc("logn",rep,exp=TRUE)[2]
    sdb2 <- get_par_loc("logsdb",rep,exp=TRUE)[2]^2
    Fmsy <- utils::tail(get_par_loc("logFmsy",rep,exp=TRUE),1)[2]
    logFs <- get_par_loc("logFs",rep)
    fff <- if (min(rep$inp$dtc) < 1) exp(annual_avg(rep$inp$time, logFs[,2])$annvec) else exp(logFs[,2])
    calc_EBinf(K, n, utils::tail(fff,1), Fmsy, sdb2)
  }
  build_msy_ellipse_abs <- function(rep) {
    tvgflag <- rep$inp$timevaryinggrowth | rep$inp$logmcovflag
    ok <- rep$opt$convergence == 0 && !tvgflag && requireNamespace("ellipse", quietly = TRUE)
    if (ok) {
      idxB <- utils::tail(which(names(rep$value) == "logBmsy"), 1)
      idxF <- utils::tail(which(names(rep$value) == "logFmsy"), 1)
      idxS <- utils::tail(which(names(rep$value) == "logBmsyPluslogFmsy"), 1)
      if (!length(idxB) || !length(idxF)) ok <- FALSE
    }
    if (ok) {
      muB <- rep$value[idxB]; sdB <- rep$sd[idxB]
      muF <- rep$value[idxF]; sdF <- rep$sd[idxF]
      if (length(idxS) && is.finite(rep$sd[idxS])) {
        varSum <- rep$sd[idxS]^2
        covBF  <- (varSum - sdB^2 - sdF^2) / 2
      } else covBF <- 0
      rho <- covBF/(sdB*sdF); rho <- pmin(pmax(rho, -0.99), 0.99)
      Elog <- ellipse::ellipse(rho, scale = c(sdB, sdF), centre = c(muB, muF), npoints = 300)
      data.frame(x = exp(Elog[,1]), y = exp(Elog[,2]))
    } else {
      data.frame(x = get_par_loc("logBmsy", rep, exp=TRUE)[2],
                 y = get_par_loc("logFmsy", rep, exp=TRUE)[2])
    }
  }

  # ----- derive core quantities -----
  check_rep(rep)
  inp <- rep$inp
  tvgflag <- isTRUE(inp$timevaryinggrowth) | isTRUE(inp$logmcovflag)
  if (tvgflag) rel.axes <- TRUE

  Bmsyall <- get_par_loc("logBmsy", rep, exp=TRUE, CI = CI)
  Fmsyall <- get_par_loc("logFmsy", rep, exp=TRUE, CI = CI)
  Bmsy <- utils::tail(Bmsyall, 1); Fmsy <- utils::tail(Fmsyall, 1)

  if (rel.axes) {
    ext <- FALSE
    bscal <- Bmsy[2]; fscal <- 1
    xlab_expr <- expression(B[t]/B[MSY])
    ylab_expr <- expression(F[t]/F[MSY])
  } else {
    bscal <- 1; fscal <- 1
    xlab_expr <- expression(B[t]); xlab_expr <- add_catchunit_loc(xlab_expr, inp$catchunit)
    ylab_expr <- expression(F[t])
  }

  Bp      <- get_par_loc("logBp",   rep, exp=TRUE, CI = CI)
  Best    <- get_par_loc("logB",    rep, exp=TRUE, CI = CI)
  logBest <- get_par_loc("logB",    rep, CI = CI)
  if (tvgflag) {
    Fest  <- get_par_loc("logFFmsy", rep, exp=TRUE, CI = CI)
    fscal <- 1; Fmsy <- c(1, 1)
  } else {
    Fest  <- get_par_loc("logFs", rep, exp=TRUE, CI = CI)
  }
  logFest <- get_par_loc("logFs", rep, CI = CI)
  ns <- nrow(Best)

  ell_abs <- build_msy_ellipse_abs(rep)

  if (min(inp$dtc) < 1) {
    alb <- annual_avg(inp$time, logBest[, 2])
    alf <- annual_avg(inp$time, logFest[, 2])
    bbb <- exp(alb$annvec)/bscal
    fff <- exp(alf$annvec)/fscal
    fbtime <- alb$anntime
  } else {
    fff <- Fest[inp$indest, 2]/fscal
    bbb <- Best[inp$indest, 2]/bscal
    fbtime <- inp$time[inp$indest]
  }

  Fl    <- utils::tail(unname(fff), 1)
  EBinf <- get_EBinf_loc(rep)/bscal

  # ---- limits + small padding ----
  if (is.null(xlim)) {
    if (min(inp$dtc) < 1) {
      xlim <- range(c(exp(alb$annvec), ell_abs$x, EBinf)/bscal, na.rm=TRUE)
    } else {
      xlim <- range(c(ell_abs$x, Best[,2], EBinf)/bscal, na.rm=TRUE)
    }
    xlim[2] <- min(c(xlim[2], 8*Bmsy[2]/bscal), 2.2*max(bbb), na.rm=TRUE)
    xlim[2] <- max(c(xlim[2], Bmsy[2]/bscal), na.rm=TRUE)
  }
  if (is.null(ylim)) {
    if (min(inp$dtc) < 1) {
      ylim <- range(c(exp(alf$annvec)/fscal, ell_abs$y/fscal), na.rm=TRUE)
    } else {
      ylim <- range(c(ell_abs$y/fscal, Fest[,2]/fscal), na.rm=TRUE)
    }
    ylim[2] <- min(c(ylim[2], 8*Fmsy[2]/fscal), 2.2*max(fff), na.rm=TRUE)
    ylim[2] <- max(c(ylim[2], Fmsy[2]/fscal), na.rm=TRUE)
    if ("man" %in% names(rep)) ylim <- range(ylim, 0)
  }

  logminval <- 1e-4
  if (isTRUE(logax)) {
    if (xlim[1] < logminval) xlim[1] <- logminval
    if (ylim[1] < logminval) ylim[1] <- logminval
  }

  pad_frac <- 0.02
  xpad <- pad_frac * diff(xlim)
  ypad <- pad_frac * diff(ylim)
  xlim <- c(if (isTRUE(logax)) max(xlim[1], logminval) else xlim[1] - xpad, xlim[2] + xpad)
  ylim <- c(ylim[1] - ypad, ylim[2] + ypad)
  if (!isTRUE(logax)) if (ylim[1] >= 0) ylim[1] <- -max(1e-6, 0.02 * (ylim[2] - 0))

  # ---- quadrant rectangles ----
  green_col  <- grDevices::rgb(0.5, 0.8, 0.4, 1)
  yellow_col <- grDevices::rgb(1, 0.925, 0.55, 1)
  red_col    <- grDevices::rgb(1, 0.188, 0.188, 1)

  bxref <- utils::tail(Bmsy, 1)[2]
  fyref <- utils::tail(Fmsy, 1)[2]

  df_green_rect   <- data.frame(xmin =  bxref, xmax =  Inf, ymin = -Inf, ymax =  fyref)
  df_yellowL_rect <- data.frame(xmin = -Inf,  xmax =  bxref, ymin = -Inf, ymax =  fyref)
  df_yellowT_rect <- data.frame(xmin =  bxref, xmax =  Inf, ymin =  fyref, ymax =  Inf)
  df_red_rect     <- data.frame(xmin = -Inf,  xmax =  bxref, ymin =  fyref, ymax =  Inf)

  # ---- path data ----
  df_traj <- data.frame(x = bbb, y = fff, ord = seq_along(bbb))
  df_pred <- df_EBseg <- NULL
  if (!(min(inp$dtc) < 1) && !("man" %in% names(rep))) {
    xpred <- get_par_loc("logB",  rep, exp=TRUE)[inp$indpred,2] / bscal
    ypred <- get_par_loc("logFs", rep, exp=TRUE)[inp$indpred,2] / fscal
    df_pred <- data.frame(x = xpred, y = ypred)
    df_EBseg <- data.frame(x = utils::tail(xpred,1), xend = EBinf, y = utils::tail(ypred,1), yend = utils::tail(ypred,1))
  }

  df_first <- data.frame(x = df_traj$x[1], y = df_traj$y[1], lab = round(fbtime[1], 2))
  df_last  <- data.frame(x = utils::tail(df_traj$x,1), y = utils::tail(df_traj$y,1), lab = round(utils::tail(fbtime,1), 2))

  # Scenarios (management path, intermediate pre-management path, and vertical jump)
  df_man <- df_man_int <- df_man_jump <- NULL
  if ("man" %in% names(rep) && length(rep$man)) {
    nman <- length(rep$man)
    leg_man <- names(rep$man); if (is.null(leg_man) || !length(leg_man)) leg_man <- paste0("Scenario ", seq_len(nman))
    lst_path <- vector("list", nman); lst_ipath <- vector("list", nman); lst_jump <- vector("list", nman)
    for (i in seq_len(nman)) {
      rp <- rep$man[[i]]
      estB <- get_par_loc("logB",  rp, exp=TRUE)[,2]
      estF <- get_par_loc("logFs", rp, exp=TRUE)[,2]
      tt <- rp$inp$time; t0 <- rp$inp$maninterval[1]; tE <- rp$inp$maneval; lastobs <- rp$inp$timerangeObs[2]
      if (min(inp$dtc) < 1) {
        alb2 <- list(anntime = unique(floor(tt)), annvec = tapply(log(estB), floor(tt), mean))
        alf2 <- list(anntime = unique(floor(tt)), annvec = tapply(log(estF), floor(tt), mean))
        bx <- exp(unname(alb2$annvec))/bscal; ttA <- as.numeric(names(alb2$annvec))
        fy <- exp(unname(alf2$annvec))/fscal
        keep_man <- (ttA >= t0) & (ttA <= tE); keep_int <- (ttA >= lastobs) & (ttA < t0); pre_mask <- ttA < t0
      } else {
        bx <- estB/bscal; fy <- estF/fscal; ttA <- tt
        keep_man <- (tt >= t0) & (tt <= tE); keep_int <- (tt >= lastobs) & (tt < t0); pre_mask <- tt < t0
      }
      n_man <- sum(keep_man, na.rm=TRUE); n_int <- sum(keep_int, na.rm=TRUE)
      lst_path[[i]]  <- data.frame(x = bx[keep_man],  y = fy[keep_man],  scenario = if (n_man) rep(leg_man[i], n_man) else character(0))
      lst_ipath[[i]] <- data.frame(x = bx[keep_int], y = fy[keep_int], scenario = if (n_int) rep(leg_man[i], n_int) else character(0))

      y0 <- if (any(pre_mask, na.rm=TRUE)) utils::tail(fy[pre_mask],1) else NA_real_
      y1 <- if (n_man > 0) fy[keep_man][1L] else NA_real_
      x1 <- if (n_man > 0) bx[keep_man][1L] else NA_real_
      x0 <- if (any(pre_mask, na.rm=TRUE)) utils::tail(bx[pre_mask],1) else NA_real_
      xv <- if (is.finite(x1)) x1 else x0
      if (is.finite(xv) && is.finite(y0) && is.finite(y1) && abs(y1-y0) > 0) {
        lst_jump[[i]] <- data.frame(x = xv, xend = xv, y = y0, yend = y1, scenario = leg_man[i])
      } else {
        lst_jump[[i]] <- data.frame(x = numeric(0), xend = numeric(0), y = numeric(0), yend = numeric(0), scenario = character(0))
      }
    }
    df_man      <- do.call(rbind, lst_path)
    df_man_int  <- do.call(rbind, lst_ipath)
    df_man_jump <- do.call(rbind, lst_jump)
  }

  # ---- build plot ----
  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = df_green_rect,   ggplot2::aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill = green_col,  color = NA) +
    ggplot2::geom_rect(data = df_yellowL_rect, ggplot2::aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill = yellow_col, color = NA) +
    ggplot2::geom_rect(data = df_yellowT_rect, ggplot2::aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill = yellow_col, color = NA) +
    ggplot2::geom_rect(data = df_red_rect,     ggplot2::aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill = red_col,    color = NA)

  if (!isTRUE(logax) && all(is.finite(xlim)) && 0 >= xlim[1] && 0 <= xlim[2]) {
    p <- p + ggplot2::geom_vline(xintercept = 0, color = "darkred", linetype = 2, linewidth = 0.4, na.rm = TRUE)
  }

  if (nrow(ell_abs) > 1) {
    p <- p + ggplot2::geom_polygon(
      data = data.frame(x = ell_abs$x/bscal, y = ell_abs$y/fscal),
      ggplot2::aes(x=x, y=y),
      fill = grDevices::rgb(0.827,0.827,0.827,0.7),
      color = grDevices::rgb(0.5,0.5,0.5,0.7),
      linewidth = 0.3
    )
  } else {
    p <- p + ggplot2::geom_point(
      data = data.frame(x = ell_abs$x/bscal, y = ell_abs$y/fscal),
      ggplot2::aes(x=x, y=y),
      color = "gray40", size = 2
    )
  }

  p <- p + ggplot2::geom_path(data = df_traj, ggplot2::aes(x=x, y=y),
                              color = grDevices::rgb(0,0,1,0.8), linewidth = 0.5)

  if (!is.null(df_pred)) {
    p <- p + ggplot2::geom_path(
      data = df_pred, ggplot2::aes(x = x, y = y),
      linetype = "33", color = grDevices::adjustcolor("blue", 0.8)
    )
    p <- p + ggplot2::geom_segment(
      data = df_EBseg, ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
      color = "blue", linetype = "33", linewidth = 1.0
    )
    p <- p + ggplot2::geom_point(
      data = data.frame(x = df_EBseg$xend, y = df_EBseg$yend),
      ggplot2::aes(x = x, y = y),
      shape = 23, fill = "gold", color = "black", size = 3, stroke = 0.5
    )
  }

  if (!is.null(df_man_int)  && nrow(df_man_int))
    p <- p + ggplot2::geom_path(data = df_man_int,  ggplot2::aes(x=x, y=y, color=scenario),
                                linetype = "33", linewidth = 0.9, show.legend = man.legend)
  if (!is.null(df_man_jump) && nrow(df_man_jump))
    p <- p + ggplot2::geom_segment(data = df_man_jump, ggplot2::aes(x=x, xend=xend, y=y, yend=yend, color=scenario),
                                   linewidth = 0.9, show.legend = man.legend)
  if (!is.null(df_man)      && nrow(df_man))
    p <- p + ggplot2::geom_path(data = df_man,      ggplot2::aes(x=x, y=y, color=scenario),
                                linewidth = 0.9, show.legend = man.legend)

  p <- p +
    ggplot2::geom_point(data = df_first, ggplot2::aes(x=x,y=y), shape=21, fill="white", color="black", size=2.5, stroke=0.5) +
    ggplot2::geom_text (data = df_first, ggplot2::aes(x=x,y=y,label=lab), vjust=-0.8, size=3) +
    ggplot2::geom_point(data = df_last,  ggplot2::aes(x=x,y=y), shape=22, fill="white", color="black", size=2.5, stroke=0.5) +
    ggplot2::geom_text (data = df_last,  ggplot2::aes(x=x,y=y,label=lab), vjust=-0.8, size=3)

  # Secondary axes (relative), with explicit breaks tied to final limits
  lab_rel_fmt <- function(vals) formatC(vals, format="f", digits=1)
  sec_factor_B  <- (utils::tail(Bmsy,1)[2] / bscal)
  sec_breaks_B  <- pretty(xlim / sec_factor_B, n = 6)
  sec_breaks_B  <- sec_breaks_B[is.finite(sec_breaks_B) & sec_breaks_B >= 0]
  if (!length(sec_breaks_B)) sec_breaks_B <- c(0.5, 1, 1.5, 2)

  if (isTRUE(logax)) {
    p <- p + ggplot2::scale_x_log10(
      limits = xlim,
      name   = if (is.null(xlabel)) xlab_expr else xlabel,
      labels = function(x) formatC(x, format="f", digits=0),
      expand = c(0,0),
      sec.axis = if (!rel.axes && isTRUE(ext)) ggplot2::sec_axis(
        ~ . / sec_factor_B,
        name   = expression(B[t]/B[MSY]),
        breaks = sec_breaks_B[sec_breaks_B > 0],
        labels = lab_rel_fmt
      ) else ggplot2::waiver()
    )
    p <- p + ggplot2::scale_y_log10(
      limits = ylim, name = ylab_expr, labels = function(x) formatC(x, format="f", digits=1), expand = c(0,0),
      sec.axis = if (!rel.axes && isTRUE(ext)) ggplot2::sec_axis(
        ~ . / (utils::tail(Fmsy,1)[2]/fscal),
        name = expression(F[t]/F[MSY]),
        breaks = function(lims) 10^pretty(log10(lims/(utils::tail(Fmsy,1)[2]/fscal))),
        labels = lab_rel_fmt
      ) else ggplot2::waiver()
    )
  } else {
    p <- p + ggplot2::scale_x_continuous(
      limits = xlim,
      name   = if (is.null(xlabel)) xlab_expr else xlabel,
      labels = function(x) formatC(x, format="f", digits=0),
      breaks = function(lims) pretty(lims, n = 6),
      expand = c(0,0),
      sec.axis = if (!rel.axes && isTRUE(ext)) ggplot2::sec_axis(
        ~ . / sec_factor_B,
        name   = expression(B[t]/B[MSY]),
        breaks = sec_breaks_B,
        labels = lab_rel_fmt
      ) else ggplot2::waiver()
    )
    p <- p + ggplot2::scale_y_continuous(
      limits = ylim, name = ylab_expr, labels = function(x) formatC(x, format="f", digits=1), expand = c(0,0),
      sec.axis = if (!rel.axes && isTRUE(ext)) ggplot2::sec_axis(
        ~ . / (utils::tail(Fmsy,1)[2]/fscal),
        name = expression(F[t]/F[MSY]),
        breaks = function(lims) pretty(lims/(utils::tail(Fmsy,1)[2]/fscal), n = 6),
        labels = lab_rel_fmt
      ) else ggplot2::waiver()
    )
  }

  # Scenario color scale for legend (if scenarios exist)
  if (!is.null(df_man) || !is.null(df_man_int) || !is.null(df_man_jump)) {
    if (is.null(scenario_colors)) {
      uniq_sc <- unique(c(
        if (!is.null(df_man)) df_man$scenario,
        if (!is.null(df_man_int)) df_man_int$scenario,
        if (!is.null(df_man_jump)) df_man_jump$scenario
      ))
      base_cols <- c('darkmagenta','cyan3','darkgreen','coral1','black','magenta','gold','green','cadetblue3','chocolate3','darkolivegreen3','cyan','darkred')
      colv <- rep(base_cols, length.out = length(uniq_sc)); names(colv) <- uniq_sc
    } else {
      colv <- scenario_colors
    }
    p <- p + ggplot2::scale_color_manual(values = colv, name = "Scenario",
                                         guide = if (man.legend) "legend" else "none")
  }

  p <- p + ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color="gray35", fill=NA, linewidth=2),
      legend.background = ggplot2::element_rect(fill = grDevices::adjustcolor("white", 0.8), color = NA),
      legend.position = if (man.legend) "top" else "none",
      legend.direction = "horizontal",
      axis.text.x      = ggplot2::element_text(size = 10, face = "bold"),
      axis.text.y      = ggplot2::element_text(size = 10, face = "bold"),
      axis.title.x     = ggplot2::element_text(size = 14, face = "bold"),
      axis.title.y     = ggplot2::element_text(size = 14, face = "bold"),
      axis.text.x.top  = ggplot2::element_text(size = 10, face = "bold"),
      axis.title.x.top = ggplot2::element_text(size = 14, face = "bold"),
      plot.margin  = ggplot2::margin(4,4,4,4)
    )

  # Mini-legend: E(B∞)
  if (isTRUE(plot.legend)) {
    px  <- xlim[2] - 0.02 * diff(xlim)
    py  <- ylim[2] - 0.04 * diff(ylim)
    gap <- 0.08 * diff(xlim)
    p <- p + ggplot2::annotate("point", x = px - gap, y = py,
                               shape = 23, size = 2.2, fill = "gold", color = "black", stroke = 0.5) +
      ggplot2::annotate("text", x = px, y = py,
                        label = "E(B[infinity])", parse = TRUE,
                        hjust = 1, vjust = 0.5, fontface = "bold", size = 5)
  }

  if (rep$opt$convergence != 0) {
    p <- p + ggplot2::annotate("text",
                               x = xlim[1] + 0.01*diff(xlim),
                               y = ylim[2] - 0.02*diff(ylim),
                               label = "!", color = "black", size = 5, fontface = "bold")
  }
  if (!is.null(stamp) && nzchar(stamp)) {
    p <- p + ggplot2::annotate("text",
                               x = xlim[2] - 0.01*diff(xlim),
                               y = ylim[1] + 0.02*diff(ylim),
                               label = stamp, hjust = 1, vjust = 0, size = 3)
  }
  p
}

# ---- NSE notes for R CMD check ----
# Put this in any file under R/ (here is fine).
utils::globalVariables(c(
  "time","lwr","upr","est","scenario","catch","catch_type",
  "index","obsrel","x","y","xend","yend","xmin","xmax","ymin","ymax","lab"
))
