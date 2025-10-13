#' Catch panel (ggplot2) with MSY band and optional scenario overlays
#'
#' @description
#' Plots observed and predicted catch with an MSY reference band.
#' Works with a fitted `spictcls` (with/without `$man`) or a raw `rep$man` list.
#' Colours, ordering, and styling are preserved.
#'
#' @param rep_man A fitted SPiCT report (`spictcls`) or a `rep$man` list.
#' @param scenario_color Optional named colour vector for scenarios; defaults to man_cols().
#' @param show_CIs Logical; draw base Cpred CI edges (≤ t0). Default TRUE.
#' @param CI Confidence level in (0,1). Default 0.95.
#' @param show_legend Logical; show legend when scenarios exist. Default FALSE.
#' @return ggplot object
#' @export
my_plot_manage_catch_panel <- function(rep_man,
                                            scenario_color = NULL,
                                            show_CIs = TRUE,
                                            CI = 0.95,
                                            show_legend = FALSE) {
  # Accept fitted spictcls or rep$man list
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

  dat <- elu2_prepare_manage_panel_data(rep_man, CI = CI)

  # Canonical scenario order + extras
  sc_order <- c("currentCatch","currentF","Fmsy","noF",
                "reduceF25","increaseF25","msyHockeyStick","ices")
  sc_present <- if (("man" %in% names(rep_man)) && length(rep_man$man)) names(rep_man$man) else character(0)
  sc_core   <- intersect(sc_order, sc_present)
  sc_other  <- setdiff(sc_present, sc_core)
  scs_final <- c(sc_core, sort(sc_other))
  has_man   <- length(scs_final) > 0

  if (!is.null(dat$catch_pred) && nrow(dat$catch_pred)) {
    dat$catch_pred$scenario <- factor(dat$catch_pred$scenario, levels = scs_final, ordered = TRUE)
  }

  # Stable named palette
  if (is.null(scenario_color)) {
    cols <- if (length(scs_final)) { x <- man_cols(length(scs_final)); names(x) <- scs_final; x } else NULL
  } else {
    cols <- scenario_color
    if (length(scs_final)) {
      if (is.null(names(cols))) names(cols) <- scs_final
      cols <- cols[match(scs_final, names(cols))]
      if (anyNA(cols)) {
        fill <- man_cols(length(scs_final))
        cols[is.na(cols)] <- fill[is.na(cols)]
      }
      names(cols) <- scs_final
    }
  }

  spict_blue_mean    <- "blue"
  spict_blue_ci_line <- grDevices::rgb(0, 0, 1, 0.20)

  # Management start (t0) for base catch mask
  t0 <- dat$t0

  # Base Cpred up to t0 (if available)
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

  # MSY band (time-varying if available)
  tvgflag <- .sfb(rep_man$inp, "timevaryinggrowth") || .sfb(rep_man$inp, "logmcovflag")
  msy_df <- NULL
  if (isTRUE(tvgflag)) {
    MSYmat <- try(get.par("logMSYvec", rep_man, exp = TRUE, CI = CI), silent = TRUE)
    if (!inherits(MSYmat, "try-error")) {
      t_msy  <- suppressWarnings(as.numeric(rownames(MSYmat)))
      if (!length(t_msy) || any(!is.finite(t_msy))) t_msy <- .sfn(rep_man$inp, "time", numeric(0))
      msy_df <- data.frame(time = t_msy, lwr = MSYmat[,1], est = MSYmat[,2], upr = MSYmat[,3])
    }
  }
  if (is.null(msy_df)) {
    MSYone <- try(get.par("logMSY", rep_man, exp = TRUE, CI = CI), silent = TRUE)
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
  if (nrow(msy_df)) {
    msy_df$lwr[msy_df$lwr < 0] <- 0
    msy_df$est[msy_df$est < 0] <- 0
    msy_df$upr[msy_df$upr < 0] <- 0
  }

  darker_color <- function(col, f = 0.6) {
    m <- grDevices::col2rgb(col) / 255
    grDevices::rgb(m[1]*f, m[2]*f, m[3]*f)
  }
  dot_blue <- darker_color(spict_blue_mean, 0.6)

  # Last observed catch start (for fitted-only dotted styling)
  obs_end <- if (!is.null(rep_man$inp$timeC) && length(rep_man$inp$timeC))
    utils::tail(rep_man$inp$timeC, 1) else NA_real_

  cu <- .sf(rep_man$inp, "catchunit", "")
  ylab <- if (nzchar(as.character(cu))) paste0("Catch (", cu, ")") else "Catch"

  pC <- ggplot2::ggplot() +
    # MSY ribbon + mean (back)
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
    # Vlines UNDER scenario overlays
    { if (has_man)
      add_management_vlines_catch_good(rep_man, color = "grey75", linetype = "solid",
                                       linewidth = 0.4, lineend = "butt")
      else
        add_management_vlines_catch_good(rep_man, color = "grey75", linetype = "solid",
                                         linewidth = 0.4, lineend = "butt")[1]
    } +
    # Scenario predicted catch (drawn above vlines)
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
    # Base predicted catch ≤ t0: CI edges (dashed) + mean (solid/dotted split at obs_end)
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
    {
      if (length(scs_final)) {
        ggplot2::scale_color_manual(
          values  = cols,
          limits  = scs_final,
          drop    = FALSE,
          guide   = if (show_legend) "legend" else "none",
          na.translate = FALSE
        )
      } else {
        ggplot2::scale_color_discrete(guide = "none")
      }
    } +
    ggplot2::scale_fill_manual(values = if (is.null(cols)) NA else cols, guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.04))) +
    theme_minimal_compact2_good_local()

  pC
}
