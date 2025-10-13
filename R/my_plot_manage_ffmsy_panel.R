#' F/Fmsy panel (ggplot2) with optional scenario overlays
#'
#' @description
#' Plots the historical F/F_MSY trajectory with optional confidence ribbons.
#' Accepts either a full fitted `spictcls` (with or without `$man`) OR a raw
#' `rep$man` list of scenarios. Colours and order are kept unchanged.
#'
#' Behavior (fitted-only, no scenarios):
#' - Base line (F/Fmsy) is BLUE and SOLID up to last observed time.
#' - From the vertical line to the right, BLUE becomes DOTTED.
#' - CI ribbon/lines are continuous across the vertical line (no gap).
#'
#' @param rep_man A fitted SPiCT report (`spictcls`) or a `rep$man` list.
#' @param scenario_color Optional named character vector of colours to use; if `NULL`, uses man_cols().
#' @param show_CIs Logical; if `TRUE`, draws base-fit CI ribbon & bounds. Default `TRUE`.
#' @param CI Confidence level in (0, 1). Default `0.95`.
#' @param show_legend Logical; show scenario legend when scenarios exist. Default `FALSE`.
#' @return A `ggplot` object.
#' @export
my_plot_manage_ffmsy_panel <- function(rep_man,
                                             scenario_color = NULL,
                                             show_CIs = TRUE,
                                             CI = 0.95,
                                             show_legend = FALSE) {
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

  sc_order <- c("currentCatch","currentF","Fmsy","noF",
                "reduceF25","increaseF25","msyHockeyStick","ices")
  sc_present <- if (("man" %in% names(rep_man)) && length(rep_man$man)) names(rep_man$man) else character(0)
  sc_core   <- intersect(sc_order, sc_present)
  sc_other  <- setdiff(sc_present, sc_core)
  scs_final <- c(sc_core, sort(sc_other))

  if (!is.null(dat$ffmsy) && nrow(dat$ffmsy)) {
    dat$ffmsy$scenario <- factor(dat$ffmsy$scenario, levels = scs_final, ordered = TRUE)
  }

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
  spict_blue_ci_fill <- grDevices::rgb(0, 0, 1, 0.10)

  base_pre <- get_base_BB_FF_pre(rep_man, CI = CI)
  has_man  <- length(scs_final) > 0

  pF <- ggplot2::ggplot() +
    { if (isTRUE(show_CIs) && NROW(base_pre$FF)) ggplot2::geom_ribbon(
      data = base_pre$FF,
      ggplot2::aes(x = time, ymin = lwr, ymax = upr),
      inherit.aes = FALSE, fill = spict_blue_ci_fill, color = NA
    ) } +
    { if (isTRUE(show_CIs) && NROW(base_pre$FF)) list(
      ggplot2::geom_line(data = base_pre$FF, ggplot2::aes(x = time, y = lwr),
                         inherit.aes = FALSE, color = spict_blue_ci_line, linewidth = 0.7),
      ggplot2::geom_line(data = base_pre$FF, ggplot2::aes(x = time, y = upr),
                         inherit.aes = FALSE, color = spict_blue_ci_line, linewidth = 0.7)
    ) } +
    { if (NROW(base_pre$FF)) ggplot2::geom_line(
      data = base_pre$FF, ggplot2::aes(x = time, y = est),
      inherit.aes = FALSE, color = spict_blue_mean, linewidth = 0.8
    ) } +
    ggplot2::geom_hline(yintercept = 1, linetype = "solid") +
    # --- MOVED EARLIER: vertical line UNDER scenario overlays ---
    { if (has_man) {
      add_management_vlines_BF_good(rep_man, color = "grey75",
                                    linetype = "solid", linewidth = 0.4, lineend = "butt")
    } else if (is.finite(base_pre$t_last_obs)) {
      ggplot2::geom_vline(xintercept = base_pre$t_last_obs,
                          color = "grey75", linetype = "solid", linewidth = 0.4)
    } else NULL } +
    # --- Scenario overlays now draw ON TOP of the vertical line ---
    ggplot2::geom_line(
      data = dat$ffmsy,
      ggplot2::aes(x = time, y = est, color = scenario),
      linewidth = 0.8, na.rm = TRUE
    ) +
    ggplot2::labs(title="Relative fishing mortality", x = "Year", y = expression(bold(F/F[MSY]))) +
    {
      if (length(scs_final)) {
        ggplot2::scale_color_manual(
          values  = cols, limits  = scs_final, drop = FALSE,
          guide   = if (show_legend) "legend" else "none",
          na.translate = FALSE
        )
      } else {
        ggplot2::scale_color_discrete(guide = "none")
      }
    } +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.02))) +
    ggplot2::geom_line(
      data = base_pre$FF,
      mapping = ggplot2::aes(x = time, y = est, group = 1),
      inherit.aes = FALSE,
      colour = I(spict_blue_mean),
      linewidth = 0.9,
      show.legend = FALSE
    ) +
    theme_minimal_compact2_good_local()

  # Fitted-only: dotted post-obs with continuous CI (unchanged)
  if (!has_man) {
    FF_all <- try(get.par("logFFmsy", rep_man, exp = TRUE, CI = CI), silent = TRUE)
    if (!inherits(FF_all, "try-error")) {
      tt <- suppressWarnings(as.numeric(rownames(FF_all)))
      if (!length(tt) || any(!is.finite(tt))) {
        tt_full <- .sfn(rep_man$inp, "time", numeric(0))
        tt <- if (length(tt_full) >= nrow(FF_all)) tt_full[seq_len(nrow(FF_all))] else rep_len(tt_full, nrow(FF_all))
      }
      df_all <- data.frame(time = tt, lwr = FF_all[,1], est = FF_all[,2], upr = FF_all[,3])

      tlo <- base_pre$t_last_obs
      if (is.finite(tlo)) {
        df_post <- df_all[df_all$time >= tlo, , drop = FALSE]

        if (nrow(df_post)) {
          if (isTRUE(show_CIs)) {
            pF <- pF +
              ggplot2::geom_ribbon(
                data = transform(df_post, ymin = lwr, ymax = upr),
                ggplot2::aes(x = time, ymin = ymin, ymax = ymax),
                inherit.aes = FALSE, fill = spict_blue_ci_fill, color = NA
              ) +
              ggplot2::geom_line(
                data = df_post, ggplot2::aes(x = time, y = lwr),
                inherit.aes = FALSE, color = spict_blue_ci_line, linewidth = 0.7
              ) +
              ggplot2::geom_line(
                data = df_post, ggplot2::aes(x = time, y = upr),
                inherit.aes = FALSE, color = spict_blue_ci_line, linewidth = 0.7
              )
          }
          pF <- pF +
            ggplot2::geom_line(
              data = df_post, ggplot2::aes(x = time, y = est),
              linetype = "dotted", color = spict_blue_mean, linewidth = 0.8
            )
        }
      }
    }
  }

  pF
}
