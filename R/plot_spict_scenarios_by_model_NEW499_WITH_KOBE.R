#' Add an Original_kobe_all_in_gg()-style mini-legend to a Kobe plot
#'
#' @description
#' Inserts a compact, **Original_kobe_all_in_gg()**-style mini-legend in the
#' **top-right** of a Kobe plot built with **ggplot2**:
#' a small gold diamond followed by a wider horizontal gap and a larger math
#' label `bold(E(B[infinity]))`. The legend is **suppressed** for seasonal models
#' (when the minimum of `inp$dtc` is `< 1`), matching the Original behavior.
#'
#' @details
#' The function reads the plotted **panel ranges** from the ggplot build object
#' and places the legend relative to the current axes (robust to ggplot2
#' version differences via several fallbacks). If the plot is seasonal
#' (i.e., `min(rep_like$inp$dtc) < 1`), the input plot is returned **unchanged**.
#'
#' @param p A `ggplot` object representing the Kobe panel.
#' @param rep_like A fitted SPiCT-like list (e.g., `spictcls`) that contains
#'   `inp$dtc` to decide whether seasonal time steps are present. Only
#'   `rep_like$inp$dtc` is inspected; the rest of the object is not modified.
#'
#' @return A `ggplot` object, with the mini-legend added (or unchanged for
#'   seasonal inputs).
#' @seealso [kobe_all_in_one_gg()]
#' @examples
#' \dontrun{
#'   # Suppose `p_kobe` is a ggplot Kobe panel and `rep` is a spictcls-like object
#'   p_kobe2 <- .add_kobe_original_mini_legend(p_kobe, rep_like = rep)
#'   print(p_kobe2)
#' }
#'
#' @keywords internal
#' @noRd
#' @importFrom ggplot2 ggplot_build annotate
# --- Helper (local): Original_kobe_all_in_gg()-style mini-legend -------------
# Adds a smaller gold diamond, more horizontal gap, and a larger math label
# "bold(E(B[infinity]))" in the top-right. Suppresses on seasonal time steps.
.add_kobe_original_mini_legend <- function(p, rep_like) {
  # If seasonal (min dtc < 1), Original suppressed this legend
  if (!is.null(rep_like$inp$dtc) && length(rep_like$inp$dtc) &&
      min(rep_like$inp$dtc, na.rm = TRUE) < 1) {
    return(p)
  }

  gb <- ggplot2::ggplot_build(p)
  # Extract axis limits robustly across ggplot2 versions
  get_range <- function(pp, axis = c("x","y")) {
    axis <- match.arg(axis)
    tryCatch({
      rng <- pp[[axis]]$range$range
      if (is.null(rng)) stop("no range slot")
      rng
    }, error = function(e) {
      # fallbacks seen in different ggplot2 versions
      if (!is.null(pp[[paste0(axis,".range")]])) {
        pp[[paste0(axis,".range")]]
      } else if (!is.null(pp[[paste0(axis,".range_c")]]) &&
                 !is.null(pp[[paste0(axis,".range_c")]]$range)) {
        pp[[paste0(axis,".range_c")]]$range
      } else {
        # last resort: limits from scales on panel layer data
        r <- range(unlist(lapply(gb$data, `[[`, if (axis == "x") "x" else "y")),
                   na.rm = TRUE)
        if (!all(is.finite(r))) c(0, 1) else r
      }
    })
  }

  pp <- gb$layout$panel_params[[1]]
  xr <- get_range(pp, "x")
  yr <- get_range(pp, "y")

  px  <- xr[2] - 0.02 * diff(xr)
  py  <- yr[2] - 0.04 * diff(yr)
  gap <- 0.12 * diff(xr)  # "more gap" like Original_kobe_all_in_gg()

  p +
    ggplot2::annotate("point", x = px - gap, y = py,
                      shape = 23, size = 2.5, fill = "gold",
                      color = "black", stroke = 0.6) +
    ggplot2::annotate("text",  x = px,       y = py,
                      label = "bold(E(B[infinity]))", parse = TRUE,
                      hjust = 1, vjust = 0.5, size = 4)
}


#' Plot SPiCT Model Scenarios (3×2) with Kobe panel (OSA path; no scenario overlays)
#'
#' @description
#' Produces a 3×2 **patchwork** of SPiCT time series and a Kobe panel for a
#' **named list** of fitted models (e.g., scenarios). Panels are:
#' **Biomass (with observed indices scaled by \eqn{\hat q})**, **B/B\eqn{_{\mathrm{MSY}}}**,
#' **Catch** (SPiCT last-observation convention), **F**, **F/F\eqn{_{\mathrm{MSY}}}**, and **Kobe**.
#' Confidence intervals are drawn as **dotted bounds** (if enabled). Thin **solid grey**
#' vertical lines mark the end of observed data (non-catch panels) and the SPiCT
#' **catch** convention on the catch panel.
#'
#' @details
#' **Kobe panel** is built via `kobe_all_in_one_gg.EE()` using the **first** model in
#' `models`. Any `$man` element is removed to **suppress scenario overlays** so the
#' **OSA** path and **E(B∞)** segment are visible; the **final-year white square**
#' is preserved. The panel’s mini-legend is replaced with an
#' **Original_kobe_all_in_gg()**-style legend using
#' `.add_kobe_original_mini_legend()`.
#'
#' Colors:
#' - If `scenario_colors` is `NULL`, a stable default palette is selected. For
#'   non-matching model names, provide your own `scenario_colors`.
#'
#' Catch panel:
#' - If `extract_catch_data` is provided, it must return a data frame with the
#'   columns `time, catch, lwr, upr, catch_type` (values `"Observed"`/`"Predicted"`),
#'   and `scenario` (the model name).
#'
#' @param models A **named list** of fitted SPiCT objects (`spictcls`), one per
#'   scenario/model. Names are used in legends and color mapping.
#' @param production_fun Deprecated; kept for backward compatibility (ignored).
#' @param extract_catch_data Optional function to extract **observed** and
#'   **predicted** catch time series for each model. Must return columns
#'   `time, catch, lwr, upr, catch_type, scenario`.
#' @param scenario_colors Optional **named** vector of colors keyed by model names.
#'   If `NULL`, an internal default palette is used (requires names to match).
#' @param return_patchwork Logical; if `TRUE` (default) returns a **patchwork**
#'   object. If `FALSE`, returns a **list** with elements
#'   `biomass, bbmsy, catch, f, ffmsy, kobe`.
#' @param lindwd Kept for compatibility (unused).
#' @param show_CIs Logical; if `TRUE` (default) draws **dotted** CI bounds.
#'
#' @return Either a **patchwork** object (default) or a named **list** of six
#' `ggplot` panels (when `return_patchwork = FALSE`).
#'
#' @section Vertical reference lines:
#' - Non-catch panels: an overall end-of-observation vertical line is computed
#'   across models from `inp$timerangeObs` (or fallbacks).
#' - Catch panel: an end-of-observed-catch line follows SPiCT’s convention,
#'   using the last **integer** time (`dtc < 1`) if seasonal.
#'
#' @examples
#' \dontrun{
#'   # models <- list(S1P = m_S1P, S1S = m_S1S, S1F = m_S1F)  # each a spictcls
#'   p <- plot_spict_scenarios_by_model_NEW499_WITH_KOBE(
#'          models,
#'          extract_catch_data = elu_extract_catch_data_spict,
#'          show_CIs = TRUE
#'        )
#'   p  # patchwork
#'
#'   # Or retrieve the individual panels:
#'   panels <- plot_spict_scenarios_by_model_NEW499_WITH_KOBE(
#'               models,
#'               extract_catch_data = elu_extract_catch_data_spict,
#'               return_patchwork = FALSE
#'             )
#'   panels$kobe
#' }
#'
#' @seealso [kobe_all_in_one_gg()], \code{.add_kobe_original_mini_legend()}
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_vline geom_hline
#' @importFrom ggplot2 scale_color_manual labs theme guide_legend guides theme_minimal
#' @importFrom dplyr bind_rows filter
#' @importFrom patchwork plot_spacer plot_layout
#' @importFrom grid unit
plot_spict_scenarios_by_model_NEW499_WITH_KOBE <- function(models,
                                                           production_fun = NULL,   # kept for compat; unused
                                                           extract_catch_data = NULL,
                                                           scenario_colors = NULL,
                                                           return_patchwork = TRUE,
                                                           lindwd = 0.8,            # kept for compat; unused
                                                           show_CIs = TRUE) {
  # deps
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(grid)

  stopifnot(length(models) >= 1)
  model_names <- names(models)
  if (is.null(model_names) || any(!nzchar(model_names))) {
    stop("`models` must be a *named* list of spictcls objects.")
  }

  # --- vertical-line style (thin solid grey) ---
  vline_col  <- "grey50"
  vline_size <- 0.2

  # --- End-of-observation for non-catch panels ---
  get_obs_end_generic <- function(m) {
    tr <- m$inp$timerangeObs
    if (!is.null(tr) && length(tr) >= 2 && is.finite(tr[2])) return(tr[2])
    if (!is.null(m$inp$timeC) && length(m$inp$timeC))       return(tail(m$inp$timeC, 1))
    if (!is.null(m$inp$timeI) && length(m$inp$timeI) &&
        length(m$inp$timeI[[1]]) > 0)                       return(tail(m$inp$timeI[[1]], 1))
    if (!is.null(m$inp$time)  && length(m$inp$time))        return(max(m$inp$time, na.rm = TRUE))
    NA_real_
  }
  cand_noncatch <- vapply(models, get_obs_end_generic, numeric(1))
  obs_end_overall <- if (all(!is.finite(cand_noncatch))) NA_real_ else suppressWarnings(max(cand_noncatch, na.rm = TRUE))

  # --- SPiCT-like end of observed catch (seasonal uses last integer timeC) ---
  get_catch_obs_end <- function(m) {
    inp <- m$inp
    if (is.null(inp$timeC) || !length(inp$timeC)) return(NA_real_)
    if (!is.null(inp$dtc) && length(inp$dtc) && min(inp$dtc, na.rm = TRUE) < 1) {
      ix <- which((inp$timeC %% 1) == 0)
      if (length(ix)) return(tail(inp$timeC[ix], 1))
      return(tail(inp$timeC, 1))
    } else {
      return(tail(inp$timeC, 1))
    }
  }
  co_ends <- vapply(models, get_catch_obs_end, numeric(1))
  catch_obs_end_global <- if (all(!is.finite(co_ends))) NA_real_ else suppressWarnings(max(co_ends, na.rm = TRUE))

  # --- color mapping ---
  if (is.null(scenario_colors)) {
    if (length(models) == 6) {
      model_colors <- setNames(c("blue","#FF7F00","#A65628","#D7191C","#2B83BA","#4DAF4A"), model_names)
    } else {
      default_scenario_colors <- c(
        "S1P"="#1b9e77","S1S"="#d95f02","S1F"="#7570b3",
        "S2P"="#e7298a","S2S"="#66a61e","S2F"="#e6ab02",
        "S3P"="#666666","S3S"="#1f78b4","S3F"="#b2df8a",
        "S4P"="#a6761d","S4S"="#fb9a99","S4F"="#8dd3c7",
        "S5P"="#66c2a5","S5S"="#fc8d62","S5F"="#8da0cb",
        "S6P"="#e78ac3","S6S"="#a6d854","S6F"="#ffd92f",
        "S7P"="#e5c494","S7S"="#b3b3b3","S7F"="#1b7837",
        "S8P"="#762a83","S8S"="#af8dc3","S8F"="#7fbf7b"
      )
      missing_colors <- setdiff(model_names, names(default_scenario_colors))
      if (length(missing_colors) > 0) {
        stop("Some model names have no defined colors: ", paste(missing_colors, collapse = ", "))
      }
      model_colors <- default_scenario_colors[model_names]
    }
  } else {
    model_colors <- scenario_colors[model_names]
  }

  # --- themes ---
  theme_minimal_compact2 <- function(base_size = 10, base_family = "") {
    theme_minimal(base_size = base_size, base_family = base_family) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        axis.text  = element_text(size = 10, face = "bold"),
        legend.position = c(0.85, 0.98),
        legend.justification = c("right","top"),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.title = element_blank(),
        legend.text  = element_text(size = 9, face = "bold"),
        legend.key.size = unit(0.7, "lines"),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, colour = "grey35", linewidth = 2),
        axis.ticks = element_line(linewidth = 0.5, color = "grey35"),
        axis.ticks.length = unit(3, "pt"),
        strip.background = element_rect(fill = "grey35", color = "grey35", linewidth = 0.5),
        strip.text = element_text(face = "bold", size = rel(1)),
        text = element_text(face = "bold", size = 10),
        plot.margin = margin(3, 3, 3, 3)
      )
  }

  # --- series extractor from get.par() ---
  get_series <- function(parname) {
    do.call(rbind, lapply(model_names, function(mod) {
      par <- get.par(parname, models[[mod]], exp = TRUE)
      df  <- as.data.frame(par)
      df$time <- as.numeric(rownames(par))
      colnames(df)[colnames(df) == "ll"] <- "lwr"
      colnames(df)[colnames(df) == "ul"] <- "upr"
      df <- df[, c("time","lwr","est","upr","sd","cv")]
      df$model <- mod
      df
    }))
  }

  # --- biomass panel with observed indices (scaled by q-hat) ---
  make_biomass_plot <- function(df) {
    df$model <- factor(df$model, levels = model_names)
    p <- ggplot(df, aes(x = time, y = est, color = model)) +
      { if (show_CIs) geom_line(aes(y = lwr), linetype = "dotted", linewidth = 0.6) } +
      { if (show_CIs) geom_line(aes(y = upr), linetype = "dotted", linewidth = 0.6) } +
      geom_line(linewidth = 0.8) +
      { if (is.finite(obs_end_overall))
        geom_vline(xintercept = obs_end_overall, color = vline_col, linewidth = vline_size) } +
      scale_color_manual(values = model_colors) +
      labs(x = "Year", y = "Biomass (tons)") +
      theme_minimal_compact2() +
      guides(color = guide_legend(override.aes = list(linewidth = 1.8)))

    # overlay observed indices from the first model (scaled by q-hat)
    first_model <- models[[1]]
    qest <- get.par("logq", first_model, exp = TRUE)
    inp  <- first_model$inp
    if (length(inp$timeI) >= 1) {
      obs1 <- data.frame(time = inp$timeI[[1]],
                         obs  = inp$obsI[[1]] / qest[inp$mapq[1], 2])
      p <- p + geom_point(data = obs1, aes(x = time, y = obs),
                          color = "blue", shape = 16, size = 2, inherit.aes = FALSE)
    }
    if (length(inp$timeI) >= 2) {
      obs2 <- data.frame(time = inp$timeI[[2]],
                         obs  = inp$obsI[[2]] / qest[inp$mapq[2], 2])
      p <- p + geom_point(data = obs2, aes(x = time, y = obs),
                          shape = 22, color = "black", fill = "green", size = 2,
                          stroke = 0.5, inherit.aes = FALSE)
    }
    p
  }

  # --- generic series panel builder ---
  make_plot <- function(df, ylab_expr, hline = NULL) {
    df$model <- factor(df$model, levels = model_names)
    p <- ggplot(df, aes(x = time, y = est, color = model)) +
      { if (show_CIs) geom_line(aes(y = lwr), linetype = "dotted", linewidth = 0.6) } +
      { if (show_CIs) geom_line(aes(y = upr), linetype = "dotted", linewidth = 0.6) } +
      geom_line(linewidth = 0.8) +
      { if (is.finite(obs_end_overall))
        geom_vline(xintercept = obs_end_overall, color = vline_col, linewidth = vline_size) } +
      scale_color_manual(values = model_colors) +
      labs(x = "Year", y = ylab_expr) +
      theme_minimal_compact2() +
      guides(color = guide_legend(override.aes = list(linewidth = 1.8)))
    if (!is.null(hline)) {
      p <- p + geom_hline(yintercept = hline, linetype = "dashed",
                          color = "black", linewidth = 0.8)
    }
    p
  }

  # --- assemble panels (except catch & kobe for now) ---
  plots <- list(
    biomass = make_biomass_plot(get_series("logB")),
    bbmsy   = make_plot(
      get_series("logBBmsy"),
      expression(bold(B / B[italic(MSY)])),
      hline = 1
    ),
    ffmsy   = make_plot(
      get_series("logFFmsy"),
      expression(bold(F / F[italic(MSY)])),
      hline = 1
    ),
    f       = make_plot(get_series("logF"), "Fishing mortality")
  )

  # --- Catch panel (via user extractor, SPiCT convention for vertical line) ---
  if (!is.null(extract_catch_data)) {
    catch_all <- bind_rows(lapply(model_names, function(mod) {
      extract_catch_data(models[[mod]], scenario_name = mod)
    }))
    predicted <- catch_all %>% filter(catch_type == "Predicted")
    observed  <- catch_all %>% filter(catch_type == "Observed")
    predicted$model <- factor(predicted$scenario, levels = model_names)
    observed$model  <- factor(observed$scenario,  levels = model_names)

    plots$catch <- ggplot() +
      { if (show_CIs) geom_line(data = predicted, aes(x = time, y = lwr, color = model),
                                linetype = "dotted", linewidth = 0.6) } +
      { if (show_CIs) geom_line(data = predicted, aes(x = time, y = upr, color = model),
                                linetype = "dotted", linewidth = 0.6) } +
      geom_line(data = predicted, aes(x = time, y = catch, color = model), linewidth = 0.8) +
      geom_point(data = observed,  aes(x = time, y = catch), color = "black", size = 1.3) +
      { if (is.finite(catch_obs_end_global))
        geom_vline(xintercept = catch_obs_end_global, color = vline_col, linewidth = vline_size) } +
      scale_color_manual(values = model_colors) +
      labs(x = "Year", y = "Catch (tons)") +
      theme_minimal_compact2() +
      guides(color = guide_legend(override.aes = list(linewidth = 2)))
  } else {
    plots$catch <- ggplot() + labs(title = "No catch data provided") + theme_void()
  }

  # --- Kobe panel (replaces Production), with Original-style mini-legend -------
  # Use the first model; drop scenarios so OSA path + EB∞ segment are drawn
  rep_kobe <- models[[ model_names[1] ]]
  if ("man" %in% names(rep_kobe)) rep_kobe$man <- NULL

  # Build Kobe from your helper; turn off its own mini-legend and re-add Original style
  p_kobe <- kobe_all_in_one_gg.EE(
    rep                = rep_kobe,
    logax              = FALSE,
    plot.legend        = FALSE,  # we'll add the Original-style legend next
    man.legend         = FALSE,  # no scenario legend
    ext                = TRUE,
    rel.axes           = FALSE,  # auto-switches internally for time-varying growth
    CI                 = 0.95,
    scenario_colors    = NULL,
    origin_offset_frac = 0.04
  )
  plots$kobe <- .add_kobe_original_mini_legend(p_kobe, rep_like = rep_kobe)

  # --- Layout & return ---
  if (return_patchwork) {
    layout <- (plots$biomass | plots$bbmsy | plots$catch) /
      patchwork::plot_spacer() /
      (plots$f | plots$ffmsy | plots$kobe)

    layout + plot_layout(heights = c(1, 0.05, 1)) &
      theme(plot.margin = margin(4, 4, 4, 4))
  } else {
    plots
  }
}
