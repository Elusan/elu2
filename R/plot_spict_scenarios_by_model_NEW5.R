#' Plot SPiCT Model Scenarios (6 panels) with SPiCT-style colored CIs & borders
#'
#' @description
#' 3×2 grid for multiple SPiCT models (e.g., S1P, S1F, S1S). Each model gets:
#' - estimate line in its own color,
#' - a semi-transparent CI ribbon in the same color,
#' - thin CI band borders (upper/lower) in the same color with reduced alpha,
#' mirroring the visual style used in `plot2_elu2_gg_good_AJUSTED.EE()`.
#' Thin solid grey vertical lines mark end-of-observations (SPiCT convention).
#'
#' @inheritParams plot_spict_scenarios_by_model_NEW5
#' @export
plot_spict_scenarios_by_model_NEW5 <- function(models,
                                               production_fun = NULL,
                                               extract_catch_data = NULL,
                                               scenario_colors = NULL,
                                               return_patchwork = TRUE,
                                               lindwd = 1,
                                               show_CIs = TRUE) {
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(grid)

  # --- hardening: dependencies & sandbox global theme/options ------------------
  if (!is.list(models) || length(models) == 0) {
    stop("`models` must be a non-empty named list of SPiCT fits.")
  }
  if (is.null(names(models)) || any(!nzchar(names(models)))) {
    stop("`models` must be a named list; all elements need non-empty names.")
  }
  if (!exists("get.par", mode = "function")) {
    stop("`get.par()` not found. Please load SPiCT (or your wrapper) before calling.")
  }
  # Avoid surprises from user-level theme_set()/options(): restore afterwards
  old_theme <- ggplot2::theme_get()
  on.exit(ggplot2::theme_set(old_theme), add = TRUE)
  ggplot2::theme_set(ggplot2::theme_gray())

  # Silence harmless dplyr informs (does not change results)
  old_opt <- options(dplyr.summarise.inform = FALSE)
  on.exit(options(old_opt), add = TRUE)

  # Lock dplyr verbs to avoid masked versions from other packages
  bind_rows  <- dplyr::bind_rows
  filter     <- dplyr::filter
  slice_max  <- dplyr::slice_max
  group_by   <- dplyr::group_by
  ungroup    <- dplyr::ungroup
  `%>%`      <- dplyr::`%>%`

  model_names <- names(models)

  # ---- Visual knobs (match the "good" style) ----------------------------------
  CI_FILL_ALPHA  <- 0.12   # ribbon opacity (fill)
  CI_EDGE_ALPHA  <- 0.25   # band-border opacity (lwr/upr lines)
  CI_EDGE_WIDTH  <- 0.6    # band-border line width
  EST_LINE_WIDTH <- 1    # main estimate line width
  VLINE_COL      <- "grey50"
  VLINE_SIZE     <- 0.2

  # ---- Theme (compact, crisp border, tidy legend) -----------------------------
  theme_minimal_compact2_good <- function(base_size = 10, base_family = "") {
    theme_minimal(base_size = base_size, base_family = base_family) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        axis.text  = element_text(size = 12, face = "bold"),
        legend.position = c(0.75, 0.98),
        legend.justification = c("left", "top"),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.title = element_blank(),
        legend.text  = element_text(size = 10),
        legend.key.size = unit(1, "lines"),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, colour = "grey45", linewidth = 2.5),
        axis.ticks = element_line(linewidth = 0.5, color = "grey45"),
        axis.ticks.length = unit(3, "pt"),
        strip.background = element_rect(fill = "grey45", color = "grey45", linewidth = 0.5),
        strip.text = element_text(face = "bold", size = rel(1)),
        text = element_text(face = "bold", size = 10),
        plot.margin = margin(3, 3, 3, 3)
      )
  }

  # Colors (same logic as your version)
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
    # --- hardening: ensure no missing/NA colors when a custom vector is supplied
    if (any(is.na(model_colors))) {
      bad <- model_names[is.na(model_colors)]
      stop("`scenario_colors` is missing named entries for: ", paste(bad, collapse = ", "))
    }
    # Validate all supplied colors are valid
    tryCatch(vapply(model_colors, grDevices::col2rgb, integer(3L)), error = function(e) {
      stop("`scenario_colors` contains an invalid color name or code: ", conditionMessage(e))
    })
  }

  # Helper: alpha blend for a named color vector
  mk_alpha <- function(cols, a) {
    vapply(cols, function(cl) grDevices::adjustcolor(cl, alpha.f = a), character(1))
  }
  fill_cols <- mk_alpha(model_colors, CI_FILL_ALPHA)   # for ribbons
  # edges share the same hue but lighter alpha via layer parameter

  # ---- SPiCT-style end-of-observation markers ---------------------------------
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
  catch_obs_end_global <- suppressWarnings(max(vapply(models, get_catch_obs_end, numeric(1)), na.rm = TRUE))

  # ---- Series extractor (ll/ul -> lwr/upr) ------------------------------------
  get_series <- function(parname) {
    do.call(rbind, lapply(model_names, function(mod) {
      par <- get.par(parname, models[[mod]], exp = TRUE)
      if (is.null(par)) stop("`get.par()` returned NULL for '", parname, "' in model '", mod, "'.")
      df <- as.data.frame(par)
      if (!all(c("ll","est","ul") %in% colnames(df))) {
        stop("`get.par('", parname, ")` must return columns ll, est, ul for model '", mod, "'.")
      }
      df$time <- as.numeric(rownames(par))
      colnames(df)[colnames(df) == "ll"] <- "lwr"
      colnames(df)[colnames(df) == "ul"] <- "upr"
      df <- df[, c("time", "lwr", "est", "upr", "sd", "cv")]
      df$model <- mod
      df
    }))
  }

  # ---- Index dot helper --------------------------------------------------------
  darker_color <- function(col, f = 0.6) {
    m <- grDevices::col2rgb(col)/255
    grDevices::rgb(m[1]*f, m[2]*f, m[3]*f)
  }
  dot_blue <- darker_color("blue", 0.6)

  # ---- Biomass (absolute) with indices ----------------------------------------
  make_biomass_plot <- function(df) {
    df$model <- factor(df$model, levels = model_names)
    p <- ggplot(df, aes(x = time)) +
      { if (show_CIs) list(
        geom_ribbon(aes(ymin = lwr, ymax = upr, fill = model), color = NA, show.legend = FALSE),
        geom_line(aes(y = lwr, color = model), linewidth = CI_EDGE_WIDTH, alpha = CI_EDGE_ALPHA, show.legend = FALSE),
        geom_line(aes(y = upr, color = model), linewidth = CI_EDGE_WIDTH, alpha = CI_EDGE_ALPHA, show.legend = FALSE)
      ) } +
      geom_line(aes(y = est, color = model), linewidth = EST_LINE_WIDTH) +
      { if (is.finite(obs_end_overall))
        geom_vline(xintercept = obs_end_overall, color = VLINE_COL, linewidth = VLINE_SIZE) } +
      scale_color_manual(values = model_colors) +
      scale_fill_manual(values = fill_cols, guide = "none") +
      labs(x = "Year", y = "Biomass (tons)") +
      theme_minimal_compact2_good() +
      guides(color = guide_legend(override.aes = list(linewidth = 1.8)))

    # Observed indices (SPiCT style)
    first_model <- models[[1]]
    qest <- get.par("logq", first_model, exp = TRUE)
    inp  <- first_model$inp
    if (!is.null(inp$timeI) && length(inp$timeI) >= 1) {
      obs1 <- data.frame(time = inp$timeI[[1]], obs = inp$obsI[[1]] / qest[inp$mapq[1], 2])
      p <- p + geom_point(data = obs1, aes(x = time, y = obs),
                          color = dot_blue, shape = 16, size = 2, inherit.aes = FALSE)
    }
    if (!is.null(inp$timeI) && length(inp$timeI) >= 2) {
      obs2 <- data.frame(time = inp$timeI[[2]], obs = inp$obsI[[2]] / qest[inp$mapq[2], 2])
      p <- p + geom_point(data = obs2, aes(x = time, y = obs),
                          shape = 22, color = "black", fill = "green",
                          size = 2, stroke = 0.5, inherit.aes = FALSE)
    }
    p
  }

  # ---- B/Bmsy with indices (relative) -----------------------------------------
  make_bbmsy_plot <- function(df) {
    df$model <- factor(df$model, levels = model_names)
    p <- ggplot(df, aes(x = time)) +
      { if (show_CIs) list(
        geom_ribbon(aes(ymin = lwr, ymax = upr, fill = model), color = NA, show.legend = FALSE),
        geom_line(aes(y = lwr, color = model), linewidth = CI_EDGE_WIDTH, alpha = CI_EDGE_ALPHA, show.legend = FALSE),
        geom_line(aes(y = upr, color = model), linewidth = CI_EDGE_WIDTH, alpha = CI_EDGE_ALPHA, show.legend = FALSE)
      ) } +
      geom_line(aes(y = est, color = model), linewidth = EST_LINE_WIDTH) +
      { if (is.finite(obs_end_overall))
        geom_vline(xintercept = obs_end_overall, color = VLINE_COL, linewidth = VLINE_SIZE) } +
      geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
      scale_color_manual(values = model_colors) +
      scale_fill_manual(values = fill_cols, guide = "none") +
      labs(x = "Year", y = expression(bold(B/B[MSY]))) +
      theme_minimal_compact2_good() +
      guides(color = guide_legend(override.aes = list(linewidth = 1.8)))

    # Relative indices
    first_model <- models[[1]]
    qest  <- get.par("logq", first_model, exp = TRUE)
    Bmsy  <- get.par("logBmsy", first_model, exp = TRUE)
    Bmsy2 <- if (is.null(nrow(Bmsy))) Bmsy[2] else Bmsy[1, 2]
    inp   <- first_model$inp
    if (!is.null(inp$timeI) && length(inp$timeI) >= 1) {
      obs1_rel <- data.frame(time = inp$timeI[[1]], obs = (inp$obsI[[1]] / qest[inp$mapq[1], 2]) / Bmsy2)
      p <- p + geom_point(data = obs1_rel, aes(x = time, y = obs),
                          color = dot_blue, shape = 16, size = 2, inherit.aes = FALSE)
    }
    if (!is.null(inp$timeI) && length(inp$timeI) >= 2) {
      obs2_rel <- data.frame(time = inp$timeI[[2]], obs = (inp$obsI[[2]] / qest[inp$mapq[2], 2]) / Bmsy2)
      p <- p + geom_point(data = obs2_rel, aes(x = time, y = obs),
                          shape = 22, color = "black", fill = "green",
                          size = 2, stroke = 0.5, inherit.aes = FALSE)
    }
    p
  }

  # ---- Generic maker for F and F/Fmsy -----------------------------------------
  make_plot <- function(df, ylab_expr, hline = NULL) {
    df$model <- factor(df$model, levels = model_names)
    p <- ggplot(df, aes(x = time)) +
      { if (show_CIs) list(
        geom_ribbon(aes(ymin = lwr, ymax = upr, fill = model), color = NA, show.legend = FALSE),
        geom_line(aes(y = lwr, color = model), linewidth = CI_EDGE_WIDTH, alpha = CI_EDGE_ALPHA, show.legend = FALSE),
        geom_line(aes(y = upr, color = model), linewidth = CI_EDGE_WIDTH, alpha = CI_EDGE_ALPHA, show.legend = FALSE)
      ) } +
      geom_line(aes(y = est, color = model), linewidth = EST_LINE_WIDTH) +
      { if (is.finite(obs_end_overall))
        geom_vline(xintercept = obs_end_overall, color = VLINE_COL, linewidth = VLINE_SIZE) } +
      scale_color_manual(values = model_colors) +
      scale_fill_manual(values = fill_cols, guide = "none") +
      labs(x = "Year", y = ylab_expr) +
      theme_minimal_compact2_good() +
      guides(color = guide_legend(override.aes = list(linewidth = 1.8)))
    if (!is.null(hline)) {
      p <- p + geom_hline(yintercept = hline, linetype = "dashed", color = "black", linewidth = 0.8)
    }
    p
  }

  plots <- list(
    biomass = make_biomass_plot(get_series("logB")),
    bbmsy   = make_bbmsy_plot(get_series("logBBmsy")),
    ffmsy   = make_plot(get_series("logFFmsy"), expression(bold(F/F[MSY])), hline = 1),
    f       = make_plot(get_series("logF"), "Fishing mortality")
  )

  # ---- Catch panel -------------------------------------------------------------
  if (!is.null(extract_catch_data)) {
    # --- hardening: check that the user supplied a function
    if (!is.function(extract_catch_data)) {
      stop("`extract_catch_data` must be a function that returns a data.frame.")
    }
    catch_all <- bind_rows(lapply(model_names, function(mod) {
      extract_catch_data(models[[mod]], scenario_name = mod)
    }))
    # --- hardening: columns contract
    req_cols <- c("time","lwr","upr","catch","catch_type","scenario")
    missing_cols <- setdiff(req_cols, names(catch_all))
    if (length(missing_cols)) {
      stop("`extract_catch_data()` must return columns: ",
           paste(req_cols, collapse = ", "), ". Missing: ",
           paste(missing_cols, collapse = ", "))
    }

    catch_all$scenario <- as.character(catch_all$scenario)
    predicted <- catch_all %>% filter(catch_type == "Predicted")
    observed  <- catch_all %>% filter(catch_type == "Observed")
    predicted$model <- factor(predicted$scenario, levels = model_names)
    observed$model  <- factor(observed$scenario,  levels = model_names)

    plots$catch <- ggplot() +
      { if (show_CIs) list(
        geom_ribbon(data = predicted,
                    aes(x = time, ymin = lwr, ymax = upr, fill = model),
                    color = NA, show.legend = FALSE),
        geom_line(data = predicted,
                  aes(x = time, y = lwr, color = model),
                  linewidth = CI_EDGE_WIDTH, alpha = CI_EDGE_ALPHA, show.legend = FALSE),
        geom_line(data = predicted,
                  aes(x = time, y = upr, color = model),
                  linewidth = CI_EDGE_WIDTH, alpha = CI_EDGE_ALPHA, show.legend = FALSE)
      ) } +
      geom_line(data = predicted, aes(x = time, y = catch, color = model), linewidth = EST_LINE_WIDTH) +
      geom_point(data = observed,  aes(x = time, y = catch), color = "black", size = 1.3) +
      { if (is.finite(catch_obs_end_global))
        geom_vline(xintercept = catch_obs_end_global, color = VLINE_COL, linewidth = VLINE_SIZE) } +
      scale_color_manual(values = model_colors) +
      scale_fill_manual(values = fill_cols, guide = "none") +
      labs(x = "Year", y = "Catch (tons)") +
      theme_minimal_compact2_good() +
      guides(color = guide_legend(override.aes = list(linewidth = 2)))
  } else {
    plots$catch <- ggplot() + labs(title = "No catch data provided") + theme_void()
  }

  # ---- Production (unchanged by the CI style; no end-of-obs line) -------------
  if (!is.null(production_fun)) {
    if (!is.function(production_fun)) {
      stop("`production_fun` must be a function that returns a data.frame.")
    }
    prod_df <- bind_rows(lapply(model_names, function(mod) {
      production_fun(models[[mod]], model_name = mod)
    }))
    # --- hardening: columns contract for production curve
    reqp <- c("B_K","Production","Model")
    missp <- setdiff(reqp, names(prod_df))
    if (length(missp)) {
      stop("`production_fun()` must return columns: ",
           paste(reqp, collapse = ", "), ". Missing: ",
           paste(missp, collapse = ", "))
    }

    prod_df$Model <- factor(prod_df$Model, levels = model_names)
    max_pts <- prod_df %>% group_by(Model) %>% slice_max(Production, n = 1) %>% ungroup()
    plots$production <- ggplot(prod_df, aes(x = B_K, y = Production, color = Model)) +
      geom_line(size = lindwd, linetype = "dashed") +
      geom_point(data = max_pts, aes(shape = Model), size = 3) +
      scale_color_manual(values = model_colors) +
      scale_shape_manual(values = 15 + seq_along(model_names)) +
      labs(x = expression(bold(B/K)), y = "Production") +
      theme_minimal_compact2_good()
  } else {
    plots$production <- ggplot() + labs(title = "No production data provided") + theme_void()
  }

  # ---- Layout ------------------------------------------------------------------
  if (return_patchwork) {
    layout <- (plots$biomass | plots$bbmsy | plots$catch) /
      patchwork::plot_spacer() /
      (plots$f | plots$ffmsy | plots$production)
    layout + plot_layout(heights = c(1, 0.05, 1)) &
      theme(plot.margin = margin(4, 4, 4, 4))
  } else {
    plots
  }
}
