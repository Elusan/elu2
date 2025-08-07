#' Plot SPiCT Model Scenarios with Optional Confidence Intervals
#'
#' This function generates a 3x2 grid of plots (biomass, B/Bmsy, catch, F, F/Fmsy, production)
#' for a set of SPiCT models. It supports coloring by scenario or model name, with automatic
#' handling of confidence intervals and color assignment.
#'
#' @param models Named list of SPiCT model fits (e.g., list(S1P = ..., S1F = ..., S1S = ...)).
#' @param production_fun Optional function to extract production data for each model.
#' @param extract_catch_data Optional function to extract catch data for each model.
#' @param scenario_colors Optional named vector of colors for each model name. If NULL, a default is used.
#' @param return_patchwork Logical. If TRUE (default), returns a patchwork grid of plots.
#' @param lindwd Numeric. Line width for production plot dashed line.
#' @param show_CIs Logical. If TRUE (default), includes confidence interval ribbons.
#'
#' @return A patchwork plot object (or a list of plots if return_patchwork = FALSE).
#' @export
plot_spict_scenarios_by_model_NEW2_legend <- function(models,
                                                      production_fun = NULL,
                                                      extract_catch_data = NULL,
                                                      scenario_colors = NULL,
                                                      return_patchwork = TRUE,
                                                      lindwd = 0.8,
                                                      show_CIs = TRUE) {
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(grid)

  model_names <- names(models)

  # Helper to extract model family from name (e.g., S1P -> Pella)
  get_base_model <- function(modname) {
    if (grepl("P", modname, fixed = TRUE)) return("Pella")
    if (grepl("F", modname, fixed = TRUE)) return("Fox")
    if (grepl("S", modname, fixed = TRUE)) return("Schaefer")
    return("Unknown")
  }

  # Assign model colors based on number of models
  # Assign model colors based on number of models
  if (is.null(scenario_colors)) {
    if (length(models) == 3) {
      # Fixed colors for 3 canonical models: Pella, Fox, Schaefer
      model_colors <- setNames(c("#D7191C", "#2B83BA", "#1A9641"), model_names)
    } else {
      # Use your predefined full scenario color set
      default_scenario_colors <- c(
        "S1P" = "#1b9e77", "S1S" = "#d95f02", "S1F" = "#7570b3",
        "S2P" = "#e7298a", "S2S" = "#66a61e", "S2F" = "#e6ab02",
        "S3P" = "#666666", "S3S" = "#1f78b4", "S3F" = "#b2df8a",
        "S4P" = "#a6761d", "S4S" = "#fb9a99", "S4F" = "#8dd3c7",
        "S5P" = "#66c2a5", "S5S" = "#fc8d62", "S5F" = "#8da0cb",
        "S6P" = "#e78ac3", "S6S" = "#a6d854", "S6F" = "#ffd92f",
        "S7P" = "#e5c494", "S7S" = "#b3b3b3", "S7F" = "#1b7837",
        "S8P" = "#762a83", "S8S" = "#af8dc3", "S8F" = "#7fbf7b"
      )

      # Ensure that model names are present in the palette
      missing_colors <- setdiff(model_names, names(default_scenario_colors))
      if (length(missing_colors) > 0) {
        stop("Some model names have no defined colors: ", paste(missing_colors, collapse = ", "))
      }

      model_colors <- default_scenario_colors[model_names]
    }
  } else {
    # User provided scenario_colors
    model_colors <- scenario_colors[model_names]
  }

  # Define compact theme
  theme_minimal_compact2 <- function(base_size = 10, base_family = "") {
    theme_minimal(base_size = base_size, base_family = base_family) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        axis.text = element_text(size = 10, face = "bold"),
        #legend.position = c(0.98, 0.98),
        #legend.justification = c("right", "top"),

        legend.position = "right",
        legend.justification = "center",
        legend.box = "vertical",
        legend.box.just = "left",

        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        legend.key.size = unit(1, "lines"),
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

  # Extract series like logB, logF etc.
  get_series <- function(parname) {
    do.call(rbind, lapply(model_names, function(mod) {
      par <- get.par(parname, models[[mod]], exp = TRUE)
      df <- as.data.frame(par)
      df$time <- as.numeric(rownames(par))
      colnames(df)[colnames(df) == "ll"] <- "lwr"
      colnames(df)[colnames(df) == "ul"] <- "upr"
      df <- df[, c("time", "lwr", "est", "upr", "sd", "cv")]
      df$model <- mod
      return(df)
    }))
  }

  # Create plot with optional ribbons
  make_plot <- function(df, ylab_expr, hline = NULL) {
    df$model <- factor(df$model, levels = model_names)
    p <- ggplot(df, aes(x = time, y = est, color = model, fill = model)) +
      { if (show_CIs) geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.22, color = NA, show.legend = FALSE) } +
      geom_line(linewidth = 0.8) +
      scale_color_manual(values = model_colors) +
      scale_fill_manual(values = model_colors) +
      labs(x = "Year", y = ylab_expr) +
      theme_minimal_compact2() +
      guides(color = guide_legend(override.aes = list(linewidth = 1.8, ncol = 2)))
    if (!is.null(hline)) {
      p <- p + geom_hline(yintercept = hline, linetype = "dashed", color = "black", linewidth = 0.8)
    }
    return(p)
  }

  plots <- list(
    biomass = make_plot(get_series("logB"), "Biomass (tons)"),
    bbmsy   = make_plot(get_series("logBBmsy"), expression(bold(B/B[MSY])), hline = 1),
    ffmsy   = make_plot(get_series("logFFmsy"), expression(bold(F/F[MSY])), hline = 1),
    f       = make_plot(get_series("logF"), "Fishing mortality")
  )

  # Catch plot if available
  if (!is.null(extract_catch_data)) {
    catch_all <- bind_rows(lapply(model_names, function(mod) {
      extract_catch_data(models[[mod]], scenario_name = mod)
    }))
    catch_all$scenario <- as.character(catch_all$scenario)

    predicted <- catch_all %>% filter(catch_type == "Predicted")
    observed  <- catch_all %>% filter(catch_type == "Observed")

    predicted$model <- factor(predicted$scenario, levels = model_names)
    observed$model  <- factor(observed$scenario,  levels = model_names)

    plots$catch <- ggplot() +
      { if (show_CIs) geom_ribbon(data = predicted, aes(x = time, ymin = lwr, ymax = upr, fill = model), alpha = 0.22, show.legend = FALSE) } +
      geom_line(data = predicted, aes(x = time, y = catch, color = model), size = 0.8) +
      geom_point(data = observed, aes(x = time, y = catch), color = "black", size = 1.3) +
      scale_color_manual(values = model_colors) +
      labs(x = "Year", y = "Catch (tons)") +
      theme_minimal_compact2() +
      guides(color = guide_legend(override.aes = list(linewidth = 2, ncol = 2)))
  } else {
    plots$catch <- ggplot() + labs(title = "No catch data provided") + theme_void()
  }

  # Production plot if available
  if (!is.null(production_fun)) {
    prod_df <- bind_rows(lapply(model_names, function(mod) {
      production_fun(models[[mod]], model_name = mod)
    }))
    prod_df$Model <- factor(prod_df$Model, levels = model_names)
    max_pts <- prod_df %>% group_by(Model) %>% slice_max(Production, n = 1) %>% ungroup()
    plots$production <- ggplot(prod_df, aes(x = B_K, y = Production, color = Model)) +
      geom_line(size = lindwd, linetype = "dashed") +
      geom_point(data = max_pts, aes(shape = Model), size = 3) +
      scale_color_manual(values = model_colors) +
      scale_shape_manual(values = 15 + seq_along(model_names)) +
      labs(x = expression(bold(B/K)), y = "Production") +
      theme_minimal_compact2()
  } else {
    plots$production <- ggplot() + labs(title = "No production data provided") + theme_void()
  }

  # Combine with patchwork
  if (return_patchwork) {
    layout <- (plots$biomass | plots$bbmsy | plots$catch) /
      patchwork::plot_spacer() /
      (plots$f | plots$ffmsy | plots$production)
    # return(layout +
    #          plot_layout(heights = c(1, 0.05, 1)) &
    #          theme(plot.margin = margin(4, 4, 4, 4)))
    return(layout +
             plot_layout(heights = c(1, 0.05, 1), guides = "collect") &
             theme(
               plot.margin = margin(4, 4, 4, 4),
               legend.position = "right",
               legend.justification = "center",
               legend.box = "vertical",
               legend.box.just = "left",
               legend.title = element_blank(),
               legend.text = element_text(size = 10, face = "bold"),
               legend.key.size = unit(0.8, "lines")
             ) &
             guides(color = guide_legend(ncol = 2)))

  } else {
    return(plots)
  }
}
