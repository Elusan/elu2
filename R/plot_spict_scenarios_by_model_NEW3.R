#' Plot SPiCT Model Scenarios with Observed Indices and Optional Confidence Intervals
#'
#' This function generates a 3x2 patchwork grid of SPiCT model outputs: biomass,
#' B/Bmsy, catch, fishing mortality (F), F/Fmsy, and production.
#' It extends the behaviour of `plot_spict_scenarios_by_model_NEW2()` by
#' adding observed indices (scaled by estimated catchability) to the biomass
#' panel, mimicking the plotting style of \code{\link[spict]{plotspict.biomass}}.
#'
#' @param models Named list of fitted SPiCT model objects (e.g., \code{list(S1P = ..., S1F = ..., S1S = ...)}).
#'   All models should have the same structure for \code{$inp$timeI} and \code{$inp$obsI}.
#' @param production_fun Optional function to extract production curve data for each model.
#'   It should return a data frame with columns \code{B_K}, \code{Production}, and \code{Model}.
#' @param extract_catch_data Optional function to extract observed and predicted catch data
#'   for each model. It should return a data frame with at least:
#'   \code{time}, \code{catch}, \code{lwr}, \code{upr}, \code{catch_type}, and \code{scenario}.
#' @param scenario_colors Optional named vector of colors for each model name. If \code{NULL},
#'   defaults are assigned based on number of models and naming convention.
#' @param return_patchwork Logical; if \code{TRUE} (default) return a patchwork plot object.
#'   If \code{FALSE}, returns a list of individual ggplot objects.
#' @param lindwd Numeric; line width for dashed production curve lines.
#' @param show_CIs Logical; if \code{TRUE} (default) include confidence intervals as shaded ribbons.
#'
#' @details
#' The biomass plot overlays observed index values after scaling by estimated
#' catchability (\eqn{q}) from the model:
#' \deqn{\mathrm{obsI}_{scaled} = \mathrm{obsI} / \hat{q}}
#' Index 1 is plotted as solid blue points, index 2 as green-filled squares
#' with black borders.
#'
#' Colors for model lines are assigned automatically based on either a fixed palette
#' (if \code{length(models) == 6}) or a predefined scenario palette for up to 8 scenarios
#' (\code{S1}–\code{S8} × Pella/Fox/Schaefer models).
#'
#' @return
#' If \code{return_patchwork = TRUE}, returns a patchwork object that can be printed or
#' saved with \code{\link[ggplot2]{ggsave}}.
#' If \code{return_patchwork = FALSE}, returns a named list of individual ggplot objects:
#' \code{$biomass}, \code{$bbmsy}, \code{$catch}, \code{$f}, \code{$ffmsy}, \code{$production}.
#'
#' @examples
#' \dontrun{
#' # Fit models first
#' models <- list(
#'   S1P = fit.elu2(inp_Pella),
#'   S1F = fit.elu2(inp_Fox),
#'   S1S = fit.elu2(inp_Schaefer)
#' )
#'
#' # Plot all panels with observed indices on biomass plot
#' plot_spict_scenarios_by_model_NEW3(models)
#' }
#'
#' @seealso
#' \code{\link[spict]{plotspict.biomass}}, \code{\link[patchwork]{plot_layout}},
#' \code{\link[ggplot2]{geom_line}}, \code{\link[ggplot2]{geom_point}}
#'
#' @export
plot_spict_scenarios_by_model_NEW3 <- function(models,
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

  get_base_model <- function(modname) {
    if (grepl("P", modname, fixed = TRUE)) return("Pella")
    if (grepl("F", modname, fixed = TRUE)) return("Fox")
    if (grepl("S", modname, fixed = TRUE)) return("Schaefer")
    return("Unknown")
  }

  if (is.null(scenario_colors)) {
    if (length(models) == 6) {
      model_colors <- setNames(c("blue","#FF7F00","#A65628", "#D7191C", "#2B83BA","#4DAF4A"), model_names)
    } else {
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
      missing_colors <- setdiff(model_names, names(default_scenario_colors))
      if (length(missing_colors) > 0) {
        stop("Some model names have no defined colors: ", paste(missing_colors, collapse = ", "))
      }
      model_colors <- default_scenario_colors[model_names]
    }
  } else {
    model_colors <- scenario_colors[model_names]
  }

  theme_minimal_compact2 <- function(base_size = 10, base_family = "") {
    theme_minimal(base_size = base_size, base_family = base_family) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        axis.text = element_text(size = 10, face = "bold"),
        legend.position = c(0.98, 0.98),
        legend.justification = c("right", "top"),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 9, face = "bold"),
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

  make_biomass_plot <- function(df) {
    df$model <- factor(df$model, levels = model_names)
    p <- ggplot(df, aes(x = time, y = est, color = model, fill = model)) +
      { if (show_CIs) geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.06, color = NA, show.legend = FALSE) } +
      geom_line(linewidth = 0.8) +
      scale_color_manual(values = model_colors) +
      scale_fill_manual(values = model_colors) +
      labs(x = "Year", y = "Biomass (tons)") +
      theme_minimal_compact2() +
      guides(color = guide_legend(override.aes = list(linewidth = 1.8)))

    # Add observed indices like SPiCT
    first_model <- models[[1]]  # assume all have same inp$timeI/obsI structure
    qest <- get.par("logq", first_model, exp = TRUE)
    inp <- first_model$inp

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
                          shape = 22, color = "black", fill = "green", size = 2, stroke = 0.5, inherit.aes = FALSE)
    }
    return(p)
  }

  make_plot <- function(df, ylab_expr, hline = NULL) {
    df$model <- factor(df$model, levels = model_names)
    p <- ggplot(df, aes(x = time, y = est, color = model, fill = model)) +
      { if (show_CIs) geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.06, color = NA, show.legend = FALSE) } +
      geom_line(linewidth = 0.8) +
      scale_color_manual(values = model_colors) +
      scale_fill_manual(values = model_colors) +
      labs(x = "Year", y = ylab_expr) +
      theme_minimal_compact2() +
      guides(color = guide_legend(override.aes = list(linewidth = 1.8)))
    if (!is.null(hline)) {
      p <- p + geom_hline(yintercept = hline, linetype = "dashed", color = "black", linewidth = 0.8)
    }
    return(p)
  }

  plots <- list(
    biomass = make_biomass_plot(get_series("logB")),
    bbmsy   = make_plot(get_series("logBBmsy"), expression(bold(B/B[MSY])), hline = 1),
    ffmsy   = make_plot(get_series("logFFmsy"), expression(bold(F/F[MSY])), hline = 1),
    f       = make_plot(get_series("logF"), "Fishing mortality")
  )

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
      { if (show_CIs) geom_ribbon(data = predicted, aes(x = time, ymin = lwr, ymax = upr, fill = model), alpha = 0.06, show.legend = FALSE) } +
      geom_line(data = predicted, aes(x = time, y = catch, color = model), linewidth = 0.8) +
      geom_point(data = observed, aes(x = time, y = catch), color = "black", size = 1.3) +
      scale_color_manual(values = model_colors) +
      labs(x = "Year", y = "Catch (tons)") +
      theme_minimal_compact2() +
      guides(color = guide_legend(override.aes = list(linewidth = 2)))
  } else {
    plots$catch <- ggplot() + labs(title = "No catch data provided") + theme_void()
  }

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

  if (return_patchwork) {
    layout <- (plots$biomass | plots$bbmsy | plots$catch) /
      patchwork::plot_spacer() /
      (plots$f | plots$ffmsy | plots$production)
    return(layout +
             plot_layout(heights = c(1, 0.05, 1)) &
             theme(plot.margin = margin(4, 4, 4, 4)))
  } else {
    return(plots)
  }
}
