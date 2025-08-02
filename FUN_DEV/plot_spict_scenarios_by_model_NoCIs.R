#' Plot Biomass, Fishing Mortality, Catch, and Production for Multiple SPiCT/ELU Model Scenarios (No Confidence Intervals)
#'
#' Generates a grid of time series plots (biomass, B/Bmsy, F/Fmsy, fishing mortality, catch, and production) for a list of fitted stock assessment models, with each scenario colored and labeled. Confidence intervals are omitted—only point estimates are shown. Returns either a combined patchwork plot or a list of individual plots.
#'
#' @param models A named list of fitted model objects. The names should indicate scenario/model codes (e.g., "S1F", "S1P", "S1S").
#' @param production_fun Optional. A function that takes a model object and model name, and returns a data frame for production plotting. See \code{get_production_data()}.
#' @param extract_catch_data Optional. A function that takes a model object and scenario name, and returns a data frame with observed and predicted catch (without CIs).
#' @param return_patchwork Logical. If \code{TRUE} (default), returns a patchwork grid of all plots. If \code{FALSE}, returns a named list of individual plots.
#' @param lindwd Numeric. Line width for plotted lines (default: 1.4).
#'
#' @details
#' This function is similar to \code{plot_spict_scenarios_by_model()}, but only plots median or mean estimates (no ribbons for confidence intervals). Each panel uses a custom compact ggplot2 theme. Catch and production panels require user-supplied functions for extracting and formatting data.
#'
#' Time series panels include:
#' \itemize{
#'   \item Biomass (tons)
#'   \item Biomass relative to Bmsy (\eqn{B/B_{MSY}})
#'   \item Fishing mortality relative to Fmsy (\eqn{F/F_{MSY}})
#'   \item Fishing mortality
#'   \item Catch (observed and predicted, if available)
#'   \item Production (if available)
#' }
#'
#' If \code{return_patchwork = TRUE}, plots are arranged in a 2x3 grid with spacers for balance.
#'
#' @return
#' Either a patchwork plot object (default), or a named list of ggplot2 plot objects.
#'
#' @seealso \code{\link{plot_spict_scenarios_by_model}}, \code{\link{get_production_data}}, \code{\link{elu_extract_catch_data}}
#'
#' @examples
#' \dontrun{
#' models <- list(
#'   S1F = fit_fox,
#'   S1P = fit_pella,
#'   S1S = fit_schaefer
#' )
#' plot_spict_scenarios_by_model_NoCIs(models,
#'   production_fun = get_production_data,
#'   extract_catch_data = elu_extract_catch_data
#' )
#' }
#' @export
plot_spict_scenarios_by_model_NoCIs <- function(models,
                                                production_fun = NULL,
                                                extract_catch_data = NULL,
                                                return_patchwork = TRUE,
                                                lindwd = 1.2) {
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(grid)

  model_colors <- c(
    "S1F" = "#2B83BA",
    "S1P" = "#D7191C",
    "S1S" = "#1A9641"
  )

  get_base_model <- function(modname) {
    if (modname == "S1P") return("S1P")
    if (modname == "S1S") return("S1S")
    if (modname == "S1F") return("S1F")
    return(NA)
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
        legend.text = element_text(size = 10, face = "bold"),
        legend.key.size = unit(0.8, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, colour = "grey45", linewidth = 2),
        axis.ticks = element_line(linewidth = 0.8, color = "grey45"),
        axis.ticks.length = unit(3, "pt"),
        strip.background = element_rect(fill = "grey45", color = "grey45", linewidth = 0.5),
        strip.text = element_text(face = "bold", size = rel(1)),
        text = element_text(face = "bold", size = 10),
        plot.margin = margin(3, 3, 3, 3)
      )
  }

  get_series <- function(parname) {
    do.call(rbind, lapply(names(models), function(mod) {
      par <- get.par(parname, models[[mod]], exp = TRUE)
      df <- as.data.frame(par)
      df$time <- as.numeric(rownames(par))
      df$model <- get_base_model(mod)
      # Only keep time and estimate for line
      df <- df[, c("time", "est", "model")]
      return(df)
    }))
  }

  make_plot <- function(df, ylab_expr, hline = NULL) {
    df$model <- factor(df$model, levels = names(model_colors))
    p <- ggplot(df, aes(x = time, y = est, color = model)) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(values = model_colors, drop = FALSE) +
      labs(x = "Year", y = ylab_expr, color = "Model") +
      theme_minimal_compact2() +
      guides(color = guide_legend(override.aes = list(linewidth = 2)))
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

  # Catch plot (no CI)
  if (!is.null(extract_catch_data)) {
    catch_all <- dplyr::bind_rows(
      lapply(names(models), function(mod) {
        extract_catch_data(models[[mod]], scenario_name = mod)
      })
    )
    predicted <- filter(catch_all, catch_type == "Predicted")
    observed  <- filter(catch_all, catch_type == "Observed")
    predicted$scenario <- as.character(predicted$scenario)
    predicted$model <- vapply(predicted$scenario, get_base_model, character(1))
    observed$model  <- vapply(observed$scenario,  get_base_model, character(1))
    predicted$model <- factor(predicted$model, levels = names(model_colors))
    observed$model  <- factor(observed$model,  levels = names(model_colors))

    plots$catch <- ggplot() +
      geom_line(data = predicted, aes(x = time, y = catch, color = model), size = lindwd) +
      geom_point(data = observed, aes(x = time, y = catch), color = "black", size = 1.2) +
      scale_color_manual(values = model_colors, drop = FALSE) +
      labs(x = "Year", y = "Catch (tons)", color = "Model") +
      theme_minimal_compact2() +
      guides(color = guide_legend(override.aes = list(linewidth = 2)))
  } else {
    plots$catch <- ggplot() + labs(title = "No extract_catch_data function provided") + theme_void()
  }

  # Production plot (by model)
  if (!is.null(production_fun)) {
    prod_df <- do.call(rbind, lapply(names(models), function(mod) {
      production_fun(models[[mod]], model_name = get_base_model(mod))
    }))
    prod_df$Model <- factor(prod_df$Model, levels = names(model_colors))
    max_points <- prod_df %>%
      group_by(Model) %>%
      slice_max(Production, n = 1) %>%
      ungroup()
    plots$production <- ggplot(prod_df, aes(x = B_K, y = Production, color = Model)) +
      geom_line(size = lindwd, linetype = "dashed") +
      geom_point(data = max_points, aes(shape = Model), size = 3) +
      scale_shape_manual(values = 15 + seq_along(model_colors)) +
      scale_color_manual(values = model_colors) +
      labs(x = expression(bold(B/K)), y = "Production", color = "Model", shape = "Model") +
      theme_minimal_compact2()
  } else {
    plots$production <- ggplot() + labs(title = "No production_fun provided") + theme_void()
  }

  if (return_patchwork) {
    spacer <- patchwork::plot_spacer() + theme_void()
    layout <- (plots$biomass | plots$bbmsy | plots$catch) /
      spacer /
      (plots$f | plots$ffmsy | plots$production)
    return(layout +
             plot_layout(heights = c(1, 0.05, 1)) &
             theme(plot.margin = margin(4, 4, 4, 4)))
  } else {
    return(plots)
  }
}
