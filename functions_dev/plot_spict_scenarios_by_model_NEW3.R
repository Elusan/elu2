#' Plot SPiCT Model Scenarios by Model Name with Optional Confidence Intervals
#'
#' This function generates a 3x2 grid of plots (biomass, B/Bmsy, F, F/Fmsy, catch, and production)
#' for a given list of SPiCT model fits. Each model should be named like 'S1F', 'S1P', etc.
#' The legend inside each panel shows the full model names (e.g., 'S1F', 'S1P', 'S1S'), while the colors
#' remain consistent across scenarios by model type (F = blue, P = red, S = green).
#'
#' @param models A named list of SPiCT model fits (e.g., list(S1F = ..., S1P = ..., S1S = ...)).
#'               Each model name should include both scenario and model type (e.g., 'S1F').
#' @param production_fun Optional. A function that extracts production data from each model.
#'                       Should return a data frame with columns \code{B_K}, \code{Production}, and \code{Model}.
#' @param extract_catch_data Optional. A function that extracts observed and predicted catch data from each model.
#'                           Should return a data frame with columns \code{time}, \code{catch}, \code{lwr}, \code{upr},
#'                           \code{catch_type}, and \code{model}.
#' @param return_patchwork Logical. If \code{TRUE} (default), returns a combined patchwork plot.
#'                         If \code{FALSE}, returns a named list of individual ggplot objects.
#' @param lindwd Line width for plot lines (default: 0.5).
#' @param show_CIs Logical. If \code{TRUE} (default), adds 95\% confidence ribbons around estimates.
#'
#' @return Either a patchwork object (default) or a list of ggplot2 plots (if \code{return_patchwork = FALSE}).
#'
#' @details
#' The function automatically extracts model types ('F', 'P', 'S') from the last character of the model names
#' and assigns consistent colors across panels. Legends appear inside each panel at the top-right.
#'
#' @examples
#' \dontrun{
#' # Example with 3 models:
#' models <- list(
#'   S1F = fit.elu2(...),
#'   S1P = fit.elu2(...),
#'   S1S = fit.elu2(...)
#' )
#' plot_spict_scenarios_by_model_NEW2(
#'   models = models,
#'   production_fun = get_production_data,
#'   extract_catch_data = elu_extract_catch_data
#' )
#' }
#'
#' @export
plot_spict_scenarios_by_model_NEW3 <- function(models,
                                               production_fun = NULL,
                                               extract_catch_data = NULL,
                                               return_patchwork = TRUE,
                                               lindwd = 0.5,
                                               show_CIs = TRUE) {
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(grid)

  model_names <- names(models)

  # Extract model type character (F/P/S)
  get_model_code <- function(model_name) substr(model_name, nchar(model_name), nchar(model_name))

  # Define model_colors based on last character
  model_types <- sapply(model_names, get_model_code)
  model_colors <- setNames(
    c("F" = "#2B83BA", "S" = "#1A9641", "P" = "#D7191C")[model_types],
    model_names
  )


  # Theme
  theme_minimal_compact2 <- function(base_size = 12, base_family = "") {
    theme_minimal(base_size = base_size, base_family = base_family) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        axis.text = element_text(size = 12, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        legend.key.size = unit(0.8, "lines"),
        legend.box = "horizontal",
        legend.position = c(0.98, 0.98),
        legend.justification = c("right", "top"),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, colour = "grey35", linewidth = 1.5),
        axis.ticks = element_line(linewidth = 0.6, color = "grey35"),
        axis.ticks.length = unit(3, "pt"),
        strip.background = element_rect(fill = "grey35", color = "grey35"),
        strip.text = element_text(face = "bold", size = rel(1)),
        plot.margin = margin(4, 4, 4, 4) #,
        #theme(legend.position = c(0.99, 0.99), legend.justification = c("right", "top"))

      )
  }

  get_series <- function(parname) {
    do.call(rbind, lapply(model_names, function(mod) {
      par <- get.par(parname, models[[mod]], exp = TRUE)
      df <- as.data.frame(par)
      df$time <- as.numeric(rownames(par))
      df$model <- mod
      colnames(df)[colnames(df) == "ll"] <- "lwr"
      colnames(df)[colnames(df) == "ul"] <- "upr"
      df[, c("time", "lwr", "est", "upr", "sd", "cv", "model")]
    }))
  }

  make_plot <- function(df, ylab_expr, hline = NULL) {
    df$model <- factor(df$model, levels = model_names)
    p <- ggplot(df, aes(x = time, y = est, color = model, fill = model)) +
      geom_line(linewidth = lindwd) +
      scale_color_manual(values = model_colors) +
      scale_fill_manual(values = model_colors) +
      labs(x = "Year", y = ylab_expr) +
      theme_minimal_compact2() +
      guides(color = guide_legend(override.aes = list(linewidth = 1.2)))

    if (show_CIs) {
      p <- p + geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, color = NA, show.legend = FALSE)
    }

    if (!is.null(hline)) {
      p <- p + geom_hline(yintercept = hline, linetype = "dashed", color = "grey30", linewidth = 0.6)
    }

    return(p)
  }

  plots <- list(
    biomass = make_plot(get_series("logB"), "Biomass (tons)"),
    bbmsy   = make_plot(get_series("logBBmsy"), expression(bold(B/B[MSY])), hline = 1),
    ffmsy   = make_plot(get_series("logFFmsy"), expression(bold(F/F[MSY])), hline = 1),
    f       = make_plot(get_series("logF"), "Fishing mortality")
  )

  # Catch
  if (!is.null(extract_catch_data)) {
    catch_df <- bind_rows(lapply(model_names, function(mod) {
      df <- extract_catch_data(models[[mod]], scenario_name = mod)
      df$model <- mod
      df
    }))

    predicted <- filter(catch_df, catch_type == "Predicted")
    observed  <- filter(catch_df, catch_type == "Observed")

    predicted$model <- factor(predicted$model, levels = model_names)
    observed$model  <- factor(observed$model,  levels = model_names)

    # Catch panel
    plots$catch <- ggplot() +
      { if (show_CIs)
        geom_ribbon(
          data = predicted,
          aes(x = time, ymin = lwr, ymax = upr, fill = model),
          alpha = 0.22,
          show.legend = FALSE
        )
      } +
      geom_line(
        data = predicted,
        aes(x = time, y = catch, color = model),
        linewidth = lindwd
      ) +
      geom_point(
        data = observed,
        aes(x = time, y = catch),
        color = "black", size = 1.3,
        show.legend = FALSE
      ) +
      scale_color_manual(values = model_colors) +
      scale_fill_manual(values = model_colors) +
      labs(x = "Year", y = "Catch (tons)") +
      theme_minimal_compact2() +
      guides(color = guide_legend(override.aes = list(linewidth = 1.2)))
  } else {
    plots$catch <- ggplot() + labs(title = "No catch data provided") + theme_void()
  }
  # Production
  if (!is.null(production_fun)) {
    prod_df <- bind_rows(lapply(model_names, function(mod) {
      df <- production_fun(models[[mod]], model_name = mod)
      df$Model <- mod
      return(df)
    }))
    prod_df$Model <- factor(prod_df$Model, levels = model_names)
    max_pts <- prod_df %>% group_by(Model) %>% slice_max(Production, n = 1) %>% ungroup()

    plots$production <- ggplot(prod_df, aes(x = B_K, y = Production, color = Model)) +
      geom_line(linewidth = lindwd, linetype = "dashed") +
      geom_point(data = max_pts, aes(shape = Model), size = 2) +
      scale_color_manual(values = model_colors) +
      scale_shape_manual(values = 15 + seq_along(model_names)) +
      labs(x = expression(bold(B/K)), y = "Production") +
      theme_minimal_compact2()
  } else {
    plots$production <- ggplot() + labs(title = "No production data provided") + theme_void()
  }

  # Patchwork
  if (return_patchwork) {
    layout <- (plots$biomass | plots$bbmsy | plots$catch) /
      patchwork::plot_spacer() /
      (plots$f | plots$ffmsy | plots$production)

    return(layout +
             plot_layout(heights = c(1, 0.05, 1)))
  } else {
    return(plots)
  }
}
