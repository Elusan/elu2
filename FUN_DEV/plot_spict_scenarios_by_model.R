#' Plot SPiCT Scenarios by Model (Fox, Schaefer, Pella)
#'
#' Generates a grid of plots summarizing key model outputs (biomass, B/Bmsy, F/Fmsy, F, catch, and production curves) for multiple SPiCT model fits (e.g., Fox, Schaefer, Pella-Tomlinson), with optional catch and production extraction functions. Returns a patchwork layout or a list of ggplot objects.
#'
#' @param models A named list of fitted SPiCT model outputs, typically with names such as "S1F", "S1P", "S1S".
#' @param production_fun Optional. A function to extract production curve data for each model. Should return a data frame with columns including `B_K`, `Production`, and `Model`.
#' @param extract_catch_data Optional. A function to extract observed and predicted catch (with uncertainty) for each model. Should return a data frame with columns including `time`, `catch`, `lwr`, `upr`, `catch_type`, and `scenario`.
#' @param return_patchwork Logical. If `TRUE` (default), returns a patchwork grid of plots. If `FALSE`, returns a named list of ggplot objects.
#' @param lindwd Numeric. Line width for main plot lines (default: 0.8).
#'
#' @details
#' This function visualizes and compares the outputs from multiple surplus production models fitted using SPiCT or compatible workflows. It supports custom color schemes and theming, and overlays confidence ribbons for uncertainty.
#'
#' - **Biomass plot:** Time series of estimated biomass with uncertainty by model.
#' - **B/Bmsy and F/Fmsy plots:** Time series with reference line at 1.
#' - **Catch plot:** Observed (black points) and predicted (colored lines and ribbons) catch.
#' - **Production plot:** Equilibrium production curve with maxima indicated per model.
#'
#' If `production_fun` or `extract_catch_data` are not provided, the corresponding plots will show placeholder messages.
#'
#' @return
#' Either a patchwork object (default), or a named list of ggplot objects (`biomass`, `bbmsy`, `ffmsy`, `f`, `catch`, `production`).
#'
#' @seealso [patchwork::plot_layout()], [ggplot2::ggplot()]
#'
#' @examples
#' \dontrun{
#'   # Example usage with three fitted models
#'   plots <- plot_spict_scenarios_by_model(
#'     models = list(S1F = fox_fit, S1P = pella_fit, S1S = schaefer_fit),
#'     production_fun = extract_production,
#'     extract_catch_data = extract_catch,
#'     return_patchwork = TRUE
#'   )
#'   plots # patchwork grid
#' }
#' @export
plot_spict_scenarios_by_model <- function(models,
                                          production_fun = NULL,
                                          extract_catch_data = NULL,
                                          return_patchwork = TRUE,
                                          lindwd = 0.7) {
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(grid)

  # model_colors <- c(
  #   "Fox" = "#2B83BA",
  #   "Pella" = "#D7191C",
  #   "Schaefer" = "#1A9641"
  # )

  if (!"scenario" %in% colnames(predicted)) {
    predicted$scenario <- names(models)[i]
  } else {
    predicted$scenario <- as.character(predicted$scenario)
  }






  model_colors <- c(
    "S1F" = "#2B83BA",
    "S1P" = "#D7191C",
    "S1S" = "#1A9641"
  )


  # get_base_model <- function(modname) {
  #   suffix <- substr(modname, nchar(modname), nchar(modname))
  #   if (suffix == "P") return("Pella")
  #   if (suffix == "S") return("Schaefer")
  #   if (suffix == "F") return("Fox")
  #   return(NA)
  # }

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
    do.call(rbind, lapply(names(models), function(mod) {
      par <- get.par(parname, models[[mod]], exp = TRUE)
      df <- as.data.frame(par)
      df$time <- as.numeric(rownames(par))
      colnames(df)[colnames(df) == "ll"] <- "lwr"
      colnames(df)[colnames(df) == "ul"] <- "upr"
      df <- df[, c("time", "lwr", "est", "upr", "sd", "cv")]
      df$model <- get_base_model(mod)
      return(df)
    }))
  }

  make_plot <- function(df, ylab_expr, hline = NULL) {
    df$model <- factor(df$model, levels = names(model_colors))
    p <- ggplot(df, aes(x = time, y = est, color = model, fill = model)) +
      geom_ribbon(aes(ymin = lwr, ymax = upr, fill = model), alpha = 0.22, color = NA, show.legend = FALSE) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(values = model_colors, drop = FALSE) +
      scale_fill_manual(values = model_colors, drop = FALSE) +
      labs(x = "Year", y = ylab_expr, color = "Model", fill = "Model") +
      theme_minimal_compact2() +
      guides(color = guide_legend(override.aes = list(linewidth = 1.8)))
    if (!is.null(hline)) {
      p <- p + geom_hline(yintercept = hline, linetype = "dashed", color = "black", linewidth = 0.8)
    }
    return(p)
  }


  # make_plot <- function(df, ylab_expr, hline = NULL) {
  #   df$model <- factor(df$model, levels = names(model_colors))
  #   p <- ggplot(df, aes(x = time, y = est, color = model, fill = model)) +
  #     geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.22, color = NA) +
  #     geom_line(linewidth = lindwd) +
  #     scale_color_manual(values = model_colors, drop = FALSE) +
  #     scale_fill_manual(values = model_colors, drop = FALSE) +
  #     labs(x = "Year", y = ylab_expr, color = "Model", fill = "Model") +
  #     theme_minimal_compact()
  #   if (!is.null(hline)) {
  #     p <- p + geom_hline(yintercept = hline, linetype = "dashed", color = "black", linewidth = 0.8)
  #   }
  #   return(p)
  # }

  plots <- list(
    biomass = make_plot(get_series("logB"), "Biomass (tons)"),
    bbmsy   = make_plot(get_series("logBBmsy"), expression(bold(B/B[MSY])), hline = 1),
    ffmsy   = make_plot(get_series("logFFmsy"), expression(bold(F/F[MSY])), hline = 1),
    f       = make_plot(get_series("logF"), "Fishing mortality")
  )

  # Catch plot
  if (!is.null(extract_catch_data)) {
    catch_all <- dplyr::bind_rows(
      lapply(names(models), function(mod) {
        extract_catch_data(models[[mod]], scenario_name = mod)
      })
    )
    predicted <- filter(catch_all, catch_type == "Predicted")
    observed  <- filter(catch_all, catch_type == "Observed")
    # --- FIX: Vectorize get_base_model over the scenario column ---
    predicted$model <- vapply(predicted$scenario, get_base_model, character(1))
    observed$model  <- vapply(observed$scenario,  get_base_model, character(1))
    predicted$model <- factor(predicted$model, levels = names(model_colors))
    observed$model  <- factor(observed$model,  levels = names(model_colors))

    # plots$catch <- ggplot() +
    #   geom_ribbon(data = predicted, aes(x = time, ymin = lwr, ymax = upr, fill = model), alpha = 0.22) +
    #   geom_line(data = predicted, aes(x = time, y = catch, color = model), size = lindwd) +
    #   geom_point(data = observed, aes(x = time, y = catch), color = "black", size = 1.2) +
    #   scale_color_manual(values = model_colors, drop = FALSE) +
    #   scale_fill_manual(values = model_colors, drop = FALSE) +
    #   labs(x = "Year", y = "Catch (tons)", color = "Model", fill = "Model") +
    #   theme_minimal_compact()

    plots$catch <- ggplot() +
      geom_ribbon(data = predicted, aes(x = time, ymin = lwr, ymax = upr, fill = model), alpha = 0.22, show.legend = FALSE) +
      geom_line(data = predicted, aes(x = time, y = catch, color = model), size = 1.4) +
      #geom_ribbon(data = predicted, aes(x = time, ymin = lwr, ymax = upr, fill = model), alpha = 0.22, show.legend = FALSE) +
      #geom_line(data = predicted, aes(x = time, y = catch, color = model), size = 1.4) +
      geom_point(data = observed, aes(x = time, y = catch), color = "black", size = 1.2) +
      scale_color_manual(values = model_colors, drop = FALSE) +
      labs(x = "Year", y = "Catch (tons)", color = "Model") +
      theme_minimal_compact2() +
      guides(fill = "none", color = guide_legend(override.aes = list(linewidth = 2)))
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
