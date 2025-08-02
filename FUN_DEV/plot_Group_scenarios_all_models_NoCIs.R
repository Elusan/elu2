#' Plot SPiCT Model Outputs (No Confidence Intervals) for Multiple Scenarios
#'
#' Generates a 2x3 grid of plots (biomass, B/Bmsy, F/Fmsy, fishing mortality, catch, production)
#' from a named list of SPiCT/ELU model fits. Each scenario is uniquely colored.
#' Only point estimates are plotted (no uncertainty ribbons).
#'
#' @param models A named list of fitted SPiCT model objects.
#'   Each name should uniquely identify a model-scenario combination (e.g., `"S1F"`, `"S2P"`).
#' @param production_fun Optional. A function of the form \code{function(model, model_name)}
#'   that returns a data frame with columns for production plotting (e.g., `B_K`, `Production`, `Model`).
#' @param extract_catch_data Optional. A function of the form \code{function(model, scenario_name)}
#'   that returns a data frame with columns: `time`, `catch`, `catch_type` (`"Observed"` or `"Predicted"`),
#'   and `scenario`.
#' @param return_patchwork Logical; if \code{TRUE} (default), returns a combined patchwork plot;
#'   if \code{FALSE}, returns a named list of individual ggplot2 plot objects.
#' @param lindwd Numeric. Line width used for time series lines. Default is \code{1.4}.
#' @param scenario_colors Optional named vector assigning unique colors to each scenario name.
#'   If not provided, colors are generated using \code{scales::hue_pal()}.
#'
#' @details
#' This function is a generalized version of \code{plot_spict_scenarios_by_model()} that
#' supports any number of scenario-model combinations (e.g., `"S1F"`, `"S2P"`, `"S3S"`, etc.).
#'
#' It creates six panels:
#' \itemize{
#'   \item Biomass (\code{logB})
#'   \item Relative biomass (\eqn{B/B_{MSY}})
#'   \item Relative fishing mortality (\eqn{F/F_{MSY}})
#'   \item Fishing mortality (\code{logF})
#'   \item Catch (observed + predicted, if provided)
#'   \item Production (if provided)
#' }
#'
#' @return
#' A patchwork plot (if \code{return_patchwork = TRUE}) or a named list of ggplot2 plots.
#'
#' @seealso \code{\link{get_par}}, \code{\link{get_production_data}}, \code{\link{elu_extract_catch_data}}
#'
#' @examples
#' \dontrun{
#' models <- list(
#'   S1P = S1_Pella, S1S = S1_Schaefer, S1F = S1_Fox,
#'   S2P = S2_Pella, S2S = S2_Schaefer, S2F = S2_Fox
#' )
#' scenario_colors <- c(
#'   "S1P" = "#1b9e77", "S1S" = "#d95f02", "S1F" = "#7570b3",
#'   "S2P" = "#e7298a", "S2S" = "#66a61e", "S2F" = "#e6ab02"
#' )
#' plot_spict_scenarios_by_model_NoCIs(
#'   models = models,
#'   production_fun = get_production_data,
#'   extract_catch_data = elu_extract_catch_data,
#'   scenario_colors = scenario_colors
#' )
#' }
#' @export
plot_Group_scenarios_all_models_NoCIs <- function(models,
                                                  production_fun = NULL,
                                                  extract_catch_data = NULL,
                                                  return_patchwork = TRUE,
                                                  lindwd = 0.7,
                                                  scenario_colors = NULL,
                                                  save_path = NULL,
                                                  save_dir = NULL,        # <--- NEW
                                                  save_width = 20,
                                                  save_height = 10,
                                                  save_dpi = 300) {

  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(grid)

  # Use provided colors or default fallback
  if (is.null(scenario_colors)) {
    all_models <- names(models)
    default_palette <- scales::hue_pal()(length(all_models))
    names(default_palette) <- all_models
    scenario_colors <- default_palette
  }

  valid_shapes <- rep(0:25, length.out = length(scenario_colors))

  # Compact theme

  theme_minimal_compactSPiCT <- function(base_size = 18, base_family = "") {
    ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
      theme(
        plot.title = ggtext::element_markdown(face = "bold", size = 18, hjust = 0.5),
        axis.title = element_text(face = "bold", size = 18),
        axis.text = element_text(size = 18, face = "bold"),
        axis.title.y.right = element_text(color = "#DC143C", size = 18, face = "bold"),
        axis.text.y.right  = element_text(color = "#DC143C", size = 18, face = "bold"),
        legend.position = c(0.99, 0.98),
        legend.justification = c("right", "top"),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.spacing.x = unit(0.5, "pt"),
        legend.spacing.y = unit(0.7, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-5, -5, -5, -5),
        legend.key = element_rect(fill = "white", color = NA),
        legend.key.height = unit(0.2, "line"),
        legend.key.width = unit(0.8, "line"),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, colour = "grey45", linewidth = 2.2),
        axis.ticks = element_line(linewidth = 0.5, color = "grey45"),
        axis.ticks.length = unit(2, "pt"),
        text = element_text(face = "bold", size = 18),
        plot.margin = margin(4, 4, 4, 4)
      )
  }



  get_series <- function(parname) {
    do.call(rbind, lapply(names(models), function(modname) {
      par <- get.par(parname, models[[modname]], exp = TRUE)
      df <- as.data.frame(par)
      df$time <- as.numeric(rownames(par))
      df$model <- modname
      df[, c("time", "est", "model")]
    }))
  }

  make_plot <- function(df, ylab_expr, hline = NULL) {
    df$model <- factor(df$model, levels = names(scenario_colors))
    gg <- ggplot(df, aes(x = time, y = est, color = model)) +
      geom_line(linewidth = lindwd) +
      scale_color_manual(values = scenario_colors, drop = FALSE) +
      labs(x = "Year", y = ylab_expr, color = "Model") +
      theme_minimal_compactSPiCT() +
      guides(color = guide_legend(override.aes = list(linewidth = 2)))
    if (!is.null(hline)) {
      gg <- gg + geom_hline(yintercept = hline, linetype = "dashed", color = "black", linewidth = 0.8)
    }
    return(gg)
  }

  plots <- list(
    biomass = make_plot(get_series("logB"), "Biomass (tons)"),
    bbmsy   = make_plot(get_series("logBBmsy"), expression(bold(B/B[MSY])), hline = 1),
    ffmsy   = make_plot(get_series("logFFmsy"), expression(bold(F/F[MSY])), hline = 1),
    f       = make_plot(get_series("logF"), "Fishing mortality")
  )

  # Catch panel
  if (!is.null(extract_catch_data)) {
    catch_all <- dplyr::bind_rows(
      lapply(names(models), function(mod) {
        extract_catch_data(models[[mod]], scenario_name = mod)
      })
    )
    catch_all$model <- catch_all$scenario
    predicted <- filter(catch_all, catch_type == "Predicted")
    observed  <- filter(catch_all, catch_type == "Observed")
    predicted$model <- factor(predicted$model, levels = names(scenario_colors))
    observed$model  <- factor(observed$model,  levels = names(scenario_colors))

    plots$catch <- ggplot() +
      geom_line(data = predicted, aes(x = time, y = catch, color = model), size = lindwd) +
      geom_point(data = observed, aes(x = time, y = catch), color = "black", size = 1.2) +
      scale_color_manual(values = scenario_colors, drop = FALSE) +
      labs(x = "Year", y = "Catch (tons)", color = "Model") +
      theme_minimal_compactSPiCT() +
      guides(color = guide_legend(override.aes = list(linewidth = 2)))
  } else {
    plots$catch <- ggplot() + labs(title = "No extract_catch_data function provided") + theme_void()
  }


  # Production panel
  if (!is.null(production_fun)) {
    prod_df <- do.call(rbind, lapply(names(models), function(modname) {
      production_fun(models[[modname]], model_name = modname)
    }))
    prod_df$Model <- factor(prod_df$Model, levels = names(scenario_colors))
    max_points <- prod_df %>%
      group_by(Model) %>%
      slice_max(Production, n = 1) %>%
      ungroup()
    plots$production <- ggplot(prod_df, aes(x = B_K, y = Production, color = Model)) +
      geom_line(size = lindwd, linetype = "dashed") +
      geom_point(data = max_points, aes(shape = Model), size = 3)+
      scale_shape_manual(values = valid_shapes)+
      scale_color_manual(values = scenario_colors) +
      labs(x = expression(bold(B/K)), y = "Production", color = "Model", shape = "Model") +
      theme_minimal_compactSPiCT()
  } else {
    plots$production <- ggplot() + labs(title = "No production_fun provided") + theme_void()
  }

  if (return_patchwork) {
    final_plot <- (plots$biomass | plots$bbmsy | plots$catch) /
      patchwork::plot_spacer() /
      (plots$f | plots$ffmsy | plots$production) +
      plot_layout(heights = c(1, 0.05, 1)) &
      theme(plot.margin = margin(4, 4, 4, 4))

    # Save combined plot
    if (!is.null(save_path)) {
      ggsave(save_path, plot = final_plot,
             width = save_width, height = save_height, dpi = save_dpi)
    }

    if (!is.null(save_dir)) {
      if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
      for (plot_name in names(plots)) {
        plot_file <- file.path(save_dir, paste0(plot_name, ".png"))
        p_clean <- plots[[plot_name]] +
          theme_minimal_compactSPiCT() +
          theme(
            panel.background = element_rect(fill = "white", color = NA),
            plot.margin = margin(4, 4, 4, 4)
          )
        ggsave(plot_file, p_clean, width = save_width, height = save_height / 3, dpi = save_dpi)
      }
    }



    return(final_plot)
  } else {
    return(plots)
  }

}
