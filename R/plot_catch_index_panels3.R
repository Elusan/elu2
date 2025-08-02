#' @title Plot Catch and Index Time Series Panels
#' @description
#' Generates vertically stacked panel plots showing observed catch (as histograms)
#' and standardized index values (as lines/points) for each index in a provided data list.
#' Each plot panel includes a markdown title with sample size summaries and index year range,
#' consistent y-axis scaling, and a shared x-axis spanning all years in the catch time series.
#'
#' @param data_list A named list with elements:
#' \describe{
#'   \item{timeC}{Numeric vector of years for observed catch}
#'   \item{obsC}{Numeric vector of observed catch values}
#'   \item{timeI}{A list of numeric vectors, each containing years for a standardized index}
#'   \item{obsI}{A list of numeric vectors, each containing values for the corresponding index}
#' }
#'
#' @param output_dir Character string specifying the folder where the PNG plot will be saved.
#'   Defaults to `"FIG"`. If the directory does not exist, it will be created.
#'
#' @return A `patchwork` object representing the combined plot. Also saves the plot to a PNG file.
#'
#' @details
#' For each index, the function:
#' \itemize{
#'   \item Restricts the display of catch bars to the range of years covered by that index
#'   \item Plots index values as red lines and points (on a scaled right y-axis)
#'   \item Sets a common y-axis range and scaling factor across all panels
#'   \item Adds rich plot titles using `ggtext::element_markdown()` with sample sizes and year ranges
#' }
#' Panel spacing and legend positioning are customized for publication-quality output.
#'
#' @import ggplot2
#' @import patchwork
#' @import ggtext
#' @export
plot_catch_index_panels3 <- function(data_list, output_dir = "FIG") {
  stopifnot(all(c("timeC", "obsC", "timeI", "obsI") %in% names(data_list)))

  library(ggplot2)
  library(patchwork)
  library(ggtext)

  # Custom compact theme with markdown for titles
  theme_minimal_compactCI <- function(base_size = 18, base_family = "") {
    theme_minimal(base_size = base_size, base_family = base_family) +
      theme(
        plot.title = element_markdown(face = "bold", size = 18, hjust = 0.5),
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
        plot.margin = margin(0, 0, 0, 0) #margin(top, right, bottom, left)
      )
  }

  # Shared global range and scale
  global_x_range <- range(data_list$timeC, na.rm = TRUE)
  global_max_catch <- max(data_list$obsC, na.rm = TRUE)
  global_max_index <- max(unlist(data_list$obsI), na.rm = TRUE)
  global_scale_factor <- global_max_catch / global_max_index

  # Plot builder
  make_panel_plot <- function(yearC, obsC, yearI, obsI, index_number) {
    all_years <- seq(global_x_range[1], global_x_range[2])

    df_catch <- data.frame(Year = yearC, Catch = obsC)
    df_catch_full <- merge(data.frame(Year = all_years), df_catch, by = "Year", all.x = TRUE)
    catch_range <- range(yearI, na.rm = TRUE)
    df_catch_full$Catch[!(df_catch_full$Year >= catch_range[1] & df_catch_full$Year <= catch_range[2])] <- NA

    df_index <- data.frame(Year = yearI, Index = obsI)
    df_index_full <- merge(data.frame(Year = all_years), df_index, by = "Year", all.x = TRUE)
    df_plot <- merge(df_catch_full, df_index_full, by = "Year", all = TRUE)

    nC <- sum(!is.na(df_plot$Catch))
    nI <- sum(!is.na(df_plot$Index))
    year_min <- min(yearI, na.rm = TRUE)
    year_max <- max(yearI, na.rm = TRUE)

    title_html <- sprintf(
      "NobsC = <b>%d</b> — <span style='color:#DC143C;'>Index %d = <b>%d</b> (%d–%d)</span>",
      nC, index_number, nI, year_min, year_max
    )

    ggplot(df_plot, aes(x = Year)) +
      geom_bar(data = subset(df_plot, !is.na(Catch)),
               aes(y = Catch, fill = "Catch"), stat = "identity", alpha = 0.7) +
      geom_line(data = subset(df_plot, !is.na(Index)),
                aes(y = Index * global_scale_factor), color = "#fcaeae", size = 0.9) +
      geom_point(data = subset(df_plot, !is.na(Index)),
                 aes(y = Index * global_scale_factor, color = "Index"), size = 4, shape = 19) +
      scale_y_continuous(
        name = "Catch (tons)",
        sec.axis = sec_axis(~ . / global_scale_factor, name = "Standardized Index"),
        limits = c(0, global_max_catch)
      ) +
      scale_fill_manual(values = c("Catch" = "steelblue")) +
      scale_color_manual(values = c("Index" = "#DC143C")) +
      coord_cartesian(xlim = global_x_range) +
      theme_minimal_compactCI() +
      labs(x = "Year", title = title_html) +
      guides(fill = guide_legend(order = 1), color = guide_legend(order = 2))
  }

  # Generate all plots
  plots <- lapply(seq_along(data_list$timeI), function(i) {
    make_panel_plot(
      yearC = data_list$timeC,
      obsC  = data_list$obsC,
      yearI = data_list$timeI[[i]],
      obsI  = data_list$obsI[[i]],
      index_number = i
    )
  })

  # Combine and save
  plots <- lapply(plots, function(p) p + theme(plot.margin = margin(10, 10, 30, 10)))

  g_all <- wrap_plots(plots, ncol = 1)

  g_all <- wrap_plots(plots, ncol = 1)


  if (!dir.exists(output_dir)) dir.create(output_dir)
  ggsave(file.path(output_dir, "plots_CATCH_INDEX_dynamic.png"),
         g_all, width = 14.5, height = 5 * length(plots), dpi = 300)

  return(g_all)
}
