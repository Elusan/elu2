#' Plot Catch and Standardized Indices (Long and Short) as Three Panels
#'
#' This function generates a 3-panel plot combining catch and standardized indices
#' for: (1) long delta-GLM index (1992–2016), (2) short delta-GLM index (2017–2022),
#' and (3) sdmTMB index (2017–2022). All required data are available internally
#' via the package as `Data_glm_short`, `Data_glm_full`, and `Data_sdm`.
#'
#' @return A patchwork plot object combining three ggplot2 panels, and saves the
#' figure to "FIG/plots_CATCH_INDEX_2x2.png"
#' @export
plot_catch_index_panels <- function() {
  stopifnot(exists("Data_glm_short"), exists("Data_sdm"))

  library(ggplot2)
  library(grid)
  library(scales)
  library(patchwork)

  theme_minimal_compact <- function(base_size = 8, base_family = "") {
    theme_minimal(base_size = base_size, base_family = base_family) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 14, face = "bold"),
        axis.title.y.right = element_text(color = "#DC143C", size = 14, face = "bold"),
        axis.text.y.right  = element_text(color = "#DC143C", size = 14, face = "bold"),
        legend.position = c(0.99, 0.98),
        legend.justification = c("right", "top"),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.key.size = unit(1, "lines"),
        legend.spacing.y = unit(0, "pt"),
        legend.spacing.x = unit(0, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, colour = "grey45", linewidth = 2.5),
        axis.ticks = element_line(linewidth = 0.5, color = "grey45"),
        axis.ticks.length = unit(2, "pt"),
        strip.background = element_rect(fill = "grey45", color = "grey45", linewidth = 0.8),
        strip.text = element_text(face = "bold", size = rel(1)),
        text = element_text(face = "bold", size = 8),
        plot.margin = margin(0, 0, 0, 0)
      )
  }

  make_panel_plot <- function(Data, year_filter, annotation_label, annotation_year, is_glm = TRUE) {
    max_catch <- max(Data$Catch, na.rm = TRUE)
    max_index <- max(if (is_glm) Data$Index_stand else Data$Index, na.rm = TRUE)
    scale_factor <- max_catch / max_index
    x_min <- min(Data$Year, na.rm = TRUE)
    x_max <- max(Data$Year, na.rm = TRUE)

    # Prepare plotting columns
    Data$Catch_plot <- ifelse(year_filter(Data$Year), Data$Catch, NA)
    Data$Index_plot <- if (is_glm) {
      ifelse(year_filter(Data$Year), Data$Index_stand, NA)
    } else {
      ifelse(year_filter(Data$Year), Data$Index, NA)
    }

    # Filtered subsets for plotting (no NA)
    Catch_data <- subset(Data, !is.na(Catch_plot))
    Index_data <- subset(Data, !is.na(Index_plot))

    ggplot() +
      geom_bar(data = Catch_data, aes(x = Year, y = Catch_plot, fill = "Catch"), stat = "identity", alpha = 0.7) +
      geom_line(data = Index_data, aes(x = Year, y = Index_plot * scale_factor), color = "#fcaeae", size = 0.8) +
      geom_point(data = Index_data, aes(x = Year, y = Index_plot * scale_factor, color = "Index"),
                 size = 1.2, shape = 19, stroke = 1.2) +
      geom_vline(xintercept = 2016.5, linetype = "dashed", color = "black", linewidth = 0.4) +
      annotate("text", x = annotation_year, y = max_catch * 0.7,
               label = annotation_label, size = 3, fontface = "bold") +
      scale_y_continuous(
        name = "Catch (tons)",
        sec.axis = sec_axis(~ . / scale_factor, name = "Standardized Index"),
        limits = c(0, max_catch)
      ) +
      scale_fill_manual(values = c("Catch" = "steelblue")) +
      scale_color_manual(values = c("Index" = "#DC143C")) +
      coord_cartesian(xlim = c(x_min, x_max)) +
      guides(fill = guide_legend(order = 1), color = guide_legend(order = 2)) +
      theme_minimal_compact(base_size = 8) +
      theme(
        axis.title.y.left = element_text(color = "steelblue", face = "bold", size = 12),
        axis.title.y.right = element_text(color = "#DC143C", face = "bold", size = 12),
        axis.text.y.left = element_text(color = "steelblue", face = "bold", size = 12),
        axis.text.y.right = element_text(color = "#DC143C", face = "bold", size = 12),
        axis.text.x = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 12),
        legend.position = c(0.99, 0.97),
        legend.justification = c("right", "top"),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.text = element_text(size = 8),
        legend.key = element_rect(fill = "white", color = NA),
        legend.key.height = unit(0.2, "line"),
        legend.key.width = unit(0.8, "line"),
        panel.border = element_rect(fill = NA, color = "grey45", linewidth = 1.8),
        axis.ticks = element_line(linewidth = 0.8, color = "grey45"),
        axis.ticks.length = unit(3, "pt"),
        text = element_text(face = "bold")
      ) +
      labs(x = "Year")
  }

  # ======================= PLOT PANELS =======================
  Q1 <- make_panel_plot(Data_glm_short, function(y) y <= 2016,
                        "Long delta-GLM \nIndex (1992–2016)", 2004)

  Q2 <- make_panel_plot(Data_glm_full, function(y) y >= 2017,
                        "Short delta-GLM \nIndex (2017–2022)", 2019.5)

  Q3 <- make_panel_plot(Data_sdm, function(y) y >= 2017,
                        "Short sdmTMB \nIndex (2017–2022)", 2019.5, is_glm = FALSE)

  g2 <- Q1 / Q2 / Q3

  if (!dir.exists("FIG")) dir.create("FIG")
  ggsave("FIG/plots_CATCH_INDEX_2x2.png", g2, width = 12, height = 12, dpi = 300)

  return(g2)
}
