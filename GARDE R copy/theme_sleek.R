#' A sleek, light-based theme for ggplot2
#' @param base_size Base font size
#' @param base_family Base font family
#' @export
theme_sleek <- function(base_size = 11, base_family = "") {
  half_line <- base_size / 2
  ggplot2::theme_light(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.major    = ggplot2::element_blank(),
      panel.grid.minor    = ggplot2::element_blank(),
      axis.ticks.length   = grid::unit(half_line / 2.2, "pt"),
      strip.background    = ggplot2::element_rect(fill = NA, colour = NA),
      strip.text.x        = ggplot2::element_text(colour = "grey20", face = "bold", size = ggplot2::rel(1.2)),
      strip.text.y        = ggplot2::element_text(colour = "grey20", face = "bold", size = ggplot2::rel(1.2)),
      axis.text           = ggplot2::element_text(colour = "grey15", size = 10),
      axis.title          = ggplot2::element_text(colour = "grey15", size = 10),
      panel.border        = ggplot2::element_rect(fill = NA, colour = "grey28", linewidth = 1.8),
      legend.key.size     = grid::unit(0.8, "lines"),
      legend.key.width    = grid::unit(0.8, "lines"),
      legend.key.height   = grid::unit(0.8, "lines"),
      legend.spacing.x    = grid::unit(0,   "pt"),
      legend.box.spacing  = grid::unit(0,   "pt"),
      legend.margin       = ggplot2::margin(0, 0, 0, 0),
      legend.title        = ggplot2::element_blank(),
      legend.text         = ggplot2::element_text(size = 10),
      legend.key          = ggplot2::element_rect(colour = NA, fill = NA),
      legend.background   = ggplot2::element_rect(colour = NA, fill = NA),
      plot.title          = ggplot2::element_text(colour = "grey20", size = 18, face = "bold"),
      plot.subtitle       = ggplot2::element_text(colour = "grey20", size = 18, face = "bold")
    )
}
