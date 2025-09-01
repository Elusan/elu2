# R/theme_minimal_compact2_hind.R

#' Minimal, compact ggplot2 theme for SPiCT hindcast panels
#'
#' A clean, compact variant of [ggplot2::theme_minimal()] tailored for
#' SPiCT-style hindcast plots. It removes facet strip backgrounds *and* strip
#' labels (the shaded header bar and any numeric strip text), keeps a crisp
#' panel border, and places a compact legend inside the panel.
#'
#' @param base_size Numeric. Base font size (points). Default `10`.
#' @param base_family Character. Base font family. Default `""`.
#'
#' @details
#' Key choices:
#' * `strip.background = element_blank()` and `strip.text = element_blank()` remove the shaded bar and its label.
#' * Strong `panel.border` while `panel.grid` is removed.
#' * Legend is inside the panel (`c(0.98, 0.98)`) with a compact white box.
#'
#' @return A [ggplot2::theme()] object to be added to a plot.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) +
#'   geom_point() +
#'   facet_wrap(~ cyl) +
#'   theme_minimal_compact2_hind()
#' }
#'
#' @export
#' @import ggplot2
#' @importFrom grid unit
theme_minimal_compact2_hind <- function(base_size = 10, base_family = "") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(face = "bold", size = 10),
      axis.text  = ggplot2::element_text(size = 10, face = "bold"),
      # Legend: keep in-panel but make the *frame* smaller/tighter
      legend.position = c(0.98, 0.98),
      legend.justification = c("right", "top"),
      legend.background = ggplot2::element_rect(fill = "white", color = "grey50", linewidth = 0.4),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      legend.margin = ggplot2::margin(1, 1, 1, 1),
      legend.spacing.x = grid::unit(0, "pt"),
      legend.spacing.y = grid::unit(0.2, "pt"),
      legend.title = ggplot2::element_blank(),
      legend.text  = ggplot2::element_text(size = 8, face = "bold"),
      legend.key.size = grid::unit(1.2, "lines"),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(fill = NA, colour = "grey35", linewidth = 0.85),
      axis.ticks = ggplot2::element_line(linewidth = 0.7, color = "grey35"),
      axis.ticks.length = grid::unit(3, "pt"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(3, 3, 3, 3)
    )
}
