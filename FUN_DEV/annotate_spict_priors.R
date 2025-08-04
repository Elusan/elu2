#' @title Annotate a SPiCT priors plot with a model ID
#'
#' @param p         A ggplot object returned by plotspict.priors.ggplot98()
#' @param model_id  Character ID to inset (e.g. "S8P", "S8F", "S8S")
#' @param size      Font size (in pts) for the label. Default 10.
#' @param color     Text color. Default "grey30".
#' @return          The same ggplot object with the model ID drawn in bold
#' @importFrom ggplot2 annotation_custom
#' @importFrom grid textGrob unit gpar
#' @export
annotate_spict_priors <- function(p, model_id, size = 10, color = "grey30") {
  if (!requireNamespace("grid", quietly = TRUE)) {
    stop("Please install the grid package (should come with base R).")
  }
  p +
    annotation_custom(
      grob = grid::textGrob(
        label = model_id,
        x     = grid::unit(0.02, "npc"),
        y     = grid::unit(0.98, "npc"),
        just  = c("left", "top"),
        gp    = grid::gpar(fontface = "bold", col = color, fontsize = size)
      ),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    )
}
