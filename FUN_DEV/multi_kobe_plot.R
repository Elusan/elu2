#' Plot a Grid of Kobe Plots for Multiple Models/Scenarios
#'
#' Generates a grid of Kobe plots (biomass vs. fishing mortality relative to MSY reference points) for a set of models, with each plot labeled by its corresponding model. Optionally, saves the grid as a PDF file.
#'
#' @param models A named list of model fit objects. Each should be accepted by \code{elu_plot_kobe()}.
#' @param model_labels Optional. Character vector of labels to use for each model panel. If NULL, uses \code{names(models)}.
#' @param scenario_label Optional. Main grid title placed above all panels.
#' @param ncol Integer. Number of columns in the grid. If NULL, uses one column per model.
#' @param width Numeric. Width (in inches) of each individual Kobe plot panel. Default is 6.
#' @param height Numeric. Height (in inches) of each individual Kobe plot panel. Default is 5.
#' @param save_pdf Optional file path (character). If provided, saves the output grid to this PDF file.
#' @param text_cex Numeric. Text size (cex) for model label inside each plot panel. Default is 1.5.
#' @param ... Additional arguments passed to \code{elu_plot_kobe()}.
#'
#' @details
#' For each model in the input list, this function generates a Kobe plot by calling \code{elu_plot_kobe()} (which must plot to the current device), writes the plot to a temporary PNG, and reads it back as a grob for flexible grid arrangement. Each panel is annotated with its model label in the top-left corner. Optionally, a scenario title can be placed above the grid. If \code{save_pdf} is specified, the grid is saved to the provided PDF path.
#'
#' @return
#' The function is called for its side effect of plotting the grid of Kobe plots to the current graphics device (or to a PDF). Invisibly returns \code{NULL}.
#'
#' @seealso \code{\link{elu_plot_kobe}}
#'
#' @examples
#' \dontrun{
#' multi_kobe_plot(
#'   models = list(Schaefer = fit_S, Fox = fit_F, Pella = fit_P),
#'   model_labels = c("Schaefer", "Fox", "Pella"),
#'   scenario_label = "Scenario 1",
#'   ncol = 3,
#'   save_pdf = "kobe_grid.pdf"
#' )
#' }
#' @export
multi_kobe_plot <- function(
    models,
    model_labels = NULL,        # Named vector for panel model names
    scenario_label = NULL,      # Main grid title
    ncol = NULL,
    width = 6, height = 5,
    save_pdf = NULL,
    text_cex = 1.5,             # Text size for model label inside plot
    ...
) {
  require(gridExtra)
  require(grid)
  require(ggplot2)

  if (is.null(model_labels)) model_labels <- names(models)
  if (is.null(ncol)) ncol <- length(models)

  # Plot a Kobe plot, then annotate with model name
  get_kobe_grob <- function(model, label) {
    temp_file <- tempfile(fileext = ".png")
    png(temp_file, width = width*100, height = height*100, res = 120)
    elu_plot_kobe(model, ...)
    # Add the model label INSIDE the plot (top left)
    usr <- par("usr")
    xleft <- usr[1] + 0.04 * (usr[2] - usr[1])
    ytop  <- usr[4] - 0.07 * (usr[4] - usr[3])
    text(x = xleft, y = ytop, label = label, adj = c(0, 1), font = 2, cex = text_cex, col = "#0B3C5D")
    dev.off()
    grob <- rasterGrob(png::readPNG(temp_file), interpolate = TRUE)
    unlink(temp_file)
    grob
  }

  plots <- mapply(
    FUN = get_kobe_grob,
    model = models,
    label = model_labels,
    SIMPLIFY = FALSE
  )

  # Function to arrange plots
  grid_fun <- function() {
    if (!is.null(scenario_label)) {
      grid.arrange(
        arrangeGrob(grobs = plots, ncol = ncol),
        top = textGrob(scenario_label, gp = gpar(fontsize = 18, fontface = "bold"))
      )
    } else {
      grid.arrange(grobs = plots, ncol = ncol)
    }
  }

  if (!is.null(save_pdf)) {
    pdf(save_pdf, width = ncol * width, height = height * ceiling(length(plots)/ncol))
    grid_fun()
    dev.off()
    message("Kobe plot grid saved to: ", save_pdf)
  } else {
    grid_fun()
  }
}
