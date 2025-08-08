#' Plot and Export a Grid of Kobe Plots for Multiple Models/Scenarios
#'
#' Generates a visually enhanced grid of Kobe plots (biomass vs. fishing mortality relative to MSY reference points) for a set of models or scenarios. For each model, an individual Kobe plot is saved as a high-resolution PNG in a dedicated directory. The full grid can also be saved as a PDF and/or PNG file. All plots are annotated with scenario/model labels for clear identification and publication-quality visuals.
#'
#' @param models A named list of model fit objects. Each object should be accepted by \code{elu_plot_kobe()}.
#' @param model_labels Optional. Character vector of labels to use for each model panel. If \code{NULL}, uses \code{names(models)}.
#' @param scenario_label Optional. Title placed above the Kobe grid (e.g., scenario name).
#' @param ncol Integer. Number of columns in the grid. If \code{NULL}, uses one column per model.
#' @param width Numeric. Width (in inches) of each individual Kobe plot panel. Default is 8.
#' @param height Numeric. Height (in inches) of each individual Kobe plot panel. Default is 5.
#' @param save_pdf Optional. File path (character) to save the Kobe plot grid as a PDF. If \code{NULL}, grid is not saved as PDF.
#' @param save_png Optional. File path (character) to save the Kobe plot grid as a PNG. If \code{NULL}, grid is not saved as PNG.
#' @param save_panel_pngs Logical. If \code{TRUE} (default), saves each individual Kobe plot panel as a PNG file in \code{panel_dir}.
#' @param panel_dir Character. Directory to save individual Kobe panel PNGs (created if it does not exist). Default is \code{"kobe_panels"}.
#' @param text_cex Numeric. Text size (cex) for model label inside each plot panel. Default is 2 (large and clear).
#' @param ... Additional arguments passed to \code{elu_plot_kobe()} for plot customization.
#'
#' @details
#' For each model in the input list, this function:
#' \itemize{
#'   \item Creates a Kobe plot using \code{elu_plot_kobe()}, saves it as a high-resolution PNG (panel), and reads it back as a grob for flexible grid arrangement.
#'   \item Annotates each plot with the corresponding model label (top-left, bold, colored).
#'   \item Optionally saves the entire grid of Kobe plots as a PDF and/or PNG, with a scenario label at the top.
#' }
#' The function ensures crisp, publication-quality images suitable for use in reports, presentations, or manuscripts.
#'
#' @return
#' The function is called for its side effects (plotting and saving image files). Invisibly returns \code{NULL}.
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
#'   save_pdf = "kobe_grid.pdf",
#'   save_png = "kobe_grid.png",
#'   save_panel_pngs = TRUE,
#'   panel_dir = "KOBE_PNG_SCENARIO1",
#'   text_cex = 2
#' )
#' }
#'
#' @export
multi_kobe_plot2 <- function(
    models,
    model_labels = NULL,
    scenario_label = NULL,
    ncol = NULL,
    width = 8, height = 5,
    save_pdf = NULL,
    save_png = NULL,           # New: Save grid as PNG
    save_panel_pngs = TRUE,    # New: Save each panel as PNG?
    panel_dir = "kobe_panels", # New: Output dir for panels
    text_cex = 2,              # Larger for clarity
    ...
) {
  require(gridExtra)
  require(grid)
  require(ggplot2)
  require(png)

  if (is.null(model_labels)) model_labels <- names(models)
  if (is.null(ncol)) ncol <- length(models)

  # Create directory for individual panel PNGs
  if (save_panel_pngs) {
    if (!dir.exists(panel_dir)) dir.create(panel_dir, recursive = TRUE)
  }

  # Function to create, save, and return a grob for each model
  get_kobe_grob <- function(model, label, model_name) {
    temp_file <- tempfile(fileext = ".png")
    out_file  <- file.path(panel_dir, paste0(model_name, "_kobe.png"))

    # Create Kobe plot in temp file (for grid) and permanent file (for user)
    for (f in c(temp_file, if (save_panel_pngs) out_file)) {
      png(f, width = width*150, height = height*150, res = 150, bg = "transparent")
      par(mar = c(4.5, 4.5, 1.5, 1)) # Clean, tight margins
      elu_plot_kobe(model, ...)
      usr <- par("usr")
      xleft <- usr[1] + 0.04 * (usr[2] - usr[1])
      ytop  <- usr[4] - 0.07 * (usr[4] - usr[3])
      text(
        x = xleft, y = ytop, label = label, adj = c(0, 1),
        font = 2, cex = text_cex, col = "#0B3C5D"
      )
      dev.off()
    }
    grob <- rasterGrob(png::readPNG(temp_file), interpolate = TRUE)
    unlink(temp_file)
    grob
  }

  # Build grobs, save individual PNGs as a side effect
  plots <- mapply(
    FUN = get_kobe_grob,
    model = models,
    label = model_labels,
    model_name = names(models),
    SIMPLIFY = FALSE
  )

  # Arrange grid function
  grid_fun <- function() {
    if (!is.null(scenario_label)) {
      grid.arrange(
        arrangeGrob(grobs = plots, ncol = ncol),
        top = textGrob(scenario_label, gp = gpar(fontsize = 22, fontface = "bold", col = "#1a1a1a"))
      )
    } else {
      grid.arrange(grobs = plots, ncol = ncol)
    }
  }

  # Save grid as PDF
  if (!is.null(save_pdf)) {
    pdf(save_pdf, width = ncol * width, height = height * ceiling(length(plots)/ncol))
    grid_fun()
    dev.off()
    message("Kobe plot grid saved to: ", save_pdf)
  }
  # Save grid as PNG
  if (!is.null(save_png)) {
    png(save_png, width = ncol * width * 150, height = height * ceiling(length(plots)/ncol) * 150, res = 150)
    grid_fun()
    dev.off()
    message("Kobe plot grid saved to: ", save_png)
  }

  # Also plot to active device
  grid_fun()

  invisible(NULL)
}
