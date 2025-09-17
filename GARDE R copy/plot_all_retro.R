#' Generate a set of peel colors
#'
#' Produces a vector of peel colors based on the JABBA color scheme, including
#' fully opaque, semi-transparent (50%), and faint (20%) versions of each color
#' (excluding black). This function is intended for consistent coloring of peels
#' in retrospective plots.
#'
#' @return A character vector of hex color codes with transparency adjustments.
#'         The vector includes the base colors and their adjusted versions.
#'
#' @examples
#' peel_colors <- cols()
#' print(peel_colors)
#'
#' @export
cols <- function() {
  cs <- c("#000000", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
  c(cs, adjustcolor(cs[-1], 0.5), adjustcolor(cs[-1], 0.2))
}


#' Batch plot retrospective diagnostics for multiple SPiCT models
#'
#' Loops over a named list of fitted SPiCT models with retro components and generates
#' retrospective plots using `retro.rho.ndiaye2()`. Saves each plot as a high-resolution
#' PNG file. Optionally allows custom peel colors or palette choice.
#'
#' @param models_retro A named list of fitted SPiCT models containing retrospective output (`$retro`).
#' @param output_dir Directory where the plots should be saved. Created if it doesn't exist.
#'                   Default is `"FIG/"`.
#' @param peel_colors Optional vector of peel colors. If NULL, defaults to `cols()` output.
#' @param palette Character string specifying the palette name, currently unused in this function,
#'                but available for future customization.
#' @param width Width of the saved plots in inches. Default is 6.
#' @param height Height of the saved plots in inches. Default is 4.
#' @param dpi Resolution in dots per inch for the saved PNG. Default is 400.
#' @param prefix Filename prefix used when saving plots. Default is `"retrospective_"`.
#'
#' @details
#' This function is designed to automate the plotting of retrospective diagnostics
#' (B, F, B/Bmsy, F/Fmsy) for multiple SPiCT model fits. Each plot is generated using
#' a call to \code{retro.rho.ndiaye2()}, and saved with a filename that includes the
#' model's name.
#'
#' The function checks for the presence of `$retro` in each model and skips models
#' that are missing or have NULL retro outputs. Errors in plotting individual models
#' are caught and reported without stopping the loop.
#'
#' @return Invisibly returns \code{NULL}. Side effect: saves .png plots to \code{output_dir}.
#'
#' @seealso \code{\link{retro.rho.ndiaye2}}, \code{\link{ggsave}}, \code{\link{cols}}
#'
#' @examples
#' \dontrun{
#' models <- list(S1F = model_fox, S1P = model_pella, S1S = model_schaefer)
#' plot_all_retro(models, output_dir = "FIG/")
#' }
#'
#' @export
plot_all_retro <- function(models_retro,
                           output_dir = "FIG/",
                           peel_colors = NULL,  # <-- NEW
                           palette = "JABBA",
                           width = 6,
                           height = 4,
                           dpi = 400,
                           prefix = "retrospective_") {

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  mycols <- cols()

  for (model_name in names(models_retro)) {
    cat("Plotting retrospective for:", model_name, "\n")
    retro_obj <- models_retro[[model_name]]

    # Check if 'retro' exists and is not NULL
    if (is.null(retro_obj$retro)) {
      cat("  -> Skipping", model_name, ": 'retro' object is missing or NULL.\n")
      next
    }

    # Try generating plot, skip on error
    g <- tryCatch({
      retro.rho.ndiaye2(retro_obj, peel_colors  = mycols)
    }, error = function(e) {
      cat("  -> Error plotting", model_name, ":", e$message, "\n")
      NULL
    })

    if (!is.null(g)) {
      outpath <- file.path(output_dir, paste0(prefix, model_name, ".png"))
      ggsave(
        filename = outpath,
        plot = g,
        width = width,
        height = height,
        dpi = dpi
      )
      cat("  -> Saved plot for", model_name, "to", outpath, "\n")
    }
  }
  cat("All done.\n")
}
