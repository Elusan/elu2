#' Plot retrospective diagnostics grid for multiple SPiCT models
#'
#' Generates a grid of retrospective plots (B, F, B/Bmsy, F/Fmsy) for each model
#' in a named list. Each column represents one model, and each row a metric.
#' The function also annotates Mohn's rho for B/Bmsy and F/Fmsy. Optionally saves
#' the output to file.
#'
#' @param models A named list of fitted SPiCT models, typically including retro results.
#'               Example: \code{list(S1F = model_fox, S1P = model_pella, S1S = model_schaefer)}
#' @param scenario_name A character string used as the title of the plot. Also used
#'                      to generate the output filename if \code{filename} is not provided.
#' @param peel_colors Optional character vector of colors to assign to each peel.
#'                    If NULL, a default palette is used (JABBA or MetBrewer).
#' @param palette A character string specifying the color palette to use for peels.
#'                Options include "JABBA" or any valid palette from \code{MetBrewer::list_met_palettes()}.
#' @param width Numeric. Width (in inches) of the output PNG file if saved.
#' @param height Numeric. Height (in inches) of the output PNG file if saved.
#' @param dpi Resolution of the saved figure in dots per inch.
#' @param output_dir Optional character string specifying the directory where the figure
#'                   should be saved. If NULL, the plot is not saved.
#' @param filename Optional filename (with .png extension) for saving the figure. If NULL,
#'                 a name will be auto-generated from \code{scenario_name}.
#'
#' @details
#' The function extracts B, F, B/Bmsy, and F/Fmsy trajectories for each model and
#' peel, computes Mohn’s rho values for B/Bmsy and F/Fmsy, and assembles them into
#' a clean grid using the \pkg{patchwork} package. Confidence intervals are shown
#' with shaded bands. The peel legend is embedded as a subtitle using colored dots.
#'
#' If \code{output_dir} is provided, the figure is saved as a high-resolution PNG.
#' This function is useful for summarizing and comparing retrospective performance
#' across multiple production model fits (e.g., Fox, Schaefer, Pella).
#'
#' @return A \code{patchwork} plot object displaying the full grid of retrospective plots.
#'         Invisibly returns the plot for further customization or saving.
#'
#' @import ggplot2
#' @import patchwork
#' @importFrom MetBrewer met.brewer
#' @importFrom purrr imap_dfr
#' @importFrom ggtext element_markdown
#' @importFrom grid unit
#'
#' @seealso \code{\link{retro}}, \code{\link{mohns_rho}}, \code{\link{plot_retro_grid.ELU5}}
#'
#' @examples
#' \dontrun{
#' models <- list(S1F = model_fox, S1P = model_pella, S1S = model_schaefer)
#' plot_retro_grid.ELU4(models, scenario_name = "Scenario 1", output_dir = "FIG")
#' }
#'
#' @export
plot_retro_grid.G4 <- function(
    models, # named list, e.g., S1F, S1P, S1S
    scenario_name = "Scenario 1",
    peel_colors = NULL,
    palette = "JABBA",
    width = 16, height = 10, dpi = 400,
    output_dir = NULL, filename = NULL
) {
  require(patchwork)
  require(ggtext)
  require(MetBrewer)
  require(dplyr)
  require(grid)

  # Color logic
  mycols <- if (is.null(peel_colors)) cols() else peel_colors
  model_order <- names(models)

  # =========== Panel builder function ===========
  get_retro_panels <- function(model_fit, model_nm, peel_colors) {
    runs <- model_fit$retro
    if (is.null(runs)) stop("Model '", model_nm, "' is missing retro.")
    names(runs)[1] <- "All"
    n_peels <- length(runs)
    peel_names <- names(runs)
    if (is.null(peel_colors)) {
      if (toupper(palette) == "JABBA") {
        peel_colors <- rep(jabba_colors, length.out = n_peels)
      } else if (palette %in% MetBrewer::list_met_palettes()) {
        peel_colors <- MetBrewer::met.brewer(palette, n_peels)
      } else {
        peel_colors <- scales::hue_pal()(n_peels)
      }
    } else {
      peel_colors <- rep(peel_colors, length.out = n_peels)
    }
    names(peel_colors) <- peel_names
    # Mohn's rho
    rho <- tryCatch(mohns_rho(model_fit, what = c("FFmsy","BBmsy")) %>% round(3), error = function(e) rep(NA,2))
    lab_BB <- paste0("Mohn*\"'s\"~rho[B/B[MSY]]==", rho["BBmsy"])
    lab_FF <- paste0("Mohn*\"'s\"~rho[F/F[MSY]]==", rho["FFmsy"])
    # Extract
    extract_df <- function(fit, peel, par_name) {
      idx  <- fit$inp$indest
      pars <- get.par(par_name, fit, exp = TRUE, CI = 0.95)[idx, 1:3]
      tibble(
        peel     = peel,
        time     = fit$inp$time[idx],
        estimate = pars[,2],
        lower    = pars[,1],
        upper    = pars[,3]
      )
    }
    df_B    <- imap_dfr(runs, ~ extract_df(.x, .y, "logB"))
    df_F    <- imap_dfr(runs, ~ extract_df(.x, .y, "logFnotS"))
    df_Bmsy <- imap_dfr(runs, ~ extract_df(.x, .y, "logBBmsy"))
    df_Fmsy <- imap_dfr(runs, ~ extract_df(.x, .y, "logFFmsynotS"))

    # Peel order

    # Peel names
    peel_names <- names(runs)

    # Separate special case "All"
    all_present <- "All" %in% peel_names
    other_peels <- setdiff(peel_names, "All")

    # Identify numeric peels
    peel_nums <- suppressWarnings(as.integer(other_peels))

    # Separate valid and non-numeric peels
    valid_peels   <- other_peels[!is.na(peel_nums)]
    invalid_peels <- other_peels[is.na(peel_nums)]

    # Order numeric peels from latest to earliest
    ordered_numeric_peels <- valid_peels[order(as.integer(valid_peels), decreasing = TRUE)]

    # Final order: All -> numeric peels -> non-numeric leftovers
    plot_peel_levels <- c(if (all_present) "All", ordered_numeric_peels, invalid_peels)




    df_B$peel    <- factor(df_B$peel, levels = plot_peel_levels)
    df_F$peel    <- factor(df_F$peel, levels = plot_peel_levels)
    df_Bmsy$peel <- factor(df_Bmsy$peel, levels = plot_peel_levels)
    df_Fmsy$peel <- factor(df_Fmsy$peel, levels = plot_peel_levels)
    ci_gray <- "#B0B0B0"
    # Panel builder
    make_panel <- function(df, ylab, rho_text = NULL) {
      ggplot(df, aes(time, estimate, group = peel)) +
        geom_ribbon(aes(ymin = lower, ymax = upper), fill = ci_gray, alpha = 0.25) +
        geom_line(aes(color = peel), size = 0.7) +
        labs(x = NULL, y = ylab) +
        theme_minimal_compact() +
        scale_color_manual(values = peel_colors, guide = "none") +
        theme(
          axis.title.y = element_text(
            margin = margin(r = 0),
            face = "bold", size = 11
          )
        ) +
        {if (!is.null(rho_text)) annotate("text", x = Inf, y = Inf, label = rho_text, parse = TRUE, hjust = 1.1, vjust = 1.5, size = 3.5) else NULL}
    }
    list(
      B    = make_panel(df_B,    expression(bold(B[t]))),
      F    = make_panel(df_F,    expression(bold(F[t]))),
      BBmsy= make_panel(df_Bmsy, expression(bold(B[t]/B[MSY])), lab_BB),
      FFmsy= make_panel(df_Fmsy, expression(bold(F[t]/F[MSY])), lab_FF)
    )
  }

  all_panels <- lapply(model_order, function(nm) {
    get_retro_panels(models[[nm]], nm, mycols)
  })

  # Header
  make_model_header_plot <- function(model_code) {
    ggplot() +
      theme_void() +
      annotate(
        "text",
        x = 0.5, y = 0.5,
        label = model_code,
        size = 6, fontface = "bold",
        hjust = 0.5, vjust = 0.5
      ) +
      theme(plot.margin = margin(0, 0, 0, 0, "pt"))
  }
  patchwork_cols <- list()
  for (j in seq_along(model_order)) {
    model_code  <- model_order[j]
    header_plot <- make_model_header_plot(model_code)
    panels <- all_panels[[j]]
    patchwork_cols[[j]] <- wrap_plots(
      header_plot,
      panels$B,
      panels$F,
      panels$BBmsy,
      panels$FFmsy,
      ncol = 1,
      heights = c(0.16, 1, 1, 1, 1)
    )
  }
  p_grid <- patchwork_cols[[1]]
  for (j in 2:length(patchwork_cols)) {
    p_grid <- p_grid | patchwork_cols[[j]]
  }

  # ---- Subtitle (peel dots) - use peels/colors from first model ----
  first_fit <- models[[model_order[1]]]
  runs <- first_fit$retro
  names(runs)[1] <- "All"
  peel_names <- names(runs)
  peel_colors_used <- rep(mycols, length.out = length(peel_names))
  subtitle_peels <- make_peel_subtitle(peel_names, peel_colors_used)

  # ---- Compose main plot with scenario title and subtitle ----
  final <- p_grid +
    plot_layout(guides = "collect") &
    theme(
      plot.title = ggtext::element_markdown(hjust = 0.5, face = "bold", size = 17, margin=margin(b=25)),
      plot.subtitle = ggtext::element_markdown(hjust = 0.5, size = 17, face = "bold", margin=margin(b=25)),
      plot.background = element_rect(fill = "white", color = NA)
    ) &
    plot_annotation(
      title = scenario_name,
      subtitle = subtitle_peels
    )

  # Save if requested
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    #if (is.null(filename)) filename <- paste0("scenario_", gsub(" ", "_", scenario_name), "_retro_grid.png")
    if (is.null(filename)) {
      # Remove spaces and special characters for file safety, add prefix
      scen_lab <- gsub("[^A-Za-z0-9]", "_", scenario_name)
      filename <- paste0("retrospective_", scen_lab, ".png")
    }
    outpath <- file.path(output_dir, filename)
    ggsave(outpath, plot = final, width = width, height = height, dpi = dpi)
    message("Saved to: ", outpath)
  }
  return(final)
}
