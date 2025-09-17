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
plot_retro_grid.G4_B_et_F2_Corrected <- function(
    models, # named list, e.g., S1F, S1P, S1S
    scenario_name = "Scenario 1",
    peel_colors = NULL,
    palette = "JABBA",
    width = 20, height = 8, dpi = 400,
    output_dir = NULL, filename = NULL
) {
  require(patchwork)
  require(ggtext)
  require(MetBrewer)
  require(dplyr)
  require(grid)

  .safe_exists_fun <- function(name) exists(name, mode = "function", inherits = TRUE)
  .safe_get_fun   <- function(name) get(name, inherits = TRUE)

  .JABBA_FIXED <- c(
    "#000000", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
    "#66A61E", "#E6AB02", "#A6761D", "#666666", "#1F78B4"
  )

  .cols_fallback <- function(nms) {
    cols <- rep(.JABBA_FIXED, length.out = length(nms))
    names(cols) <- nms
    cols
  }

  .apply_theme <- local({
    if (.safe_exists_fun("theme_minimal_compact")) {
      .safe_get_fun("theme_minimal_compact")
    } else {
      function(base_size = 8, base_family = "") {
        ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 10),
            axis.title = ggplot2::element_text(face = "bold", size = 6),
            axis.text = ggplot2::element_text(size = 6, face = "bold"),
            legend.position = "bottom",
            legend.title = ggplot2::element_blank(),
            legend.text = ggplot2::element_text(size = 10, face = "bold"),
            legend.key.size = grid::unit(0.6, "lines"),
            legend.spacing.y = grid::unit(0, "pt"),
            legend.spacing.x = grid::unit(0, "pt"),
            legend.margin = ggplot2::margin(0, 0, 0, 0),
            legend.box.margin = ggplot2::margin(0, 0, 0, 0),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            panel.background = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(fill = NA, colour = "grey45", linewidth = 1.5),
            axis.ticks = ggplot2::element_line(linewidth = 0.5, color = "grey45"),
            axis.ticks.length = grid::unit(3, "pt"),
            strip.background = ggplot2::element_rect(fill = "grey45", color = "grey45", linewidth = 0.5),
            strip.text = ggplot2::element_text(face = "bold", size = ggplot2::rel(1)),
            text = ggplot2::element_text(face = "bold", size = 10),
            plot.margin = ggplot2::margin(2, 2, 2, 2)
          )
      }
    }
  })

  .get_par <- function(...) {
    if (.safe_exists_fun("get.par")) {
      .safe_get_fun("get.par")(...)
    } else if (requireNamespace("spict", quietly = TRUE)) {
      spict::get.par(...)
    } else {
      stop("get.par() not found (neither your version nor spict::get.par()).")
    }
  }

  .mohns_rho <- function(...) {
    if (.safe_exists_fun("mohns_rho")) {
      .safe_get_fun("mohns_rho")(...)
    } else {
      c(FFmsy = NA_real_, BBmsy = NA_real_, B = NA_real_, F = NA_real_)
    }
  }

  .make_peel_subtitle <- local({
    if (.safe_exists_fun("make_peel_subtitle")) {
      .safe_get_fun("make_peel_subtitle")
    } else {
      function(peel_names, peel_cols) {
        paste(
          mapply(function(nm, col) {
            paste0("<span style='color:", col, "'>\u25CF</span>&nbsp;<b>", nm, "</b>")
          }, peel_names, peel_cols),
          collapse = "&nbsp;&nbsp;&nbsp;"
        )
      }
    }
  })

  mycols <- if (is.null(peel_colors)) {
    if (.safe_exists_fun("cols")) .safe_get_fun("cols") else NULL
  } else peel_colors

  model_order <- names(models)

  get_retro_panels <- function(model_fit, model_nm, peel_colors_input) {
    runs <- model_fit$retro
    if (is.null(runs)) stop("Model '", model_nm, "' is missing retro.")
    names(runs)[1] <- "All"

    peel_names_raw <- names(runs)
    n_peels <- length(runs)

    peel_colors_final <- NULL
    if (!is.null(peel_colors_input)) {
      peel_colors_final <- rep(peel_colors_input, length.out = n_peels)
    } else if (toupper(palette) == "JABBA" && exists("jabba_colors", inherits = TRUE)) {
      peel_colors_final <- rep(get("jabba_colors", inherits = TRUE), length.out = n_peels)
    } else if (toupper(palette) == "JABBA") {
      peel_colors_final <- rep(.JABBA_FIXED, length.out = n_peels)
    } else if (!is.null(mycols)) {
      peel_colors_final <- rep(mycols, length.out = n_peels)
    } else if (!is.null(palette) &&
               isTRUE(palette %in% tryCatch(MetBrewer::list_met_palettes(), error = function(e) character()))) {
      peel_colors_final <- MetBrewer::met.brewer(palette, n_peels)
    } else {
      peel_colors_final <- scales::hue_pal()(n_peels)
    }
    names(peel_colors_final) <- peel_names_raw

    rho_vals <- tryCatch(
      round(.mohns_rho(model_fit, what = c("FFmsy","BBmsy","B","F")), 3),
      error = function(e) c(FFmsy = NA_real_, BBmsy = NA_real_, B = NA_real_, F = NA_real_)
    )

    # *** FIXED: use italic(MSY) instead of \emph{MSY} ***
    val_BB <- if (is.na(rho_vals["BBmsy"])) "NA" else sprintf("%.3f", rho_vals["BBmsy"])
    val_FF <- if (is.na(rho_vals["FFmsy"])) "NA" else sprintf("%.3f", rho_vals["FFmsy"])
    val_B  <- if (is.na(rho_vals["B"]))      "NA" else sprintf("%.3f", rho_vals["B"])
    val_F  <- if (is.na(rho_vals["F"]))      "NA" else sprintf("%.3f", rho_vals["F"])

    lab_BB <- paste0("Mohn*\"'s\"~rho[B/B[italic(MSY)]]==", val_BB)
    lab_FF <- paste0("Mohn*\"'s\"~rho[F/F[italic(MSY)]]==", val_FF)
    lab_B  <- paste0("Mohn*\"'s\"~rho[B]==",                 val_B)
    lab_F  <- paste0("Mohn*\"'s\"~rho[F]==",                 val_F)

    extract_df <- function(fit, peel, par_name) {
      idx  <- fit$inp$indest
      pars <- .get_par(par_name, fit, exp = TRUE, CI = 0.95)[idx, 1:3]
      tibble::tibble(
        peel     = peel,
        time     = fit$inp$time[idx],
        estimate = pars[, 2],
        lower    = pars[, 1],
        upper    = pars[, 3]
      )
    }

    df_B    <- purrr::imap_dfr(runs, ~ extract_df(.x, .y, "logB"))
    df_F    <- purrr::imap_dfr(runs, ~ extract_df(.x, .y, "logFnotS"))
    df_Bmsy <- purrr::imap_dfr(runs, ~ extract_df(.x, .y, "logBBmsy"))
    df_Fmsy <- purrr::imap_dfr(runs, ~ extract_df(.x, .y, "logFFmsynotS"))

    peel_names <- unique(df_B$peel)
    all_present <- "All" %in% peel_names
    other_peels <- setdiff(peel_names, "All")
    peel_nums <- suppressWarnings(as.integer(other_peels))
    valid_peels   <- other_peels[!is.na(peel_nums)]
    invalid_peels <- other_peels[is.na(peel_nums)]
    ordered_numeric_peels <- valid_peels[order(as.integer(valid_peels), decreasing = TRUE)]
    plot_peel_levels <- c(if (all_present) "All", ordered_numeric_peels, invalid_peels)

    df_B$peel    <- factor(df_B$peel,    levels = plot_peel_levels)
    df_F$peel    <- factor(df_F$peel,    levels = plot_peel_levels)
    df_Bmsy$peel <- factor(df_Bmsy$peel, levels = plot_peel_levels)
    df_Fmsy$peel <- factor(df_Fmsy$peel, levels = plot_peel_levels)

    if (!all(plot_peel_levels %in% names(peel_colors_final))) {
      peel_colors_final <- .cols_fallback(plot_peel_levels)
    } else {
      peel_colors_final <- peel_colors_final[plot_peel_levels]
    }

    ci_gray <- "#B0B0B0"

    make_panel <- function(df, ylab, rho_text = NULL) {
      ggplot2::ggplot(df, ggplot2::aes(time, estimate, group = peel)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = ci_gray, alpha = 0.25) +
        ggplot2::geom_line(ggplot2::aes(color = peel), linewidth = 0.5) +
        ggplot2::labs(x = NULL, y = ylab) +
        .apply_theme() +
        ggplot2::scale_color_manual(values = peel_colors_final, guide = "none") +
        ggplot2::theme(
          axis.title.y = ggplot2::element_text(
            margin = ggplot2::margin(r = 0),
            face = "bold", size = 8
          )
        ) +
        {
          if (!is.null(rho_text)) {
            ggplot2::annotate(
              "text", x = Inf, y = Inf, label = rho_text,
              parse = TRUE, hjust = 1.1, vjust = 1.5, size = 3
            )
          } else NULL
        }
    }

    list(
      B     = make_panel(df_B,    expression(bold(B[t])),                        rho_text = lab_B),
      F     = make_panel(df_F,    expression(bold(F[t])),                        rho_text = lab_F),
      BBmsy = make_panel(df_Bmsy, expression(bold(B[t]/B[italic(MSY)])),         rho_text = lab_BB),
      FFmsy = make_panel(df_Fmsy, expression(bold(F[t]/F[italic(MSY)])),         rho_text = lab_FF)
    )
  }

  all_panels <- lapply(model_order, function(nm) {
    get_retro_panels(models[[nm]], nm, mycols)
  })

  make_model_header_plot <- function(model_code) {
    ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::annotate(
        "text",
        x = 0.5, y = 0.5,
        label = model_code,
        size = 5, fontface = "bold",
        hjust = 0.5, vjust = 0.5
      ) +
      ggplot2::theme(plot.margin = ggplot2::margin(2.5, 0, 0, 0, "pt"))
  }

  patchwork_cols <- list()
  for (j in seq_along(model_order)) {
    model_code  <- model_order[j]
    header_plot <- make_model_header_plot(model_code)
    panels <- all_panels[[j]]
    patchwork_cols[[j]] <- patchwork::wrap_plots(
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
  if (length(patchwork_cols) >= 2L) {
    for (j in 2:length(patchwork_cols)) {
      p_grid <- p_grid | patchwork_cols[[j]]
    }
  }

  first_fit <- models[[model_order[1]]]
  runs <- first_fit$retro
  names(runs)[1] <- "All"
  peel_names_first <- names(runs)

  if (is.null(peel_colors)) {
    if (.safe_exists_fun("cols")) {
      base_cols <- .safe_get_fun("cols")
      peel_colors_used <- rep(base_cols, length.out = length(peel_names_first))
    } else if (toupper(palette) == "JABBA" && exists("jabba_colors", inherits = TRUE)) {
      peel_colors_used <- rep(get("jabba_colors", inherits = TRUE), length.out = length(peel_names_first))
    } else if (toupper(palette) == "JABBA") {
      peel_colors_used <- rep(.JABBA_FIXED, length.out = length(peel_names_first))
    } else if (!is.null(palette) &&
               isTRUE(palette %in% tryCatch(MetBrewer::list_met_palettes(), error = function(e) character()))) {
      peel_colors_used <- MetBrewer::met.brewer(palette, length(peel_names_first))
    } else {
      peel_colors_used <- scales::hue_pal()(length(peel_names_first))
    }
  } else {
    peel_colors_used <- rep(peel_colors, length.out = length(peel_names_first))
  }
  names(peel_colors_used) <- peel_names_first

  subtitle_peels <- .make_peel_subtitle(peel_names_first, peel_colors_used)

  final <- p_grid +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(
      plot.title = ggtext::element_markdown(hjust = 0.5, face = "bold", size = 17, margin = ggplot2::margin(b = 25)),
      plot.subtitle = ggtext::element_markdown(hjust = 0.5, size = 17, face = "bold", margin = ggplot2::margin(b = 25)),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    ) &
    patchwork::plot_annotation(
      title = scenario_name,
      subtitle = subtitle_peels
    )

  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    if (is.null(filename)) {
      scen_lab <- gsub("[^A-Za-z0-9]", "_", scenario_name)
      filename <- paste0("retrospective_", scen_lab, ".png")
    }
    outpath <- file.path(output_dir, filename)
    ggplot2::ggsave(outpath, plot = final, width = width, height = height, dpi = dpi)
    message("Saved to: ", outpath)
  }
  return(final)
}
