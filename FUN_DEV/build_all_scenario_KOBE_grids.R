kobe_safe <- function(rep_obj,
                      logax = FALSE, plot.legend = TRUE, man.legend = TRUE,
                      ext = TRUE, rel.axes = FALSE,
                      xlim = NULL, ylim = NULL,
                      labpos = c(1, 1), xlabel = NULL, stamp = NULL,
                      verbose = TRUE, CI = 0.95) {
  res <- try(
    {
      p <- kobe_all_in_one_gg(rep = rep_obj,
                              logax = logax, plot.legend = plot.legend, man.legend = man.legend,
                              ext = ext, rel.axes = rel.axes,
                              xlim = xlim, ylim = ylim,
                              labpos = labpos, xlabel = xlabel, stamp = stamp,
                              verbose = verbose, CI = CI)
      # Force it to be visible so patchwork & ggsave can use it
      p
    },
    silent = TRUE
  )
  if (inherits(res, "try-error") || is.null(res)) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
    p_err <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0,
                        label = "Plot error", size = 4, colour = "red", fontface = "bold") +
      ggplot2::theme_void()
    return(p_err)
  }
  res
}


# Add a model ID label at top-left inside the plotting area
# add_id_label <- function(p, id_text) {
#   if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
#   p + ggplot2::annotate("text",
#                         x = -Inf, y = Inf, label = id_text,
#                         hjust = -0.2, vjust = 2,  # nudged inward
#                         fontface = "bold", size = 6) +
#       ggplot2::coord_cartesian(clip = "off") +
#       ggplot2::theme(plot.margin = ggplot2::margin(6, 6, 6, 6))
# }
# Add a model ID label at top-left inside the plotting area
add_id_label <- function(p, id_text, id_size = 4) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
  p + ggplot2::annotate(
    "text",
    x = -Inf, y = Inf, label = id_text,
    hjust = -0.2, vjust = 2,         # keep your current nudge
    fontface = "bold", size = id_size
  ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(plot.margin = ggplot2::margin(6, 6, 6, 6))
}


# Internal: pick an element by suffix (.SDM/.GLM) and model letter (P/S/F)
.pick_model <- function(lst, model_letter, data_type) {
  # model_letter in {"P","S","F"}, data_type in {"SDM","GLM"}
  # Example names like "S1P.GLM", "S1S.SDM", etc.
  nms <- names(lst)
  target <- paste0(model_letter, ".", data_type)
  idx <- grep(paste0(model_letter, "\\.", data_type, "$"), nms)
  if (length(idx) == 0) return(NULL)
  lst[[idx[1]]]
}


#' Build Kobe plot grids for all scenarios
#'
#' @description
#' Generates and optionally saves Kobe plot grids for every scenario in a list
#' of fitted models. For each scenario, this function calls
#' [elu2::build_scenario_KOBE_grid()] to create the Kobe grid (e.g., "3x2"
#' layout) using `kobe_all_in_one_gg()` for each submodel in the scenario.
#'
#' The result is a named list of `ggplot`/patchwork objects, one per scenario.
#' Optionally, PNG files are saved to a specified directory.
#'
#' @details
#' This is a **batch runner** that wraps around
#' [elu2::build_scenario_KOBE_grid()]. It is intended for workflows where
#' multiple management scenarios (e.g., S1–S16) have been fitted with multiple
#' SPiCT-based models (e.g., Pella, Schaefer, Fox).
#'
#' The `all_models` argument should be a **named list** where each element name
#' is a scenario ID (e.g., `"S1"`, `"S2"`), and each element is itself a list of
#' fitted SPiCT model objects for that scenario.
#'
#' @param all_models Named list of fitted model lists, indexed by scenario name.
#'   Each scenario element should be a named list of fitted `spictcls` objects
#'   (result of `elu2::fit.elu2()` or compatible).
#' @param layout Character string, `"3x2"` (default) or `"3x1"`, controlling
#'   the patchwork grid arrangement for each scenario’s Kobe plots.
#' @param save Logical; if `TRUE` (default), saves each grid as a PNG to `output_dir`.
#' @param output_dir Directory path where PNGs will be saved when `save = TRUE`.
#'   If `NULL` (default), files are saved in the current working directory.
#' @param width,height Numeric; width/height in pixels for saved PNGs when
#'   `save = TRUE`. Defaults are `1500` × `800` (approx. 12.5 × 6.7 inches at 120 dpi).
#' @param dpi Numeric; resolution for saved PNGs. Default `120`.
#' @param ... Additional arguments passed to
#'   [elu2::build_scenario_KOBE_grid()] and/or `kobe_all_in_one_gg()`.
#'
#' @return
#' Named list of Kobe plot grids (`patchwork` objects), one per scenario.
#'
#' If `save = TRUE`, each scenario grid is also saved to a PNG file with the
#' filename pattern:
#' ```
#' paste0("KobeGrid_", scenario_name, "_", layout, ".png")
#' ```
#'
#' @section Workflow:
#' This function is typically run **after** all scenarios have been fitted and
#' stored in a nested list (`all_models`). For example:
#' \preformatted{
#' all_models <- list(
#'   S1 = list(Pella = fit.elu2(inp_S1P), Schaefer = fit.elu2(inp_S1S)),
#'   S2 = list(Pella = fit.elu2(inp_S2P), Schaefer = fit.elu2(inp_S2S))
#' )
#' }
#'
#' @seealso
#' [elu2::build_scenario_KOBE_grid()] for creating a Kobe grid for one scenario.
#' [elu2::kobe_all_in_one_gg()] for generating a single Kobe plot.
#'
#' @examples
#' \dontrun{
#' # Assuming all_models is a named list of scenario fits:
#' plots <- build_all_scenario_KOBE_grids(all_models, layout = "3x2", save = TRUE)
#'
#' # Access the plot for Scenario 1
#' plots$S1
#' }
#'
#' @author
#' Elhadji Ndiaye & Contributors
#'
#' @encoding UTF-8
#'
#' @import ggplot2
#' @export
build_all_scenario_KOBE_grids <- function(all_models,
                                scenario_name,
                                layout = c("3x2", "3x1"),
                                # passthrough plotting args to your Kobe function:
                                logax = FALSE, plot.legend = TRUE, man.legend = TRUE,
                                ext = TRUE, rel.axes = FALSE,
                                xlim = NULL, ylim = NULL,
                                xlabel = NULL, stamp = NULL, verbose = TRUE, CI = 0.95) {

  layout <- match.arg(layout)

  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required for the grid layout.")
  }

  if (!scenario_name %in% names(all_models)) {
    stop("Scenario '", scenario_name, "' not found in all_models.")
  }

  sc <- all_models[[scenario_name]]
  if (!is.list(sc) || length(sc) == 0) stop("Empty scenario list: ", scenario_name)

  # Prepare the six panels if layout == "3x2"
  # Row order: Pella (P), Schaefer (S), Fox (F)
  # Column order: SDM (left), GLM (right)
  if (layout == "3x2") {
    # Retrieve objects
    reps <- list(
      P_SDM = .pick_model(sc, "P", "SDM"),
      P_GLM = .pick_model(sc, "P", "GLM"),
      S_SDM = .pick_model(sc, "S", "SDM"),
      S_GLM = .pick_model(sc, "S", "GLM"),
      F_SDM = .pick_model(sc, "F", "SDM"),
      F_GLM = .pick_model(sc, "F", "GLM")
    )

    # Make plots (with IDs = names present in scenario sublist)
    get_id <- function(letter, dtype) {
      base <- paste0(scenario_name, letter, ".", dtype)
      nm <- names(sc)
      hit <- nm[grep(paste0(letter, "\\.", dtype, "$"), nm)]
      if (length(hit)) hit[1] else base
    }

    p_P_SDM <- add_id_label(
      kobe_safe(reps$P_SDM, logax, plot.legend, man.legend, ext, rel.axes,
                xlim, ylim, c(1, 1), xlabel, stamp, verbose, CI),
      get_id("P", "SDM"))
    p_P_GLM <- add_id_label(
      kobe_safe(reps$P_GLM, logax, plot.legend, man.legend, ext, rel.axes,
                xlim, ylim, c(1, 1), xlabel, stamp, verbose, CI),
      get_id("P", "GLM"))

    p_S_SDM <- add_id_label(
      kobe_safe(reps$S_SDM, logax, plot.legend, man.legend, ext, rel.axes,
                xlim, ylim, c(1, 1), xlabel, stamp, verbose, CI),
      get_id("S", "SDM"))
    p_S_GLM <- add_id_label(
      kobe_safe(reps$S_GLM, logax, plot.legend, man.legend, ext, rel.axes,
                xlim, ylim, c(1, 1), xlabel, stamp, verbose, CI),
      get_id("S", "GLM"))

    p_F_SDM <- add_id_label(
      kobe_safe(reps$F_SDM, logax, plot.legend, man.legend, ext, rel.axes,
                xlim, ylim, c(1, 1), xlabel, stamp, verbose, CI),
      get_id("F", "SDM"))
    p_F_GLM <- add_id_label(
      kobe_safe(reps$F_GLM, logax, plot.legend, man.legend, ext, rel.axes,
                xlim, ylim, c(1, 1), xlabel, stamp, verbose, CI),
      get_id("F", "GLM"))

    # Assemble with patchwork
    #g <- (p_P_SDM | p_P_GLM) /
    #(p_S_SDM | p_S_GLM) /
    #(p_F_SDM | p_F_GLM)

    # Assemble with patchwork
    g <- (p_P_SDM | p_S_SDM | p_F_SDM)/
      (p_P_GLM | p_S_GLM | p_F_GLM)


    # Add a big title for the scenario
    #g <- g + patchwork::plot_annotation(
    #title = paste0("Scenario ", scenario_name)
    #)

    return(g)
  }

  # Alternative layout: 3x1 (one panel per model using a priority order SDM > GLM)
  if (layout == "3x1") {
    pick_best <- function(letter) {
      obj <- .pick_model(sc, letter, "SDM")
      if (is.null(obj)) obj <- .pick_model(sc, letter, "GLM")
      obj
    }
    id_best  <- function(letter) {
      nm <- names(sc)
      hit <- nm[grep(paste0(letter, "\\.(SDM|GLM)$"), nm)]
      if (length(hit)) hit[1] else paste0(scenario_name, letter)
    }

    repP <- pick_best("P"); idP <- id_best("P")
    repS <- pick_best("S"); idS <- id_best("S")
    repF <- pick_best("F"); idF <- id_best("F")

    pP <- add_id_label(
      kobe_safe(repP, logax, plot.legend, man.legend, ext, rel.axes,
                xlim, ylim, c(1, 1), xlabel, stamp, verbose, CI),
      idP)
    pS <- add_id_label(
      kobe_safe(repS, logax, plot.legend, man.legend, ext, rel.axes,
                xlim, ylim, c(1, 1), xlabel, stamp, verbose, CI),
      idS)
    pF <- add_id_label(
      kobe_safe(repF, logax, plot.legend, man.legend, ext, rel.axes,
                xlim, ylim, c(1, 1), xlabel, stamp, verbose, CI),
      idF)

    g <- pP / pS / pF
    g <- g + patchwork::plot_annotation(
      title = paste0("Scenario ", scenario_name, " (3×1)")
    )
    return(g)
  }

  stop("Unknown layout.")
}

# Batch build all scenarios and (optionally) save each to a file
# ---- helper: robust saver with fallback to base png() ----
.save_plot_safely <- function(plot_obj, file_path, width, height, dpi) {
  # 1) Try ggsave (fast path)
  ok <- try({
    ggplot2::ggsave(
      filename = file_path,
      plot     = plot_obj,
      device   = "png",
      width    = width,
      height   = height,
      units    = "in",
      dpi      = dpi,
      bg       = "white"  # avoids transparent bg issues in Beamer/Word
    )
    TRUE
  }, silent = TRUE)

  if (!inherits(ok, "try-error")) return(invisible(TRUE))

  # 2) Fallback: base png() device (very reliable)
  gr <- try({
    grDevices::png(
      filename = file_path,
      width    = width,
      height   = height,
      units    = "in",
      res      = dpi,
      type     = "cairo" # better text rendering (if available)
    )
    print(plot_obj)
    grDevices::dev.off()
    TRUE
  }, silent = TRUE)

  if (inherits(gr, "try-error")) {
    warning("Failed to save plot to: ", file_path)
    return(invisible(FALSE))
  }
  invisible(TRUE)
}

# ---- REPLACE your build_all_scenario_grids() with this version ----
build_all_scenario_grids <- function(all_models,
                                     layout = "3x2",
                                     save = TRUE,
                                     out_dir = file.path("FIG", "KobePhases"),
                                     width = if (layout == "3x2") 15 else 6,
                                     height = if (layout == "3x2") 6.5 else 10,
                                     dpi = 300,
                                     ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("patchwork is required.")

  # Create FIG/KobePhases (case-sensitive filesystems matter)
  if (save && !dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  scens <- names(all_models)
  out <- vector("list", length(scens)); names(out) <- scens

  i <- 1L
  while (i <= length(scens)) {
    scen <- scens[i]

    # Build the patchwork object for the scenario
    g <- build_scenario_grid(all_models, scen, layout = layout, ...)
    out[[i]] <- g

    # Save per-scenario PNG
    if (save) {
      fpath <- file.path(out_dir, paste0("KobePhases_", scen, "_", layout, ".png"))
      .save_plot_safely(g, fpath, width = width, height = height, dpi = dpi)
    }
    i <- i + 1L
  }
  invisible(out)
}
