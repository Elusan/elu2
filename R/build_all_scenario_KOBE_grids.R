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


#' Build a Kobe plot grid for a single scenario
#'
#' @description
#' Constructs a grid of Kobe plots for one scenario from `all_models`, using
#' [kobe_all_in_one_gg()] for each model/data-source combination and arranging
#' the panels with **patchwork**.
#'
#' @details
#' Layout `"3x2"` creates six panels: rows = models (Pella, Schaefer, Fox),
#' columns = data types (SDM, GLM). Layout `"3x1"` creates three panels, one
#' per model (preferring SDM if both SDM and GLM exist).
#'
#' @param all_models Named list of scenarios; each scenario is a named list of
#'   fitted SPiCT objects (e.g., `"S1P.SDM"`, `"S1S.GLM"`).
#' @param scenario_name Character; scenario to plot (must exist in `all_models`).
#' @param layout Either `"3x2"` (default) or `"3x1"`.
#' @param logax,plot.legend,man.legend,ext,rel.axes,xlim,ylim,xlabel,stamp,verbose,CI
#'   Passed through to [kobe_all_in_one_gg()]. See its documentation for details.
#'
#' @return
#' A **patchwork** object representing the scenario grid (also a `ggplot`).
#'
#' @examples
#' \dontrun{
#' # One scenario, 3x2 grid:
#' build_scenario_KOBE_grid(all_models, "S1", layout = "3x2")
#'
#' # One scenario, 3x1 grid (prefers SDM when both exist):
#' build_scenario_KOBE_grid(all_models, "S1", layout = "3x1")
#' }
#'
#' @seealso [build_all_scenario_KOBE_grids()], [kobe_all_in_one_gg()]
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
