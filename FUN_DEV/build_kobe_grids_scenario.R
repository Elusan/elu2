#' Safely build a single Kobe plot (robust to failures)
#'
#' @description
#' Thin wrapper around [elu2::kobe_all_in_one_gg()] that traps errors and
#' returns a minimal error plot instead of failing. This is useful when you
#' batch-generate many Kobe plots and want the pipeline to continue even if
#' one fit is problematic.
#'
#' @param rep_obj A fitted model object compatible with
#'   [elu2::kobe_all_in_one_gg()] (e.g., `spictcls` or `elu2::fit.elu2` result).
#' @param logax Logical; log-scale both axes? Default `FALSE`.
#' @param plot.legend Logical; show legend for MSY / E(B[∞]) / True. Default `TRUE`.
#' @param man.legend Logical; show management-scenario legend (if present). Default `TRUE`.
#' @param ext Logical; add relative axes (B/Bmsy, F/Fmsy) on top/right when
#'   `rel.axes = FALSE`. Default `TRUE`.
#' @param rel.axes Logical; plot axes as relative units (B/Bmsy vs F/Fmsy).
#'   Default `FALSE` (absolute Ft vs biomass).
#' @param xlim,ylim Numeric length-2; plot limits. If `NULL`, computed internally.
#' @param labpos Integer length-2; kept for API parity (not used by ggplot). Default `c(1, 1)`.
#' @param xlabel Optional x-axis label override.
#' @param stamp Optional character version stamp. If `NULL`, no stamp.
#' @param verbose Logical; print diagnostics. Default `TRUE`.
#' @param CI Numeric in (0,1); confidence level for ellipses / ribbons. Default `0.95`.
#'
#' @return A `ggplot` object. If an internal error occurs, a simple red
#'   "Plot error" placeholder plot is returned instead.
#'
#' @seealso [elu2::kobe_all_in_one_gg()], [elu2::build_kobe_01()],
#'   [elu2::build_kobe_grids_scenario()]
#'
#' @examples
#' \dontrun{
#' p <- kobe_safe(S1_Pella.SDM, rel.axes = FALSE, CI = 0.95)
#' print(p)
#' }
#'
#' @import ggplot2
#' @export
kobe_safe <- function(rep_obj,
                      logax = FALSE, plot.legend = TRUE, man.legend = TRUE,
                      ext = TRUE, rel.axes = FALSE,
                      xlim = NULL, ylim = NULL,
                      labpos = c(1, 1), xlabel = NULL, stamp = NULL,
                      verbose = TRUE, CI = 0.95) {
  res <- try(
    {
      # delegate the real work to your ggplot builder
      p <- kobe_all_in_one_gg(
        rep = rep_obj,
        logax = logax, plot.legend = plot.legend, man.legend = man.legend,
        ext = ext, rel.axes = rel.axes,
        xlim = xlim, ylim = ylim,
        labpos = labpos, xlabel = xlabel, stamp = stamp,
        verbose = verbose, CI = CI
      )
      p
    },
    silent = TRUE
  )

  # If anything failed, return a minimal but valid ggplot
  if (inherits(res, "try-error") || is.null(res)) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
    p_err <- ggplot2::ggplot() +
      ggplot2::annotate(
        "text", x = 0, y = 0,
        label = "Plot error",
        size = 4, colour = "red", fontface = "bold"
      ) +
      ggplot2::xlim(-1, 1) + ggplot2::ylim(-1, 1) +
      ggplot2::theme_void()
    return(p_err)
  }
  res
}

#' Add a model ID label inside a plot (top-left)
#'
#' @description
#' Annotates a `ggplot` with a bold ID (e.g., `"S1P.SDM"`) in the top-left
#' inside the plotting area. Coordinates use `-Inf`/`Inf` with small nudges,
#' and clipping is disabled so text is always visible.
#'
#' @param p A `ggplot` object.
#' @param id_text Character; the label to draw (e.g., `"S1P.SDM"`).
#' @param id_size Numeric; text size passed to `annotate()`. Default `4`.
#'
#' @return The input `ggplot` with an added text annotation and slightly
#'   enlarged plot margins.
#'
#' @import ggplot2
#' @export
add_id_label <- function(p, id_text, id_size = 4) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
  p +
    ggplot2::annotate(
      "text",
      x = -Inf, y = Inf, label = id_text,
      hjust = -0.2, vjust = 2,     # small, empirically chosen nudges
      fontface = "bold", size = id_size
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(plot.margin = ggplot2::margin(6, 6, 6, 6))
}

# INTERNAL: pick a submodel by suffix (".SDM"/".GLM") and model letter ("P"/"S"/"F")
#'
#' @keywords internal
#' @noRd
.pick_model <- function(lst, model_letter, data_type) {
  # Example names like "S1P.GLM", "S1S.SDM", etc.
  # We find the FIRST name in `lst` that ends with "<letter>.<data_type>"
  # and return that element; otherwise NULL.
  if (!length(lst)) return(NULL)
  nms <- names(lst)
  idx <- grep(paste0(model_letter, "\\.", data_type, "$"), nms)
  if (length(idx) == 0) return(NULL)
  lst[[idx[1]]]
}

#' Build a Kobe-plot grid for one scenario
#'
#' @description
#' Create a `patchwork` grid of Kobe plots for a single scenario held inside
#' `all_models[[scenario_name]]`. Two layouts are supported:
#'
#' * `"3x2"` (default): rows = models P / S / F; columns = SDM | GLM.
#' * `"3x1"`: one column; prefers SDM if present, otherwise GLM for each model.
#'
#' Each panel is produced via [elu2::kobe_safe()], with a bold per-panel ID
#' placed at the top-left using [elu2::add_id_label()].
#'
#' @param all_models Named list of scenarios. Each element is a named list of
#'   fitted models (objects accepted by [elu2::kobe_all_in_one_gg()]) whose
#'   names end with `".SDM"` or `".GLM"` and also contain the model letter
#'   just before the dot (e.g., `"S1P.SDM"`, `"S1S.GLM"`, `"S1F.SDM"`).
#' @param scenario_name Character; name of the scenario to render (must exist
#'   in `names(all_models)`).
#' @param layout Either `"3x2"` (default) or `"3x1"`.
#' @param logax,plot.legend,man.legend,ext,rel.axes,xlim,ylim,xlabel,stamp,verbose,CI
#'   Passed through to [elu2::kobe_safe()] / [elu2::kobe_all_in_one_gg()].
#'
#' @return A `patchwork` object (a composed grid of six or three Kobe plots).
#'
#' @examples
#' \dontrun{
#' g <- build_kobe_01(all_models, "S1", layout = "3x2", CI = 0.95)
#' g
#' }
#'
#' @import ggplot2
#' @importFrom patchwork plot_annotation
#' @export
build_kobe_01 <- function(all_models,
                          scenario_name,
                          layout = c("3x2", "3x1"),
                          logax = FALSE, plot.legend = TRUE, man.legend = TRUE,
                          ext = TRUE, rel.axes = FALSE,
                          xlim = NULL, ylim = NULL,
                          xlabel = NULL, stamp = NULL, verbose = TRUE, CI = 0.95) {
  layout <- match.arg(layout)
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("patchwork is required.")
  if (!scenario_name %in% names(all_models)) {
    stop("Scenario '", scenario_name, "' not found in all_models.")
  }

  sc <- all_models[[scenario_name]]
  if (!is.list(sc) || length(sc) == 0) stop("Empty scenario list: ", scenario_name)

  # Helper: recover the *actual* panel ID from names in the scenario list.
  get_id <- function(letter, dtype) {
    nm  <- names(sc)
    hit <- nm[grep(paste0(letter, "\\.", dtype, "$"), nm)]
    if (length(hit)) hit[1] else paste0(scenario_name, letter, ".", dtype)
  }

  if (layout == "3x2") {
    # Left column = SDM; Right column = GLM; Rows = P, S, F
    reps <- list(
      P_SDM = .pick_model(sc, "P", "SDM"),
      P_GLM = .pick_model(sc, "P", "GLM"),
      S_SDM = .pick_model(sc, "S", "SDM"),
      S_GLM = .pick_model(sc, "S", "GLM"),
      F_SDM = .pick_model(sc, "F", "SDM"),
      F_GLM = .pick_model(sc, "F", "GLM")
    )

    p_P_SDM <- add_id_label(kobe_safe(reps$P_SDM, logax, plot.legend, man.legend,
                                      ext, rel.axes, xlim, ylim, c(1, 1),
                                      xlabel, stamp, verbose, CI),
                            get_id("P", "SDM"))
    p_S_SDM <- add_id_label(kobe_safe(reps$S_SDM, logax, plot.legend, man.legend,
                                      ext, rel.axes, xlim, ylim, c(1, 1),
                                      xlabel, stamp, verbose, CI),
                            get_id("S", "SDM"))
    p_F_SDM <- add_id_label(kobe_safe(reps$F_SDM, logax, plot.legend, man.legend,
                                      ext, rel.axes, xlim, ylim, c(1, 1),
                                      xlabel, stamp, verbose, CI),
                            get_id("F", "SDM"))

    p_P_GLM <- add_id_label(kobe_safe(reps$P_GLM, logax, plot.legend, man.legend,
                                      ext, rel.axes, xlim, ylim, c(1, 1),
                                      xlabel, stamp, verbose, CI),
                            get_id("P", "GLM"))
    p_S_GLM <- add_id_label(kobe_safe(reps$S_GLM, logax, plot.legend, man.legend,
                                      ext, rel.axes, xlim, ylim, c(1, 1),
                                      xlabel, stamp, verbose, CI),
                            get_id("S", "GLM"))
    p_F_GLM <- add_id_label(kobe_safe(reps$F_GLM, logax, plot.legend, man.legend,
                                      ext, rel.axes, xlim, ylim, c(1, 1),
                                      xlabel, stamp, verbose, CI),
                            get_id("F", "GLM"))

    # Final 3x2 layout: SDM column | GLM column
    g <- (p_P_SDM | p_P_GLM) /
      (p_S_SDM | p_S_GLM) /
      (p_F_SDM | p_F_GLM)

    return(g)
  }

  # layout == "3x1": prefer SDM; if missing, fall back to GLM
  pick_best <- function(letter) {
    obj <- .pick_model(sc, letter, "SDM")
    if (is.null(obj)) obj <- .pick_model(sc, letter, "GLM")
    obj
  }
  id_best <- function(letter) {
    nm  <- names(sc)
    hit <- nm[grep(paste0(letter, "\\.(SDM|GLM)$"), nm)]
    if (length(hit)) hit[1] else paste0(scenario_name, letter)
  }

  repP <- pick_best("P"); idP <- id_best("P")
  repS <- pick_best("S"); idS <- id_best("S")
  repF <- pick_best("F"); idF <- id_best("F")

  pP <- add_id_label(kobe_safe(repP, logax, plot.legend, man.legend,
                               ext, rel.axes, xlim, ylim, c(1, 1),
                               xlabel, stamp, verbose, CI), idP)
  pS <- add_id_label(kobe_safe(repS, logax, plot.legend, man.legend,
                               ext, rel.axes, xlim, ylim, c(1, 1),
                               xlabel, stamp, verbose, CI), idS)
  pF <- add_id_label(kobe_safe(repF, logax, plot.legend, man.legend,
                               ext, rel.axes, xlim, ylim, c(1, 1),
                               xlabel, stamp, verbose, CI), idF)

  g <- pP / pS / pF
  # Optionally add a title if you like:
  # g <- g + patchwork::plot_annotation(title = paste0("Scenario ", scenario_name, " (3x1)"))
  g
}

# INTERNAL: robust saver with fallback to base png()
#'
#' First tries `ggplot2::ggsave()`. If that fails, falls back to
#' `grDevices::png()` + `print(plot)` + `dev.off()`.
#'
#' @param plot_obj A `ggplot` or `patchwork` object.
#' @param file_path Destination PNG path.
#' @param width,height Numeric; inches.
#' @param dpi Numeric; resolution (dots per inch).
#'
#' @return (Invisibly) `TRUE` on success; `FALSE` with a warning if both
#'   methods fail.
#'
#' @keywords internal
#' @noRd
#' @importFrom grDevices png dev.off
.save_plot_safely <- function(plot_obj, file_path, width, height, dpi) {
  ok <- try({
    ggplot2::ggsave(
      filename = file_path,
      plot     = plot_obj,
      device   = "png",
      width    = width,
      height   = height,
      units    = "in",
      dpi      = dpi,
      bg       = "white"  # helps with Beamer/Word transparency issues
    )
    TRUE
  }, silent = TRUE)

  if (!inherits(ok, "try-error")) return(invisible(TRUE))

  gr <- try({
    grDevices::png(
      filename = file_path,
      width    = width,
      height   = height,
      units    = "in",
      res      = dpi,
      type     = "cairo"
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

#' Build Kobe-plot grids for *all* scenarios (optionally save)
#'
#' @description
#' Batch runner across every element in `all_models`. For each scenario,
#' calls [elu2::build_kobe_01()] to compose the grid, stores the result in a
#' named list, and (optionally) saves a PNG per scenario.
#'
#' @param all_models Named list of scenarios. Each scenario is a list of fitted
#'   models (see [elu2::build_kobe_01()] for naming expectations).
#' @param layout `"3x2"` (default) or `"3x1"`. Passed to [elu2::build_kobe_01()].
#' @param save Logical; if `TRUE`, save each grid as a PNG to `out_dir`. Default `TRUE`.
#' @param out_dir Output directory for PNG files (created if needed).
#'   Default `file.path("FIG", "KobePhases")`.
#' @param width,height Numeric; **inches** for saved PNGs. Defaults depend on
#'   `layout`: for `"3x2"` use `15 x 6.5`, otherwise `6 x 10`.
#' @param dpi Numeric; dots per inch for outputs. Default `300`.
#' @param ... Additional arguments forwarded to [elu2::build_kobe_01()]
#'   and ultimately to [elu2::kobe_safe()] / [elu2::kobe_all_in_one_gg()].
#'
#' @return A named list of `patchwork` grids, one per scenario (invisible if
#'   you only care about saved files).
#'
#' @examples
#' \dontrun{
#' grids <- build_kobe_grids_scenario(all_models, layout = "3x2",
#'                                    save = TRUE, out_dir = "FIG/KobePhases")
#' grids$S1  # view Scenario 1 grid
#' }
#'
#' @import ggplot2
#' @export
build_kobe_grids_scenario <- function(all_models,
                                      layout = "3x2",
                                      save = TRUE,
                                      out_dir = file.path("FIG", "KobePhases"),
                                      width = if (layout == "3x2") 15 else 6,
                                      height = if (layout == "3x2") 6.5 else 10,
                                      dpi = 300,
                                      ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("patchwork is required.")

  # Create output directory if needed (avoid warnings on case-sensitive FS)
  if (save && !dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  scens <- names(all_models)
  out <- vector("list", length(scens))
  names(out) <- scens

  i <- 1L
  while (i <= length(scens)) {
    scen <- scens[i]

    # Build the patchwork object for this scenario
    g <- build_kobe_01(all_models, scen, layout = layout, ...)

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

