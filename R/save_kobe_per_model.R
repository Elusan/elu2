#' Save individual Kobe plots for one or more scenarios
#'
#' @description
#' Iterates through the fitted models in each scenario of `all_models`,
#' builds a Kobe plot via [elu2::kobe_safe()], and saves one PNG per model
#' under `FIG/per_model/<scenario>/`. Returns the file paths invisibly.
#'
#' @param all_models Named list of scenarios. Each scenario is a named list of
#'   fitted model objects accepted by [elu2::kobe_all_in_one_gg()].
#' @param scenario_names Character vector of scenario names to process.
#'   Default: `names(all_models)` (process all).
#' @param out_dir Base directory for per-model PNGs (subdirs per scenario
#'   are created automatically). Default `file.path("FIG","per_model")`.
#' @param width,height Numeric; inches for saved PNGs. Default `6 x 5`.
#' @param dpi Numeric; dots per inch for PNGs. Default `300`.
#' @param add_id Logical; add top-left ID label using [elu2::add_id_label()].
#'   Default `TRUE`.
#' @param ... Additional arguments forwarded to [elu2::kobe_safe()]
#'   (e.g., `rel.axes`, `CI`, `logax`, etc.).
#'
#' @return Invisibly returns a named list: for each scenario, a named character
#'   vector of file paths keyed by model name.
#'
#' @examples
#' \dontrun{
#' save_kobe_per_model(all_models, scenario_names = c("S1","S2"),
#'                     rel.axes = FALSE, CI = 0.95)
#' }
#' @export
save_kobe_per_model <- function(all_models,
                                scenario_names = names(all_models),
                                out_dir = file.path("FIG", "per_model"),
                                width = 6, height = 5, dpi = 300,
                                add_id = TRUE,
                                ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
  if (is.null(scenario_names) || !length(scenario_names)) {
    stop("No scenarios to process: 'scenario_names' is empty.")
  }
  missing_sc <- setdiff(scenario_names, names(all_models))
  if (length(missing_sc)) {
    stop("These scenarios are not in 'all_models': ", paste(missing_sc, collapse = ", "))
  }

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  out <- vector("list", length(scenario_names))
  names(out) <- scenario_names

  si <- 1L
  while (si <= length(scenario_names)) {
    scen <- scenario_names[si]
    sc_dir <- file.path(out_dir, scen)
    if (!dir.exists(sc_dir)) dir.create(sc_dir, recursive = TRUE, showWarnings = FALSE)

    scen_list <- all_models[[scen]]
    if (!is.list(scen_list) || !length(scen_list)) {
      warning("Empty scenario '", scen, "'. Skipping.")
      out[[si]] <- character(0)
      si <- si + 1L
      next
    }

    fpaths <- character(length(scen_list))
    names(fpaths) <- names(scen_list)

    mi <- 1L
    while (mi <= length(scen_list)) {
      mname <- names(scen_list)[mi]
      rep_obj <- scen_list[[mi]]

      p <- elu2::kobe_safe(rep_obj, ...)   # single Kobe plot
      if (isTRUE(add_id)) {
        p <- elu2::add_id_label(p, mname, id_size = 4)
      }

      fpath <- file.path(sc_dir, paste0(mname, ".png"))
      # robust save
      ok <- try({
        ggplot2::ggsave(filename = fpath, plot = p, device = "png",
                        width = width, height = height, units = "in",
                        dpi = dpi, bg = "white")
        TRUE
      }, silent = TRUE)
      if (inherits(ok, "try-error")) {
        gr <- try({
          grDevices::png(filename = fpath, width = width, height = height,
                         units = "in", res = dpi, type = "cairo")
          print(p)
          grDevices::dev.off()
          TRUE
        }, silent = TRUE)
        if (inherits(gr, "try-error")) warning("Failed to save plot: ", fpath)
      }
      fpaths[mi] <- fpath
      mi <- mi + 1L
    }

    out[[si]] <- fpaths
    si <- si + 1L
  }

  invisible(out)
}
