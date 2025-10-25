#' Locate a hindcast extractor function
#'
#' Tries to find `extract.hindcast.info()` in the current session. It first
#' checks the global environment, then the **spict** namespace. If not found,
#' it throws a clear error explaining how to proceed.
#'
#' @return A function object corresponding to `extract.hindcast.info()`.
#' @details
#' This is a lightweight resolver so callers don’t have to know whether the
#' extractor lives in your wrapper package or in `spict`.
#'
#' @seealso spict::hindcast
#' @keywords internal
#' @noRd
.get_hc_extractor <- function() {
  if (exists("extract.hindcast.info", mode = "function")) return(get("extract.hindcast.info"))
  if (requireNamespace("spict", quietly = TRUE) &&
      exists("extract.hindcast.info", where = asNamespace("spict"), mode = "function")) {
    return(get("extract.hindcast.info", envir = asNamespace("spict")))
  }
  stop("extract.hindcast.info() not found. Load spict or your wrapper.")
}

#' Ensure a model contains hindcast results
#'
#' Returns a copy of `rep` that contains a computed `rep$hindcast`. If
#' `rep$hindcast` already exists, the input is returned unchanged. Otherwise
#' it calls `spict::hindcast()` with the supplied parameters.
#'
#' @param rep A SPiCT result list (as returned by `fit.spict()` or equivalent).
#' @param npeels Integer; number of hindcast peels (default `2L`).
#' @param peel.dtc Logical; whether to peel the data time coverage (default `FALSE`).
#' @param verbose Logical; message when hindcast is computed (default `FALSE`).
#'
#' @return The input `rep`, possibly augmented with a `hindcast` element.
#' @examples
#' \dontrun{
#' rep2 <- .ensure_hindcast(rep, npeels = 3L, peel.dtc = TRUE, verbose = TRUE)
#' }
#' @seealso spict::hindcast
#' @keywords internal
#' @noRd
.ensure_hindcast <- function(rep, npeels = 2L, peel.dtc = FALSE, verbose = FALSE) {
  if (is.null(rep) || !is.list(rep)) stop("Model is NULL or not a list.")
  if (!is.null(rep$hindcast) && length(rep$hindcast) >= 1L) return(rep)
  if (!requireNamespace("spict", quietly = TRUE)) stop("Need spict to compute hindcasts.")
  out <- try(spict::hindcast(rep, npeels = as.integer(npeels), peel.dtc = isTRUE(peel.dtc)), silent = TRUE)
  if (inherits(out, "try-error") || is.null(out$hindcast)) stop("hindcast() failed for a model.")
  if (isTRUE(verbose)) message("hindcast computed (npeels=", npeels, ", peel.dtc=", peel.dtc, ")")
  out
}

#' Preflight for hindcast plotting
#'
#' Validates that a model (1) has a hindcast and (2) the extractor returns
#' at least one index with data. It throws informative errors otherwise.
#'
#' @param rep A SPiCT result list.
#' @param CI Numeric scalar in (0,1); confidence level passed to the extractor
#'   (default `0.95`).
#' @param verbose Logical; passed to helper routines (default `FALSE`).
#'
#' @return A `rep` object guaranteed to contain valid hindcast content.
#' @examples
#' \dontrun{
#' rep_checked <- .preflight_hc(rep, CI = 0.9)
#' }
#' @seealso spict::hindcast
#' @keywords internal
#' @noRd
.preflight_hc <- function(rep, CI = 0.95, verbose = FALSE) {
  rep2 <- .ensure_hindcast(rep, verbose = verbose)
  extractor <- .get_hc_extractor()
  hcInfo <- try(extractor(rep2, CI = CI, verbose = verbose), silent = TRUE)
  if (inherits(hcInfo, "try-error") || is.null(hcInfo$index) || length(hcInfo$index) < 1)
    stop("Extractor returned no index data; check indices/iuse.")
  rep2
}

#' Hindcast grid (2×3) for a single scenario (SDM vs GLM × Pella/Schaefer/Fox)
#'
#' Builds a 2×3 panel grid for one scenario. Columns are **Pella**, **Schaefer**,
#' **Fox**; top row is **SDM**, bottom row is **GLM**. Each panel calls
#' `plotspict.hindcast_elu2_gg_exact_FOR_GRIDS()` and adds a bold tag in the
#' top-left corner (e.g., `"S1P"`). Missing models or plot errors are shown as
#' placeholder panels with an informative label. The grid title is intentionally
#' omitted (panel tags carry the identification).
#'
#' @param all_models Named list-of-lists by scenario; e.g. `all_models$S1$S1P.SDM`.
#'   Inside a scenario, expected names are: `SxP.SDM`, `SxS.SDM`, `SxF.SDM`,
#'   `SxP.GLM`, `SxS.GLM`, `SxF.GLM` (where `Sx` is the scenario name).
#' @param scenario Character; the scenario to plot (e.g., `"S1"`).
#' @param add.mase Logical; if `TRUE`, overlay MASE diagnostics when available.
#' @param CI Numeric scalar in (0,1); confidence level passed downstream.
#' @param verbose Logical; verbose messages from helpers.
#' @param show_legends Logical; keep per-panel legends (default `TRUE`).
#' @param title Ignored; kept for API compatibility (grid has no top title).
#' @param npeels Integer; hindcast peels (default `2L`).
#' @param peel.dtc Logical; data-time-coverage peeling (default `FALSE`).
#' @param save Logical; if `TRUE`, save a PNG next to returning the grid.
#' @param out_dir Output directory for PNG (created if needed).
#' @param width,height,dpi Graphics device parameters for saved PNG.
#'
#' @return A patchwork object representing the 2×3 grid. Invisibly saves a PNG
#'   when `save = TRUE`.
#'
#' @examples
#' \dontrun{
#' # Single scenario:
#' p <- plot_hindcast_grid_scenario123(
#'   all_models, scenario = "S1",
#'   npeels = 2L, peel.dtc = FALSE,
#'   add.mase = TRUE, CI = 0.95, verbose = FALSE,
#'   show_legends = TRUE, save = TRUE,
#'   out_dir = file.path("FIG", "Hindcast"),
#'   width = 14, height = 6, dpi = 300
#' )
#' }
#' @seealso spict::hindcast, plotspict.hindcast_elu2_gg_exact_FOR_GRIDS
#' @export
plot_hindcast_grid_scenario123 <- function(
    all_models,              # list: top-level names are scenarios (S1, S2, ...)
    scenario,                # e.g. "S1"
    add.mase = TRUE, CI = 0.95, verbose = FALSE,
    show_legends = TRUE,
    title = NULL,            # <- ignored now (no grid title)
    npeels = 2L, peel.dtc = FALSE,
    save = TRUE, out_dir = file.path("FIG", "Hindcast"),
    width = 14, height = 6, dpi = 300
){
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required.")
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("patchwork required.")

  if (!(scenario %in% names(all_models))) {
    stop("Scenario '", scenario, "' not found in all_models.")
  }
  scen_list <- all_models[[scenario]]
  if (!is.list(scen_list) || is.null(names(scen_list))) {
    stop("all_models[['", scenario, "']] must be a named list of models.")
  }
  # --------------------------------------------------------------------------
  # Put this inside plot_hindcast_grid_scenario123() where you define helpers
  .short_tag <- function(code) sub("\\..*$", "", code)  # "S1P.SDM" -> "S1P"

  # Use ggplot2's tag (safer than patchwork::plot_annotation for per-panel tags)
  .add_tag <- function(p, tag_txt) {
    p +
      ggplot2::labs(tag = tag_txt) +
      ggplot2::theme(
        plot.tag = ggplot2::element_text(face = "bold", size = 11),
        # (x,y) in [0,1] coordinates relative to plot area
        plot.tag.position = c(0, 1),         # left & top
        plot.margin = ggplot2::margin(t = 6, r = 6, b = 6, l = 6)
      )
  }

  # --------------------------------------------------------------------------

  .placeholder <- function(label, why = NULL) {
    lab <- if (is.null(why)) paste0(label, "\n(missing)") else paste0(label, "\n", why)
    ggplot2::ggplot() +
      ggplot2::geom_text(ggplot2::aes(0,0,label = lab), size = 4, fontface = 2, lineheight = 1.05) +
      ggplot2::xlim(-1,1) + ggplot2::ylim(-1,1) + ggplot2::theme_void()
  }

  .safe_panel <- function(mod, label) {
    if (is.null(mod)) return(.placeholder(label))
    mod2 <- try(.preflight_hc(.ensure_hindcast(mod, npeels = npeels, peel.dtc = peel.dtc, verbose = verbose),
                              CI = CI, verbose = verbose), silent = TRUE)
    if (inherits(mod2, "try-error")) return(.placeholder(label, why = sub("^.*?: ", "", as.character(mod2))))
    p <- try(plotspict.hindcast_elu2_gg_exact_FOR_GRIDS(mod2, add.mase = add.mase, CI = CI, verbose = verbose),
             silent = TRUE)
    if (inherits(p, "try-error") || is.null(p)) return(.placeholder(label, why = "plot error"))
    if (!isTRUE(show_legends)) p <- p + ggplot2::theme(legend.position = "none")
    # --- NEW: add top-left external tag with short code
    .add_tag(p, .short_tag(label))
  }

  .get <- function(code) if (code %in% names(scen_list)) scen_list[[code]] else NULL

  # Codes stored INSIDE the scenario list (exactly as you provided)
  codes_SDM <- c(P = paste0(scenario, "P.SDM"),
                 S = paste0(scenario, "S.SDM"),
                 F = paste0(scenario, "F.SDM"))
  codes_GLM <- c(P = paste0(scenario, "P.GLM"),
                 S = paste0(scenario, "S.GLM"),
                 F = paste0(scenario, "F.GLM"))

  p_SDM_P <- .safe_panel(.get(codes_SDM["P"]), codes_SDM["P"])
  p_SDM_S <- .safe_panel(.get(codes_SDM["S"]), codes_SDM["S"])
  p_SDM_F <- .safe_panel(.get(codes_SDM["F"]), codes_SDM["F"])

  p_GLM_P <- .safe_panel(.get(codes_GLM["P"]), codes_GLM["P"])
  p_GLM_S <- .safe_panel(.get(codes_GLM["S"]), codes_GLM["S"])
  p_GLM_F <- .safe_panel(.get(codes_GLM["F"]), codes_GLM["F"])

  grid <- (p_SDM_P | p_SDM_S | p_SDM_F) /
    (p_GLM_P | p_GLM_S | p_GLM_F)

  # --- CHANGED: remove grid title entirely (no S1–Hindcast, etc.)
  # grid <- grid + patchwork::plot_annotation(title = ttl)

  if (isTRUE(save)) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(file.path(out_dir, paste0("hindcast_grid_", scenario, ".png")),
                    grid, width = width, height = height, dpi = dpi)
  }
  grid
}

#' Batch hindcast grids across scenarios
#'
#' Convenience wrapper that loops over scenarios and calls
#' [plot_hindcast_grid_scenario123()] for each one. Returns a named list of
#' patchwork objects (invisible) and lets each call save its PNG if requested.
#'
#' @param all_models Named list-of-lists by scenario (see
#'   [plot_hindcast_grid_scenario123()] for expected inner names).
#' @param scenarios Optional character vector of scenario names to plot. If
#'   `NULL`, all names in `all_models` are used.
#' @param ... Passed through to [plot_hindcast_grid_scenario123()] (e.g.,
#'   `npeels`, `peel.dtc`, `CI`, `add.mase`, `save`, `out_dir`, `width`,
#'   `height`, `dpi`, etc.).
#'
#' @return A named list of patchwork objects, one per scenario (invisible).
#'
#' @examples
#' \dontrun{
#' # Batch scenarios S1–S4:
#' res <- plot_hindcast_grid_many(
#'   all_models,
#'   scenarios = paste0("S", 1:4),
#'   npeels = 2, peel.dtc = FALSE,
#'   CI = 0.95, add.mase = TRUE, verbose = FALSE,
#'   show_legends = TRUE, save = TRUE,
#'   out_dir = file.path("FIG", "Hindcast"),
#'   width = 14, height = 6, dpi = 300
#' )
#' }
#' @seealso plot_hindcast_grid_scenario123
#' @export
plot_hindcast_grid_many <- function(all_models, scenarios = NULL, ...) {
  stopifnot(is.list(all_models))
  if (is.null(scenarios)) scenarios <- names(all_models)
  scenarios <- scenarios[scenarios %in% names(all_models)]
  if (!length(scenarios)) stop("No matching scenarios in all_models.")

  out <- setNames(vector("list", length(scenarios)), scenarios)
  for (i in seq_along(scenarios)) {
    out[[i]] <- plot_hindcast_grid_scenario123(all_models, scenario = scenarios[i], ...)
  }
  invisible(out)
}
