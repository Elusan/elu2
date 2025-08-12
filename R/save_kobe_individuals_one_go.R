#' Save all individual Kobe plots under FIG/KobeIndividuals
#'
#' @description
#' For each scenario in `all_models` and for each submodel (e.g., "S1P.SDM",
#' "S1S.GLM", etc.), this function:
#' 1) Ensures the directory `FIG/KobeIndividuals/<Scenario>/` exists
#' 2) Builds the plot via [elu2::kobe_safe()] + [elu2::add_id_label()]
#' 3) Saves a PNG named `<submodel>.png` inside that scenario subfolder.
#'
#' Returns a data.frame log (scenario, submodel, file_path, ok).
#'
#' @param all_models Named list of scenarios. Each `all_models[[scenario]]` is
#'   a named list of fitted objects. Names should contain the submodel ID you
#'   want on the PNG (e.g., "S1P.SDM", "S1F.GLM").
#' @param out_root Root output folder (created if needed). Default `"FIG/KobeIndividuals"`.
#' @param per_scenario_subfolders Logical; if `TRUE` (default), save each scenario
#'   under `FIG/KobeIndividuals/<Scenario>/`. If `FALSE`, save all PNGs flat in `out_root`.
#' @param width,height PNG size in **inches**. Default `6 x 5`.
#' @param dpi PNG resolution. Default `300`.
#' @param ... Passed to [elu2::kobe_safe()] (e.g., `rel.axes = FALSE`, `CI = 0.95`, etc.).
#'
#' @return A data.frame with columns: `scenario`, `submodel`, `file_path`, `ok` (logical).
#'
#' @examples
#' \dontrun{
#' log_df <- save_kobe_individuals(all_models,
#'                                 out_root = file.path("FIG", "KobeIndividuals"),
#'                                 rel.axes = FALSE, CI = 0.95)
#' head(log_df)
#' }
#'
#' @import ggplot2
#' @export
save_kobe_individuals_one_go <- function(all_models,
                                  out_root = file.path("FIG", "KobeIndividuals"),
                                  per_scenario_subfolders = TRUE,
                                  width = 6,
                                  height = 5,
                                  dpi = 300,
                                  ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")

  # Ensure FIG and out_root exist (create recursively)
  if (!dir.exists("FIG")) dir.create("FIG", recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(out_root)) dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

  scens <- names(all_models)
  if (is.null(scens) || !length(scens)) stop("all_models must be a named list of scenarios.")

  # Preallocate a simple log
  log_list <- vector("list", 0L)

  s <- 1L
  while (s <= length(scens)) {
    scen <- scens[s]
    scen_list <- all_models[[scen]]
    if (!is.list(scen_list) || !length(scen_list)) {
      s <- s + 1L
      next
    }

    # Scenario folder
    out_dir <- if (isTRUE(per_scenario_subfolders)) {
      file.path(out_root, scen)
    } else {
      out_root
    }
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    submods <- names(scen_list)
    if (is.null(submods) || !length(submods)) {
      s <- s + 1L
      next
    }

    m <- 1L
    while (m <= length(submods)) {
      subm <- submods[m]
      rep_obj <- scen_list[[m]]

      # Fallback ID if name is missing
      id_text <- if (!is.null(subm) && nzchar(subm)) subm else paste0(scen, "_model", m)

      # Sanitize a safe filename
      fname <- gsub("[^A-Za-z0-9._-]+", "_", paste0(id_text, ".png"))
      fpath <- file.path(out_dir, fname)

      ok <- save_kobe_plot_single(rep_obj,
                                  id_text = id_text,
                                  file_path = fpath,
                                  width = width,
                                  height = height,
                                  dpi = dpi,
                                  ...)

      log_list[[length(log_list) + 1L]] <- data.frame(
        scenario  = scen,
        submodel  = id_text,
        file_path = fpath,
        ok        = isTRUE(ok),
        stringsAsFactors = FALSE
      )

      m <- m + 1L
    }

    s <- s + 1L
  }

  log_df <- do.call(rbind, log_list)
  rownames(log_df) <- NULL
  return(log_df)
}
