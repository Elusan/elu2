# R/diag_pairs_all_osar.R
# ------------------------------------------------------------------------------

#' Ensure OSAR residuals exist on a fitted object
#' @keywords internal
ensure_osar_resid <- function(fit) {
  # If it's not a list-like object, just return it unchanged
  if (is.null(fit) || !is.list(fit)) return(fit)
  # If already has OSAR, done
  if ("osar" %in% names(fit)) return(fit)
  # Try to compute OSAR; swallow errors and return original if it fails
  out <- try(calc.osa.resid(fit), silent = TRUE)
  if (!inherits(out, "try-error")) fit <- out
  fit
}

#' Locate the OSAR diagnostic plotting function
#' @keywords internal
.get_osar_diag_fun <- function() {
  if (exists("plotspict.diagnostic.osar_gg", mode = "function"))
    return(get("plotspict.diagnostic.osar_gg"))
  # Optional alternate name hook (mirror your process helper pattern)
  if (exists("plot.elu2.diagnostic.osar_gg", mode = "function"))
    return(get("plot.elu2.diagnostic.osar_gg"))
  stop("Cannot find 'plotspict.diagnostic.osar_gg()' (or 'plot.elu2.diagnostic.osar_gg()').")
}

#' Side-by-side OSAR diagnostics for SDM and GLM
#'
#' Calls your OSAR diagnostic function on each fit and arranges the
#' two outputs side-by-side. Adds a single centered red label drawn at the
#' very top of each half-plot canvas (using cowplot, as in your example).
#'
#' @export
diag_pair_side_by_side_osar <- function(fit_sdm, fit_glm,
                                        label_sdm = "SDM", label_glm = "GLM",
                                        lag.max = 4, qlegend = TRUE, plot.data = TRUE,
                                        stamp = get.version(),
                                        save = FALSE, file = NULL,
                                        width = 16, height = 10, dpi = 400) {
  diag_fun <- .get_osar_diag_fun()

  # Make sure OSAR residuals are present
  fit_sdm <- ensure_osar_resid(fit_sdm)
  fit_glm <- ensure_osar_resid(fit_glm)

  # Build each side
  p_sdm <- diag_fun(fit_sdm, lag.max = lag.max, qlegend = qlegend,
                    plot.data = plot.data, stamp = stamp)
  p_glm <- diag_fun(fit_glm, lag.max = lag.max, qlegend = qlegend,
                    plot.data = plot.data, stamp = stamp)

  # Top labels (same cowplot approach as your process-based helper)
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    stop("Package 'cowplot' is required for top labels. Please install it.")
  }

  p_sdm_lab <- cowplot::ggdraw(p_sdm) +
    cowplot::draw_label(
      label_sdm, x = 0.35, y = 0.995,
      vjust = 1, hjust = 0.5, fontface = "bold", color = "red", size = 12
    )
  p_glm_lab <- cowplot::ggdraw(p_glm) +
    cowplot::draw_label(
      label_glm, x = 0.35, y = 0.995,
      vjust = 1, hjust = 0.5, fontface = "bold", color = "red", size = 12
    )

  # Side-by-side
  paired <- p_sdm_lab | p_glm_lab

  if (save) {
    if (is.null(file)) {
      out_dir <- file.path("FIG", "Diag_OSAR_SDM_GLM")
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      file <- file.path(out_dir, paste0("DiagnosticsOSAR_", gsub("\\.SDM$","", label_sdm), "_pair.png"))
    }
    ggplot2::ggsave(filename = file, plot = paired,
                    width = width, height = height, dpi = dpi, limitsize = FALSE)
    return(invisible(normalizePath(file, winslash = "/", mustWork = FALSE)))
  }

  paired
}

#' Batch: build OSAR SDM–GLM pairs for every scenario/model in a nested list
#'
#' Expects a nested \code{all_models} like yours:
#' \preformatted{
#' all_models <- list(
#'   S1 = list(S1P.GLM=..., S1S.GLM=..., S1F.GLM=..., S1P.SDM=..., S1S.SDM=..., S1F.SDM=...),
#'   S2 = list(...), S3 = list(...), S4 = list(...)
#' )
#' }
#' For each scenario and each base key (e.g., "S1P", "S1S", "S1F"), it:
#'   1) ensures OSAR residuals exist,
#'   2) builds the side-by-side SDM/GLM OSAR diagnostic,
#'   3) saves to \code{FIG/Diag_OSAR_SDM_GLM/<Scenario>/DiagnosticsOSAR_<Scenario>_<Key>_pair.png}.
#'
#' @param all_models Nested named list exactly as above.
#' @param keys Optional subset of base keys per scenario (e.g., c("S1P","S1S","S1F")).
#'             If NULL, inferred from names present.
#' @param save,width,height,dpi Passed to ggsave via \code{diag_pair_side_by_side_osar()}.
#' @param verbose Print where each file is written.
#' @return A \code{data.frame} with columns: \code{scenario}, \code{key}, \code{path}, \code{exists}.
#' @export
make_osar_diag_pairs_for_all <- function(all_models,
                                         keys = NULL,
                                         save = TRUE,
                                         width = 18, height = 10, dpi = 400,
                                         verbose = TRUE) {
  if (!length(all_models)) stop("`all_models` is empty.")
  scn_names <- names(all_models)
  if (is.null(scn_names) || any(!nzchar(scn_names))) stop("`all_models` must be a named list of scenarios.")

  recs <- list()

  for (scn in scn_names) {
    mods <- all_models[[scn]]
    nm   <- names(mods)
    if (is.null(nm) || !length(nm)) {
      warning(sprintf("Scenario '%s' is empty; skipping.", scn))
      next
    }

    # Infer base keys inside scenario (e.g., "S1P" from "S1P.SDM" / "S1P.GLM")
    base_keys <- sort(unique(sub("\\.[^.]+$", "", nm)))
    if (!is.null(keys)) base_keys <- intersect(base_keys, keys)
    if (!length(base_keys)) {
      warning(sprintf("Scenario '%s': no base keys found; skipping.", scn))
      next
    }

    # Ensure output dir per scenario
    out_dir <- file.path("FIG", "Diag_OSAR_SDM_GLM", scn)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    out_dir_abs <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)

    for (k in base_keys) {
      lab_sdm <- paste0(k, ".SDM")
      lab_glm <- paste0(k, ".GLM")

      if (!(lab_sdm %in% nm) && !(lab_glm %in% nm)) {
        warning(sprintf("Scenario '%s': no SDM/GLM entries for key '%s'; skipping.", scn, k))
        next
      }
      if (!(lab_sdm %in% nm)) warning(sprintf("Scenario '%s': missing %s; will use placeholder.", scn, lab_sdm))
      if (!(lab_glm %in% nm)) warning(sprintf("Scenario '%s': missing %s; will use placeholder.", scn, lab_glm))

      fit_sdm <- if (lab_sdm %in% nm) mods[[lab_sdm]] else structure(list(), names = character(0))
      fit_glm <- if (lab_glm %in% nm) mods[[lab_glm]] else structure(list(), names = character(0))

      file_out <- file.path(out_dir_abs, sprintf("DiagnosticsOSAR_%s_%s_pair.png", scn, k))
      invisible(
        try(
          diag_pair_side_by_side_osar(
            fit_sdm = fit_sdm,
            fit_glm = fit_glm,
            label_sdm = lab_sdm, label_glm = lab_glm,
            save = save, file = file_out,
            width = width, height = height, dpi = dpi
          ),
          silent = TRUE
        )
      )

      ok <- !save || file.exists(file_out)
      if (verbose) message(sprintf("%s %s", if (ok) "Saved:" else "Save failed:", file_out))
      recs[[length(recs) + 1]] <- data.frame(
        scenario = scn, key = k, path = file_out, exists = ok, stringsAsFactors = FALSE
      )
    }
  }

  if (length(recs)) do.call(rbind, recs) else
    data.frame(scenario = character(0), key = character(0), path = character(0), exists = logical(0))
}

# ------------------------------------------------------------------------------

# EXAMPLE USAGE (mirrors your process-based example)
# 1) (Optional) Pre-compute OSAR residuals so subsequent plotting is fast/robust:
# for (scn in names(all_models)) {
#   for (nm in names(all_models[[scn]])) {
#     all_models[[scn]][[nm]] <- ensure_osar_resid(all_models[[scn]][[nm]])
#   }
# }
#
# 2) Build & save every OSAR pair (SDM | GLM) for S1..S4 and P/S/F keys.
# res_osar <- make_osar_diag_pairs_for_all(all_models, verbose = TRUE)
# res_osar
