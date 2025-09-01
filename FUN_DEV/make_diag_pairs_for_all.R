# R/diag_pairs_all.R

#' Ensure process residuals exist on a fitted object
#' @keywords internal
ensure_process_resid <- function(fit) {
  if (!("process.resid" %in% names(fit))) {
    out <- try(calc.process.resid(fit), silent = TRUE)
    if (!inherits(out, "try-error")) fit <- out
  }
  fit
}

#' Locate the diagnostic plotting function you already have
#' @keywords internal
.get_diag_fun <- function() {
  if (exists("plotspict.diagnostic.process_gg", mode = "function"))
    return(get("plotspict.diagnostic.process_gg"))
  if (exists("plot.elu2.diagnostic.process_gg", mode = "function"))
    return(get("plot.elu2.diagnostic.process_gg"))
  stop("Cannot find 'plotspict.diagnostic.process_gg()' (or 'plot.elu2.diagnostic.process_gg()').")
}

# ===== NEW: centered, top label (inside the panel), red, fully transparent inset =====
#' Add a single centered label at the top of a panel (transparent overlay)
#' @keywords internal
add_top_center_label <- function(p, label, size = 10, color = "red") {
  # A tiny, fully transparent overlay with one centered label near the top
  band <- ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 0.5, y = 0.97,                 # centered, safely inside the top edge
      label = label,
      hjust = 0.5, vjust = 1,
      size = size / 3,
      fontface = "bold",
      colour = color
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    ggplot2::theme_void() +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = NA, colour = NA),
      plot.background  = ggplot2::element_rect(fill = NA, colour = NA),
      plot.margin      = ggplot2::margin(0, 0, 0, 0)
    )

  # Add a very thin, transparent top band spanning the panel width
  p + patchwork::inset_element(
    band,
    left = 0.00, right = 1.00,
    bottom = 0.90, top = 0.99,
    align_to = "panel"
  )
}

# (Kept for reference but no longer used)
#' Add a small in-panel label at top-left (safe padding; transparent inset)
#' @keywords internal
add_inset_label <- function(p, label, size = 10) {
  lab <- ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 0.01, y = 0.965,
      label = label,
      hjust = 0, vjust = 1,
      size = size / 3,
      fontface = "bold"
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    ggplot2::theme_void() +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = NA, colour = NA),
      plot.background  = ggplot2::element_rect(fill = NA, colour = NA),
      plot.margin      = ggplot2::margin(0, 0, 0, 0)
    )

  p + patchwork::inset_element(
    lab,
    left = 0.012, right = 0.5,
    bottom = 0.885, top = 0.9,
    align_to = "panel"
  )
}

# R/diag_pairs_all.R

# ... keep ensure_process_resid() and .get_diag_fun() as-is ...

#' Side-by-side diagnostics for SDM and GLM
#'
#' Calls your existing diagnostic function on each fit and arranges the
#' two outputs side-by-side. Adds a single centered red label drawn at the
#' very top of each half-plot canvas (visually between the "Biomass" and
#' "Fishing mortality" columns).
#'
#' @export
diag_pair_side_by_side <- function(fit_sdm, fit_glm,
                                   label_sdm = "SDM", label_glm = "GLM",
                                   lag.max = 4, qlegend = TRUE, plot.data = TRUE,
                                   add.loess = FALSE, span = 0.75, stamp = get.version(),
                                   save = FALSE, file = NULL,
                                   width = 16, height = 10, dpi = 400) {
  diag_fun <- .get_diag_fun()

  fit_sdm <- ensure_process_resid(fit_sdm)
  fit_glm <- ensure_process_resid(fit_glm)

  p_sdm <- diag_fun(fit_sdm, lag.max = lag.max, qlegend = qlegend,
                    plot.data = plot.data, add.loess = add.loess,
                    span = span, stamp = stamp)

  p_glm <- diag_fun(fit_glm, lag.max = lag.max, qlegend = qlegend,
                    plot.data = plot.data, add.loess = add.loess,
                    span = span, stamp = stamp)

  # --- NEW: draw one top-centered label in the plot canvas (not in panels) ---
  # Require cowplot
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    stop("Package 'cowplot' is required for top labels. Please install it.")
  }

  p_sdm_lab <- cowplot::ggdraw(p_sdm) +
    cowplot::draw_label(
      label_sdm, x = 0.54, y = 0.995,            # centered, right at the top edge
      vjust = 1, hjust = 0.5,
      fontface = "bold", color = "red", size = 12
    )

  p_glm_lab <- cowplot::ggdraw(p_glm) +
    cowplot::draw_label(
      label_glm, x = 0.54, y = 0.995,
      vjust = 1, hjust = 0.5,
      fontface = "bold", color = "red", size = 12
    )

  # Place the two labeled halves side-by-side
  paired <- p_sdm_lab | p_glm_lab

  if (save) {
    if (is.null(file)) {
      out_dir <- file.path("FIG", "Diag_SDM_GLM")
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      file <- file.path(out_dir, paste0("Diagnostics_", gsub("\\.SDM$","", label_sdm), "_pair.png"))
    }
    ggplot2::ggsave(filename = file, plot = paired,
                    width = width, height = height, dpi = dpi, limitsize = FALSE)
    return(invisible(normalizePath(file, winslash = "/", mustWork = FALSE)))
  }

  paired
}



#' Batch: build SDM–GLM pairs for every scenario/model in a nested list
#'
#' Expects a nested \code{all_models} like you provided:
#' \preformatted{
#' all_models <- list(
#'   S1 = list(S1P.GLM=..., S1S.GLM=..., S1F.GLM=..., S1P.SDM=..., S1S.SDM=..., S1F.SDM=...),
#'   S2 = list(...), S3 = list(...), S4 = list(...)
#' )
#' }
#' For each scenario and each base key (e.g., S1P, S1S, S1F), it:
#' 1) ensures process residuals exist,
#' 2) builds the side-by-side SDM/GLM diagnostic,
#' 3) saves to \code{FIG/Diag_SDM_GLM/<Scenario>/Diagnostics_<Scenario>_<Key>_pair.png}.
#'
#' @param all_models Nested named list exactly as above.
#' @param keys Optional subset of base keys (e.g., c("S1P","S1S","S1F")) per scenario.
#'             If NULL, inferred from names present.
#' @param save,width,height,dpi Passed to ggsave through \code{diag_pair_side_by_side()}.
#' @param verbose Print where each file is written.
#' @return A \code{data.frame} with columns: \code{scenario}, \code{key}, \code{path}, \code{exists}.
#' @export
make_diag_pairs_for_all <- function(all_models,
                                    keys = NULL,
                                    save = TRUE,
                                    width = 16, height = 10, dpi = 400,
                                    verbose = TRUE) {
  if (!length(all_models)) stop("`all_models` is empty.")
  scn_names <- names(all_models)
  if (is.null(scn_names) || any(!nzchar(scn_names))) stop("`all_models` must be a named list of scenarios.")

  # results collector
  recs <- list()

  for (scn in scn_names) {
    mods <- all_models[[scn]]
    nm   <- names(mods)
    if (is.null(nm) || !length(nm)) {
      warning(sprintf("Scenario '%s' is empty; skipping.", scn))
      next
    }

    # infer base keys inside scenario (e.g., "S1P" from "S1P.SDM" / "S1P.GLM")
    base_keys <- sort(unique(sub("\\.[^.]+$", "", nm)))
    if (!is.null(keys)) base_keys <- intersect(base_keys, keys)
    if (!length(base_keys)) {
      warning(sprintf("Scenario '%s': no base keys found; skipping.", scn))
      next
    }

    # ensure output dir per scenario
    out_dir <- file.path("FIG", "Diag_SDM_GLM", scn)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    out_dir_abs <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)

    for (k in base_keys) {
      lab_sdm <- paste0(k, ".SDM")
      lab_glm <- paste0(k, ".GLM")

      if (!(lab_sdm %in% nm) && !(lab_glm %in% nm)) {
        warning(sprintf("Scenario '%s': no SDM/GLM entries for key '%s'; skipping.", scn, k))
        next
      }
      if (!(lab_sdm %in% nm)) {
        warning(sprintf("Scenario '%s': missing %s; will use placeholder.", scn, lab_sdm))
      }
      if (!(lab_glm %in% nm)) {
        warning(sprintf("Scenario '%s': missing %s; will use placeholder.", scn, lab_glm))
      }

      # falls back to placeholder if one side is missing
      fit_sdm <- if (lab_sdm %in% nm) mods[[lab_sdm]] else structure(list(), names = character(0))
      fit_glm <- if (lab_glm %in% nm) mods[[lab_glm]] else structure(list(), names = character(0))

      file_out <- file.path(out_dir_abs, sprintf("Diagnostics_%s_%s_pair.png", scn, k))
      invisible(
        try(
          diag_pair_side_by_side(
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
