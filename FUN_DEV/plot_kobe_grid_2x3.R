#' Unified 2×3 Kobe grid builder (single or multiple scenarios)
#'
#' @description
#' A single robust entry point that:
#' * If given a **single scenario list** (e.g., `list(S1P.SDM=..., S1S.SDM=..., ...)`),
#'   builds **one** 2×3 grid (rows = SDM/GLM; cols = P/S/F).
#' * If given a **multi-scenario list** (e.g., `list(S1=list(...), S2=list(...))`),
#'   iterates scenarios and builds **one** grid per scenario in subfolders.
#'
#' All rendering is delegated to your existing
#' \code{my_plot_kobe_all_management_scenario()} (unchanged).
#'
#' Each panel shows a small bold in-panel tag with the original list name
#' (e.g., "S1P.SDM"). Blank cells are rendered when a model is missing.
#'
#' @param x Either:
#'   - A **single scenario list** with up to six model fits named like
#'     `SxP.SDM`, `SxS.SDM`, `SxF.SDM`, `SxP.GLM`, `SxS.GLM`, `SxF.GLM`, OR
#'   - A **multi-scenario list**, where each element is such a single scenario list.
#' @param man.legend Logical; forwarded to \code{my_plot_kobe_all_management_scenario()}.
#' @param text_size Numeric; in-panel tag size.
#' @param save Logical; if TRUE, write PNG(s) to disk.
#' @param out_dir Character; parent output directory. For multi-scenario input,
#'   a subfolder per scenario is created.
#' @param width,height Device size (inches) for saving.
#' @param dpi Resolution used by \code{ggsave()}.
#' @param file_basename Optional base name (without extension). If NULL, a stable
#'   name is derived from the first entry (single scenario) or the scenario name (multi).
#' @param mode One of \code{"auto"}, \code{"single"}, \code{"multi"}.
#'   \code{"auto"} attempts to infer based on structure.
#' @param verbose Logical; print progress and skips.
#' @param ... Extra args forwarded to \code{my_plot_kobe_all_management_scenario()}.
#'
#' @return
#' * If single scenario: invisibly returns \code{list(plot = patchwork_obj, file = path_or_NULL)}.
#' * If multi-scenario: invisibly returns a named list of such lists (one per scenario).
#' @export
plot_kobe_grid_2x3 <- function(x,
                               man.legend    = FALSE,
                               text_size     = 4,
                               save          = TRUE,
                               out_dir       = "FIG/SS_PP",
                               width         = 20,
                               height        = 8,
                               dpi           = 300,
                               file_basename = NULL,
                               mode          = c("auto", "single", "multi"),
                               verbose       = TRUE,
                               ...) {
  mode <- match.arg(mode)

  # ---- deps ----
  if (!requireNamespace("ggplot2",  quietly = TRUE)) stop("Package 'ggplot2' is required.")
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("Package 'patchwork' is required.")

  # ---- layout spec (stable) ----
  data_types <- c("SDM", "GLM")  # rows
  mod_keys   <- c("P", "S", "F") # cols
  mod_names  <- c(P = "Pella", S = "Schaefer", F = "Fox")

  # ---- tiny helpers ----
  .is_named_list <- function(obj) is.list(obj) && !is.null(names(obj)) && length(obj) > 0
  .has_model_like_names <- function(nms) any(grepl("^[S]\\w*[PSF]\\.(SDM|GLM)$", nms))
  .has_any_expected <- function(obj) .is_named_list(obj) && .has_model_like_names(names(obj))

  .infer_mode <- function(obj) {
    if (!.is_named_list(obj)) return("single") # safest fallback: treat as single
    # If *all* elements look like single-scenario lists => multi
    els <- obj
    looks_single <- vapply(els, function(e) .has_any_expected(e), logical(1))
    if (length(els) >= 1 && all(looks_single, na.rm = TRUE)) return("multi")
    # If this very list itself looks like a single-scenario models_list:
    if (.has_any_expected(obj)) return("single")
    # Default to single to avoid surprising recursion
    "single"
  }

  .first_or_null <- function(x) if (length(x) >= 1) x[[1]] else NULL
  .find_entry_name <- function(models_list, key, dtype) {
    pat <- paste0(key, "\\.", dtype, "$")   # e.g., "P.SDM$"
    nm  <- names(models_list)
    hit <- nm[grep(pat, nm)]
    .first_or_null(hit)
  }

  .blank_panel <- function(tag) {
    ggplot2::ggplot() +
      ggplot2::geom_blank() +
      ggplot2::expand_limits(x = c(0, 1), y = c(0, 1)) +
      ggplot2::theme_void() +
      ggplot2::theme(
        panel.border = ggplot2::element_rect(color = "grey80", fill = NA, linewidth = 0.6),
        plot.margin  = ggplot2::margin(4, 4, 4, 4)
      ) +
      ggplot2::annotate("text",
                        x = -Inf, y =  Inf, label = tag,
                        hjust = -0.1, vjust = 2,
                        fontface = "bold", size = text_size) +
      ggplot2::coord_cartesian(clip = "off")
  }

  .add_inpanel_tag <- function(p, tag) {
    p +
      ggplot2::annotate("text",
                        x = -Inf, y =  Inf, label = tag,
                        hjust = -0.1, vjust = 2,
                        fontface = "bold", size = text_size) +
      ggplot2::coord_cartesian(clip = "off")
  }

  .derive_single_basename <- function(models_list, fallback = "Scenario") {
    nm_all <- names(models_list)
    base <- if (length(nm_all) >= 1) sub("([PSF]).*$", "", nm_all[1]) else fallback
    if (!nzchar(base)) base <- fallback
    paste0(base, "KobeGrid_2x3")
  }

  .build_single_grid <- function(models_list,
                                 man.legend, text_size, save, out_dir, width, height, dpi,
                                 file_basename, verbose, ...) {
    if (!.is_named_list(models_list)) stop("Single scenario: `models_list` must be a named list.")
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    plots <- vector("list", length(data_types) * length(mod_keys))
    idx <- 0L

    for (dtype in data_types) {
      for (key in mod_keys) {
        idx <- idx + 1L
        nm <- .find_entry_name(models_list, key, dtype)
        if (!is.null(nm) && !is.null(models_list[[nm]])) {
          rep_obj <- models_list[[nm]]
          # Delegate: DO NOT alter your underlying plotting logic
          p <- my_plot_kobe_all_management_scenario(rep = rep_obj,
                                                    man.legend = man.legend, ...)
          plots[[idx]] <- .add_inpanel_tag(p, nm)
        } else {
          plots[[idx]] <- .blank_panel(sprintf("(%s – %s)", dtype, mod_names[[key]]))
        }
      }
    }

    g_row1 <- plots[[1]] | plots[[2]] | plots[[3]]
    g_row2 <- plots[[4]] | plots[[5]] | plots[[6]]
    g <- g_row1 / g_row2

    out_path <- NULL
    if (isTRUE(save)) {
      if (is.null(file_basename) || !nzchar(file_basename)) {
        file_basename <- .derive_single_basename(models_list)
      }
      out_path <- file.path(out_dir, paste0(file_basename, ".png"))
      ggplot2::ggsave(filename = out_path, plot = g,
                      width = width, height = height, dpi = dpi, units = "in")
      if (isTRUE(verbose)) message(sprintf("✔ Saved: %s", out_path))
    }

    print(g)
    invisible(list(plot = g, file = out_path))
  }

  # ---- Decide mode ----
  mode_eff <- if (mode == "auto") .infer_mode(x) else mode

  # ---- Execute ----
  if (mode_eff == "single") {
    return(.build_single_grid(models_list   = x,
                              man.legend    = man.legend,
                              text_size     = text_size,
                              save          = save,
                              out_dir       = out_dir,
                              width         = width,
                              height        = height,
                              dpi           = dpi,
                              file_basename = file_basename,
                              verbose       = verbose,
                              ...))
  }

  # multi-scenario branch
  if (!.is_named_list(x)) stop("Multi-scenario: `x` must be a named list of scenario lists.")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  results <- vector("list", length(x))
  names(results) <- names(x)

  for (sc in names(x)) {
    models_list <- x[[sc]]
    if (!.is_named_list(models_list)) {
      if (isTRUE(verbose)) message(sprintf("[SKIP] %s: not a named list.", sc))
      results[[sc]] <- NULL
      next
    }
    if (!.has_any_expected(models_list)) {
      if (isTRUE(verbose)) message(sprintf("[SKIP] %s: no expected entries found.", sc))
      results[[sc]] <- NULL
      next
    }

    sc_out_dir <- file.path(out_dir, sc)
    if (!dir.exists(sc_out_dir)) dir.create(sc_out_dir, recursive = TRUE, showWarnings = FALSE)

    base <- if (is.null(file_basename) || !nzchar(file_basename)) paste0(sc, "_KobeGrid_2x3") else file_basename
    if (isTRUE(verbose)) message(sprintf("[PLOT] %s → %s/%s.png", sc, sc_out_dir, base))

    results[[sc]] <- .build_single_grid(models_list   = models_list,
                                        man.legend    = man.legend,
                                        text_size     = text_size,
                                        save          = save,
                                        out_dir       = sc_out_dir,
                                        width         = width,
                                        height        = height,
                                        dpi           = dpi,
                                        file_basename = base,
                                        verbose       = verbose,
                                        ...)
  }

  invisible(results)
}
