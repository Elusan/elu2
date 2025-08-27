#' @title Minimal Compact Theme for ggplot2
#'
#' @description
#' A lightweight, publication-style theme for dense multi-panel figures.
#' Removes grid lines, draws a visible panel border, and uses compact bold
#' text. Designed to pair with prior–posterior panels, but generally useful
#' for any small panels or patchwork grids.
#'
#' @param base_size Numeric. Base font size. Default: \code{8}.
#' @param base_family Character. Base font family. Default: \code{""} (use device default).
#'
#' @return A \code{ggplot2} theme object to be added with \code{+}.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(iris, aes(Sepal.Length, Sepal.Width, color = Species)) +
#'   geom_point() +
#'   theme_minimal_compact()
#' }
#'
#' @seealso \code{\link[ggplot2]{theme_minimal}}
#' @importFrom grid unit
#' @importFrom ggplot2 margin element_text element_rect element_blank theme_minimal
#' @export
theme_minimal_compact_this <- function(base_size = 8, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.title = element_text(face = "bold", size = 8),
      axis.text = element_text(size = 8, face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 10, face = "bold"),
      legend.key.size = unit(0.6, "lines"),
      legend.spacing.y = unit(0, "pt"),
      legend.spacing.x = unit(0, "pt"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, colour = "grey45", linewidth = 1.5),
      axis.ticks = element_line(linewidth = 0.5, color = "grey45"),
      axis.ticks.length = unit(3, "pt"),
      strip.background = element_rect(fill = "grey45", color = "grey45", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = rel(1)),
      text = element_text(face = "bold", size = 10),
      plot.margin = margin(2, 2, 2, 2)
    )
}

#' @title Prior–Posterior Grid (rows = models, cols = priors)
#'
#' @description
#' Builds a patchwork grid where **each row** corresponds to a model and
#' **each column** corresponds to a prior parameter (with multi-row priors
#' expanded in numeric order, e.g., \code{logsdi1}, \code{logsdi2}). The
#' function calls a panel-builder (by default \code{priors_as_named_list})
#' for each model, aligns columns across models by parameter name, and
#' optionally saves the resulting grid to disk.
#'
#' @param models Named list of fitted ELU/SPiCT model objects. Names are used
#'   to label rows and to infer a default file prefix.
#' @param priors_fun Function that takes a single model and returns a **named
#'   list of ggplot objects** (one per active prior). Default:
#'   \code{priors_as_named_list}.
#' @param width,height Numbers. Output image size (in inches) if saving. Defaults:
#'   \code{width = 12}, \code{height = 11}.
#' @param dpi Integer. Output resolution (dots per inch) if saving. Default: \code{300}.
#' @param outdir Character or \code{NULL}. Output directory to save PNG. If \code{NULL},
#'   nothing is written. Default: \code{"FIG"} (created if missing).
#' @param file_prefix Character or \code{NULL}. File name stem for the PNG. If \code{NULL}
#'   or empty, a prefix is inferred from the first model name (scenario part).
#' @param return_patchwork Logical. If \code{TRUE}, return the combined patchwork object;
#'   otherwise return \code{NULL} (invisibly). Default: \code{TRUE}.
#'
#' @return A \pkg{patchwork} object (if \code{return_patchwork = TRUE}) or
#'   \code{NULL} invisibly after saving.
#'
#' @details
#' Column order is derived from the union of prior names across all models,
#' arranged by each model's \code{$inp$priorsuseflags} order, with multi-row
#' priors expanded as \code{name1}, \code{name2}, ... in numeric order.
#' Missing cells (a model without a given prior) are filled with a spacer.
#'
#' @examples
#' \dontrun{
#' g <- elu_prior_posterior_grid_999(models_all$S1,
#'                                   priors_fun = priors_as_named_list,
#'                                   outdir = "FIG",
#'                                   file_prefix = "prior_grid_S1")
#' g
#' }
#'
#' @seealso \code{\link{priors_as_named_list}}, \code{\link{theme_minimal_compact}}
#'
#' @importFrom patchwork wrap_plots plot_spacer
#' @importFrom ggplot2 theme_void ggsave
#' @export
elu_prior_posterior_grid_999 <- function(models,
                                         priors_fun = priors_as_named_list,  # use adapter above
                                         width = 12, height = 11, dpi = 300,
                                         outdir = "FIG", file_prefix = NULL,
                                         return_patchwork = TRUE) {
  require(patchwork)
  require(ggplot2)

  # ---- filename prefix (same logic as yours) ----
  scenario_prefix <- ""
  model_names <- names(models)
  if (length(model_names) > 0) {
    scenario_prefix <- sub("([A-Za-z]+)$", "", model_names[1])
    scenario_prefix <- gsub("_+$", "", scenario_prefix)
  }
  if (is.null(file_prefix) || file_prefix == "") {
    file_prefix <- paste0("prior_posterior_grid_", scenario_prefix)
  }

  # ---- 1) Build named plot lists for each model (actual keys: logsdi1, logq2, logsdc, ...) ----
  named_plots_by_model <- lapply(seq_along(models), function(i) {
    priors_fun(models[[i]], model_id = names(models)[i], do.plot = NULL, stamp = NULL)
  })
  names(named_plots_by_model) <- names(models)

  # Union of all plot names across models
  plot_name_union <- unique(unlist(lapply(named_plots_by_model, names)))

  # ---- 2) Base prior order from priorsuseflags (keeps your intended ordering) ----
  active_base_order <- unique(unlist(lapply(models, function(rep) {
    uf <- rep$inp$priorsuseflags
    names(rep$inp$priors)[which(uf == 1)]
  })))

  base_of <- function(x) sub("\\d+$", "", x)                       # "logsdi1" -> "logsdi"
  sufnum  <- function(x) suppressWarnings(as.integer(sub(".*?(\\d+)$", "\\1", x)))  # "logsdi12" -> 12

  # Expand each base prior into its numbered variants in numeric order
  all_plot_names <- character(0)
  for (bp in active_base_order) {
    hits <- plot_name_union[base_of(plot_name_union) == bp]
    if (length(hits)) {
      # Separate numbered and unnumbered (e.g., "logsdc" is unnumbered)
      numbered   <- hits[grepl("\\d+$", hits)]
      unnumbered <- setdiff(hits, numbered)
      if (length(numbered)) {
        numbered <- numbered[order(sufnum(numbered))]
        all_plot_names <- c(all_plot_names, numbered)
      }
      if (length(unnumbered)) {
        all_plot_names <- c(all_plot_names, unnumbered)
      }
    }
  }

  # Safety: append any leftover names not in active_base_order (rare)
  remaining <- setdiff(plot_name_union, all_plot_names)
  if (length(remaining)) all_plot_names <- c(all_plot_names, remaining)

  # Edge case: nothing to plot
  if (!length(all_plot_names)) {
    message("No active priors across the supplied models. Nothing to plot.")
    if (return_patchwork) {
      return(patchwork::wrap_plots(patchwork::plot_spacer() + theme_void(), nrow = 1, ncol = 1))
    } else {
      return(invisible(NULL))
    }
  }

  n_models <- length(models)
  n_cols   <- length(all_plot_names)

  # ---- 3) Cell getter: exact name match; else spacer ----
  get_plot_or_spacer <- function(i_model, plot_name) {
    np <- named_plots_by_model[[i_model]]
    if (!is.null(np) && plot_name %in% names(np)) return(np[[plot_name]])
    patchwork::plot_spacer() + theme_void()
  }

  # ---- 4) Build grid: models as rows, expanded plot names as columns ----
  plots_grid <- lapply(seq_len(n_models), function(i) {
    lapply(all_plot_names, function(pn) get_plot_or_spacer(i, pn))
  })
  plots_vec <- unlist(plots_grid, recursive = FALSE)
  final_patchwork <- wrap_plots(plots_vec, nrow = n_models, ncol = n_cols)

  # ---- 5) Save & return ----
  if (!is.null(outdir)) {
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
    out_file <- file.path(outdir, paste0(file_prefix, ".png"))
    ggsave(out_file, final_patchwork, width = width, height = height, dpi = dpi)
  }
  if (return_patchwork) final_patchwork else invisible(NULL)
}
