# ======================================================================
# params_grid_extract.R  —  ELU2 utilities (extract + grid plot + helper)
# ======================================================================

# ----------------------------------------------------------------------
# Helper (internal): keep ONLY the requested year for ratio panels
# ----------------------------------------------------------------------

#' Keep a Single Year for Ratio Panels (Internal)
#'
#' Filters a tidy parameter table to keep **only one target year** for the
#' ratio parameters \code{"B/Bmsy"} and \code{"F/Fmsy"}. Non-ratio parameters
#' are left unchanged. Within each group
#' \code{(scenario, type, family, param)}, if \code{param} is a ratio:
#' rows with \code{year == keep_year} are retained; if multiple rows match,
#' the first is kept; if none match, the group is **dropped** (strict).
#'
#' Robustness:
#' \itemize{
#'   \item If \code{df} is \code{NULL} or has zero rows, it is returned unchanged.
#'   \item Does not modify non-ratio parameters.
#'   \item For ratio parameters with no \code{keep_year} present, the entire group is dropped.
#' }
#'
#' @param df A data frame containing at least the columns
#'   \code{scenario}, \code{type}, \code{family}, \code{param}, \code{year}.
#'   Additional columns (e.g., \code{est}, \code{lo}, \code{hi}) are preserved.
#' @param keep_year Integer scalar; the target year to keep for ratio parameters.
#'   Default \code{2022}.
#'
#' @return A data frame where ratio parameters contain only the specified
#'   \code{keep_year} (one row per ratio-group), and non-ratio parameters are
#'   unchanged. Row names are reset. If all ratio groups are dropped and no
#'   non-ratio rows remain, a 0-row data frame with the same columns is returned.
#'
#' @examples
#' \dontrun{
#' df2 <- .keep_ratio_year(param_df, keep_year = 2022)
#' }
#'
#' @keywords internal
#' @noRd
.keep_ratio_year <- function(df, keep_year = 2022) {
  if (is.null(df) || !is.data.frame(df) || !nrow(df)) return(df)

  out <- vector("list", 0L)
  keys <- paste(df$scenario, df$type, df$family, df$param, sep = "||")
  ukeys <- unique(keys)

  for (k in ukeys) {
    idx <- which(keys == k)
    block <- df[idx, , drop = FALSE]

    if (isTRUE(block$param[1] %in% c("B/Bmsy", "F/Fmsy"))) {
      keep_idx <- which(!is.na(block$year) & as.integer(block$year) == as.integer(keep_year))
      if (length(keep_idx) >= 1L) {
        block <- block[keep_idx, , drop = FALSE]
        if (nrow(block) > 1L) block <- block[1L, , drop = FALSE]
      } else {
        next  # strict: drop if requested year is absent
      }
    }
    out[[length(out) + 1L]] <- block
  }

  if (!length(out)) return(df[0, , drop = FALSE])
  res <- tryCatch(do.call(rbind, out), error = function(e) df[0, , drop = FALSE])
  rownames(res) <- NULL
  res
}


# ----------------------------------------------------------------------
# Extractor: K, r, MSY, B/Bmsy, F/Fmsy (+ CIs) from SDM + GLM SPiCT models
# ----------------------------------------------------------------------

#' Extract K, r, MSY, B/Bmsy, F/Fmsy (+ CIs) from SDM + GLM SPiCT models
#'
#' Iterates over a named list of fitted models and returns a tidy table of
#' parameter estimates and confidence intervals for the key quantities:
#' \code{K}, \code{r}, \code{MSY}, \code{B/Bmsy}, \code{F/Fmsy}.
#' The function expects your \code{make_spict_summary_df1_nopipe()} helper to be
#' available and uses it to standardize summaries across models.
#'
#' @section Expected input:
#' \itemize{
#'   \item \code{model_list}: a \strong{named} list. Names encode scenario/model, e.g. \code{"S1P"}, \code{"S1F"}, \code{"S1S"}.
#'   \item Each element is a fitted object compatible with \code{make_spict_summary_df1_nopipe()}.
#' }
#'
#' @param model_list A \strong{named} list of fitted model objects. Names are used
#'   to derive scenario/model family (e.g., \code{"S1P"} \textrightarrow{} scenario \code{"S1"}, family \code{"Pella"}).
#' @param scenarios Optional character vector of scenario labels to keep/order
#'   (e.g., \code{paste0("S", 1:7)}). If \code{NULL}, scenarios are inferred from names.
#' @param verbose Logical; if \code{TRUE} (default) print minimal progress messages.
#'
#' @return A tidy \code{data.frame} with columns including (at least):
#'   \code{scenario}, \code{type} (\code{"SDM"}/\code{"GLM"}), \code{family}
#'   (\code{"Pella"}/\code{"Fox"}/\code{"Schaefer"}), \code{param}
#'   (\code{"K"}, \code{"r"}, \code{"MSY"}, \code{"B/Bmsy"}, \code{"F/Fmsy"}),
#'   \code{year} (for ratios, may be \code{NA} otherwise), \code{est}, \code{lo}, \code{hi}.
#'
#' @details
#' \itemize{
#'   \item The function parses names like \code{"S1P"} to infer scenario \code{"S1"} and family
#'         (\code{"P"} \textrightarrow{} \code{"Pella"}, \code{"F"} \textrightarrow{} \code{"Fox"}, \code{"S"} \textrightarrow{} \code{"Schaefer"}).
#'   \item Uses your \code{make_spict_summary_df1_nopipe()} to produce standardized parameter summaries.
#'   \item Rows with non-finite estimates/intervals are dropped.
#' }
#'
#' @examples
#' \dontrun{
#' tbl <- extract_params_SDM_GLM(models_all, scenarios = paste0("S", 1:7))
#' head(tbl)
#' }
#'
#' @seealso \code{\link{plot_param_grid_3rows_by_family888}}
#' @export
extract_params_SDM_GLM <- function(model_list, scenarios = NULL, verbose = TRUE) {
  if (is.null(model_list) || !is.list(model_list) || is.null(names(model_list))) {
    stop("`model_list` must be a named list.", call. = FALSE)
  }

  family_from_code <- function(code) {
    z <- toupper(trimws(as.character(code)))
    if (nchar(z) > 1L) z <- substr(z, 1L, 1L)
    if (z == "F") "Fox" else if (z == "S") "Schaefer" else if (z == "P") "Pella" else NA_character_
  }
  parse_name <- function(nm) {
    base <- sub("\\..*$", "", nm)               # "S1F" from "S1F.GLM" etc.
    scen <- sub("^((S\\d+)).*$", "\\1", base)   # "S1"
    famc <- sub("^S\\d+([A-Za-z]).*$", "\\1", base)
    fam  <- family_from_code(famc)
    list(scenario = scen, family = fam)
  }

  out <- vector("list", length(model_list))
  nn  <- names(model_list)

  for (i in seq_along(model_list)) {
    nm <- nn[i]
    meta <- parse_name(nm)
    if (verbose) message("Summarizing: ", nm, " -> ", meta$scenario, "/", meta$family)

    if (!exists("make_spict_summary_df1_nopipe", mode = "function")) {
      stop("`make_spict_summary_df1_nopipe()` not found in scope.", call. = FALSE)
    }
    sumdf <- make_spict_summary_df1_nopipe(model_list[[i]])

    sumdf$scenario <- meta$scenario
    sumdf$family   <- meta$family
    out[[i]] <- sumdf
  }

  res <- do.call(rbind, out)
  rownames(res) <- NULL

  if (!is.null(scenarios)) {
    keep <- intersect(scenarios, unique(res$scenario))
    res  <- res[res$scenario %in% keep, , drop = FALSE]
    res$scenario <- factor(res$scenario, levels = keep)
    res <- res[order(res$scenario), , drop = FALSE]
    res$scenario <- as.character(res$scenario)
  }

  for (nm in c("est","lo","hi")) if (!nm %in% names(res)) res[[nm]] <- NA_real_
  ok <- is.finite(res$est) | (is.na(res$est) & is.na(res$lo) & is.na(res$hi))
  res <- res[ok, , drop = FALSE]

  res
}


# ----------------------------------------------------------------------
# Plot: 3 rows (families) × 5 cols (K, r, MSY, B/Bmsy, F/Fmsy) — SDM vs GLM
# ----------------------------------------------------------------------

#' Parameter Grid (SDM vs GLM): 3 Rows by Family × 5 Columns
#'
#' Produces a clean \strong{3 rows × 5 columns} grid of panels showing
#' \strong{SDM} and \strong{GLM} point estimates with confidence intervals
#' across scenarios, split by family. Rows correspond to families
#' (default: \code{"Pella"}, \code{"Fox"}, \code{"Schaefer"}) and columns to
#' \code{"K"}, \code{"r"}, \code{"MSY"}, \code{"B/Bmsy"}, \code{"F/Fmsy"}.
#'
#' Design choices:
#' \itemize{
#'   \item X-axis shows scenario indices \code{1..N} (no x-title).
#'   \item Only the \emph{rightmost} panel in each row gets a \emph{secondary y-axis title}
#'         displaying the family name (no ticks/labels on the secondary axis).
#'   \item White background, subtle panel borders, compact bold strip labels.
#'   \item For ratio panels (\code{B/Bmsy}, \code{F/Fmsy}), only a single year is kept
#'         (default \code{2022}) via \code{.keep_ratio_year()}.
#'   \item Defensive checks: required columns, packages, finite values; non-finite rows dropped.
#' }
#'
#' @section Expected input:
#' \code{param_df} must be a tidy data frame with columns:
#' \describe{
#'   \item{scenario}{Scenario ID (e.g., \code{"S1"}, \code{"S2"}, ...).}
#'   \item{type}{Model type: \code{"SDM"} or \code{"GLM"}.}
#'   \item{family}{One of \code{"Pella"}, \code{"Fox"}, \code{"Schaefer"} (others allowed).}
#'   \item{param}{One of \code{"K"}, \code{"r"}, \code{"MSY"}, \code{"B/Bmsy"}, \code{"F/Fmsy"}.}
#'   \item{year}{Numeric year for ratio parameters (may be \code{NA} for others).}
#'   \item{est, lo, hi}{Point estimate and CI bounds (numeric; non-finite dropped).}
#' }
#'
#' @param param_df Data frame as described above.
#' @param scenario_levels Optional character vector controlling scenario order on the x-axis
#'   (e.g., \code{paste0("S", 1:7)}). If \code{NULL}, order is inferred numerically from \code{scenario}.
#'   Levels not present are ignored; scenarios present but not listed are appended.
#' @param families_row Character vector giving the row order of families.
#'   Default \code{c("Pella","Fox","Schaefer")}. Families not present are dropped; present-but-unlisted
#'   families are appended at the end in alphabetical order.
#' @param params_cols Character vector giving the column order of parameters.
#'   Default \code{c("K","r","MSY","B/Bmsy","F/Fmsy")}. Parameters not present are dropped; present-but-unlisted
#'   parameters are appended at the end in alphabetical order.
#' @param keep_ratio_year Numeric year to keep for \code{B/Bmsy} and \code{F/Fmsy} panels. Default \code{2022}.
#' @param width,height,dpi Device size and resolution used only when saving to file.
#' @param output_dir Optional directory to save a PNG. If \code{NULL}, no file is saved.
#' @param filename File name for the PNG (used when \code{output_dir} is non-\code{NULL}).
#' @param ci_linewidth Numeric CI line thickness (default \code{0.45}).
#' @param ci_alpha Numeric CI opacity in \code{[0,1]} (default \code{0.95}).
#' @param point_size Numeric point size for estimates (default \code{2.4}).
#' @param dodge_width Horizontal dodge between SDM and GLM at each scenario (default \code{0.45}).
#' @param use_errorbar Logical; if \code{TRUE}, use \code{geom_errorbar()} (hairline option)
#'   instead of \code{geom_linerange()} for CIs (default \code{FALSE}).
#' @param errorbar_cap_width Width of errorbar caps (set \code{0} for hairlines; default \code{0}).
#'
#' @return A \pkg{patchwork} object representing the combined grid.
#'   If \code{output_dir} is provided, a PNG is also saved.
#'
#' @details
#' Internally, ratio panels are reduced to \code{year = keep_ratio_year} via
#' \code{.keep_ratio_year()}. Non-finite rows in \code{est/lo/hi} are dropped
#' before plotting. The rightmost column in each family row gets a secondary y-axis
#' title with the family name for visual grouping.
#'
#' @examples
#' \dontrun{
#' # Basic usage (order scenarios 1..7 on x):
#' p <- plot_param_grid_3rows_by_family888(
#'   param_df,
#'   scenario_levels = paste0("S", 1:7),
#'   keep_ratio_year = 2022
#' )
#' print(p)
#'
#' # Hairline CIs with error bars:
#' p2 <- plot_param_grid_3rows_by_family888(
#'   param_df,
#'   scenario_levels     = paste0("S", 1:7),
#'   use_errorbar        = TRUE,
#'   ci_linewidth        = 0.35,
#'   errorbar_cap_width  = 0
#' )
#' }
#'
#' @seealso \pkg{ggplot2}, \pkg{patchwork}
#' @import ggplot2
#' @importFrom patchwork wrap_plots plot_layout
#' @export
plot_param_grid_3rows_by_family888 <- function(
    param_df,
    scenario_levels = NULL,
    families_row   = c("Pella","Fox","Schaefer"),
    params_cols    = c("K","r","MSY","B/Bmsy","F/Fmsy"),
    keep_ratio_year = 2022,
    # aesthetics & export
    width = 16, height = 10, dpi = 400,
    output_dir = NULL, filename = "param_grid_3rows.png",
    # CI and points controls
    ci_linewidth = 0.45,
    ci_alpha     = 0.95,
    point_size   = 2.4,
    dodge_width  = 0.45,
    use_errorbar = FALSE,
    errorbar_cap_width = 0
) {
  if (is.null(param_df) || !is.data.frame(param_df) || !nrow(param_df)) {
    stop("`param_df` must be a non-empty data frame.", call. = FALSE)
  }
  need <- c("scenario","family","type","param","est")
  if (!all(need %in% names(param_df))) {
    stop("`param_df` must include: scenario, family, type, param, est (and lo, hi for CI; year for ratios).", call. = FALSE)
  }
  if (!requireNamespace("ggplot2",  quietly = TRUE)) stop("Package 'ggplot2' is required.", call. = FALSE)
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("Package 'patchwork' is required.", call. = FALSE)

  gg <- ggplot2

  # Ultra-clean theme aligned with your compact theme (fallbacks included)
  .clean_panel_theme <- {
    if (exists("theme_minimal_compact_this", mode = "function")) {
      theme_minimal_compact_this(base_size = 10)
    } else if (exists("theme_minimal_compact", mode = "function")) {
      theme_minimal_compact(base_size = 10)
    } else {
      ggplot2::theme_minimal(base_size = 10)
    }
  } +
    ggplot2::theme(
      plot.background   = ggplot2::element_rect(fill = "white", color = NA),
      panel.background  = ggplot2::element_rect(fill = "white", color = NA),
      panel.grid        = ggplot2::element_blank(),
      panel.border      = ggplot2::element_rect(fill = NA, colour = "grey45", linewidth = 1.6),
      strip.background  = ggplot2::element_rect(fill = "white", color = "grey80", linewidth = 0.5),
      strip.text        = ggplot2::element_text(face = "bold"),
      axis.title.y      = ggplot2::element_text(face = "bold"),
      axis.title.x      = ggplot2::element_blank(),
      axis.text.x       = ggplot2::element_text(color = "black"),
      axis.ticks.x      = ggplot2::element_line(color = "grey45"),
      legend.position   = "bottom",
      legend.title      = ggplot2::element_blank(),
      plot.margin       = ggplot2::margin(6, 6, 6, 6)
    )

  # Column orders and factorization
  df <- param_df
  df$scenario <- as.character(df$scenario)
  df$type     <- as.character(df$type)
  df$family   <- as.character(df$family)
  df$param    <- as.character(df$param)
  if ("year" %in% names(df)) df$year <- suppressWarnings(as.integer(df$year))
  for (nm in c("est","lo","hi")) if (nm %in% names(df)) df[[nm]] <- suppressWarnings(as.numeric(df[[nm]]))

  # Drop non-finite rows
  finite_ok <- is.finite(df$est)
  if ("lo" %in% names(df) && "hi" %in% names(df)) {
    finite_ok <- finite_ok & is.finite(df$lo) & is.finite(df$hi)
  }
  df <- df[finite_ok, , drop = FALSE]
  if (!nrow(df)) stop("No finite rows remain in `param_df` after cleaning.", call. = FALSE)

  # Establish parameter & family orders (drop missing, append extras)
  def_params <- params_cols
  found_params <- intersect(def_params, unique(df$param))
  extras_params <- setdiff(unique(df$param), def_params)
  params_use <- c(found_params, sort(extras_params))

  def_fams <- families_row
  found_fams <- intersect(def_fams, unique(df$family))
  extras_fams <- setdiff(unique(df$family), def_fams)
  families_use <- c(found_fams, sort(extras_fams))
  if (!length(families_use)) stop("No valid families found in `param_df`.", call. = FALSE)

  # Scenario order
  if (is.null(scenario_levels)) {
    sc_all <- unique(df$scenario)
    ord <- order(suppressWarnings(as.integer(sub("^\\D*(\\d+).*", "\\1", sc_all))))
    scenario_levels <- sc_all[ord]
  } else {
    scenario_levels <- unique(c(scenario_levels, setdiff(unique(df$scenario), scenario_levels)))
  }
  if (!length(scenario_levels)) stop("No scenario levels could be determined.", call. = FALSE)

  # Factorize with orders
  df$param    <- factor(df$param,  levels = params_use)
  df$type     <- factor(df$type,   levels = c("SDM","GLM"))
  df$family   <- factor(df$family, levels = families_use)
  df$scenario <- factor(df$scenario, levels = scenario_levels)
  df$x_id     <- as.integer(df$scenario)

  # Keep only the requested year for ratio panels
  df <- .keep_ratio_year(df, keep_year = keep_ratio_year)
  if (!nrow(df)) stop("No data remains after filtering ratio panels to requested year.", call. = FALSE)

  # y-labels (bold math for ratios)
  ylab_for <- function(param) {
    if (param == "B/Bmsy") return(expression(bold(B/B[MSY])))
    if (param == "F/Fmsy") return(expression(bold(F/F[MSY])))
    if (param == "MSY")    return("MSY")
    if (param == "K")      return("K")
    if (param == "r")      return("r")
    as.character(param)
  }
  ratio_limits <- c(0, NA)

  # Colors and dodge
  type_cols  <- c(SDM = "#1B9E77", GLM = "#D95F02")
  pos_dodge  <- ggplot2::position_dodge(width = dodge_width)

  # CI layer factory (linerange or errorbar)
  .ci_layer <- function() {
    if (isTRUE(use_errorbar)) {
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$lo, ymax = .data$hi),
        position = pos_dodge,
        width = errorbar_cap_width,
        linewidth = ci_linewidth,
        alpha = ci_alpha,
        na.rm = TRUE
      )
    } else {
      ggplot2::geom_linerange(
        ggplot2::aes(ymin = .data$lo, ymax = .data$hi),
        position = pos_dodge,
        linewidth = ci_linewidth,
        alpha = ci_alpha,
        na.rm = TRUE
      )
    }
  }

  # Build a single panel (row=family, col=param)
  make_panel <- function(dd, fam_name, one_param, show_col_title = FALSE, add_right_label = FALSE, n_x = 7) {
    if (is.null(dd) || !nrow(dd)) {
      return(ggplot2::ggplot() + ggplot2::theme_void() +
               ggplot2::annotate("text", x = 0.5, y = 0.5,
                                 label = paste(fam_name, one_param, "— no data"),
                                 size = 4))
    }
    if (!("lo" %in% names(dd))) dd$lo <- NA_real_
    if (!("hi" %in% names(dd))) dd$hi <- NA_real_

    p <- ggplot2::ggplot(dd, ggplot2::aes(x = .data$x_id, y = .data$est, color = .data$type)) +
      .ci_layer() +
      ggplot2::geom_point(position = pos_dodge, size = point_size, stroke = 0.36, na.rm = TRUE) +
      ggplot2::scale_color_manual(values = type_cols, breaks = c("SDM","GLM")) +
      ggplot2::scale_x_continuous(breaks = seq_len(n_x), labels = seq_len(n_x)) +
      ggplot2::labs(y = ylab_for(one_param)) +
      .clean_panel_theme

    if (one_param %in% c("B/Bmsy","F/Fmsy")) {
      p <- p + ggplot2::coord_cartesian(ylim = ratio_limits)
    }

    if (isTRUE(add_right_label)) {
      p <- p +
        ggplot2::scale_y_continuous(
          sec.axis = ggplot2::sec_axis(~ ., name = fam_name, breaks = NULL)
        ) +
        ggplot2::theme(axis.text.y.right  = ggplot2::element_blank(),
                       axis.ticks.y.right = ggplot2::element_blank())
    }

    if (isTRUE(show_col_title)) {
      p <- p + ggplot2::labs(title = one_param) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
    }
    p
  }

  # Build rows in families_row order
  plots_list <- list()
  n_x <- length(levels(df$scenario))
  last_col_idx <- length(params_use)

  for (i in seq_along(families_use)) {
    fam <- families_use[i]
    subfam <- df[df$family == fam & !is.na(df$param) & !is.na(df$type), , drop = FALSE]

    row_plots <- vector("list", length(params_use))
    for (j in seq_along(params_use)) {
      prm <- as.character(params_use[j])
      dd  <- subfam[subfam$param == prm, , drop = FALSE]
      row_plots[[j]] <- make_panel(
        dd, fam, prm,
        show_col_title  = (i == 1),
        add_right_label = (j == last_col_idx),
        n_x = n_x
      )
    }

    row_grid <- do.call(patchwork::wrap_plots, c(row_plots, list(ncol = length(params_use), nrow = 1)))
    plots_list[[i]] <- row_grid
  }

  # Stack rows
  grid <- plots_list[[1]]
  if (length(plots_list) >= 2L) {
    for (k in 2:length(plots_list)) grid <- grid / plots_list[[k]]
  }

  # One legend at bottom; white overall
  grid <- grid +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom",
                   plot.background = ggplot2::element_rect(fill = "white", color = NA))

  # Save if requested
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    outpath <- file.path(output_dir, filename)
    ggplot2::ggsave(outpath, grid, width = width, height = height, dpi = dpi)
    message("Saved: ", outpath)
  }

  grid
}
