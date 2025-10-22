#' Keep a Single Year for Ratio Panels (Internal, Robust)
#'
#' Filters a tidy parameter table to keep **only one target year** for the
#' ratio parameters \code{"B/Bmsy"} and \code{"F/Fmsy"} and returns the
#' reduced data frame. Non-ratio parameters are left unchanged.
#'
#' For each group \code{(scenario, type, family, param)}, if \code{param} is a
#' ratio, rows with \code{year == keep_year} are retained. If multiple rows
#' match (unexpected), the first is kept. If no row matches, the entire group
#' is dropped (strict behavior).
#'
#' Robustness:
#' \itemize{
#'   \item If \code{df} is empty or \code{NULL}, it is returned as-is.
#'   \item If any required columns are missing, a clear error is thrown.
#' }
#'
#' @param df A data frame containing the columns:
#'   \code{scenario}, \code{type}, \code{family}, \code{param}, \code{year}.
#' @param keep_year Integer scalar; the target year to keep for ratio
#'   parameters. Default \code{2022}.
#'
#' @return A data frame filtered such that ratio parameters contain only the
#'   specified \code{keep_year}. Row names are reset.
#'
#' @examples
#' \dontrun{
#' df2 <- .keep_ratio_year(param_df, keep_year = 2022)
#' }
#'
#' @keywords internal
#' @noRd
.keep_ratio_year <- function(df, keep_year = 2022) {
  # Basic validation
  if (is.null(df) || !is.data.frame(df) || !nrow(df)) return(df)

  req_cols <- c("scenario", "type", "family", "param", "year")
  miss <- setdiff(req_cols, names(df))
  if (length(miss)) {
    stop("`.keep_ratio_year()` missing required columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  out <- vector("list", 0L)
  keys <- paste(df$scenario, df$type, df$family, df$param, sep = "||")
  ukeys <- unique(keys)

  for (k in ukeys) {
    idx <- which(keys == k)
    block <- df[idx, , drop = FALSE]

    # Only enforce on ratio parameters
    if (isTRUE(block$param[1] %in% c("B/Bmsy","F/Fmsy"))) {
      keep_idx <- which(!is.na(block$year) & as.integer(block$year) == as.integer(keep_year))
      if (length(keep_idx) >= 1L) {
        block <- block[keep_idx, , drop = FALSE]
        if (nrow(block) > 1L) block <- block[1L, , drop = FALSE]
      } else {
        next  # strict: drop this group entirely if requested year is absent
      }
    }
    out[[length(out) + 1L]] <- block
  }

  if (!length(out)) return(df[0, , drop = FALSE])
  res <- tryCatch(do.call(rbind, out), error = function(e) df[0, , drop = FALSE])
  rownames(res) <- NULL
  res
}


#' Grid Plot of SDM vs GLM CIs by Family (3×5 Panels, Robust)
#'
#' Produces a clean **3 rows × 5 columns** grid of panels showing **SDM** and **GLM**
#' point estimates with confidence intervals across scenarios, split by family:
#' rows are **Pella**, **Fox**, **Schaefer** (by default), and columns are
#' **K**, **r**, **MSY**, **B/Bmsy**, **F/Fmsy**.
#'
#' Features:
#' \itemize{
#'   \item X-axis shows scenario indices \code{1..N} (no x-title).
#'   \item Only the **rightmost panel in each row** carries a **secondary y-axis title**
#'         with the **family name** (no numbers or ticks on that secondary axis).
#'   \item Panels have a clean white background with subtle borders.
#'   \item For ratio panels (\code{B/Bmsy}, \code{F/Fmsy}), only **year = 2022** is kept.
#'   \item Defensive checks for inputs, columns, and missing data; non-finite values are dropped.
#' }
#'
#' @section Expected input:
#' \code{param_df} must be a tidy data frame with columns:
#' \describe{
#'   \item{scenario}{Scenario ID (e.g., \code{"S1"}, \code{"S2"}, ...).}
#'   \item{type}{Model type: \code{"SDM"} or \code{"GLM"}.}
#'   \item{family}{One of \code{"Pella"}, \code{"Fox"}, \code{"Schaefer"}.}
#'   \item{param}{One of \code{"K"}, \code{"r"}, \code{"MSY"}, \code{"B/Bmsy"}, \code{"F/Fmsy"}.}
#'   \item{year}{Numeric year for ratio parameters (e.g., \code{2022}); may be \code{NA} for others.}
#'   \item{est, lo, hi}{Point estimate and CI bounds (numeric).}
#' }
#'
#' @param param_df Data frame as described above.
#' @param scenario_levels Optional character vector controlling scenario order on the x-axis
#'   (e.g., \code{paste0("S", 1:7)}). If \code{NULL}, order is inferred numerically from \code{scenario}.
#'   Any levels not present in the data are ignored; any scenarios present but not listed will be appended.
#' @param families_row Character vector giving the row order of families.
#'   Default \code{c("Pella","Fox","Schaefer")}. Families not present are dropped; present-but-unlisted
#'   families are appended at the end in alphabetical order.
#' @param params_cols Character vector giving the column order of parameters.
#'   Default \code{c("K","r","MSY","B/Bmsy","F/Fmsy")}. Parameters not present are dropped; present-but-unlisted
#'   parameters are appended at the end in alphabetical order.
#' @param width,height,dpi Device size and resolution used only when saving to file.
#' @param output_dir Optional directory to save a PNG. If \code{NULL}, no file is saved.
#' @param filename File name for the PNG (used when \code{output_dir} is non-\code{NULL}).
#' @param ci_linewidth Numeric CI line thickness (smaller = thinner). Default \code{0.45}.
#' @param ci_alpha Numeric CI opacity in \code{[0,1]}. Default \code{0.95}.
#' @param point_size Numeric point size for estimates. Default \code{2.4}.
#' @param dodge_width Horizontal dodge between SDM and GLM at each scenario. Default \code{0.45}.
#' @param use_errorbar Logical; if \code{TRUE}, use \code{geom_errorbar()} (hairline option)
#'   instead of \code{geom_linerange()} for CIs. Default \code{FALSE}.
#' @param errorbar_cap_width Width of errorbar caps (set \code{0} for hairlines). Default \code{0}.
#'
#' @return A \pkg{patchwork} object representing the combined grid. If \code{output_dir}
#'   is provided, a PNG is also saved.
#'
#' @details
#' Internally, the function filters ratio panels to \code{year = 2022} via an internal helper.
#' Only the rightmost panel in each family row receives a secondary y-axis with the family
#' name and no tick labels on that secondary axis. Non-finite rows in \code{est/lo/hi} are dropped.
#'
#' @examples
#' \dontrun{
#' # Assuming 'param_df' was built by your extractor:
#' p <- plot_param_grid_3rows_by_family888(
#'   param_df,
#'   scenario_levels = paste0("S", 1:7),
#'   ci_linewidth = 0.35,
#'   dodge_width  = 0.40,
#'   point_size   = 2.2
#' )
#' print(p)
#'
#' # Hairline CIs
#' p2 <- plot_param_grid_3rows_by_family888(
#'   param_df,
#'   scenario_levels = paste0("S", 1:7),
#'   use_errorbar   = TRUE,
#'   ci_linewidth   = 0.35,
#'   errorbar_cap_width = 0
#' )
#' }
#'
#' @seealso \pkg{ggplot2}, \pkg{patchwork}
#' @import ggplot2
#' @importFrom patchwork wrap_plots plot_layout
#' @export
plot_param_grid_3rows_by_family888 <- function(
    param_df,
    scenario_levels = NULL,                    # e.g., paste0("S", 1:7)
    families_row   = c("Pella","Fox","Schaefer"),
    params_cols    = c("K","r","MSY","B/Bmsy","F/Fmsy"),
    # aesthetics & export
    width = 16, height = 10, dpi = 400,
    output_dir = NULL, filename = "param_grid_3rows.png",
    # CI and points controls (for thinner CIs etc.)
    ci_linewidth = 0.45,                       # thickness of CI lines
    ci_alpha     = 0.95,                       # CI opacity
    point_size   = 2.4,                        # estimate point size
    dodge_width  = 0.45,                       # SDM/GLM horizontal spacing
    use_errorbar = FALSE,                      # TRUE -> geom_errorbar hairlines; FALSE -> geom_linerange
    errorbar_cap_width = 0                     # 0 = hairline; try 0.08 for tiny caps
) {
  # Package checks with clear errors
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required. Please install it.", call. = FALSE)
  }

  # Input checks
  if (is.null(param_df) || !is.data.frame(param_df) || !nrow(param_df)) {
    stop("`param_df` must be a non-empty data frame.", call. = FALSE)
  }

  # Required columns
  req <- c("scenario","type","family","param","year","est","lo","hi")
  miss <- setdiff(req, names(param_df))
  if (length(miss)) {
    stop("`param_df` is missing required columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  # Coerce basic types safely
  df <- param_df
  df$scenario <- as.character(df$scenario)
  df$type     <- as.character(df$type)
  df$family   <- as.character(df$family)
  df$param    <- as.character(df$param)
  df$year     <- suppressWarnings(as.integer(df$year))
  df$est      <- suppressWarnings(as.numeric(df$est))
  df$lo       <- suppressWarnings(as.numeric(df$lo))
  df$hi       <- suppressWarnings(as.numeric(df$hi))

  # Drop non-finite estimate rows (avoid drawing issues)
  finite_ok <- is.finite(df$est) & is.finite(df$lo) & is.finite(df$hi)
  if (!all(finite_ok)) {
    df <- df[finite_ok, , drop = FALSE]
  }
  if (!nrow(df)) stop("No finite rows remain in `param_df` after cleaning.", call. = FALSE)

  # Establish parameter & family orders robustly (drop missing, append extras)
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
    # Infer numeric order by stripping leading 'S'
    sc_all <- unique(df$scenario)
    ord <- order(suppressWarnings(as.integer(sub("^\\D*(\\d+).*", "\\1", sc_all))))
    scenario_levels <- sc_all[ord]
  } else {
    # Use provided order but ensure all present scenarios included
    scenario_levels <- unique(c(scenario_levels, setdiff(unique(df$scenario), scenario_levels)))
  }
  if (!length(scenario_levels)) stop("No scenario levels could be determined.", call. = FALSE)

  # Factorize with orders
  df$param  <- factor(df$param,  levels = params_use)
  df$type   <- factor(df$type,   levels = c("SDM","GLM"))
  df$family <- factor(df$family, levels = families_use)
  df$scenario <- factor(df$scenario, levels = scenario_levels)
  df$x_id <- as.integer(df$scenario)

  # Strictly keep 2022 for ratio panels
  df <- .keep_ratio_year(df, keep_year = 2022)

  # After filtering, it may become empty
  if (!nrow(df)) stop("No data remains after filtering ratio panels to year 2022.", call. = FALSE)

  # y-axis labels builder
  ylab_for <- function(param) {
    if (identical(param, "B/Bmsy")) return(ggplot2::label_parsed("bold(B/B[MSY])"))
    if (identical(param, "F/Fmsy")) return(ggplot2::label_parsed("bold(F/F[MSY])"))
    param
  }
  # For expression labels in labs(y=...), use expressions not labelers:
  ylab_expr <- function(param) {
    if (identical(param, "B/Bmsy")) return(expression(bold(B/B[MSY])))
    if (identical(param, "F/Fmsy")) return(expression(bold(F/F[MSY])))
    if (identical(param, "MSY"))    return("MSY")
    if (identical(param, "K"))      return("K")
    if (identical(param, "r"))      return("r")
    return(param)
  }
  ratio_limits <- c(0, NA)

  # Colors and positions
  type_cols  <- c(SDM = "#1B9E77", GLM = "#D95F02")
  pos_dodge  <- ggplot2::position_dodge(width = dodge_width)

  # Theme (clean white + subtle frame; keep x ticks)
  .exists_fun <- function(name) exists(name, mode = "function", inherits = TRUE)
  .get_fun    <- function(name) get(name, inherits = TRUE)
  .clean_panel_theme <- local({
    base_th <- if (.exists_fun("theme_minimal_compact")) {
      .get_fun("theme_minimal_compact")(base_size = 10)
    } else {
      ggplot2::theme_minimal(base_size = 10)
    }
    base_th + ggplot2::theme(
      plot.background   = ggplot2::element_rect(fill = "white", color = NA),
      panel.background  = ggplot2::element_rect(fill = "white", color = NA),
      panel.grid        = ggplot2::element_blank(),
      panel.border      = ggplot2::element_rect(fill = NA, colour = "grey45", linewidth = 1),
      strip.background  = ggplot2::element_rect(fill = "white", color = "grey80", linewidth = 0.5),
      strip.text        = ggplot2::element_text(face = "bold"),
      axis.title.y      = ggplot2::element_text(face = "bold"),
      axis.title.x      = ggplot2::element_blank(),   # no "Model index" title
      axis.text.x       = ggplot2::element_text(color = "black"),  # keep 1..N
      axis.ticks.x      = ggplot2::element_line(color = "grey45"),
      legend.position   = "bottom",
      legend.title      = ggplot2::element_blank(),
      plot.margin       = ggplot2::margin(6, 6, 6, 6)
    )
  })

  # CI layer factory (linerange or errorbar)
  .ci_layer <- function() {
    if (isTRUE(use_errorbar)) {
      return(ggplot2::geom_errorbar(
        ggplot2::aes(ymin = lo, ymax = hi),
        position = pos_dodge,
        width = errorbar_cap_width,          # 0 = hairlines
        linewidth = ci_linewidth,
        alpha = ci_alpha
      ))
    } else {
      return(ggplot2::geom_linerange(
        ggplot2::aes(ymin = lo, ymax = hi),
        position = pos_dodge,
        linewidth = ci_linewidth,
        alpha = ci_alpha
      ))
    }
  }

  # Build one panel (single family row + single parameter column)
  # add_right_label: TRUE only for the rightmost column in each row
  make_panel <- function(dd, fam_name, one_param, show_col_title = FALSE, add_right_label = FALSE, n_x = 7) {
    if (!nrow(dd)) {
      return(ggplot2::ggplot() + ggplot2::theme_void() +
               ggplot2::annotate("text", x = 0.5, y = 0.5,
                                 label = paste(fam_name, one_param, "— no data"),
                                 size = 4))
    }

    p <- ggplot2::ggplot(dd, ggplot2::aes(x = x_id, y = est, color = type)) +
      .ci_layer() +
      ggplot2::geom_point(position = pos_dodge, size = point_size, stroke = 0.36) +
      ggplot2::scale_color_manual(values = type_cols) +
      ggplot2::scale_x_continuous(breaks = seq_len(n_x), labels = seq_len(n_x)) +
      ggplot2::labs(y = ylab_expr(one_param)) +
      .clean_panel_theme

    if (identical(one_param, "B/Bmsy") || identical(one_param, "F/Fmsy")) {
      p <- p + ggplot2::coord_cartesian(ylim = ratio_limits)
    }

    # Secondary y-axis title ONLY on rightmost column; NO numbers/ticks on secondary axis
    if (isTRUE(add_right_label)) {
      p <- p +
        ggplot2::scale_y_continuous(
          sec.axis = ggplot2::sec_axis(~ ., name = fam_name, breaks = NULL)
        ) +
        ggplot2::theme(
          axis.text.y.right  = ggplot2::element_blank(),
          axis.ticks.y.right = ggplot2::element_blank()
        )
    }

    # Column headers only on the top row
    if (isTRUE(show_col_title)) {
      p <- p + ggplot2::labs(title = one_param) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
    }
    p
  }

  # Build row by row
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

  # Single legend at bottom; white overall
  grid <- grid +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom",
                   plot.background = ggplot2::element_rect(fill = "white", color = NA))

  # Save if requested
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }
    outpath <- file.path(output_dir, filename)
    ggplot2::ggsave(outpath, grid, width = width, height = height, dpi = dpi)
    message("Saved: ", outpath)
  }

  grid
}
