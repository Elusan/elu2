#' Plot SPiCT parameter panels (SDM vs GLM) by scenario for one model
#'
#' @description
#' Builds a facetted panel plot of key SPiCT parameters with 95\% CIs, contrasting
#' \emph{SDM} vs \emph{GLM} estimates per scenario for a single model family
#' (\code{"Pella"}, \code{"Schaefer"}, or \code{"Fox"}). Each facet is a parameter,
#' x-axis is scenario, points are estimates, and vertical line ranges are CIs.
#'
#' @param model_list A named list of fitted SPiCT objects (class \code{"spictcls"})
#'   with names of the form \code{"<scenario><model_code>.<type>"} where
#'   \code{<scenario>} is e.g. \code{"S1"}, \code{<model_code>} is \code{"P"}, \code{"S"}, or \code{"F"}
#'   (Pella/Schaefer/Fox), and \code{<type>} is \code{"SDM"} or \code{"GLM"}.
#'   For example: \code{c("S1F.SDM", "S1S.SDM", "S1P.SDM", "S1F.GLM", ...)}.
#' @param model_name Character scalar. Which model family to plot columns for:
#'   one of \code{c("Pella","Schaefer","Fox")} (matched via \code{match.arg()}).
#'   This selects the three columns consumed from the summary data frame:
#'   \code{<model_name>_estimate}, \code{<model_name>_low}, \code{<model_name>_up}.
#' @param scenario_levels Character vector of scenario labels (factor order on x-axis).
#'   Default: \code{c("S1","S2","S3","S4")}.
#' @param year Optional integer. Year used in facet labels (e.g. \eqn{B_y, F_y}).
#'   If \code{NULL}, it is inferred from the first available SPiCT report via
#'   \code{rep$inp$timerangeObs[2]} (terminal observation year).
#' @param font_size_base Numeric. Base font size for the built-in compact theme.
#'   Default: \code{11}.
#' @param save_path Optional file path. If provided and the extension is one of
#'   \code{c("pdf","png","jpeg","jpg","tiff")}, the plot is saved via
#'   \code{ggplot2::ggsave()} (20×12 inches @ 300 dpi).
#' @param dodge_width Numeric. Horizontal dodge between SDM and GLM points/intervals.
#'   Default: \code{0.36}.
#' @param x_expand_mult Numeric length-2. Expansion multiplicative padding for the
#'   discrete x-scale (left, right). Default: \code{c(0.02, 0.02)}.
#'
#' @return A \code{ggplot} object (returned invisibly). If \code{save_path} is provided
#'   with a supported extension, the plot is also written to disk as a side effect.
#'
#' @details
#' This function requires an auxiliary summariser
#' \code{make_spict_summary_df1_nopipe(models, model_names)} to be available in scope.
#' The summariser must return a data frame with:
#' \itemize{
#' \item Column \code{Parameters}: parameter labels (LaTeX/TeX strings), e.g.,
#'   \eqn{B_{y}/B_{\mathrm{MSY}}}, \eqn{B_{\mathrm{MSY}}}, \eqn{B_y}, \eqn{K}, \eqn{n},
#'   \eqn{F_{y}/F_{\mathrm{MSY}}}, \eqn{F_y}, \eqn{F_{\mathrm{MSY}}}, \eqn{\mathrm{MSY}}, \eqn{r}.
#' \item Three numeric columns per model family: \code{<model_name>_estimate},
#'   \code{<model_name>_low}, \code{<model_name>_up}.
#' }
#'
#' @section Built-in compact theme:
#' The plot uses a private, built-in theme (crisp border, bold text, inside top-right
#' legend, thicker axis ticks, blank grid, inside facet strips). It mirrors your
#' previous \code{theme_minimal_compact2_for_PANELS()} but is not exported.
#'
#' @examples
#' \dontrun{
#' model_list <- list(
#'   "S1F.SDM" = fit_fox_sdm, "S1S.SDM" = fit_schaefer_sdm, "S1P.SDM" = fit_pella_sdm,
#'   "S1F.GLM" = fit_fox_glm, "S1S.GLM" = fit_schaefer_glm, "S1P.GLM" = fit_pella_glm
#' )
#' p <- plot_param_panels_by_model(
#'   model_list,
#'   model_name = "Pella",
#'   scenario_levels = c("S1","S2","S3","S4")
#' )
#' }
#'
#' @seealso \code{make_spict_summary_df1_nopipe}
#' @export
#' @import ggplot2
#' @importFrom grid unit
#' @importFrom tools file_ext
plot_param_panels_by_model <- function(model_list,
                                       model_name = c("Pella","Schaefer","Fox"),
                                       scenario_levels = c("S1","S2","S3","S4"),
                                       year = NULL,
                                       font_size_base = 11,
                                       save_path = NULL,
                                       dodge_width = 0.36,
                                       x_expand_mult = c(0.02, 0.02)) {

  # --- Private compact theme (not exported) ---------------------------------
  .theme_minimal_compact2_for_PANELS <- function(base_size = 10, base_family = "") {
    ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = ggplot2::element_text(face = "bold", size = 12),
        axis.text  = ggplot2::element_text(size = 12, face = "bold"),
        legend.position = c(0.88, 0.98),
        legend.justification = c("left", "top"),
        legend.background = ggplot2::element_blank(),
        legend.box.background = ggplot2::element_blank(),
        legend.title = ggplot2::element_blank(),
        legend.text  = ggplot2::element_text(size = 12),
        legend.key.size = grid::unit(1, "lines"),
        panel.grid = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(fill = NA, colour = "grey35", linewidth = 2),
        axis.ticks = ggplot2::element_line(linewidth = 0.7, color = "grey35"),
        axis.ticks.length = grid::unit(3, "pt"),
        strip.background = ggplot2::element_rect(fill = "grey35", color = "grey35", linewidth = 0.5),
        strip.text = ggplot2::element_text(face = "bold", size = ggplot2::rel(1)),
        text = ggplot2::element_text(face = "bold", size = 11),
        plot.margin = ggplot2::margin(3, 3, 3, 3)
      )
  }

  # --- Guards ---------------------------------------------------------------
  model_name <- match.arg(model_name)
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
  if (!exists("make_spict_summary_df1_nopipe"))
    stop("make_spict_summary_df1_nopipe() must be available.")

  # --- Small helpers --------------------------------------------------------
  .first_spict_rep <- function(x) {
    for (i in seq_along(x)) if (inherits(x[[i]], "spictcls")) return(x[[i]])
    stop("No 'spictcls' in model_list.")
  }
  .infer_year <- function(x) {
    rep0 <- .first_spict_rep(x)
    yy <- tryCatch(rep0$inp$timerangeObs[2], error = function(e) NA_real_)
    if (!is.finite(yy)) stop("Cannot infer terminal year.")
    as.integer(floor(yy))
  }
  .tex_label_map <- function(y) c(
    BBmsy_year = sprintf("$B_{%d}/B_{\\mathrm{MSY}}$", y),
    Bmsy = "$B_{\\mathrm{MSY}}$", B_year = sprintf("$B_{%d}$", y), K = "$K$", n = "$n$",
    FFmsy_year = sprintf("$F_{%d}/F_{\\mathrm{MSY}}$", y), F_year = sprintf("$F_{%d}$", y),
    Fmsy = "$F_{\\mathrm{MSY}}$", MSY = "$\\mathrm{MSY}$", r = "$r$"
  )
  .plotmath_string_map <- function(y) {
    tex <- .tex_label_map(y)
    stats::setNames(
      c(sprintf("B[%d]/B[MSY]", y), "B[MSY]", sprintf("B[%d]", y), "K", "n",
        sprintf("F[%d]/F[MSY]", y), sprintf("F[%d]", y), "F[MSY]", "MSY", "r"),
      tex
    )
  }
  .get_models <- function(scen, type) list(
    Fox = model_list[[paste0(scen, "F.", type)]],
    Schaefer = model_list[[paste0(scen, "S.", type)]],
    Pella = model_list[[paste0(scen, "P.", type)]]
  )
  .model_column_names <- function(model_name) paste0(model_name, c("_estimate", "_low", "_up"))

  # --- Year and labels ------------------------------------------------------
  if (is.null(year)) year <- .infer_year(model_list)
  tex_map <- .tex_label_map(year)
  target_order_keys <- c("BBmsy_year","Bmsy","B_year","K","n","FFmsy_year","F_year","Fmsy","MSY","r")
  target_labels_tex <- unname(tex_map[target_order_keys])

  est_col <- .model_column_names(model_name)[1]
  lo_col  <- .model_column_names(model_name)[2]
  up_col  <- .model_column_names(model_name)[3]

  # --- Build long data ------------------------------------------------------
  all_rows <- list(); k <- 0L
  for (scen in scenario_levels) {
    mods_SDM <- .get_models(scen, "SDM")
    df_SDM <- make_spict_summary_df1_nopipe(mods_SDM, names(mods_SDM))
    sub_SDM <- df_SDM[!is.na(match(df_SDM$Parameters, target_labels_tex)),
                      c("Parameters", est_col, lo_col, up_col), drop = FALSE]
    if (nrow(sub_SDM)) for (i in seq_len(nrow(sub_SDM))) {
      k <- k + 1L; all_rows[[k]] <- data.frame(
        scenario = scen, dtype = "SDM", parameter_tex = as.character(sub_SDM$Parameters[i]),
        est = as.numeric(sub_SDM[i, est_col]), lwr = as.numeric(sub_SDM[i, lo_col]),
        upr = as.numeric(sub_SDM[i, up_col]), stringsAsFactors = FALSE
      )
    }
    mods_GLM <- .get_models(scen, "GLM")
    df_GLM <- make_spict_summary_df1_nopipe(mods_GLM, names(mods_GLM))
    sub_GLM <- df_GLM[!is.na(match(df_GLM$Parameters, target_labels_tex)),
                      c("Parameters", est_col, lo_col, up_col), drop = FALSE]
    if (nrow(sub_GLM)) for (i in seq_len(nrow(sub_GLM))) {
      k <- k + 1L; all_rows[[k]] <- data.frame(
        scenario = scen, dtype = "GLM", parameter_tex = as.character(sub_GLM$Parameters[i]),
        est = as.numeric(sub_GLM[i, est_col]), lwr = as.numeric(sub_GLM[i, lo_col]),
        upr = as.numeric(sub_GLM[i, up_col]), stringsAsFactors = FALSE
      )
    }
  }
  if (!length(all_rows)) stop("No data found for the requested model / scenarios.")
  dat <- do.call(rbind, all_rows)

  dat$scenario      <- factor(dat$scenario, levels = scenario_levels, ordered = TRUE)
  dat$dtype         <- factor(dat$dtype, levels = c("SDM","GLM"), ordered = TRUE)
  dat$parameter_tex <- factor(dat$parameter_tex, levels = target_labels_tex, ordered = TRUE)

  keep <- is.finite(dat$est) & is.finite(dat$lwr) & is.finite(dat$upr)
  dat  <- dat[keep, , drop = FALSE]
  if (!nrow(dat)) stop("All rows had missing estimate/CI.")

  # --- Plot -----------------------------------------------------------------
  facet_map_strings <- .plotmath_string_map(year)
  pd <- ggplot2::position_dodge(width = dodge_width)

  gg <- ggplot2::ggplot(dat, ggplot2::aes(x = scenario, group = dtype)) +
    ggplot2::geom_linerange(
      ggplot2::aes(ymin = lwr, ymax = upr, linetype = dtype),
      position = pd, linewidth = 0.6
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = est, shape = dtype),
      position = pd, size = 2.5, colour = "red"
    ) +
    ggplot2::facet_wrap(
      ~ parameter_tex,
      scales = "free_y",
      ncol = 2,
      labeller = ggplot2::as_labeller(facet_map_strings, ggplot2::label_parsed)
    ) +
    ggplot2::labs(
      x = NULL, y = NULL,
      title = paste0(model_name, " model — Parameters with 95% CI (SDM vs GLM, year ", year, ")"),
      subtitle = paste("Scenarios:", paste(scenario_levels, collapse = ", "))
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.14))) +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = x_expand_mult)) +
    ggplot2::scale_linetype_manual(values = c(SDM = "solid", GLM = "dashed")) +
    ggplot2::scale_shape_manual(values = c(SDM = 16, GLM = 17)) +
    ggplot2::guides(linetype = ggplot2::guide_legend(title = NULL),
                    shape    = ggplot2::guide_legend(title = NULL)) +
    .theme_minimal_compact2_for_PANELS(base_size = font_size_base) +
    ggplot2::theme(
      strip.placement = "inside",
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.25),
                                         margin = ggplot2::margin(t = 2, r = 3, b = 2, l = 3)),
      plot.subtitle = ggplot2::element_text(hjust = 0.5)
    )

  # --- Save (optional) ------------------------------------------------------
  if (!is.null(save_path)) {
    ext <- tolower(tools::file_ext(save_path))
    if (ext %in% c("pdf","png","jpeg","jpg","tiff")) {
      ggplot2::ggsave(filename = save_path, plot = gg, width = 20, height = 12, units = "in", dpi = 300)
    } else {
      warning("Unknown file extension for save_path; skipping save.")
    }
  }
  invisible(gg)
}
