#' JABBA Color Set for Peel Plots
#'
#' Predefined color palette inspired by the JABBA package, suitable for use in
#' retrospective "peel" plots or any situation requiring distinct colors.
#'
#' @format A character vector of 8 hexadecimal color codes.
#' @details
#' The colors are: red, blue, green, purple, orange, brown, pink, and grey.
#'
#' @examples
#' barplot(rep(1, 8), col = jabba_colors)
#'
#' @export
jabba_colors <- c(
  "#e7298a",  # Red
  "#1f78b4",  # Blue
  "#1b9e77",  # Green
  "#984EA3",  # Purple
  "#d95f02",  # Orange
  "#A65628",  # Brown
  "#F781BF",  # Pink
  "#8dd3c7"   # Grey
)

#' Create a Peel Plot Subtitle with Colored Symbols
#'
#' Generates an HTML-formatted subtitle for retrospective (peel) plots,
#' displaying colored dots and bold peel names for each peel.
#'
#' @param peel_names Character vector of peel labels (e.g., c("All", "-1", "-2")).
#' @param peel_colors Character vector of color codes, one for each peel.
#'
#' @return A single character string containing HTML for use as a plot subtitle.
#'
#' @details
#' Designed for use with \code{ggtext::element_markdown()} in ggplot2 plot
#' subtitles, so that each peel's name is shown in its associated color.
#'
#' @examples
#' subtitle <- make_peel_subtitle(c("All", "-1", "-2"), c("black", "red", "blue"))
#' @export
make_peel_subtitle <- function(peel_names, peel_colors) {
  paste0(
    "Peels: ",
    paste(
      mapply(function(nm, col) {
        sprintf("<span style='color:%s'>&#9679;</span> <b>%s</b>", col, nm)
      }, peel_names, peel_colors, SIMPLIFY = TRUE),
      collapse = "   "
    )
  )
}

#' Retrospective Analysis Plot for Stock Assessment Models
#'
#' Plots results of retrospective analysis ("peel" plots) for a fitted model,
#' displaying time series of biomass, fishing mortality, and their ratios to
#' reference points, with optional display of Mohn's rho values.
#'
#' @param rep A fitted model object containing a \code{$retro} element
#'   (as produced by \code{retro()}).
#' @param add_mohn Logical; whether to add Mohn's rho statistics to B/Bmsy and
#'   F/Fmsy panels. Default: \code{TRUE}.
#' @param CI Confidence interval width (default 0.95).
#' @param peel_colors Optional vector of color codes for each peel. If
#'   \code{NULL}, uses the JABBA color set or MetBrewer/RColorBrewer palettes.
#' @param palette String; color palette to use ("JABBA" for default, or any
#'   valid \code{MetBrewer} or \code{RColorBrewer} palette name).
#'
#' @details
#' The function extracts time series from all model peels, aligns them for
#' comparison, and plots them using consistent coloring and ordering. It annotates
#' panels with Mohn's rho statistics if requested. Subtitles show colored peel
#' names using \code{ggtext::element_markdown()} for clarity.
#'
#' @return A \code{patchwork} composite ggplot2 object, suitable for display or
#'   saving with \code{ggsave()}.
#'
#' @seealso \code{retro()}, \code{mohns_rho()}
#'
#' @examples
#' \dontrun{
#' g <- retro.rho.ndiaye2(model_with_retro)
#' print(g)
#' }
#' @export
retro.rho.ndiaye2 <- function(rep, add_mohn = TRUE, CI = 0.95,
                              peel_colors = NULL, palette = "Hokusai1") {

  runs <- rep$retro
  if (is.null(runs)) stop("rep$retro is NULL: you must provide a fitted retrospective model list")
  names(runs)[1] <- "All"
  n_peels    <- length(runs)
  peel_names <- names(runs)

  ## ----- Color logic (kept as-is; robust fallbacks) --------------------------
  if (!is.null(peel_colors)) {
    peel_colors <- rep(peel_colors, length.out = n_peels)
  } else if (identical(toupper(palette), "JABBA")) {
    peel_colors <- rep(jabba_colors, length.out = n_peels)
  } else if (requireNamespace("MetBrewer", quietly = TRUE) &&
             palette %in% MetBrewer::list_met_palettes()) {
    peel_colors <- MetBrewer::met.brewer(palette, n_peels)
  } else if (requireNamespace("RColorBrewer", quietly = TRUE) &&
             palette %in% rownames(RColorBrewer::brewer.pal.info)) {
    peel_colors <- RColorBrewer::brewer.pal(max(n_peels, 3L), palette)[seq_len(n_peels)]
  } else {
    warning("Unknown palette, using ggplot2 default.")
    peel_colors <- scales::hue_pal()(n_peels)
  }
  names(peel_colors) <- peel_names
  ## --------------------------------------------------------------------------

  ## ----- Mohn's rho labels (plotmath-friendly) -------------------------------
  lab_BB <- lab_FF <- NULL
  if (add_mohn) {
    rho <- round(mohns_rho(rep, what = c("FFmsy", "BBmsy")), 3)
    # Use italic(MSY) so parse=TRUE works (plotmath does not know \emph{}):
    lab_BB <- paste0("Mohn*\"'s\"~rho[B/B[italic(MSY)]]==", rho["BBmsy"])
    lab_FF <- paste0("Mohn*\"'s\"~rho[F/F[italic(MSY)]]==", rho["FFmsy"])
  }
  ## --------------------------------------------------------------------------

  ## ----- Extract time series per peel ---------------------------------------
  extract_df <- function(fit, peel, par_name) {
    idx  <- fit$inp$indest
    pars <- get.par(par_name, fit, exp = TRUE, CI = CI)[idx, 1:3]
    tibble::tibble(
      peel     = peel,
      time     = fit$inp$time[idx],
      estimate = pars[, 2],
      lower    = pars[, 1],
      upper    = pars[, 3]
    )
  }
  df_B    <- purrr::imap_dfr(runs, ~ extract_df(.x, .y, "logB"))
  df_F    <- purrr::imap_dfr(runs, ~ extract_df(.x, .y, "logFnotS"))
  df_Bmsy <- purrr::imap_dfr(runs, ~ extract_df(.x, .y, "logBBmsy"))
  df_Fmsy <- purrr::imap_dfr(runs, ~ extract_df(.x, .y, "logFFmsynotS"))

  ## ----- Reorder peels for plotting (All first, then -1, -2, ...) -----------
  other_peels      <- setdiff(names(runs), "All")
  peel_nums        <- suppressWarnings(as.integer(other_peels))
  valid_peels      <- other_peels[!is.na(peel_nums)]
  ordered_peels    <- valid_peels[order(as.integer(valid_peels), decreasing = TRUE)]
  plot_peel_levels <- c("All", ordered_peels)

  df_B$peel    <- factor(df_B$peel,    levels = plot_peel_levels)
  df_F$peel    <- factor(df_F$peel,    levels = plot_peel_levels)
  df_Bmsy$peel <- factor(df_Bmsy$peel, levels = plot_peel_levels)
  df_Fmsy$peel <- factor(df_Fmsy$peel, levels = plot_peel_levels)

  ## ----- Panel builder (labels unchanged in spirit; italic(MSY) for plotmath)
  make_panel <- function(df, ylab, rho_text = NULL) {
    ci_gray <- "#B0B0B0"
    ggplot2::ggplot(df, ggplot2::aes(time, estimate, group = peel)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                           fill = ci_gray, alpha = 0.25) +
      ggplot2::geom_line(ggplot2::aes(color = peel), linewidth = 0.6) +
      ggplot2::labs(x = NULL, y = ylab) +
      theme_minimal_compact() +
      ggplot2::scale_color_manual(values = peel_colors, guide = "none") +
      {
        if (!is.null(rho_text)) {
          ggplot2::annotate(
            "text", x = Inf, y = Inf, label = rho_text,
            parse = TRUE, hjust = 1.1, vjust = 1.5, size = 2
          )
        } else NULL
      }
  }

  p1 <- make_panel(df_B,    expression(bold(B[t])))
  p2 <- make_panel(df_F,    expression(bold(F[t])))
  p3 <- make_panel(df_Bmsy, expression(bold(B[t]/B[italic(MSY)])), lab_BB)
  p4 <- make_panel(df_Fmsy, expression(bold(F[t]/F[italic(MSY)])), lab_FF)

  ## ----- Colored-dot subtitle ------------------------------------------------
  subtitle_peels <- make_peel_subtitle(peel_names, peel_colors)

  ## ----- Assemble (patchwork) -----------------------------------------------
  (p1 + p2) / (p3 + p4) +
    patchwork::plot_layout(guides = "collect", heights = c(1, 1)) &
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = ggtext::element_markdown(hjust = 0.5, size = 11, face = "bold")
    ) &
    patchwork::plot_annotation(
      title = "Retrospective Analysis",
      subtitle = subtitle_peels
    )
}
