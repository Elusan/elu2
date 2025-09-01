# R/osar_diag_plots.R

#' Compact diagnostic theme (no frills, bold panels, no legends)
#'
#' A minimal, high-contrast theme tailored for OSAR/diagnostic panels:
#' bold titles/axes, visible panel borders, no legends, and no grid lines.
#'
#' @param base_size Base font size. Defaults to `10`.
#' @param base_family Base font family. Defaults to `""` (use device default).
#'
#' @return A \pkg{ggplot2} theme object.
#'
#' @examples
#' if (interactive()) {
#'   library(ggplot2)
#'   ggplot(mtcars, aes(wt, mpg)) +
#'     geom_point() +
#'     theme_minimal_compact2_diagOSA()
#' }
#'
#' @export
#' @import ggplot2
#' @importFrom grid unit
theme_minimal_compact2_diagOSA <- function(base_size = 10, base_family = "") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 11),
      axis.title = ggplot2::element_text(face = "bold", size = 10),
      axis.text  = ggplot2::element_text(size = 8, face = "bold"),
      legend.position = "none",
      legend.background = ggplot2::element_blank(),
      legend.box.background = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank(),
      legend.text  = ggplot2::element_text(size = 12),
      legend.key.size = grid::unit(0.8, "lines"),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(fill = NA, colour = "grey45", linewidth = 1),
      axis.ticks = ggplot2::element_line(linewidth = 0.6, color = "grey45"),
      axis.ticks.length = grid::unit(2, "pt"),
      strip.background = ggplot2::element_rect(fill = "grey45", color = "grey45", linewidth = 0.5),
      strip.text = ggplot2::element_text(face = ggplot2::rel(1), size = ggplot2::rel(1)),
      text = ggplot2::element_text(face = "bold", size = 10),
      plot.margin = ggplot2::margin(1, 1, 1, 1)
    )
}


#' OSAR diagnostics (ggplot2) faithfully matching SPiCT base panels
#'
#' Builds a multi-row, multi-column grid of diagnostic panels (optionally
#' including raw data at the top) for **Catch**, optional **Effort**, and each
#' **Index** series found in a SPiCT `rep` object. Rows show, in order:
#' data (optional), OSA residuals, residual ACF, and QQ ("Sample quantiles")
#' panels. Axes, p-value coloring, and QQ scaling are aligned with SPiCT’s base
#' plotting style.
#'
#' **Notes**
#' - Expects `rep$osar` and `rep$diagn` produced by `calc.osa.resid(rep)` and
#'   related SPiCT routines.
#' - X-axis year breaks are integer (annual) if the span is short; otherwise
#'   they fall every 5 years for readability.
#' - QQ panels emulate base `qqnorm/qqline` with ~4% padding
#'   (as in `par(xaxs = "r", yaxs = "r")`).
#'
#' @param rep A fitted SPiCT result list (with `$inp`, `$osar`, `$diagn`).
#' @param lag.max Maximum lag to display in ACF panels. Default `4`.
#' @param qlegend Kept for API parity; legends are not drawn. Default `TRUE`.
#' @param plot.data If `TRUE`, include the top row with log-data by series.
#'   Default `TRUE`.
#' @param mfcol Kept for API parity; not used (layout via \pkg{patchwork}).
#'   Default `FALSE`.
#' @param stamp Caption text for the full figure (e.g., a version string).
#'   Default `get.version()`.
#'
#' @return A \pkg{patchwork} / \pkg{ggplot2} object combining all panels.
#'
#' @examples
#' \dontrun{
#'   # Assuming `rep` is a SPiCT report list with OSA residuals:
#'   p <- plotspict.diagnostic.osar_gg(rep)
#'   p
#' }
#'
#' @export
#' @import ggplot2
#' @importFrom grDevices colorRamp rgb
#' @importFrom grid unit
#' @importFrom stats acf qnorm ppoints quantile na.omit
#' @importFrom patchwork plot_spacer plot_layout plot_annotation
plotspict.diagnostic.osar_gg <- function(rep,
                                         lag.max = 4,
                                         qlegend = TRUE,   # API parity; no legends drawn
                                         plot.data = TRUE,
                                         mfcol = FALSE,    # API parity; not used (patchwork)
                                         stamp = get.version()) {
  check.rep(rep)
  if (!("osar" %in% names(rep)))
    stop("rep$osar not found. Run calc.osa.resid(rep) first.")

  inp <- rep$inp
  thm <- theme_minimal_compact2_diagOSA()

  # ----------------------- helpers ------------------------------------------
  dynamic_year_breaks <- function(xmin, xmax, threshold_years = 8L) {
    if (!is.finite(xmin) || !is.finite(xmax) || xmin == xmax) return(floor(xmin))
    lo <- floor(xmin); hi <- ceiling(xmax)
    nyrs <- hi - lo + 1L
    if (nyrs < threshold_years) {
      return(seq(lo, hi, by = 1L))
    } else {
      start <- lo + ((5L - (lo %% 5L)) %% 5L)
      brks <- seq(start, hi, by = 5L)
      if (!length(brks)) brks <- seq(lo, hi, by = 1L)
      return(brks)
    }
  }
  year_labels_int <- function(x) as.integer(round(x))

  season_bins <- function(tt, ncols = 13) {
    brks <- seq(0, 1, length = ncols)
    ind  <- cut(tt %% 1, brks, right = FALSE, labels = FALSE)
    ind[is.na(ind)] <- 1L
    as.integer(ind)
  }
  season_palette <- function(ncols = 13) {
    cr <- grDevices::colorRamp(c('blue', 'green', 'gold', 'red', 'blue'))
    cols <- cr(seq(0, 1, length.out = ncols))
    grDevices::rgb(cols[,1]/255, cols[,2]/255, cols[,3]/255)
  }

  add_sig_label <- function(g, xr, yr, txt) {
    if (!nzchar(txt)) return(g)
    dx <- xr[2] - xr[1]; dy <- yr[2] - yr[1]
    x_pad <- max(0.03 * dx, 0.15)
    y_pad <- max(0.06 * dy, 0.05)
    g + ggplot2::annotate("text",
                          x = xr[2] - x_pad, y = yr[2] - y_pad,
                          label = txt, hjust = 1, vjust = 1,
                          color = "red", size = 3)
  }

  # ----------------------- SPiCT-consistent Shapiro p-values -----------------
  shapiroC_raw <- as.list(rep$diagn)$shapiroC.p
  shapiroE_raw <- if (inp$nobsE > 0) as.list(rep$diagn)$shapiroE.p else NA_real_

  inds <- grep("shapiroI", names(rep$diagn))
  nms  <- names(rep$diagn)[inds]
  nos  <- as.numeric(unlist(regmatches(nms, gregexpr("[0-9]+", nms))))
  shapiroI_map <- list()
  if (length(nms)) {
    for (k in seq_along(nms)) {
      shapiroI_map[[as.character(nos[k])]] <- as.list(rep$diagn)[[nms[k]]]
    }
  }

  # ----------------------- builders -----------------------------------------
  ylab_data_text <- function(main_label) {
    if (identical(main_label, "Catch"))  return("log catch data")
    if (identical(main_label, "Effort")) return("log effort data")
    if (grepl("^Index", main_label)) {
      j <- sub("^Index\\s+", "", main_label)
      return(paste0("log index ", j, " data"))
    }
    "log data"
  }

  make_data_panel <- function(time, obs, main, vline_at) {
    ncols <- 13
    pal   <- season_palette(ncols)
    bins  <- season_bins(time, ncols)
    df <- data.frame(time = time, y = log(obs), bin = bins)
    xmin <- min(df$time, na.rm = TRUE); xmax <- max(df$time, na.rm = TRUE)

    ggplot2::ggplot(df, ggplot2::aes(time, y)) +
      ggplot2::geom_line(linewidth = 0.5, color = "lightgray", na.rm = TRUE) +
      ggplot2::geom_point(ggplot2::aes(fill = factor(bin)),
                          shape = 21, stroke = 0.3, size = 1.2, color = "black", na.rm = TRUE) +
      ggplot2::geom_vline(xintercept = vline_at, linetype = 3, color = "gray") +
      ggplot2::labs(x = "Time", y = ylab_data_text(main), title = main) +
      ggplot2::scale_x_continuous(
        breaks = dynamic_year_breaks(xmin, xmax), labels = year_labels_int
      ) +
      ggplot2::scale_fill_manual(values = pal, guide = "none") +
      thm
  }

  ylab_resid_text <- function(nm) {
    if (nm == "C") return("Catch OSA residuals")
    if (nm == "E") return("Effort OSA residuals")
    paste0("Index ", sub("^I", "", nm), " OSA residuals")
  }

  make_resid_panel <- function(time, res, panel_key, pval, title_col, vline_at) {
    df   <- data.frame(time = time, res = res)
    dfna <- df[is.na(df$res), , drop = FALSE]
    xmin <- min(df$time, na.rm = TRUE); xmax <- max(df$time, na.rm = TRUE)

    ncols <- 13
    pal   <- season_palette(ncols)
    bins  <- season_bins(df$time, ncols)
    df$bin <- bins

    base <- ggplot2::ggplot(df, ggplot2::aes(time, res)) +
      ggplot2::geom_hline(yintercept = 0, linetype = 3) +
      ggplot2::geom_line(linewidth = 0.7, color = "lightgray", na.rm = TRUE) +
      ggplot2::geom_point(ggplot2::aes(fill = factor(bin)),
                          shape = 21, stroke = 0.3, size = 1.5, color = "black", na.rm = TRUE) +
      ggplot2::scale_fill_manual(values = pal, guide = "none")

    base +
      ggplot2::geom_vline(xintercept = vline_at, linetype = 3, color = "gray") +
      ggplot2::geom_text(data = dfna, ggplot2::aes(time, y = 0, label = "NA"),
                         inherit.aes = FALSE, size = 3, vjust = -0.4) +
      ggplot2::labs(x = "Time", y = ylab_resid_text(panel_key),
                    title = paste0("Bias p-val: ", round(pval, 4))) +
      ggplot2::scale_x_continuous(
        breaks = dynamic_year_breaks(xmin, xmax), labels = year_labels_int
      ) +
      thm +
      ggplot2::theme(plot.title = ggplot2::element_text(color = title_col))
  }

  acf_df <- function(x, lag.max) {
    x <- stats::na.omit(x)
    out <- stats::acf(x, lag.max = lag.max, plot = FALSE, na.action = na.omit)
    data.frame(lag = as.numeric(out$lag[, , 1]),
               acf = as.numeric(out$acf[, , 1]),
               n   = length(x))
  }

  make_acf_panel <- function(x, ylab, LBox_p) {
    dacf <- acf_df(x, lag.max)
    ci <- 1.96 / sqrt(unique(dacf$n))
    pos <- dacf$lag > 0
    sig <- which(abs(dacf$acf[pos]) > ci)
    sig_txt <- if (length(sig)) paste0("lag.signf: ", paste0(dacf$lag[pos][sig], collapse = ",")) else ""
    col_title <- ifelse(LBox_p < 0.05, "red", "forestgreen")

    g <- ggplot2::ggplot(dacf, ggplot2::aes(lag, acf)) +
      ggplot2::geom_hline(yintercept = 0) +
      ggplot2::geom_hline(yintercept = c(-ci, ci), linetype = 2, color = "blue") +
      ggplot2::geom_segment(ggplot2::aes(xend = lag, y = 0, yend = acf), linewidth = 0.5) +
      ggplot2::labs(x = "Lag", y = ylab, title = paste0("LBox p-val: ", round(LBox_p, 4))) +
      thm +
      ggplot2::theme(plot.title = ggplot2::element_text(color = col_title))

    xr <- range(dacf$lag, na.rm = TRUE)
    yr <- range(c(dacf$acf, -ci, ci), na.rm = TRUE)
    add_sig_label(g, xr, yr, sig_txt)
  }

  # ===== Pixel-parity QQ panel with base-R style 4% padding (xaxs/yaxs = "r") =====
  make_qq_panel <- function(x, shap_p_raw) {
    x_ok <- stats::na.omit(x)
    n <- length(x_ok)

    title_col <- if (!is.na(shap_p_raw) && shap_p_raw < 0.05) "red" else "forestgreen"
    title_txt <- if (is.na(shap_p_raw)) "Shapiro p-val: NA"
    else sprintf("Shapiro p-val: %.4f", round(shap_p_raw, 4))

    if (n < 2L) {
      return(
        ggplot2::ggplot() +
          ggplot2::labs(x = "Theoretical Quantiles", y = "Sample Quantiles", title = title_txt) +
          thm +
          ggplot2::theme(
            plot.title   = ggplot2::element_text(color = title_col),
            panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.8)
          )
      )
    }

    # Base R qqnorm ingredients
    z  <- stats::qnorm(stats::ppoints(n))  # theoretical
    xs <- sort(x_ok)                       # sample

    # Base R qqline via quartiles
    yq <- stats::quantile(xs, probs = c(0.25, 0.75), names = FALSE)
    xq <- stats::qnorm(c(0.25, 0.75))
    slope <- (yq[2] - yq[1]) / (xq[2] - xq[1])
    intercept <- yq[1] - slope * xq[1]

    # Exact data ranges (qqnorm reference)
    xlim <- range(z,  finite = TRUE)
    ylim <- range(xs, finite = TRUE)

    df <- data.frame(theoretical = z, sample = xs)

    ggplot2::ggplot(df, ggplot2::aes(theoretical, sample)) +
      ggplot2::geom_abline(intercept = intercept, slope = slope, linewidth = 0.7, color = "black") +
      ggplot2::geom_point(shape = 1, size = 1.5, na.rm = TRUE) +  # open circles, small size
      ggplot2::labs(x = "Theoretical Quantiles", y = "Sample Quantiles", title = title_txt) +
      # Emulate par(xaxs="r", yaxs="r") => ~4% padding each side
      ggplot2::scale_x_continuous(limits = xlim, expand = ggplot2::expansion(mult = 0.04)) +
      ggplot2::scale_y_continuous(limits = ylim, expand = ggplot2::expansion(mult = 0.04)) +
      thm +
      ggplot2::theme(
        plot.title   = ggplot2::element_text(color = title_col),
        panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.8)
      )
  }

  # ----------------------- gather series -------------------------------------
  series <- list()

  # Catch
  series[["C"]] <- list(
    data_time   = inp$timeC,
    data_obs    = inp$obsC,
    osar_time   = rep$osar$timeC,
    osar_res    = rep$osar$logCpres,
    bias_p      = as.list(rep$diagn)$biasC.p,
    LBox_p      = as.list(rep$diagn)$LBoxC.p,
    shap_raw    = shapiroC_raw,
    main_label  = "Catch"
  )

  # Effort (optional)
  if (inp$nobsE > 0) {
    series[["E"]] <- list(
      data_time   = inp$timeE,
      data_obs    = inp$obsE,
      osar_time   = rep$osar$timeE,
      osar_res    = rep$osar$logEpres,
      bias_p      = as.list(rep$diagn)$biasE.p,
      LBox_p      = as.list(rep$diagn)$LBoxE.p,
      shap_raw    = shapiroE_raw,
      main_label  = "Effort"
    )
  }

  # Indices (I1..In)
  if (inp$nindex > 0) {
    for (j in seq_len(inp$nindex)) {
      shap_raw_j <- NA_real_
      key <- as.character(j)
      if (key %in% names(shapiroI_map)) shap_raw_j <- shapiroI_map[[key]]
      if (is.na(shap_raw_j)) {
        nm_fallback <- paste0("shapiroI", j, ".p")
        if (nm_fallback %in% names(rep$diagn)) shap_raw_j <- as.list(rep$diagn)[[nm_fallback]]
      }

      series[[paste0("I", j)]] <- list(
        data_time   = inp$timeI[[j]],
        data_obs    = inp$obsI[[j]],
        osar_time   = rep$osar$timeI[[j]],
        osar_res    = rep$osar$logIpres[[j]],
        bias_p      = as.list(rep$diagn)[[paste0("biasI", j, ".p")]],
        LBox_p      = as.list(rep$diagn)[[paste0("LBoxI", j, ".p")]],
        shap_raw    = shap_raw_j,
        main_label  = paste("Index", j)
      )
    }
  }

  # ----------------------- build panels per column ---------------------------
  cols_data <- list(); cols_res <- list(); cols_acf <- list(); cols_qq <- list()

  for (nm in names(series)) {
    s <- series[[nm]]
    col_res <- if (!is.na(s$bias_p) && s$bias_p < 0.05) "red" else "forestgreen"

    if (plot.data) {
      cols_data[[nm]] <- make_data_panel(time = s$data_time, obs = s$data_obs,
                                         main = s$main_label,
                                         vline_at = s$osar_time[1])
    }

    cols_res[[nm]] <- make_resid_panel(time = s$osar_time, res = s$osar_res,
                                       panel_key = nm,
                                       pval = s$bias_p,
                                       title_col = col_res,
                                       vline_at = s$osar_time[1])

    cols_acf[[nm]] <- make_acf_panel(x = s$osar_res,
                                     ylab = switch(nm,
                                                   C = "Catch ACF",
                                                   E = "Effort ACF",
                                                   paste0("Index ", sub("^I", "", nm), " ACF")),
                                     LBox_p = s$LBox_p)

    cols_qq[[nm]]  <- make_qq_panel(x = s$osar_res, shap_p_raw = s$shap_raw)
  }

  # ----------------------- assemble grid (columns by series) -----------------
  gap <- patchwork::plot_spacer()
  row_gap_rel <- 0.14

  make_row <- function(lst) {
    if (length(lst) == 1) return(lst[[1]])
    Reduce(`+`, lst)
  }

  if (plot.data) {
    top    <- make_row(cols_data)
    mid1   <- make_row(cols_res)
    mid2   <- make_row(cols_acf)
    bot    <- make_row(cols_qq)
    grid   <- top / gap / mid1 / gap / mid2 / gap / bot
    grid   <- grid + patchwork::plot_layout(heights = c(1, row_gap_rel, 1, row_gap_rel, 1, row_gap_rel, 1))
  } else {
    mid1   <- make_row(cols_res)
    mid2   <- make_row(cols_acf)
    bot    <- make_row(cols_qq)
    grid   <- mid1 / gap / mid2 / gap / bot
    grid   <- grid + patchwork::plot_layout(heights = c(1, row_gap_rel, 1, row_gap_rel, 1))
  }

  grid <- grid +
    patchwork::plot_annotation(caption = stamp) &
    ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 1, size = 10))

  return(grid)
}
