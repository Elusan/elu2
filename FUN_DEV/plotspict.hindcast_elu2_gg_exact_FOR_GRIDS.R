# R/plotspict_hindcast_elu2_gg_exact_FOR_GRIDS.R

#' SPiCT hindcast plot (elu2 exact version; grid-friendly)
#'
#' Produce SPiCT-style hindcast plots with ggplot2 for one or more index series.
#' Draws the reference fit (ribbon + line), colored peel prediction lines for
#' each hindcast peel, and a short **dashed** grey tail segment for the final
#' prediction step of each peel. Observations are shown as large circles
#' (white pre-hindcast; colored during hindcast). Legend is forced to **two
#' columns**: column 1 lists *Ref* and peel labels; column 2 shows only *obs*
#' and *pred* entries.
#'
#' @param rep A fitted SPiCT result object **with** `$hindcast` attached and
#'   consumable by `extract.hindcast.info()`.
#' @param add.mase Logical; if `TRUE`, annotate panels with MASE values. Default `TRUE`.
#' @param CI Confidence level for ribbons (e.g., `0.95`). Default `0.95`.
#' @param verbose Logical; handed to the extractor. Default `TRUE`.
#' @param xlim,ylim Optional numeric length-2 vectors for axis limits.
#' @param xlab,ylab Axis labels (kept for signature compatibility). The plot
#'   sets `"Year"` and `"log(Index)"` for SPiCT parity when `ylab` is `NA`.
#' @param plot.log Logical; if `FALSE`, exponentiate `obs/pred/lc/uc`. Default `TRUE`.
#' @param legend.title Optional legend title; default `NULL` (no title).
#' @param legend.pos,legend.ncol,asp Ignored (kept for compatibility).
#' @param stamp Character stamp (e.g., `get.version()`), currently not drawn
#'   (reserved for compatibility).
#'
#' @details
#' * Final segment of each peel is `linetype = "dashed"`, `color = "grey30"`.
#' * Legend uses **two columns**: *Ref* + peels (col 1), *obs*, *pred* (col 2).
#' * Facet strip background/text removed via [theme_minimal_compact2_hind()].
#' * The *two stray circles at (0,0)* problem is avoided by using `NA` dummies
#'   with `na.rm = TRUE` to create legend keys without drawing points.
#'
#' Requires a helper `extract.hindcast.info()` available in the search path
#' (e.g., from **spict** or your own wrapper in **elu2**).
#'
#' @return A [ggplot2::ggplot] object (faceted if multiple indices).
#'
#' @seealso [theme_minimal_compact2_hind()], [plot_hindcast_grid_scenario()]
#'
#' @export
#' @import ggplot2
#' @importFrom stats na.omit
plotspict.hindcast_elu2_gg_exact_FOR_GRIDS <- function(
    rep,
    add.mase = TRUE, CI = 0.95, verbose = TRUE,
    xlim = NULL, ylim = NULL,
    xlab = "Year", ylab = NA,
    plot.log = TRUE,
    legend.title = NULL,
    legend.pos = "topright",
    legend.ncol = 1,
    asp = 2,
    stamp = get.version()
) {
  hcInfo    <- extract.hindcast.info(rep, CI = CI, verbose = verbose)
  inpin     <- rep$inp
  hindcast  <- rep$hindcast
  npeels    <- length(hindcast) - 1L
  peels     <- seq_len(npeels)

  if (max(unlist(inpin$timeI)) == inpin$lastCatchObs) peeling <- 0:(npeels-1L) else peeling <- 1:npeels
  peel.dtc <- ifelse(length(which(abs(diff(
    sapply(hindcast[-1], function(x) max(c(x$inp$timeC + x$inp$dtc,
                                           x$inp$timeE + x$inp$dte))))) < 1)) >= 2, 1, 0)

  if (peel.dtc) {
    tmp <- seq(floor(inpin$timerange[1]), ceiling(inpin$timerange[2] + 100), tail(inpin$dtc, 1))
    hindcastTimes <- tmp[as.integer(cut(inpin$timerangeObs[2], tmp, right = FALSE)) + 1] -
      cumsum(rev(inpin$dtc))[peeling]
  } else {
    hindcastTimes <- ceiling(inpin$timerangeObs[2]) - peeling
  }

  conv <- sapply(hindcast[-1], function(x) x$opt$convergence)
  conv <- rev(ifelse(conv == 0, TRUE, FALSE))

  # palette (SPiCT-like)
  spict_cols_fun <- if (exists("cols", mode = "function")) get("cols", mode = "function") else {
    function() {
      cs <- c("#000000", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
              "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
      c(cs, grDevices::adjustcolor(cs[-1], 0.5), grDevices::adjustcolor(cs[-1], 0.2))
    }
  }
  full_cols <- spict_cols_fun()
  cols <- full_cols[-1][1:npeels]

  peel_index_for_time <- function(t) {
    if (!length(hindcastTimes) || !is.finite(t)) return(NA_integer_)
    w <- which(hindcastTimes <= t)
    if (length(w)) w[1] else NA_integer_
  }
  peel_lab <- paste0(hindcastTimes)

  safe_rbind <- function(lst) {
    lst <- Filter(function(x) !is.null(x) && is.data.frame(x) && nrow(x) > 0, lst)
    if (!length(lst)) return(data.frame())
    do.call(rbind, lst)
  }

  nind.val <- length(hcInfo$index)
  if (nind.val < 1) stop("No valid index in hindcast info (nothing to plot).")

  make_index_df <- function(ii) {
    dat  <- hcInfo$index[[ii]]$dat
    peel0 <- min(dat$peel, na.rm = TRUE)
    dat0  <- dat[dat$peel == peel0, , drop = FALSE]
    if (!plot.log) {
      tr <- c("obs","pred","lc","uc"); dat[,tr] <- lapply(dat[,tr], exp); dat0[,tr] <- lapply(dat0[,tr], exp)
    }

    by_peel <- split(dat, dat$peel)
    peel_draw <- lapply(seq_along(by_peel), function(k) {
      dd <- by_peel[[k]]; this_peel <- unique(dd$peel)
      if (length(this_peel) != 1 || this_peel == 0) return(NULL)
      if (!isTRUE(conv[this_peel])) return(NULL)
      if (!isTRUE(tail(dd$iuse, 1) == 0)) return(NULL)
      dd
    })
    peel_draw_df <- Filter(function(dd) is.data.frame(dd) && nrow(dd) > 0, peel_draw)

    # peel lines (colored by last time of the peel)
    df_peels <- safe_rbind(peel_draw)
    if (nrow(df_peels)) {
      df_peels <- do.call(rbind, lapply(split(df_peels, df_peels$peel), function(dd) {
        if (!nrow(dd)) return(dd)
        pid <- peel_index_for_time(max(dd$time, na.rm = TRUE))
        dd$.peel_lab <- factor(ifelse(is.na(pid), NA, peel_lab[pid]), levels = peel_lab)
        dd
      }))
    }

    # short dashed grey tail segment (last 2 points)
    seg_list <- lapply(peel_draw_df, function(dd) {
      if (!is.data.frame(dd) || nrow(dd) < 2) return(NULL)
      nd <- nrow(dd); pid <- peel_index_for_time(dd$time[nd])
      data.frame(x = dd$time[nd-1], xend = dd$time[nd],
                 y = dd$pred[nd-1], yend = dd$pred[nd],
                 peel = unique(dd$peel),
                 .peel_lab = factor(ifelse(is.na(pid), NA, peel_lab[pid]), levels = peel_lab))
    })
    df_segments <- safe_rbind(seg_list)

    # observations split
    t0        <- suppressWarnings(min(hindcastTimes, na.rm = TRUE))
    obs_pre   <- dat0[dat0$time <  t0, , drop = FALSE]
    obs_hind  <- dat0[dat0$time >= t0, , drop = FALSE]
    if (nrow(obs_hind)) {
      times <- if ("time" %in% names(obs_hind)) obs_hind$time else numeric(0)
      peel_ids <- if (length(times)) vapply(times, peel_index_for_time, integer(1)) else rep(NA_integer_, nrow(obs_hind))
      if (length(peel_ids) != nrow(obs_hind)) peel_ids <- rep(NA_integer_, nrow(obs_hind))
      obs_hind$.peel_id  <- peel_ids
      ok <- !is.na(peel_ids) & peel_ids >= 1 & peel_ids <= length(conv)
      obs_hind$.is_conv  <- FALSE; if (any(ok)) obs_hind$.is_conv[ok] <- conv[peel_ids[ok]]
      obs_hind$.peel_lab <- factor(ifelse(is.na(peel_ids), NA, peel_lab[peel_ids]), levels = peel_lab)
    }

    # end points per peel
    end_list <- lapply(peel_draw_df, function(dd) if (nrow(dd)) dd[nrow(dd), , drop = FALSE] else NULL)
    df_pred_end <- safe_rbind(end_list)
    if (nrow(df_pred_end)) {
      peel_ids <- vapply(df_pred_end$time, peel_index_for_time, integer(1))
      df_pred_end$.peel_lab <- factor(ifelse(is.na(peel_ids), NA, peel_lab[peel_ids]), levels = peel_lab)
    }

    dat$Index <- dat0$Index <- ii
    if (nrow(df_peels))    df_peels$Index    <- ii
    if (nrow(df_segments)) df_segments$Index <- ii
    if (nrow(obs_pre))     obs_pre$Index     <- ii
    if (nrow(obs_hind))    obs_hind$Index    <- ii
    if (nrow(df_pred_end)) df_pred_end$Index <- ii

    list(dat0 = dat0, peels = df_peels, segments = df_segments,
         obs_pre = obs_pre, obs_hind = obs_hind, pred_end = df_pred_end)
  }

  idx_list   <- lapply(seq_len(nind.val), make_index_df)
  bind0      <- safe_rbind(lapply(idx_list, `[[`, "dat0"))
  bind_peels <- safe_rbind(lapply(idx_list, `[[`, "peels"))
  bind_segs  <- safe_rbind(lapply(idx_list, `[[`, "segments"))
  bind_pre   <- safe_rbind(lapply(idx_list, `[[`, "obs_pre"))
  bind_hind  <- safe_rbind(lapply(idx_list, `[[`, "obs_hind"))
  bind_end   <- safe_rbind(lapply(idx_list, `[[`, "pred_end"))

  used_peel_levels <- character(0)
  if (nrow(bind_hind)) {
    used_peel_levels <- as.character(stats::na.omit(unique(bind_hind$.peel_lab[bind_hind$.is_conv %in% TRUE])))
  }
  if (!length(used_peel_levels) && nrow(bind_peels)) {
    used_peel_levels <- as.character(stats::na.omit(unique(bind_peels$.peel_lab)))
  }
  used_peel_levels <- used_peel_levels[order(match(used_peel_levels, peel_lab))]
  peel_color_values <- setNames(cols[match(used_peel_levels, peel_lab)], used_peel_levels)

  p <- ggplot2::ggplot()

  if (nrow(bind0)) {
    p <- p +
      ggplot2::geom_ribbon(data = bind0,
                           ggplot2::aes(x = time, ymin = lc, ymax = uc),
                           inherit.aes = FALSE, fill = "grey80") +
      ggplot2::geom_line(data = bind0,
                         ggplot2::aes(x = time, y = pred, color = "Ref"),
                         inherit.aes = FALSE, linewidth = 0.8)
  }

  if (nrow(bind_peels)) {
    p <- p +
      ggplot2::geom_line(data = bind_peels,
                         ggplot2::aes(x = time, y = pred,
                                      group = interaction(Index, peel),
                                      color = .peel_lab),
                         linewidth = 0.8)
  }

  if (nrow(bind_segs)) {
    p <- p +
      ggplot2::geom_segment(
        data = bind_segs,
        ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
        linewidth = 0.8, linetype = "dashed", color = "grey30"
      )
  }

  size_obs  <- 2.8
  size_pred <- 1.2

  if (nrow(bind_pre)) {
    p <- p +
      ggplot2::geom_point(data = bind_pre,
                          ggplot2::aes(x = time, y = obs),
                          shape = 21, fill = "white", color = "black",
                          stroke = 0.5, size = size_obs)
  }
  if (nrow(bind_hind)) {
    if (any(bind_hind$.is_conv %in% TRUE)) {
      p <- p +
        ggplot2::geom_point(data = subset(bind_hind, .is_conv %in% TRUE),
                            ggplot2::aes(x = time, y = obs, fill = .peel_lab),
                            shape = 21, color = "black", stroke = 0.5, size = size_obs,
                            show.legend = FALSE)
    }
    if (any(!(bind_hind$.is_conv %in% TRUE))) {
      p <- p +
        ggplot2::geom_point(data = subset(bind_hind, !(.is_conv %in% TRUE)),
                            ggplot2::aes(x = time, y = obs),
                            shape = 21, fill = "white", color = "black",
                            stroke = 0.5, size = size_obs,
                            show.legend = FALSE)
    }
  }
  if (nrow(bind_end)) {
    p <- p +
      ggplot2::geom_point(data = bind_end,
                          ggplot2::aes(x = time, y = pred, fill = .peel_lab),
                          shape = 21, color = "black", stroke = 0.5, size = size_pred,
                          show.legend = FALSE)
  }

  # Legend (2 columns; col1 = Ref + peels, col2 = obs/pred) without drawing stray points
  base_breaks <- c("Ref", used_peel_levels)
  A <- length(base_breaks)
  n_rows <- max(A, 2L)
  color_breaks <- c(base_breaks, "obs", "pred")
  color_values <- c("Ref" = "#000000",
                    peel_color_values,
                    "obs" = "#000000", "pred" = "#000000")

  dummy_obs  <- data.frame(x = NA_real_, y = NA_real_, key = "obs")
  dummy_pred <- data.frame(x = NA_real_, y = NA_real_, key = "pred")

  p <- p +
    ggplot2::geom_point(
      data = dummy_obs,
      ggplot2::aes(x = x, y = y, color = "obs"),
      shape = 21, fill = "white", stroke = 0.5, size = size_obs,
      inherit.aes = FALSE, na.rm = TRUE
    ) +
    ggplot2::geom_point(
      data = dummy_pred,
      ggplot2::aes(x = x, y = y, color = "pred"),
      shape = 21, fill = "white", stroke = 0.5, size = size_pred,
      inherit.aes = FALSE, na.rm = TRUE
    ) +
    ggplot2::scale_color_manual(
      limits = color_breaks,
      values = color_values,
      drop   = FALSE,
      name   = if (is.null(legend.title)) "" else legend.title,
      guide  = ggplot2::guide_legend(
        ncol = 2,
        nrow = n_rows,
        byrow = FALSE,
        keywidth  = grid::unit(0.8, "lines"),
        keyheight = grid::unit(0.6, "lines"),
        label.hjust = 0,
        override.aes = list(
          linetype  = c(rep(1, A), NA, NA),
          shape     = c(rep(NA, A), 21, 21),
          fill      = c(rep(NA, A), "white", "white"),
          size      = c(rep(0.9, A), 1.8, 1.4),
          linewidth = c(rep(0.9, A), NA, NA)
        )
      )
    ) +
    ggplot2::scale_fill_manual(values = setNames(cols, peel_lab), guide = "none") +
    ggplot2::facet_wrap(~ Index, scales = "free_x") +
    theme_minimal_compact2_hind() +
    ggplot2::theme(strip.text = ggplot2::element_blank()) +
    ggplot2::labs(x = "Year", y = if (is.na(ylab)) "log(Index)" else ylab)

  if (!is.null(xlim) && !is.list(xlim) && length(xlim) == 2) p <- p + ggplot2::coord_cartesian(xlim = xlim)
  if (!is.null(ylim) && !is.list(ylim) && length(ylim) == 2) p <- p + ggplot2::coord_cartesian(ylim = ylim)

  if (isTRUE(add.mase) && nrow(hcInfo$mase) && nrow(bind0)) {
    place_df <- do.call(rbind, lapply(split(bind0, bind0$Index), function(d) {
      if (!nrow(d)) return(NULL)
      xr <- range(d$time, na.rm = TRUE); yr <- range(c(d$lc, d$uc), na.rm = TRUE)
      data.frame(Index = unique(d$Index)[1], x = mean(xr), y = yr[2] - 0.02*diff(yr))
    }))
    if (!is.null(place_df) && nrow(place_df)) {
      mase_df <- data.frame(Index = seq_len(nind.val),
                            label = paste0("Index ", seq_len(nind.val), ": MASE = ",
                                           signif(hcInfo$mase[, "MASE"], 3)))
      mase_df <- merge(mase_df, place_df, by = "Index", all.x = TRUE)
      p <- p + ggplot2::geom_text(data = mase_df,
                                  ggplot2::aes(x = x, y = y, label = label),
                                  inherit.aes = FALSE, fontface = 2, size = 3,
                                  hjust = 0.5, vjust = 1)
    }
  }

  p
}
