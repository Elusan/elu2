#' Reusable colour palette for management scenarios
#'
#' @description
#' Returns a repeatable vector of colours to map across an arbitrary number
#' of management scenarios. Colours recycle to the requested length.
#'
#' @param n Integer. Number of colours required.
#' @return A character vector of length `n` with colour names.
#' @keywords internal
#' @noRd
man_cols <- function(n) {
  base <- c('darkmagenta','cyan3','darkgreen','coral1','black',
            'magenta','gold','green','cadetblue3','chocolate3',
            'darkolivegreen3','cyan','darkred')
  rep(base, length.out = n)
}

# ---- small safe-field helpers (internal) ----
.sf  <- function(x, nm, default = NULL) if (!is.null(x) && nm %in% names(x)) x[[nm]] else default
.sfb <- function(x, nm) isTRUE(.sf(x, nm, FALSE))
.sfi <- function(x, nm, default = integer(0)) .sf(x, nm, default)
.sfn <- function(x, nm, default = numeric(0)) .sf(x, nm, default)

# ---- wrapper: accept a fitted spictcls OR a rep$man list ----------------------
.as_spict_like <- function(x) {
  if (inherits(x, "spictcls")) return(x)
  if (is.list(x) && length(x) && all(vapply(x, function(z) inherits(z, "spictcls"), logical(1)))) {
    first <- x[[1]]
    inp_fallback <- if (!is.null(first$inp)) first$inp else list()
    out <- list(inp = inp_fallback, man = x)
    class(out) <- "spictcls"
    return(out)
  }
  if (is.list(x) && !length(x)) {
    out <- list(inp = list(), man = list()); class(out) <- "spictcls"; return(out)
  }
  stop("Input must be a fitted spictcls, or a list like rep$man containing spictcls objects.")
}

#' Observed indices on the relative B/Bmsy scale
#' @keywords internal
#' @noRd
obs_indices_rel_df <- function(rep_main) {
  stopifnot(inherits(rep_main, "spictcls"))
  qest <- try(get.par("logq", rep_main, exp = TRUE), silent = TRUE)
  Bmsy <- try(get.par("logBmsy", rep_main, exp = TRUE), silent = TRUE)
  if (inherits(qest, "try-error") || inherits(Bmsy, "try-error")) return(NULL)

  Bmsy2 <- if (is.null(dim(Bmsy))) as.numeric(Bmsy[2]) else as.numeric(utils::tail(Bmsy[, 2], 1))
  if (!is.finite(Bmsy2) || Bmsy2 <= 0) return(NULL)

  inp  <- rep_main$inp
  mapq <- .sfi(inp, "mapq", integer(0))

  .one_index <- function(i) {
    ti <- .sf(inp, "timeI", list())[[i]]
    oi <- .sf(inp, "obsI",  list())[[i]]
    if (is.null(ti) || is.null(oi) || length(oi) != length(ti)) return(NULL)

    qi_idx <- if (length(mapq) >= i) mapq[i] else i
    if (is.null(dim(qest))) {
      qhat <- as.numeric(qest[2])
    } else {
      if (!is.finite(qi_idx) || qi_idx < 1 || qi_idx > nrow(qest)) return(NULL)
      qhat <- as.numeric(qest[qi_idx, 2])
    }
    if (!is.finite(qhat) || qhat <= 0) return(NULL)

    data.frame(index  = paste0("Index", i),
               time   = ti,
               obsrel = (oi / qhat) / Bmsy2,
               row.names = NULL)
  }

  out <- Filter(NROW, list(.one_index(1L), .one_index(2L)))
  if (!length(out)) NULL else do.call(rbind, out)
}

#' Historical B/Bmsy and F/Fmsy up to the last observed time
#' @keywords internal
#' @noRd
get_base_BB_FF_pre <- function(rep_man, CI = 0.95) {
  stopifnot(inherits(rep_man, "spictcls"))
  t_last_obs <- .sfn(rep_man$inp, "timerangeObs", c(NA_real_, NA_real_))[2]
  tt <- .sfn(rep_man$inp, "time", numeric(0))
  if (!length(tt)) return(list(t_last_obs = NA_real_, BB = data.frame(), FF = data.frame()))

  keep <- which(tt <= t_last_obs)

  BB <- try(get.par("logBBmsy", rep_man, exp = TRUE), silent = TRUE)
  FF <- try(get.par("logFFmsy", rep_man, exp = TRUE), silent = TRUE)
  if (inherits(BB, "try-error") || inherits(FF, "try-error")) {
    return(list(t_last_obs = t_last_obs, BB = data.frame(), FF = data.frame()))
  }

  if (NROW(BB) > length(tt)) BB <- BB[seq_along(tt), , drop = FALSE]
  if (NROW(FF) > length(tt)) FF <- FF[seq_along(tt), , drop = FALSE]

  BB <- BB[keep, , drop = FALSE]
  FF <- FF[keep, , drop = FALSE]

  list(
    t_last_obs = t_last_obs,
    BB = if (length(keep)) data.frame(time = tt[keep], lwr = BB[,1], est = BB[,2], upr = BB[,3], row.names = NULL) else data.frame(),
    FF = if (length(keep)) data.frame(time = tt[keep], lwr = FF[,1], est = FF[,2], upr = FF[,3], row.names = NULL) else data.frame()
  )
}

#' Prepare data frames for management panels (fit or manage)
#' @keywords internal
#' @noRd
elu2_prepare_manage_panel_data <- function(rep_man, CI = 0.95) {
  stopifnot(inherits(rep_man, "spictcls"))

  has_man <- ("man" %in% names(rep_man)) && length(rep_man$man)

  co_df <- if (!is.null(rep_man$inp$timeC) && !is.null(rep_man$inp$obsC)) {
    data.frame(time = rep_man$inp$timeC,
               catch = rep_man$inp$obsC,
               catch_type = "Observed",
               row.names = NULL)
  } else NULL

  if (has_man) {
    scenarios <- names(rep_man$man)
    ns <- length(scenarios)

    bb_list <- vector("list", ns)
    ff_list <- vector("list", ns)
    cp_list <- vector("list", ns)
    eval_list <- vector("list", ns)

    for (i in seq_len(ns)) {
      sc <- scenarios[i]
      rp <- rep_man$man[[sc]]

      BB <- try(get.par("logBBmsy", rp, exp = TRUE), silent = TRUE)
      FF <- try(get.par("logFFmsy", rp, exp = TRUE), silent = TRUE)
      if (inherits(BB, "try-error") || inherits(FF, "try-error")) next

      rnt <- suppressWarnings(as.numeric(rownames(BB)))
      if (!length(rnt) || any(!is.finite(rnt))) {
        ti <- .sfn(rp$inp, "time", numeric(0))
        rnt <- if (length(ti) >= NROW(BB)) ti[seq_len(NROW(BB))] else seq_len(NROW(BB))
      }

      bb_list[[i]] <- data.frame(time = rnt, lwr = BB[,1], est = BB[,2], upr = BB[,3],
                                 scenario = sc, row.names = NULL)
      ff_list[[i]] <- data.frame(time = rnt, lwr = FF[,1], est = FF[,2], upr = FF[,3],
                                 scenario = sc, row.names = NULL)

      CP  <- try(get.par("logCpred", rp, exp = TRUE), silent = TRUE)
      if (!inherits(CP, "try-error")) {
        tcp <- .sfn(rp$inp, "timeCpred", numeric(0))
        if (!length(tcp)) tcp <- seq_len(NROW(CP))
        tcp <- tcp[seq_len(min(length(tcp), NROW(CP)))]
        cp_list[[i]] <- data.frame(time = tcp, lwr = CP[seq_along(tcp),1], catch = CP[seq_along(tcp),2], upr = CP[seq_along(tcp),3],
                                   scenario = sc, catch_type = "Predicted", row.names = NULL)
      }

      ti  <- .sfn(rp$inp, "time", numeric(0))
      tE  <- .sfn(rp$inp, "maneval", NA_real_)
      if (length(ti) && is.finite(tE)) {
        idx <- which.min(abs(ti - tE))
        if (idx >= 1 && idx <= NROW(FF)) {
          eval_list[[i]] <- data.frame(scenario = sc, maneval = tE,
                                       BB_est = BB[idx,2], FF_est = FF[idx,2],
                                       row.names = NULL)
        }
      }
    }

    t0 <- if (!is.null(rep_man$inp$maninterval) && length(rep_man$inp$maninterval) >= 1L)
      rep_man$inp$maninterval[1] else NA_real_

    list(
      bbmsy      = if (length(Filter(NROW, bb_list))) do.call(rbind, bb_list) else data.frame(time=numeric(0), lwr=numeric(0), est=numeric(0), upr=numeric(0), scenario=character(0)),
      ffmsy      = if (length(Filter(NROW, ff_list))) do.call(rbind, ff_list) else data.frame(time=numeric(0), lwr=numeric(0), est=numeric(0), upr=numeric(0), scenario=character(0)),
      catch_pred = if (length(Filter(NROW, cp_list))) do.call(rbind, cp_list) else data.frame(time=numeric(0), lwr=numeric(0), catch=numeric(0), upr=numeric(0), scenario=character(0), catch_type=character(0)),
      catch_obs  = co_df,
      eval_pts   = if (length(Filter(NROW, eval_list))) do.call(rbind, eval_list) else NULL,
      obsI_rel   = obs_indices_rel_df(rep_man),
      t0         = t0
    )

  } else {
    list(
      bbmsy      = data.frame(time = numeric(0), lwr = numeric(0), est = numeric(0),
                              upr = numeric(0), scenario = character(0)),
      ffmsy      = data.frame(time = numeric(0), lwr = numeric(0), est = numeric(0),
                              upr = numeric(0), scenario = character(0)),
      catch_pred = data.frame(time = numeric(0), lwr = numeric(0), catch = numeric(0),
                              upr = numeric(0), scenario = character(0), catch_type = character(0)),
      catch_obs  = co_df,
      eval_pts   = NULL,
      obsI_rel   = obs_indices_rel_df(rep_man),
      t0         = NA_real_
    )
  }
}

#' Robust time step for catch series
#' @keywords internal
#' @noRd
.safe_dt_step <- function(rep) {
  if (!is.null(rep$inp$timeC) && length(rep$inp$timeC) > 1L) {
    dt <- stats::median(diff(rep$inp$timeC))
    if (is.finite(dt) && dt > 0) return(dt)
  }
  if (!is.null(rep$inp$dtc)) {
    dt <- suppressWarnings(min(rep$inp$dtc, na.rm = TRUE))
    if (is.finite(dt) && dt > 0) return(dt)
  }
  1
}

#' Management interval vertical lines for B/Bmsy and F/Fmsy panels
#' @keywords internal
#' @noRd
add_management_vlines_BF_good <- function(rep,
                                          color = "grey30",
                                          linetype = "dashed",
                                          linewidth = 0.6,
                                          lineend = "butt") {
  mi <- .sfn(rep$inp, "maninterval", numeric(0))
  if (!length(mi) || !is.finite(mi[1])) return(list())

  left  <- mi[1]
  right <- if (length(mi) >= 2L && is.finite(mi[2])) mi[2] else NA_real_

  out <- list(ggplot2::geom_vline(
    xintercept = left, color = color, linetype = linetype,
    linewidth = linewidth, lineend = lineend
  ))
  if (is.finite(right)) {
    out[[length(out) + 1L]] <- ggplot2::geom_vline(
      xintercept = right, color = color, linetype = linetype,
      linewidth = linewidth, lineend = lineend
    )
  }
  out
}

#' Management vertical lines for Catch panel (not used here but kept)
#' @keywords internal
#' @noRd
add_management_vlines_catch_good <- function(rep,
                                             color = "grey30",
                                             linetype = "dashed",
                                             linewidth = 0.6,
                                             lineend = "butt") {
  mi <- .sfn(rep$inp, "maninterval", numeric(0))
  last_obs_start <- if (!is.null(rep$inp$timeC) && length(rep$inp$timeC)) {
    utils::tail(rep$inp$timeC, 1)
  } else NA_real_

  if (!length(mi) || !is.finite(mi[1])) {
    if (is.finite(last_obs_start)) {
      return(list(ggplot2::geom_vline(
        xintercept = last_obs_start,
        color = color, linetype = linetype,
        linewidth = linewidth, lineend = lineend
      )))
    } else {
      return(list())
    }
  }

  old_left <- mi[1]
  dt_step  <- .safe_dt_step(rep)
  new_left  <- if (is.finite(last_obs_start)) last_obs_start else (old_left - dt_step)
  new_right <- old_left

  list(
    ggplot2::geom_vline(xintercept = new_left,
                        color = color, linetype = linetype,
                        linewidth = linewidth, lineend = lineend),
    ggplot2::geom_vline(xintercept = new_right,
                        color = color, linetype = linetype,
                        linewidth = linewidth, lineend = lineend)
  )
}

#' B/Bmsy panel (ggplot2) with optional scenario overlays
#'
#' B/Bmsy panel (ggplot2) with optional scenario overlays
#'
#' @export
my_plot_manage_bbmsy_panel <-  function(rep_man,
                                              scenario_color = NULL,
                                              show_CIs = TRUE,
                                              CI = 0.95,
                                              show_legend = TRUE) {
  rep_man <- .as_spict_like(rep_man)
  stopifnot(inherits(rep_man, "spictcls"))

  theme_minimal_compact2_good_local <- function(base_size = 10, base_family = "") {
    ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12),
        axis.title = ggplot2::element_text(face = "bold", size = 12),
        axis.text  = ggplot2::element_text(size = 10, face = "bold"),
        legend.position = c(0.4, 0.98),
        legend.justification = c("left", "top"),
        legend.background = ggplot2::element_blank(),
        legend.box.background = ggplot2::element_blank(),
        legend.title = ggplot2::element_blank(),
        legend.text  = ggplot2::element_text(size = 10),
        legend.key.size = grid::unit(1, "lines"),
        panel.grid = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(fill = NA, colour = "grey45", linewidth = 1.2),
        axis.ticks = ggplot2::element_line(linewidth = 0.5, color = "grey45"),
        axis.ticks.length = grid::unit(3, "pt"),
        strip.background = ggplot2::element_rect(fill = "grey45", color = "grey45", linewidth = 0.5),
        strip.text = ggplot2::element_text(face = "bold", size = ggplot2::rel(1)),
        text = ggplot2::element_text(face = "bold", size = 10),
        plot.margin = ggplot2::margin(3, 3, 3, 3)
      )
  }

  dat <- elu2_prepare_manage_panel_data(rep_man, CI = CI)

  sc_order <- c("currentCatch","currentF","Fmsy","noF",
                "reduceF25","increaseF25","msyHockeyStick","ices")
  label_map <- c(
    currentCatch   = "1. Keep current catch",
    currentF       = "2. Keep current F",
    Fmsy           = "3. Fish at Fmsy",
    noF            = "4. No fishing",
    reduceF25      = "5. Reduce F by 25%",
    increaseF25    = "6. Increase F by 25%",
    msyHockeyStick = "7. MSY hockey-stick rule",
    ices           = "8. ICES advice rule"
  )

  sc_present <- if (("man" %in% names(rep_man)) && length(rep_man$man)) names(rep_man$man) else character(0)
  sc_core   <- intersect(sc_order, sc_present)
  sc_other  <- setdiff(sc_present, sc_core)
  scs_final <- c(sc_core, sort(sc_other))

  if (!is.null(dat$bbmsy) && nrow(dat$bbmsy)) {
    dat$bbmsy$scenario <- factor(dat$bbmsy$scenario, levels = scs_final, ordered = TRUE)
  }

  if (is.null(scenario_color)) {
    cols <- if (length(scs_final)) { x <- man_cols(length(scs_final)); names(x) <- scs_final; x } else NULL
  } else {
    cols <- scenario_color
    if (length(scs_final)) {
      if (is.null(names(cols))) names(cols) <- scs_final
      cols <- cols[match(scs_final, names(cols))]
      if (anyNA(cols)) {
        fill <- man_cols(length(scs_final))
        cols[is.na(cols)] <- fill[is.na(cols)]
      }
      names(cols) <- scs_final
    }
  }

  lab_final <- if (length(scs_final)) ifelse(scs_final %in% names(label_map),
                                             unname(label_map[scs_final]),
                                             scs_final) else character(0)

  spict_blue_mean    <- "blue"
  spict_blue_ci_line <- grDevices::rgb(0, 0, 1, 0.20)
  spict_blue_ci_fill <- grDevices::rgb(0, 0, 1, 0.10)
  darker_color <- function(col, f = 0.8) {
    m <- grDevices::col2rgb(col) / 255
    grDevices::rgb(m[1]*f, m[2]*f, m[3]*f)
  }
  dot_blue <- darker_color(spict_blue_mean, 0.6)

  base_pre <- get_base_BB_FF_pre(rep_man, CI = CI)

  pB <- ggplot2::ggplot() +
    { if (show_CIs && NROW(base_pre$BB)) ggplot2::geom_ribbon(
      data = base_pre$BB,
      ggplot2::aes(x = time, ymin = lwr, ymax = upr),
      inherit.aes = FALSE, fill = spict_blue_ci_fill, color = NA
    ) } +
    { if (show_CIs && NROW(base_pre$BB)) list(
      ggplot2::geom_line(data = base_pre$BB, ggplot2::aes(x = time, y = lwr),
                         inherit.aes = FALSE, color = spict_blue_ci_line, linewidth = 0.6),
      ggplot2::geom_line(data = base_pre$BB, ggplot2::aes(x = time, y = upr),
                         inherit.aes = FALSE, color = spict_blue_ci_line, linewidth = 0.6)
    ) } +
    { if (NROW(base_pre$BB)) ggplot2::geom_line(
      data = base_pre$BB,
      ggplot2::aes(x = time, y = est),
      inherit.aes = FALSE,
      color = spict_blue_mean,
      linewidth = 0.8
    ) } +
    ggplot2::geom_line(
      data = dat$bbmsy,
      ggplot2::aes(x = time, y = est, color = scenario),
      linewidth = 0.8, na.rm = TRUE
    ) +
    { if (!is.null(dat$obsI_rel)) {
      list(
        ggplot2::geom_point(
          data = dat$obsI_rel[dat$obsI_rel$index == "Index1", , drop = FALSE],
          ggplot2::aes(x = time, y = obsrel), inherit.aes = FALSE,
          shape = 16, size = 3, color = dot_blue
        ),
        ggplot2::geom_point(
          data = dat$obsI_rel[dat$obsI_rel$index == "Index2", , drop = FALSE],
          ggplot2::aes(x = time, y = obsrel), inherit.aes = FALSE,
          shape = 22, size = 3, fill = "green", color = "black", stroke = 0.8
        )
      )
    } } +
    ggplot2::geom_hline(yintercept = 1, linetype = "solid") +
    { if (length(scs_final)) {
      add_management_vlines_BF_good(rep_man, color = "grey75",
                                    linetype = "solid", linewidth = 0.4, lineend = "butt")
    } else if (is.finite(base_pre$t_last_obs)) {
      ggplot2::geom_vline(xintercept = base_pre$t_last_obs,
                          color = "grey75", linetype = "solid", linewidth = 0.4)
    } else NULL } +
    ggplot2::labs(title= "Relative biomass", x = "Year", y = expression(bold(B/B[MSY]))) +
    {
      if (length(scs_final)) {
        ggplot2::scale_color_manual(
          values  = cols, limits = scs_final, drop = FALSE,
          labels  = lab_final, name = NULL,
          guide   = if (show_legend) "legend" else "none",
          na.translate = FALSE
        )
      } else {
        ggplot2::scale_color_discrete(guide = "none")
      }
    } +
    ggplot2::scale_fill_manual(values = if (is.null(cols)) NA else cols, guide = "none") +
    ggplot2::geom_line(
      data = base_pre$BB,
      mapping = ggplot2::aes(x = time, y = est, group = 1),
      inherit.aes = FALSE,
      colour = I(spict_blue_mean),
      linewidth = 0.9,
      show.legend = FALSE
    ) +
    theme_minimal_compact2_good_local()

  # ---------- Fitted-only mode: keep continuity across the vertical line ----------
  if (!length(scs_final)) {
    BB_all <- try(get.par("logBBmsy", rep_man, exp = TRUE, CI = CI), silent = TRUE)
    if (!inherits(BB_all, "try-error")) {
      tt <- suppressWarnings(as.numeric(rownames(BB_all)))
      if (!length(tt) || any(!is.finite(tt))) {
        tt_full <- .sfn(rep_man$inp, "time", numeric(0))
        tt <- if (length(tt_full) >= nrow(BB_all)) tt_full[seq_len(nrow(BB_all))] else rep_len(tt_full, nrow(BB_all))
      }
      df_all <- data.frame(time = tt, lwr = BB_all[,1], est = BB_all[,2], upr = BB_all[,3])

      tlo <- base_pre$t_last_obs
      if (is.finite(tlo)) {
        # CHANGED: start post at >= t_last_obs so ribbons/lines touch the vline
        df_post <- df_all[df_all$time >= tlo, , drop = FALSE]

        if (nrow(df_post)) {
          if (isTRUE(show_CIs)) {
            pB <- pB +
              ggplot2::geom_ribbon(
                data = transform(df_post, ymin = lwr, ymax = upr),
                ggplot2::aes(x = time, ymin = ymin, ymax = ymax),
                inherit.aes = FALSE, fill = spict_blue_ci_fill, color = NA
              ) +
              ggplot2::geom_line(
                data = df_post, ggplot2::aes(x = time, y = lwr),
                inherit.aes = FALSE, color = spict_blue_ci_line, linewidth = 0.6
              ) +
              ggplot2::geom_line(
                data = df_post, ggplot2::aes(x = time, y = upr),
                inherit.aes = FALSE, color = spict_blue_ci_line, linewidth = 0.6
              )
          }
          pB <- pB +
            ggplot2::geom_line(
              data = df_post, ggplot2::aes(x = time, y = est),
              linetype = "dotted", color = spict_blue_mean, linewidth = 0.8
            )
        }
      }
    }
  }

  pB
}
