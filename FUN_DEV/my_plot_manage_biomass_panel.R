# ---------- Internal helpers (no export) --------------------------------------

#' Reusable colour palette for management scenarios
#' @keywords internal
#' @noRd
man_cols <- function(n) {
  base <- c('darkmagenta','cyan3','darkgreen','coral1','black',
            'magenta','gold','green','cadetblue3','chocolate3',
            'darkolivegreen3','cyan','darkred')
  rep(base, length.out = n)
}

# Tiny safe accessors
.sf  <- function(x, nm, default = NULL) if (!is.null(x) && nm %in% names(x)) x[[nm]] else default
.sfb <- function(x, nm) isTRUE(.sf(x, nm, FALSE))
.sfi <- function(x, nm, default = integer(0)) .sf(x, nm, default)
.sfn <- function(x, nm, default = numeric(0)) .sf(x, nm, default)

# Prefer rownames(time) else fallback to inp$time
# par_matrix is result from spict::get.par(...)
.spict_time_from_par <- function(model, par_matrix) {
  rn <- suppressWarnings(as.numeric(rownames(par_matrix)))
  if (length(rn) == nrow(par_matrix) && any(is.finite(rn))) return(rn)
  t_full <- as.numeric(.sfn(model$inp, "time", numeric(0)))
  if (length(t_full) >= nrow(par_matrix)) return(t_full[seq_len(nrow(par_matrix))])
  if (length(t_full)) return(rep_len(t_full, nrow(par_matrix)))
  seq_len(nrow(par_matrix)) # last resort
}

# Overall “end of observation” time used when $man is absent
.spict_obs_end_overall <- function(model) {
  inp <- model$inp
  tr  <- .sfn(inp, "timerangeObs", c(NA_real_, NA_real_))
  if (length(tr) >= 2 && is.finite(tr[2])) return(tr[2])
  if (!is.null(inp$timeC) && length(inp$timeC))       return(utils::tail(inp$timeC, 1))
  if (!is.null(inp$timeI) && length(inp$timeI) && length(inp$timeI[[1]]) > 0)
    return(utils::tail(inp$timeI[[1]], 1))
  if (!is.null(inp$time) && length(inp$time))         return(max(inp$time, na.rm = TRUE))
  NA_real_
}

# Build a ribbon dataframe from lwr/upr vectors
.spict_make_ribbon_df <- function(time, lwr, upr) {
  ok <- is.finite(time) & is.finite(lwr) & is.finite(upr)
  data.frame(time = time[ok], ymin = pmax(lwr[ok], 0), ymax = pmax(upr[ok], 0))
}

# q-scaled index points for first two indices (optional)
.spict_index_points <- function(model) {
  out <- list()
  inp <- model$inp
  if (!is.null(inp$timeI) && length(inp$timeI) >= 1 && !is.null(inp$obsI)) {
    qest <- try(elu2::get.par("logq", model, exp = TRUE), silent = TRUE)
    if (!inherits(qest, "try-error")) {
      for (i in seq_along(inp$timeI)) {
        ti <- inp$timeI[[i]]; oi <- inp$obsI[[i]]
        if (is.null(ti) || is.null(oi) || !length(ti) || length(oi) != length(ti)) next
        qrow <- if (!is.null(inp$mapq) && length(inp$mapq) >= i) inp$mapq[i] else i
        qhat <- if (is.null(dim(qest))) as.numeric(qest[2]) else if (qrow %in% seq_len(nrow(qest))) qest[qrow, 2] else NA_real_
        if (!is.finite(qhat) || qhat <= 0) next
        out[[length(out) + 1L]] <- data.frame(time = ti, obs = oi / qhat, idx = i)
      }
    }
  }
  out
}

# Management interval vertical lines for B/Bmsy & biomass panels
# Returns a list of ggplot2 layers
add_management_vlines_BF_good <- function(rep,
                                          color = "grey30",
                                          linetype = "dashed",
                                          linewidth = 0.6,
                                          lineend = "butt") {
  mi <- .sfn(rep$inp, "maninterval", numeric(0))
  if (!length(mi) || !is.finite(mi[1])) return(list())
  left  <- mi[1]
  right <- if (length(mi) >= 2L && is.finite(mi[2])) mi[2] else NA_real_
  out <- list(ggplot2::geom_vline(xintercept = left,  color = color, linetype = linetype,
                                  linewidth = linewidth, lineend = lineend))
  if (is.finite(right)) {
    out[[length(out) + 1L]] <- ggplot2::geom_vline(xintercept = right, color = color, linetype = linetype,
                                                   linewidth = linewidth, lineend = lineend)
  }
  out
}

# Compact theme (kept consistent with your other panels)
.theme_minimal_compact2 <- function(base_size = 10, base_family = "") {
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

# Safe Bmsy(t) band builder with fallbacks (uses logBmsy)
# Returns list(band_df, line_df, Bmsy_est)
.build_Bmsy_band <- function(rep_obj, CI = 0.95) {
  # Try to use get.manmax()/get.msyvec() if available; otherwise fallback
  repmax <- if (exists("get.manmax", mode = "function")) try(get.manmax(rep_obj), silent = TRUE) else rep_obj
  if (inherits(repmax, "try-error") || is.null(repmax$inp)) repmax <- rep_obj

  Bmsy_all <- try(elu2::get.par("logBmsy", repmax, exp = TRUE, CI = CI), silent = TRUE)
  tt       <- .sfn(repmax$inp, "time", numeric(0))
  if (inherits(Bmsy_all, "try-error") || !length(tt)) {
    return(list(
      band_df  = data.frame(time = numeric(0), ymin = numeric(0), ymax = numeric(0)),
      line_df  = data.frame(time = numeric(0), y = numeric(0)),
      Bmsy_est = NA_real_
    ))
  }

  # get.msyvec() fallback: constant band from the last Bmsy estimate if needed
  msyvec_fun <- if (exists("get.msyvec", mode = "function")) get.msyvec else NULL

  if (is.function(msyvec_fun)) {
    mv <- try(msyvec_fun(repmax$inp, Bmsy_all), silent = TRUE)
    if (!inherits(mv, "try-error") && all(c("ll","msy","ul") %in% names(mv))) {
      band <- data.frame(time = tt, ymin = pmax(mv$ll, 0), ymax = pmax(mv$ul, 0))
      line <- data.frame(time = tt, y = pmax(mv$msy, 0))
      Bhat <- if (is.null(dim(Bmsy_all))) as.numeric(Bmsy_all[2]) else as.numeric(utils::tail(Bmsy_all[, 2], 1))
      return(list(band_df = band, line_df = line, Bmsy_est = Bhat))
    }
  }

  # Fallback band: repeat final Bmsy +/- CI across time
  Bhat <- if (is.null(dim(Bmsy_all))) as.numeric(Bmsy_all[2]) else as.numeric(utils::tail(Bmsy_all[, 2], 1))
  Bll  <- if (is.null(dim(Bmsy_all))) as.numeric(Bmsy_all[1]) else as.numeric(utils::tail(Bmsy_all[, 1], 1))
  Bul  <- if (is.null(dim(Bmsy_all))) as.numeric(Bmsy_all[3]) else as.numeric(utils::tail(Bmsy_all[, 3], 1))
  band <- data.frame(time = tt, ymin = pmax(Bll, 0), ymax = pmax(Bul, 0))
  line <- data.frame(time = tt, y = pmax(Bhat, 0))
  list(band_df = band, line_df = line, Bmsy_est = Bhat)
}

# ---------- Main function (exported) ------------------------------------------

#' Biomass panel (absolute B[t]) with Bmsy band, mapped B/Bmsy ribbon, scenarios
#'
#' @description
#' Draws an absolute-biomass panel: time-varying \eqn{B_{MSY}} band (ribbon+mean),
#' a \eqn{B/B_{MSY}} CI ribbon mapped back to biomass, optional scenario biomass
#' paths, vertical markers (management interval or last observation), and the
#' in-sample vs one-ahead styling consistent with your other panels.
#'
#' @param rep_man Fitted SPiCT report object (`spictcls`), optionally with `$man`.
#' @param scenario_color Named vector of colours for scenarios (optional).
#' @param show_CIs Logical; draw CI ribbons/edges. Default `TRUE`.
#' @param CI Confidence level in (0, 1). Default `0.95`.
#' @param show_legend Logical; show scenario legend if scenarios exist. Default `TRUE`.
#'
#' @return A `ggplot` object.
#' @export
my_plot_manage_biomass_panel <- function(rep_man,
                                         scenario_color = NULL,
                                         show_CIs = TRUE,
                                         CI = 0.95,
                                         show_legend = TRUE) {
  stopifnot(inherits(rep_man, "spictcls"))

  # Canonical scenario order & labels
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

  # Scenario colors
  if (is.null(scenario_color)) {
    cols <- if (length(scs_final)) { x <- man_cols(length(scs_final)); names(x) <- scs_final; x } else NULL
  } else {
    cols <- scenario_color
    if (length(scs_final)) {
      if (is.null(names(cols))) names(cols) <- scs_final
      cols <- cols[match(scs_final, names(cols))]; names(cols) <- scs_final
    }
  }
  lab_final <- if (length(scs_final)) ifelse(scs_final %in% names(label_map),
                                             unname(label_map[scs_final]), scs_final) else character(0)

  # Base B[t] and B/Bmsy
  Best <- elu2::get.par("logB",     rep_man, exp = TRUE, CI = CI)
  BB   <- elu2::get.par("logBBmsy", rep_man, exp = TRUE, CI = CI)

  t_B <- .spict_time_from_par(rep_man, Best)
  df_B <- data.frame(time = t_B, lwr = Best[,1], est = Best[,2], upr = Best[,3])

  inp    <- rep_man$inp
  indest <- .sfi(inp, "indest", integer(0))
  indpred<- .sfi(inp, "indpred", integer(0))

  # In-sample ends at last element of 'indest'; prediction uses 'indpred'
  ind_in <- if (length(indest)) indest[-length(indest)] else integer(0)
  ind_pr <- if (!length(scs_final)) indpred else integer(0)

  df_B_in <- if (length(ind_in)) df_B[ind_in, , drop = FALSE] else df_B[0, ]
  df_B_pr <- if (length(ind_pr)) df_B[ind_pr, , drop = FALSE] else df_B[0, ]

  # Bmsy ribbon/line (with robust fallback)
  ms <- .build_Bmsy_band(rep_man, CI = CI)
  df_Bmsy_band <- ms$band_df
  df_Bmsy_line <- ms$line_df
  Bmsy_est     <- ms$Bmsy_est

  # B/Bmsy-on-B ribbon
  t_BB <- .spict_time_from_par(rep_man, BB)
  df_BB_rib <- if (is.finite(Bmsy_est)) {
    .spict_make_ribbon_df(time = t_BB, lwr = BB[,1] * Bmsy_est, upr = BB[,3] * Bmsy_est)
  } else {
    data.frame(time = numeric(0), ymin = numeric(0), ymax = numeric(0))
  }

  # Scenario biomass paths (managed case) — drawn before base mean
  df_scen_B <- NULL
  if (length(scs_final)) {
    lst <- vector("list", length(scs_final))
    for (i in seq_along(scs_final)) {
      sc <- scs_final[i]
      rp <- rep_man$man[[sc]]
      Bi <- elu2::get.par("logB", rp, exp = TRUE, CI = CI)
      ti <- .spict_time_from_par(rp, Bi)
      lst[[i]] <- data.frame(time = ti, est = Bi[,2], scenario = sc)
    }
    df_scen_B <- do.call(rbind, lst)
    df_scen_B$scenario <- factor(df_scen_B$scenario, levels = scs_final)
  }

  # Colors
  spict_blue_mean    <- "#0000FF"
  spict_blue_ci_line <- grDevices::rgb(0, 0, 1, 0.20)
  spict_blue_ci_fill <- grDevices::rgb(0, 0, 1, 0.10)

  # ----- Build plot -----
  p <- ggplot2::ggplot()

  # Bmsy band + mean
  if (isTRUE(show_CIs) && nrow(df_Bmsy_band)) {
    p <- p + ggplot2::geom_ribbon(
      data = df_Bmsy_band,
      ggplot2::aes(x = time, ymin = ymin, ymax = ymax),
      fill = "lightgray", colour = NA
    )
  }
  if (nrow(df_Bmsy_line)) {
    p <- p + ggplot2::geom_line(
      data = df_Bmsy_line, ggplot2::aes(x = time, y = y),
      color = "black", linewidth = 0.7
    )
  }

  # B/Bmsy-on-B ribbon + edges
  if (isTRUE(show_CIs) && nrow(df_BB_rib)) {
    p <- p +
      ggplot2::geom_ribbon(
        data = df_BB_rib,
        ggplot2::aes(x = time, ymin = ymin, ymax = ymax),
        fill = spict_blue_ci_fill, colour = NA
      ) +
      ggplot2::geom_line(
        data = transform(df_BB_rib, y = ymin),
        ggplot2::aes(x = time, y = y),
        color = spict_blue_ci_line, linewidth = 0.6, linetype = "solid"
      ) +
      ggplot2::geom_line(
        data = transform(df_BB_rib, y = ymax),
        ggplot2::aes(x = time, y = y),
        color = spict_blue_ci_line, linewidth = 0.6, linetype = "solid"
      )
  }

  # V-lines (management interval for managed; last obs otherwise)
  if (length(scs_final)) {
    p <- p + add_management_vlines_BF_good(rep_man, color = "grey75",
                                           linetype = "solid", linewidth = 0.4, lineend = "butt")
  } else {
    obs_end <- .spict_obs_end_overall(rep_man)
    if (is.finite(obs_end)) {
      p <- p + ggplot2::geom_vline(xintercept = obs_end, color = "grey75",
                                   linetype = "solid", linewidth = 0.4)
    }
  }

  # Scenario biomass (under the blue mean so mean stays visible)
  if (!is.null(df_scen_B) && nrow(df_scen_B)) {
    p <- p + ggplot2::geom_line(
      data = df_scen_B, ggplot2::aes(x = time, y = est, color = scenario),
      linewidth = 0.8, na.rm = TRUE
    )
  }

  # In-sample Bt: dashed CI edges + solid mean in blue
  if (nrow(df_B_in)) {
    if (isTRUE(show_CIs)) {
      p <- p +
        ggplot2::geom_line(
          data = df_B_in, ggplot2::aes(x = time, y = lwr),
          linetype = "dashed", linewidth = 0.6,
          colour = spict_blue_mean, inherit.aes = FALSE
        ) +
        ggplot2::geom_line(
          data = df_B_in, ggplot2::aes(x = time, y = upr),
          linetype = "dashed", linewidth = 0.6,
          colour = spict_blue_mean, inherit.aes = FALSE
        )
    }
    p <- p + ggplot2::geom_line(
      data = df_B_in, ggplot2::aes(x = time, y = est),
      linewidth = 0.9, colour = spict_blue_mean, inherit.aes = FALSE
    )
  }

  # One-ahead (raw fit only): dotted mean + dashed CI
  if (!length(scs_final) && nrow(df_B_pr)) {
    if (isTRUE(show_CIs)) {
      p <- p +
        ggplot2::geom_line(
          data = df_B_pr, ggplot2::aes(x = time, y = lwr),
          linetype = "dashed", linewidth = 0.6, colour = spict_blue_mean
        ) +
        ggplot2::geom_line(
          data = df_B_pr, ggplot2::aes(x = time, y = upr),
          linetype = "dashed", linewidth = 0.6, colour = spict_blue_mean
        )
    }
    p <- p + ggplot2::geom_line(
      data = df_B_pr, ggplot2::aes(x = time, y = est),
      linetype = "dotted", linewidth = 0.8, colour = spict_blue_mean
    )
  }

  # Optional q-scaled index points
  idx <- .spict_index_points(rep_man)
  if (length(idx) >= 1) p <- p + ggplot2::geom_point(
    data = idx[[1]], ggplot2::aes(x = time, y = obs),
    color = "blue", shape = 16, size = 2, inherit.aes = FALSE
  )
  if (length(idx) >= 2) p <- p + ggplot2::geom_point(
    data = idx[[2]], ggplot2::aes(x = time, y = obs),
    shape = 22, color = "black", fill = "green", size = 2, stroke = 0.5, inherit.aes = FALSE
  )

  # Labels, scaling, legend, theme
  p +
    ggplot2::labs(title = "Absolute biomass", x = "Year", y = expression(B[t])) +
    { if (length(scs_final)) ggplot2::scale_color_manual(
      values = cols, breaks = scs_final, labels = lab_final, name = NULL,
      guide = if (show_legend) "legend" else "none"
    ) else ggplot2::scale_color_discrete(guide = "none") } +
    ggplot2::scale_fill_manual(values = if (is.null(cols)) NA else cols, guide = "none") +
    ggplot2::scale_y_continuous(
      sec.axis = if (is.finite(Bmsy_est)) ggplot2::sec_axis(~ . / Bmsy_est, name = expression(B[t]/B[MSY]))
      else ggplot2::waiver()
    ) +
    .theme_minimal_compact2()
}
