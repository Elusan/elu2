#' Kobe panel (complete, ggplot2) — OSA on, MSY ellipse, scenarios, E(B∞) mini-legend
#'
#' @description
#' Draws a Kobe-style phase plot for a fitted SPiCT model. Supports:
#' \itemize{
#'   \item Relative or absolute axes (B\eqn{_t}/B\eqn{_{MSY}} vs F\eqn{_t}/F\eqn{_{MSY}} or B\eqn{_t} vs F\eqn{_t})
#'   \item MSY uncertainty ellipse (when available)
#'   \item Historical trajectory, predicted segment, and E(B\eqn{\infty}) marker
#'   \item Optional management scenario trajectories and legend
#'   \item OSA/time-varying growth handling (forces relative axes)
#' }
#'
#' @param model Fitted SPiCT object (`spictcls`) with `inp$reportmode == 0`.
#' @param logax Logical; use log10 scales on both axes. Default `FALSE`.
#' @param plot.legend Logical; show the small top-right legend for **E(B\[infinity])**. Default `TRUE`.
#' @param man.legend Logical; if scenarios exist in `model$man`, show their legend. Default `TRUE`.
#' @param ext Logical; add secondary relative axes when `rel.axes = FALSE`. Default `TRUE`.
#' @param rel.axes Logical; plot as B\eqn{_t}/B\eqn{_{MSY}} vs F\eqn{_t}/F\eqn{_{MSY}}. Forces `ext = FALSE`. Default `FALSE`.
#' @param xlim,ylim Optional numeric length-2 axis limits (computed if `NULL`).
#' @param labpos Reserved; integer length-2 offsets for first/last-year labels. Default `c(1, 1)`.
#' @param xlabel Optional x-axis label override.
#' @param stamp Optional small text stamped at lower-right (e.g., version tag).
#' @param verbose Logical; emit extraction warnings. Default `TRUE`.
#' @param CI Confidence level for parameter extraction (0,1). Default `0.95`.
#' @param print_it Logical; if `TRUE`, print the plot and return it. Default `FALSE`.
#'
#' @return A `ggplot` object.
#'
#' @details
#' This function avoids depending on SPiCT helper functions by inlining
#' minimal equivalents (see internal helpers below). If the package
#' \pkg{ellipse} is installed, an MSY ellipse is drawn; otherwise a point.
#'
#' @examples
#' \dontrun{
#'   rep <- fit.spict(inp)
#'   p <- plot_elu2_panel_kobe(rep, rel.axes = TRUE)
#'   print(p)
#' }
#'
#' @import ggplot2
#' @importFrom scales pretty_breaks
#' @importFrom stats qnorm cov2cor
#' @importFrom grDevices adjustcolor
#' @export
plot_elu2_panel_kobe <- function(model,
                                 logax       = FALSE,
                                 plot.legend = TRUE,
                                 man.legend  = TRUE,
                                 ext         = TRUE,
                                 rel.axes    = FALSE,
                                 xlim        = NULL,
                                 ylim        = NULL,
                                 labpos      = c(1, 1),
                                 xlabel      = NULL,
                                 stamp       = NULL,
                                 verbose     = TRUE,
                                 CI          = 0.95,
                                 print_it    = FALSE) {

  # ----------------------- Internal helpers (file-local) ---------------------- #

  #' @keywords internal
  #' @noRd
  .elu2_check_rep <- function(rep, reportmode0 = TRUE) {
    if (!inherits(rep, "spictcls") || !"opt" %in% names(rep)) {
      stop("The argument 'model' must be a fitted spict object (fit.spict()).")
    } else if (reportmode0 && rep$inp$reportmode != 0) {
      stop("All states must be reported! Set 'inp$reportmode <- 0' and refit.")
    }
    TRUE
  }

  #' Extract parameter (random/fixed/SDreport/optim fallback)
  #' @keywords internal
  #' @noRd
  .elu2_get_par <- function(parname, rep, exp = FALSE, CI = 0.95) {
    if (CI <= 0 || CI >= 1) stop("CI must be in (0,1).")
    z <- stats::qnorm(CI + (1 - CI) / 2)

    ind_ran <- which(names(rep$par.random) == parname)
    ind_fix <- which(names(rep$par.fixed)  == parname)
    ind_sdr <- which(names(rep$value)      == parname)
    ind_opt <- which(names(rep$opt$par)    == parname)

    est <- ll <- ul <- sdv <- NULL

    if (length(ind_ran)) {
      est <- rep$par.random[ind_ran]
      sdv <- sqrt(rep$diag.cov.random[ind_ran])
      ll  <- est - z * sdv; ul <- est + z * sdv
    }
    if (length(ind_fix)) {
      est <- rep$par.fixed[ind_fix]
      cf  <- rep$cov.fixed
      sdv <- if (is.null(cf)) rep(NA_real_, length(est)) else sqrt(diag(cf))[ind_fix]
      ll  <- est - z * sdv; ul <- est + z * sdv
    }
    if (length(ind_sdr)) {
      est <- rep$value[ind_sdr]
      sdv <- rep$sd[ind_sdr]
      ll  <- est - z * sdv; ul <- est + z * sdv
    }

    if (length(est) == 0) {
      if (length(ind_opt)) {
        est <- rep$opt$par[ind_opt]
        sdv <- ll <- ul <- rep(NA_real_, length(est))
      } else if ("phases" %in% names(rep$inp) &&
                 parname %in% names(rep$inp$phases) &&
                 rep$inp$phases[[parname]] == -1) {
        est <- rep$inp$parlist[[parname]]
        sdv <- rep(0, length(est)); ll <- est; ul <- est
      } else if (!is.na(parname) && identical(parname, "P")) {
        B <- .elu2_get_par("logB",     rep, exp = TRUE, CI = CI)
        C <- .elu2_get_par("logCpred", rep, exp = TRUE, CI = CI)
        ic <- rep$inp$ic; nc <- rep$inp$nc
        B0 <- B[ic, 2]; B1 <- B[ic + nc, 2]
        T0 <- rep$inp$time[ic]; T1 <- rep$inp$time[ic + nc]
        est <- (B1 - B0 + C[, 2]) / (T1 - T0)
        sdv <- ll <- ul <- rep(NA_real_, length(est))
      } else {
        if (verbose) warning("get_par: could not extract '", parname, "'. Returning NA.")
        est <- sdv <- ll <- ul <- NA_real_
      }
    }

    n <- length(est); if (length(sdv) == 0) sdv <- rep(NA_real_, n)
    if (length(ll) == 0) ll <- rep(NA_real_, n)
    if (length(ul) == 0) ul <- rep(NA_real_, n)

    if (isTRUE(exp)) {
      cv <- ifelse(is.finite(sdv), sqrt(exp(sdv^2) - 1), NA_real_)
      ll <- exp(ll); ul <- exp(ul); est <- exp(est)
      ul[is.infinite(ul)] <- exp(705)
    } else {
      cv <- sdv / est
    }
    out <- cbind(ll, est, ul, sdv, cv)
    if (parname %in% c("logB","logF","logBBmsy","logFFmsy"))
      rownames(out) <- rep$inp$time
    out
  }

  #' @keywords internal
  #' @noRd
  .elu2_add_catchunit <- function(lab, cu) {
    cu <- as.character(cu)
    if (nzchar(cu)) eval(bquote(.(lab[[1]]) * ',' ~ .(cu))) else lab
  }

  #' MSY ellipse in (logBmsy, logFmsy) plane
  #' @keywords internal
  #' @noRd
  .elu2_make_rp_ellipse <- function(rep) {
    inds <- c(max(which(names(rep$value) == "logBmsy")),
              max(which(names(rep$value) == "logFmsy")))
    sds  <- rep$sd[inds]
    s_sum <- rep$sd[which(names(rep$value) == "logBmsyPluslogFmsy")]
    cova <- (s_sum^2 - sds[1]^2 - sds[2]^2) / 2
    covBF <- matrix(c(sds[1]^2, cova, cova, sds[2]^2), 2, 2, byrow = TRUE)
    parBF <- rep$value[inds]
    if (requireNamespace("ellipse", quietly = TRUE)) {
      ellipse::ellipse(stats::cov2cor(covBF)[1, 2],
                       scale  = sqrt(diag(covBF)),
                       centre = parBF, npoints = 300)
    } else matrix(c(parBF[1], parBF[2]), ncol = 2,
                  dimnames = list(NULL, c("x","y")))
  }

  #' Annual mean helper (for sub-annual dt)
  #' @keywords internal
  #' @noRd
  .elu2_annual_avg <- function(intime, vec, type = "mean") {
    fun <- match.fun(type)
    anntime   <- unique(floor(intime))
    floortime <- floor(intime)
    nstepvec  <- vapply(anntime, function(a) sum(a == floortime), numeric(1))
    anntime   <- anntime[which(nstepvec == max(nstepvec))]
    annvec    <- vapply(anntime, function(a) fun(vec[which(a == floortime)]), numeric(1))
    list(anntime = anntime, annvec = annvec)
  }

  #' E(B_inf) for given last F and uncertainty correction
  #' @keywords internal
  #' @noRd
  .elu2_calc_EBinf <- function(K, n, Fl, Fmsy, sdb2) {
    base <- 1 - (n - 1)/n * (Fl/Fmsy)
    base <- max(0, base)
    corr <- 1 - n/2 / (1 - (1 - n*Fmsy + (n - 1)*Fl))
    EB   <- K * (base)^(1/(n - 1)) * (1 - corr*sdb2)
    max(0, EB)
  }

  #' Wrapper to compute E(B_inf) from report
  #' @keywords internal
  #' @noRd
  .elu2_get_EBinf <- function(rep, CI = 0.95) {
    K    <- .elu2_get_par("logK", rep, exp = TRUE, CI = CI)[2]
    n    <- .elu2_get_par("logn", rep, exp = TRUE, CI = CI)[2]
    sdb2 <- .elu2_get_par("logsdb", rep, exp = TRUE, CI = CI)[2]^2
    Fmsy <- tail(.elu2_get_par("logFmsy", rep, exp = TRUE, CI = CI), 1)[2]
    logFs <- .elu2_get_par("logFs", rep, CI = CI)
    if (min(rep$inp$dtc) < 1) {
      alf <- .elu2_annual_avg(rep$inp$time, logFs[, 2])
      fff <- exp(alf$annvec)
    } else fff <- exp(logFs[, 2])
    Fl <- tail(unname(fff), 1)
    .elu2_calc_EBinf(K, n, Fl, Fmsy, sdb2)
  }

  #' Scenario color set
  #' @keywords internal
  #' @noRd
  .elu2_man_cols <- function() {
    colvec <- c('darkmagenta','cyan3','darkgreen','coral1','black','magenta',
                'gold','green','cadetblue3','chocolate3','darkolivegreen3',
                'cyan','darkred')
    rep(colvec, 3)
  }

  #' Format helpers
  #' @keywords internal
  #' @noRd
  .elu2_fmt1 <- function(x) { x <- as.numeric(x); ifelse(is.finite(x),
                                                         formatC(x, format = "f", digits = 1),
                                                         NA_character_) }
  #' @keywords internal
  #' @noRd
  .elu2_lab_rel_1 <- function(vals) {
    out <- .elu2_fmt1(vals); out[abs(vals - 1) < 1e-9] <- "1"; out
  }

  #' Clamp helper
  #' @keywords internal
  #' @noRd
  .elu2_clamp <- function(val, lo, hi) pmin(pmax(val, lo), hi)

  # ------------------ Validate & derive quantities --------------- #
  rep <- model
  .elu2_check_rep(rep)

  tvgflag <- rep$inp$timevaryinggrowth | rep$inp$logmcovflag
  if (tvgflag) rel.axes <- TRUE

  Bmsy_all <- .elu2_get_par("logBmsy", rep, exp = TRUE, CI = CI)
  Fmsy_all <- .elu2_get_par("logFmsy", rep, exp = TRUE, CI = CI)
  Bmsy <- tail(Bmsy_all, 1L); Fmsy <- tail(Fmsy_all, 1L)

  if (rel.axes) {
    ext <- FALSE
    bscal <- Bmsy[2]; fscal <- Fmsy[2]
    xlab_expr <- expression(B[t]/B[MSY]); ylab_expr <- expression(F[t]/F[MSY])
  } else {
    bscal <- 1; fscal <- 1
    xlab_expr <- .elu2_add_catchunit(expression(B[t]), rep$inp$catchunit)
    ylab_expr <- expression(F[t])
  }

  Best   <- .elu2_get_par("logB", rep, exp = TRUE, CI = CI)
  logBest <- .elu2_get_par("logB", rep, CI = CI)
  if (tvgflag) {
    Fest <- .elu2_get_par("logFFmsy", rep, exp = TRUE, CI = CI); fscal <- 1; Fmsy <- c(1,1)
  } else {
    Fest <- .elu2_get_par("logFs",     rep, exp = TRUE, CI = CI)
  }
  logFest <- .elu2_get_par("logFs", rep, CI = CI)

  cl <- try(.elu2_make_rp_ellipse(rep), silent = TRUE)
  if (inherits(cl,"try-error")) cl <- matrix(c(log(Bmsy[2]), log(Fmsy[2])), ncol = 2)

  if (min(rep$inp$dtc) < 1) {
    alb <- .elu2_annual_avg(rep$inp$time, logBest[,2])
    alf <- .elu2_annual_avg(rep$inp$time, logFest[,2])
    bbb <- exp(alb$annvec) / bscal
    fff <- exp(alf$annvec) / fscal
    fbtime <- alb$anntime
  } else {
    bbb <- Best[rep$inp$indest, 2] / bscal
    fff <- Fest[rep$inp$indest, 2] / fscal
    fbtime <- rep$inp$time[rep$inp$indest]
  }

  Fl <- tail(unname(fff), 1)
  Bl <- tail(unname(bbb), 1)
  EBinf <- .elu2_get_EBinf(rep, CI = CI) / bscal

  if (is.null(xlim)) {
    xlim <- range(c(exp(cl[,1]), Best[,2], EBinf) / bscal, na.rm = TRUE)
    if (min(rep$inp$dtc) < 1)
      xlim <- range(c(exp(alb$annvec), exp(cl[,1]), EBinf) / bscal, na.rm = TRUE)
    xlim[2] <- min(c(xlim[2], 8*Bmsy[2]/bscal), 2.2*max(bbb), na.rm = TRUE)
    xlim[2] <- max(c(xlim[2], Bmsy[2]/bscal), na.rm = TRUE)
  }
  if (is.null(ylim)) {
    ylim <- range(c(exp(cl[,2]), Fest[,2]) / fscal, na.rm = TRUE)
    if (min(rep$inp$dtc) < 1)
      ylim <- range(c(exp(logFest[,2]) / fscal, exp(cl[,2]) / fscal), na.rm = TRUE)
    ylim[2] <- min(c(ylim[2], 8*Fmsy[2]/fscal), 2.2*max(fff), na.rm = TRUE)
    ylim[2] <- max(c(ylim[2], Fmsy[2]/fscal), na.rm = TRUE)
    if ("man" %in% names(rep)) ylim <- range(ylim, 0)
  }

  logminval <- 1e-4
  if (isTRUE(logax)) {
    if (xlim[1] < logminval) xlim[1] <- logminval
    if (ylim[1] < logminval) ylim[1] <- logminval
  }

  # Slight padding so fills/paths don't hug the frame
  pad_frac <- 0.02; xpad <- pad_frac * diff(xlim); ypad <- pad_frac * diff(ylim)
  xlim <- c(max(if (isTRUE(logax)) logminval else 0, xlim[1] - xpad), xlim[2] + xpad)
  ylim <- c(max(if (isTRUE(logax)) logminval else 0, ylim[1] - ypad), ylim[2] + ypad)

  # ------------------------- Data frames for ggplot ------------------------- #
  xminq <- xlim[1]; xmaxq <- xlim[2]; yminq <- ylim[1]; ymaxq <- ylim[2]
  bx <- if (rel.axes) 1 else Bmsy[2]; bx <- .elu2_clamp(bx, xminq, xmaxq)
  fy <- if (rel.axes) 1 else Fmsy[2]; fy <- .elu2_clamp(fy, yminq, ymaxq)

  # Quadrants (fill panel fully)
  df_green_rect   <- data.frame(xmin = bx,   xmax =  Inf, ymin = -Inf, ymax =  fy)
  df_yellowL_rect <- data.frame(xmin = -Inf, xmax =  bx,  ymin = -Inf, ymax =  fy)
  df_yellowT_rect <- data.frame(xmin = bx,   xmax =  Inf, ymin =  fy,  ymax =  Inf)
  df_red_rect     <- data.frame(xmin = -Inf, xmax =  bx,  ymin =  fy,  ymax =  Inf)

  df_ellipse <- if (ncol(as.matrix(cl)) == 2 && nrow(as.matrix(cl)) > 1) {
    data.frame(x = exp(cl[,1]) / (if (rel.axes) Bmsy[2] else 1),
               y = exp(cl[,2]) / (if (rel.axes) Fmsy[2] else 1))
  } else {
    data.frame(x = exp(cl[1]) / (if (rel.axes) Bmsy[2] else 1),
               y = exp(cl[2]) / (if (rel.axes) Fmsy[2] else 1))
  }

  epsx <- 0.005*diff(xlim); epsy <- 0.005*diff(ylim)
  clamp_in <- function(v, lo, hi, eps) pmin(pmax(v, lo + eps), hi - eps)

  n_hist <- length(bbb)
  df_traj <- data.frame(x = clamp_in(bbb, xminq, xmaxq, epsx),
                        y = clamp_in(fff, yminq, ymaxq, epsy),
                        ord = seq_len(n_hist))

  df_pred <- NULL; df_EBseg <- NULL
  if (!(min(rep$inp$dtc) < 1) && !("man" %in% names(rep))) {
    xpred <- .elu2_get_par("logB", rep, exp = TRUE)[rep$inp$indpred, 2] / (if (rel.axes) Bmsy[2] else 1)
    ypred <- .elu2_get_par("logFs", rep, exp = TRUE)[rep$inp$indpred, 2] / (if (rel.axes) Fmsy[2] else 1)
    df_pred <- data.frame(x = xpred, y = ypred)
    Bll <- tail(xpred, 1); Fll <- tail(ypred, 1)
    EBinf <- .elu2_get_EBinf(rep, CI = CI) / (if (rel.axes) Bmsy[2] else 1)
    df_EBseg <- data.frame(x = Bll, xend = EBinf, y = Fll, yend = Fll)
  }

  df_first <- data.frame(x = bbb[1], y = fff[1], lab = round(fbtime[1], 2))
  df_last  <- data.frame(x = tail(bbb,1), y = tail(fff,1), lab = round(tail(fbtime, 1), 2))

  nr <- length(rep$inp$ini$logr)
  df_msy_prev <- df_msy_curr <- NULL
  if (nr > 1) {
    df_msy_prev <- data.frame(x = Bmsy_all[1:(nr-1),2] / (if (rel.axes) Bmsy[2] else 1),
                              y = Fmsy_all[1:(nr-1),2] / (if (rel.axes) Fmsy[2] else 1))
    df_msy_curr <- data.frame(x = Bmsy_all[nr,2] / (if (rel.axes) Bmsy[2] else 1),
                              y = Fmsy_all[nr,2] / (if (rel.axes) Fmsy[2] else 1))
  }

  df_true <- NULL
  if ("true" %in% names(rep$inp)) {
    df_true <- data.frame(x = rep$inp$true$Bmsy / (if (rel.axes) Bmsy[2] else 1),
                          y = rep$inp$true$Fmsy / (if (rel.axes) Fmsy[2] else 1))
  }

  df_man <- df_man_int <- NULL
  leg_man <- NULL
  if ("man" %in% names(rep) && length(rep$man)) {
    nman <- length(rep$man)
    leg_man <- names(rep$man); if (is.null(leg_man) || !length(leg_man)) leg_man <- paste0("Scenario ", seq_len(nman))
    df_list <- vector("list", nman); df_int_list <- vector("list", nman)
    for (i in seq_len(nman)) {
      rp <- rep$man[[i]]
      estF <- .elu2_get_par("logF", rp, exp = TRUE)[,2]
      estB <- .elu2_get_par("logB", rp, exp = TRUE)[,2]
      time <- rp$inp$time
      manint <- rp$inp$maninterval
      indmanstart <- which(time >= manint[1])
      if (length(indmanstart)) {
        maninds <- indmanstart[1]:tail(indmanstart, 1)
        df_list[[i]] <- data.frame(x = estB[maninds] / (if (rel.axes) Bmsy[2] else 1),
                                   y = estF[maninds] / (if (rel.axes) Fmsy[2] else 1),
                                   scenario = leg_man[i])
      } else df_list[[i]] <- data.frame(x = numeric(0), y = numeric(0), scenario = character(0))
      lastobs <- rp$inp$timerangeObs[2]
      ind <- which(time >= lastobs & time < manint[1])
      if (length(ind)) {
        df_int_list[[i]] <- data.frame(x = estB[ind] / (if (rel.axes) Bmsy[2] else 1),
                                       y = estF[ind] / (if (rel.axes) Fmsy[2] else 1),
                                       scenario = leg_man[i])
      } else df_int_list[[i]] <- data.frame(x = numeric(0), y = numeric(0), scenario = character(0))
    }
    df_man <- do.call(rbind, df_list)
    df_man_int <- do.call(rbind, df_int_list)
  }

  # --------------------------- Build ggplot ------------------------------ #
  p <- ggplot2::ggplot()
  p <- p + ggplot2::geom_rect(data = df_green_rect,
                              ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                              fill = grDevices::adjustcolor("darkseagreen3", alpha.f = 1), color = NA)
  p <- p + ggplot2::geom_rect(data = df_yellowL_rect,
                              ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                              fill = grDevices::adjustcolor("khaki1", alpha.f = 1), color = NA)
  p <- p + ggplot2::geom_rect(data = df_yellowT_rect,
                              ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                              fill = grDevices::adjustcolor("khaki1", alpha.f = 1), color = NA)
  p <- p + ggplot2::geom_rect(data = df_red_rect,
                              ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                              fill = grDevices::adjustcolor("indianred1", alpha.f = 1), color = NA)

  if (nrow(df_ellipse) > 1) {
    p <- p + ggplot2::geom_polygon(data = df_ellipse, ggplot2::aes(x = x, y = y),
                                   fill = grDevices::adjustcolor("gray80", alpha.f = 0.8),
                                   color = grDevices::adjustcolor("gray80", alpha.f = 0.8),
                                   linewidth = 0.3)
  } else {
    p <- p + ggplot2::geom_point(data = df_ellipse, ggplot2::aes(x = x, y = y),
                                 color = "gray40", size = 2)
  }

  p <- p + ggplot2::geom_path(data = df_traj, ggplot2::aes(x = x, y = y),
                              linewidth = 0.4, color = grDevices::adjustcolor("blue", alpha.f = 0.8))

  if (!is.null(df_pred)) {
    p <- p + ggplot2::geom_path(data = df_pred, ggplot2::aes(x = x, y = y),
                                linetype = "11", color = grDevices::adjustcolor("blue", alpha.f = 0.8))
  }
  if (!is.null(df_EBseg)) {
    p <- p + ggplot2::geom_segment(data = df_EBseg,
                                   ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
                                   color = "blue", linetype = "11", linewidth = 0.7) +
      ggplot2::geom_point(data = data.frame(x = df_EBseg$xend, y = df_EBseg$yend),
                          ggplot2::aes(x = x, y = y), shape = 23, fill = "gold",
                          color = "black", size = 2, stroke = 0.5)
  }

  if (!is.null(df_man) && nrow(df_man)) {
    p <- p + ggplot2::geom_path(data = df_man, ggplot2::aes(x = x, y = y, color = scenario),
                                linewidth = 0.7, inherit.aes = FALSE, show.legend = man.legend)
  }
  if (!is.null(df_man_int) && nrow(df_man_int)) {
    p <- p + ggplot2::geom_path(data = df_man_int, ggplot2::aes(x = x, y = y, color = scenario),
                                linewidth = 0.7, linetype = "33", inherit.aes = FALSE, show.legend = man.legend)
  }

  if (!is.null(df_msy_prev) && nrow(df_msy_prev))
    p <- p + ggplot2::geom_point(data = df_msy_prev, ggplot2::aes(x = x, y = y),
                                 shape = 3, color = "magenta", size = 2)
  if (!is.null(df_msy_curr) && nrow(df_msy_curr))
    p <- p + ggplot2::geom_point(data = df_msy_curr, ggplot2::aes(x = x, y = y),
                                 shape = 3, color = "black", size = 2)

  if (!is.null(df_true)) {
    p <- p + ggplot2::geom_point(data = df_true, ggplot2::aes(x = x, y = y),
                                 shape = 25, fill = grDevices::adjustcolor(rgb(1, 165/255, 0), alpha.f = 0.7),
                                 color = "black", size = 3)
  }

  p <- p + ggplot2::geom_point(data = df_first, ggplot2::aes(x = x, y = y),
                               shape = 21, fill = "white", color = "black", size = 2, stroke = 0.5) +
    ggplot2::geom_text(data = df_first, ggplot2::aes(x = x, y = y, label = lab),
                       vjust = -0.8, size = 3) +
    ggplot2::geom_point(data = df_last, ggplot2::aes(x = x, y = y),
                        shape = 22, fill = "white", color = "black", size = 2, stroke = 0.5) +
    ggplot2::geom_text(data = df_last, ggplot2::aes(x = x, y = y, label = lab),
                       vjust = -0.8, size = 3)

  if (!isTRUE(logax) && (0 >= xlim[1]) && (0 <= xlim[2])) {
    p <- p + ggplot2::geom_vline(xintercept = 0, color = "darkred", linetype = 2)
  }

  # Labels & scales (with optional secondary relative axes)
  xlab_final <- if (is.null(xlabel)) xlab_expr else xlabel
  ylab_final <- ylab_expr
  if (isTRUE(logax)) {
    p <- p + ggplot2::scale_x_log10(
      limits = xlim, name = xlab_final,
      labels = function(x) formatC(x, format = "f", digits = 0), expand = c(0, 0),
      sec.axis = if (!rel.axes && isTRUE(ext)) ggplot2::sec_axis(
        trans  = ~ . / (Bmsy[2] / 1),
        name   = expression(B[t]/B[MSY]),
        breaks = function(lims) { br <- scales::pretty_breaks(n = 6)(lims); sort(unique(c(br, 1))) },
        labels = .elu2_lab_rel_1
      ) else ggplot2::waiver()
    ) +
      ggplot2::scale_y_log10(
        limits = ylim, name = ylab_final, labels = .elu2_fmt1, expand = c(0, 0),
        sec.axis = if (!rel.axes && isTRUE(ext)) ggplot2::sec_axis(
          trans  = ~ . / (Fmsy[2] / 1),
          name   = expression(F[t]/F[MSY]),
          breaks = function(lims) { br <- scales::pretty_breaks(n = 6)(lims); sort(unique(c(br, 1))) },
          labels = .elu2_lab_rel_1
        ) else ggplot2::waiver()
      )
  } else {
    p <- p + ggplot2::scale_x_continuous(
      limits = xlim, name = xlab_final,
      labels = function(x) formatC(x, format = "f", digits = 0), expand = c(0, 0),
      sec.axis = if (!rel.axes && isTRUE(ext)) ggplot2::sec_axis(
        trans  = ~ . / (Bmsy[2] / 1),
        name   = expression(B[t]/B[MSY]),
        breaks = function(lims) { br <- scales::pretty_breaks(n = 6)(lims); sort(unique(c(br, 1))) },
        labels = .elu2_lab_rel_1
      ) else ggplot2::waiver()
    ) +
      ggplot2::scale_y_continuous(
        limits = ylim, name = ylab_final, labels = .elu2_fmt1, expand = c(0, 0),
        sec.axis = if (!rel.axes && isTRUE(ext)) ggplot2::sec_axis(
          trans  = ~ . / (Fmsy[2] / 1),
          name   = expression(F[t]/F[MSY]),
          breaks = function(lims) { br <- scales::pretty_breaks(n = 6)(lims); sort(unique(c(br, 1))) },
          labels = .elu2_lab_rel_1
        ) else ggplot2::waiver()
      )
  }

  if (!is.null(df_man) && nrow(df_man) && man.legend) {
    uniq_sc <- unique(df_man$scenario)
    colv <- .elu2_man_cols()[seq_along(uniq_sc)]; names(colv) <- uniq_sc
    p <- p + ggplot2::scale_color_manual(values = colv, name = "Scenario")
  } else {
    p <- p + ggplot2::scale_color_discrete(guide = "none")
  }

  p <- p + ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "gray25", fill = NA, linewidth = 2),
      legend.background = ggplot2::element_rect(fill = grDevices::adjustcolor("white", alpha.f = 0.8), color = NA),
      legend.position = if (man.legend) "top" else "none",
      legend.direction = "horizontal",
      axis.text.x  = ggplot2::element_text(size = 10, face = "bold", color = "gray25"),
      axis.text.y  = ggplot2::element_text(size = 10, face = "bold", color = "gray25"),
      axis.title.x = ggplot2::element_text(size = 14, face = "bold", color = "gray25"),
      axis.title.y = ggplot2::element_text(size = 14, face = "bold", color = "gray25"),
      plot.margin = ggplot2::margin(4, 4, 4, 4)
    )

  # Mini-legend for E(B∞)
  if (plot.legend && !(min(rep$inp$dtc) < 1)) {
    px <- xlim[2] - 0.02 * diff(xlim)
    py <- ylim[2] - 0.04 * diff(ylim)
    gap <- 0.12 * diff(xlim)
    p <- p + ggplot2::annotate("point", x = px - gap, y = py, shape = 23,
                               size = 2.5, fill = "gold", color = "black", stroke = 0.6) +
      ggplot2::annotate("text", x = px, y = py,
                        label = "bold(E(B[infinity]))", parse = TRUE,
                        hjust = 1, vjust = 0.5, size = 4)
  }

  # Non-convergence badge
  if (rep$opt$convergence != 0) {
    p <- p + ggplot2::annotate("text",
                               x = xlim[1] + 0.01 * diff(xlim),
                               y = ylim[2] + 0.02 * diff(ylim),
                               label = "▲", color = "black", size = 6) +
      ggplot2::annotate("text",
                        x = xlim[1] + 0.01 * diff(xlim),
                        y = ylim[2] + 0.02 * diff(ylim),
                        label = "!", color = "black", size = 3, vjust = 0.35)
  }

  # Optional stamp
  if (!is.null(stamp) && nzchar(stamp)) {
    p <- p + ggplot2::annotate("text",
                               x = xlim[2] - 0.01 * diff(xlim),
                               y = ylim[1] - 0.04 * diff(ylim),
                               label = stamp, hjust = 1, vjust = 0, size = 3)
  }

  if (isTRUE(print_it)) print(p)
  p
}
