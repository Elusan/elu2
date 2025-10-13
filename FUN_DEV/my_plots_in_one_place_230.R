# -------------------------------------------------------------------
# Robust plotting + helpers for SPiCT (ELU2 style)
# Keeps public API and visuals as-is; adds defensive guards throughout.
# -------------------------------------------------------------------

#' SPiCT-compatible getter alias (internal)
#' Ensures internal calls to get.par(...) resolve to get_par(...).
#' @keywords internal
#' @noRd
if (!exists("get.par", mode = "function")) {
  get.par <- function(...) get_par(...)
}

# ===================================================================
# KOBE PANEL
# ===================================================================

#' Kobe plot for fitted or managed SPiCT objects (ggplot2)
#'
#' @title Multi-scenario Kobe-style phase plot (B vs F)
#'
#' @description
#' Draws a Kobe-style phase plot of biomass on the x-axis and fishing
#' mortality on the y-axis for a fitted **SPiCT** report object (`spictcls`).
#' The plot shows:
#' \itemize{
#'   \item the four management quadrants (green/yellow/red rectangles) defined by
#'         \eqn{B_{MSY}} and \eqn{F_{MSY}},
#'   \item the recent historical trajectory \eqn{(B_t, F_t)},
#'   \item optional one-step-ahead prediction and the segment to \eqn{E(B_\infty)},
#'   \item optional management-scenario trajectories (if `rep$man` is present),
#'   \item optional MSY points when growth is time-varying and the true point if available.
#' }
#'
#' Axes can be drawn either on the absolute scale or on the relative scale
#' (i.e., \eqn{B_t/B_{MSY}} and \eqn{F_t/F_{MSY}}) when `rel.axes = TRUE`.
#' A secondary axis can be shown for relative values when `ext = TRUE` and
#' `rel.axes = FALSE`.
#'
#' @param rep A fitted SPiCT report object (class `spictcls`). May include a
#'   `rep$man` list created by [spict::manage()] containing scenario fits to plot.
#' @param logax Logical; if `TRUE`, both axes are shown on logarithmic scales.
#' @param plot.legend Logical; if `TRUE`, draws a small legend label for
#'   \eqn{E(B_\infty)} in the top-right when available (not for sub-annual `dtc`).
#' @param man.legend Logical; if `TRUE`, shows a colour legend for management
#'   scenarios when present.
#' @param ext Logical; if `TRUE` and `rel.axes = FALSE`, adds secondary axes
#'   for the relative scales (\eqn{B_t/B_{MSY}}, \eqn{F_t/F_{MSY}}).
#' @param rel.axes Logical; if `TRUE`, the primary axes are the relative scales.
#' @param xlim,ylim Optional limits.
#' @param labpos Ignored (kept for compat).
#' @param xlabel Optional x label override.
#' @param stamp Optional version stamp bottom-right.
#' @param verbose Logical; enable extraction warnings.
#' @param CI Confidence level (0,1) for intervals.
#' @return Invisibly, a `ggplot` object.
#' @export
my_plot_kobe_all_management_scenario <- function(rep,
                                                 logax = FALSE,
                                                 plot.legend = TRUE,
                                                 man.legend = TRUE,
                                                 ext = TRUE,
                                                 rel.axes = FALSE,
                                                 xlim = NULL,
                                                 ylim = NULL,
                                                 labpos = c(1, 1),
                                                 xlabel = NULL,
                                                 stamp = NULL,
                                                 verbose = TRUE,
                                                 CI = 0.95) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")

  ## ----------------------- Local helpers ---------------------- ##
  check_rep <- function(rep, reportmode0 = TRUE) {
    if (!inherits(rep, "spictcls") || !"opt" %in% names(rep)) {
      stop("The argument 'rep' must be a fitted spict object (fit.spict()).")
    } else if (reportmode0 && !is.null(rep$inp$reportmode) && rep$inp$reportmode != 0) {
      stop("All states must be reported! Set 'inp$reportmode <- 0' and refit.")
    }
    TRUE
  }

  get_par <- function(parname, rep, exp = FALSE, CI = 0.95) {
    if (CI <= 0 || CI >= 1) stop("CI must be in (0,1).")
    z <- stats::qnorm(CI + (1 - CI) / 2)

    ind_ran <- which(names(rep$par.random %||% numeric()) == parname)
    ind_fix <- which(names(rep$par.fixed  %||% numeric()) == parname)
    ind_sdr <- which(names(rep$value      %||% numeric()) == parname)
    ind_opt <- which(names(rep$opt$par    %||% numeric()) == parname)

    est <- ll <- ul <- sdv <- NULL

    if (length(ind_ran)) {
      est <- rep$par.random[ind_ran]
      sdv <- sqrt(rep$diag.cov.random[ind_ran])
      ll  <- est - z * sdv; ul <- est + z * sdv
    }
    if (length(ind_fix)) {
      est <- rep$par.fixed[ind_fix]
      cf  <- rep$cov.fixed %||% matrix(NA_real_, nrow = length(rep$par.fixed), ncol = length(rep$par.fixed))
      sdv <- sqrt(diag(cf))[ind_fix]
      ll  <- est - z * (sdv %||% 0); ul <- est + z * (sdv %||% 0)
    }
    if (length(ind_sdr)) {
      est <- rep$value[ind_sdr]; sdv <- rep$sd[ind_sdr]
      ll  <- est - z * sdv;      ul  <- est + z * sdv
    }

    if (length(est) == 0) {
      if (length(ind_opt)) {
        est <- rep$opt$par[ind_opt]
        sdv <- rep(NA_real_, length(est))
        ll  <- ul <- rep(NA_real_, length(est))
      } else if (!is.null(rep$inp$phases) &&
                 parname %in% names(rep$inp$phases) &&
                 isTRUE(rep$inp$phases[[parname]] == -1)) {
        est <- rep$inp$parlist[[parname]]
        sdv <- rep(0, length(est)); ll <- est; ul <- est
      } else if (!is.na(parname) && identical(parname, "P")) {
        B <- get_par("logB",     rep, exp = TRUE, CI = CI)
        C <- get_par("logCpred", rep, exp = TRUE, CI = CI)
        ic <- rep$inp$ic; nc <- rep$inp$nc
        B0 <- B[ic, 2]; B1 <- B[ic + nc, 2]
        T0 <- rep$inp$time[ic]; T1 <- rep$inp$time[ic + nc]
        est <- (B1 - B0 + C[, 2]) / (T1 - T0)
        sdv <- ll <- ul <- rep(NA_real_, length(est))
      } else {
        if (verbose) warning("get_par: could not extract '", parname, "'. Returning NA.")
        est <- ll <- ul <- sdv <- NA_real_
      }
    }

    n <- length(est)
    if (length(sdv) == 0) sdv <- rep(NA_real_, n)
    if (length(ll)  == 0) ll  <- rep(NA_real_, n)
    if (length(ul)  == 0) ul  <- rep(NA_real_, n)

    if (isTRUE(exp)) {
      cv <- ifelse(is.finite(sdv), sqrt(exp(sdv^2) - 1), NA_real_)
      ll <- exp(ll); ul <- exp(ul); est <- exp(est)
      ul[is.infinite(ul)] <- exp(705)
    } else {
      cv <- sdv / est
    }

    out <- cbind(ll, est, ul, sdv, cv)
    if (parname %in% c("logB", "logF", "logBBmsy", "logFFmsy")) {
      rn <- suppressWarnings(as.numeric(rep$inp$time))
      if (length(rn)) rownames(out) <- rep$inp$time
    }
    out
  }

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  add_catchunit <- function(lab, cu) {
    cu <- as.character(cu %||% "")
    if (nzchar(cu)) eval(bquote(.(lab[[1]]) * ',' ~ .(cu))) else lab
  }

  make_rp_ellipse <- function(rep) {
    idxB <- which(names(rep$value %||% numeric()) == "logBmsy")
    idxF <- which(names(rep$value %||% numeric()) == "logFmsy")
    idxS <- which(names(rep$value %||% numeric()) == "logBmsyPluslogFmsy")
    if (!length(idxB) || !length(idxF)) {
      Bmsy <- try(get_par("logBmsy", rep, exp = TRUE), silent = TRUE)
      Fmsy <- try(get_par("logFmsy", rep, exp = TRUE), silent = TRUE)
      bx <- if (!inherits(Bmsy,"try-error") && is.finite(Bmsy[2])) log(Bmsy[2]) else 0
      fy <- if (!inherits(Fmsy,"try-error") && is.finite(Fmsy[2])) log(Fmsy[2]) else 0
      return(matrix(c(bx, fy), ncol = 2, dimnames = list(NULL, c("x","y"))))
    }

    sds  <- (rep$sd %||% numeric())[c(max(idxB), max(idxF))]
    s_sum <- if (length(idxS)) (rep$sd %||% numeric())[max(idxS)] else NA_real_
    cova <- if (is.finite(s_sum)) (s_sum^2 - sds[1]^2 - sds[2]^2) / 2 else 0
    covBF <- matrix(c(sds[1]^2, cova, cova, sds[2]^2), 2, 2, byrow = TRUE)
    parBF <- (rep$value %||% numeric())[c(max(idxB), max(idxF))]

    if (requireNamespace("ellipse", quietly = TRUE)) {
      ellipse::ellipse(stats::cov2cor(covBF)[1, 2],
                       scale = sqrt(diag(covBF)),
                       centre = parBF, npoints = 300)
    } else {
      matrix(c(parBF[1], parBF[2]), ncol = 2, dimnames = list(NULL, c("x","y")))
    }
  }

  annual_avg <- function(intime, vec, type = "mean") {
    fun <- match.fun(type)
    anntime   <- unique(floor(intime))
    floortime <- floor(intime)
    nstepvec  <- vapply(anntime, function(a) sum(a == floortime), numeric(1))
    anntime   <- anntime[which(nstepvec == max(nstepvec))]
    annvec    <- vapply(anntime, function(a) fun(vec[which(a == floortime)]), numeric(1))
    list(anntime = anntime, annvec = annvec)
  }

  calc_EBinf <- function(K, n, Fl, Fmsy, sdb2) {
    base <- 1 - (n - 1) / n * (Fl / Fmsy)
    base <- max(0, base)
    corr <- 1 - n / 2 / (1 - (1 - n * Fmsy + (n - 1) * Fl))
    EB  <- K * (base)^(1 / (n - 1)) * (1 - corr * sdb2)
    max(0, EB)
  }

  get_EBinf <- function(rep) {
    K     <- get_par("logK",   rep, exp = TRUE)[2]
    n     <- get_par("logn",   rep, exp = TRUE)[2]
    sdb2  <- get_par("logsdb", rep, exp = TRUE)[2]^2
    Fmsy  <- tail(get_par("logFmsy", rep, exp = TRUE), 1)[2]
    logFs <- get_par("logFs",  rep, CI = CI)
    if (min(rep$inp$dtc) < 1) {
      alf <- annual_avg(rep$inp$time, logFs[, 2])
      fff <- exp(alf$annvec)
    } else {
      fff <- exp(logFs[, 2])
    }
    Fl <- tail(unname(fff), 1)
    calc_EBinf(K, n, Fl, Fmsy, sdb2)
  }

  man_cols_local <- function() {
    colvec <- c('darkmagenta','cyan3','darkgreen','coral1','black',
                'magenta','gold','green','cadetblue3','chocolate3',
                'darkolivegreen3','cyan','darkred')
    rep(colvec, 3)
  }

  fmt1 <- function(x) {
    x <- as.numeric(x)
    ifelse(is.finite(x), formatC(x, format = "f", digits = 1), NA_character_)
  }

  ## Secondary-axis helpers
  lab_rel_1 <- function(vals) {
    out <- fmt1(vals)
    is_one <- abs(vals - 1) < 1e-9
    out[is_one] <- "1"
    out
  }

  clamp <- function(val, lo, hi) pmin(pmax(val, lo), hi)

  `%nin%` <- function(x, y) !(x %in% y)

  ## ------------------ Validate & derive quantities --------------- ##
  check_rep(rep)

  tvgflag <- isTRUE(rep$inp$timevaryinggrowth) | isTRUE(rep$inp$logmcovflag)
  if (tvgflag) rel.axes <- TRUE

  Bmsy_all <- get_par("logBmsy", rep, exp = TRUE, CI = CI)
  Fmsy_all <- get_par("logFmsy", rep, exp = TRUE, CI = CI)
  Bmsy <- tail(Bmsy_all, 1L)
  Fmsy <- tail(Fmsy_all, 1L)

  if (rel.axes) {
    ext <- FALSE
    bscal <- Bmsy[2]; fscal <- Fmsy[2]
    xlab_expr <- expression(B[t]/B[MSY])
    ylab_expr <- expression(F[t]/F[MSY])
  } else {
    bscal <- 1; fscal <- 1
    xlab_expr <- expression(B[t]); xlab_expr <- add_catchunit(xlab_expr, rep$inp$catchunit)
    ylab_expr <- expression(F[t])
  }

  Best     <- get_par("logB",  rep, exp = TRUE, CI = CI)
  logBest  <- get_par("logB",  rep,        CI = CI)
  if (tvgflag) {
    Fest  <- get_par("logFFmsy", rep, exp = TRUE, CI = CI)
    fscal <- 1; Fmsy <- c(1, 1)
  } else {
    Fest  <- get_par("logFs",  rep, exp = TRUE, CI = CI)
  }
  logFest <- get_par("logFs",  rep,        CI = CI)

  cl <- try(make_rp_ellipse(rep), silent = TRUE)
  if (inherits(cl, "try-error")) {
    cl <- matrix(c(log(Bmsy[2] %||% 1), log(Fmsy[2] %||% 1)), ncol = 2)
  }

  if (min(rep$inp$dtc) < 1) {
    alb <- annual_avg(rep$inp$time, logBest[, 2])
    alf <- annual_avg(rep$inp$time, logFest[, 2])
    bbb <- exp(alb$annvec) / bscal
    fff <- exp(alf$annvec) / fscal
    fbtime <- alb$anntime
  } else {
    bbb <- Best[rep$inp$indest, 2] / bscal
    fff <- Fest[rep$inp$indest, 2] / fscal
    fbtime <- rep$inp$time[rep$inp$indest]
  }

  EBinf <- get_EBinf(rep) / bscal

  if (is.null(xlim)) {
    xlim <- range(c(exp(cl[, 1]), Best[, 2], EBinf) / bscal, na.rm = TRUE)
    if (min(rep$inp$dtc) < 1)
      xlim <- range(c(exp(alb$annvec), exp(cl[, 1]), EBinf) / bscal, na.rm = TRUE)
    xlim[2] <- min(c(xlim[2], 8 * (Bmsy[2] %||% 1) / bscal), 2.2 * max(bbb), na.rm = TRUE)
    xlim[2] <- max(c(xlim[2], (Bmsy[2] %||% 1) / bscal), na.rm = TRUE)
  }
  if (is.null(ylim)) {
    ylim <- range(c(exp(cl[, 2]), Fest[, 2]) / fscal, na.rm = TRUE)
    if (min(rep$inp$dtc) < 1)
      ylim <- range(c(exp(logFest[, 2]) / fscal, exp(cl[, 2]) / fscal), na.rm = TRUE)
    ylim[2] <- min(c(ylim[2], 8 * (Fmsy[2] %||% 1) / fscal), 2.2 * max(fff), na.rm = TRUE)
    ylim[2] <- max(c(ylim[2], (Fmsy[2] %||% 1) / fscal), na.rm = TRUE)
    if ("man" %in% names(rep)) ylim <- range(ylim, 0)
  }
  logminval <- 1e-4
  if (isTRUE(logax)) {
    if (xlim[1] < logminval) xlim[1] <- logminval
    if (ylim[1] < logminval) ylim[1] <- logminval
  }

  if (!isTRUE(logax)) {
    if (is.finite(ylim[2]) && (ylim[1] >= 0 || !is.finite(ylim[1]))) {
      ylim[1] <- -max(1e-6, 0.02 * (ylim[2] - 0))
    }
  }

  pad_frac <- 0.02
  xpad <- pad_frac * diff(xlim)
  ypad <- pad_frac * diff(ylim)
  xlim <- c(max(if (isTRUE(logax)) logminval else 0, xlim[1] - xpad), xlim[2] + xpad)
  ylim <- c(max(if (isTRUE(logax)) logminval else 0, ylim[1] - ypad), ylim[2] + ypad)

  # Final guard for log-scales
  if (isTRUE(logax)) {
    if (!is.finite(xlim[1]) || xlim[1] <= 0) xlim[1] <- logminval
    if (!is.finite(ylim[1]) || ylim[1] <= 0) ylim[1] <- logminval
    if (xlim[2] <= xlim[1]) xlim[2] <- xlim[1] * 10
    if (ylim[2] <= ylim[1]) ylim[2] <- ylim[1] * 10
  }

  ## ------------------------- Data frames for ggplot ------------------------- ##
  xminq <- xlim[1]; xmaxq <- xlim[2]
  yminq <- ylim[1]; ymaxq <- ylim[2]

  bx <- if (rel.axes) 1 else (Bmsy[2] %||% 1) / bscal
  bx <- clamp(bx, xminq, xmaxq)
  fy <- if (rel.axes) 1 else (Fmsy[2] %||% 1) / fscal
  fy <- clamp(fy, yminq, ymaxq)

  df_green_rect   <- data.frame(xmin = bx,   xmax =  Inf, ymin = -Inf, ymax =  fy)
  df_yellowL_rect <- data.frame(xmin = -Inf, xmax =  bx,  ymin = -Inf, ymax =  fy)
  df_yellowT_rect <- data.frame(xmin = bx,   xmax =  Inf, ymin =  fy,  ymax =  Inf)
  df_red_rect     <- data.frame(xmin = -Inf, xmax =  bx,  ymin =  fy,  ymax =  Inf)

  Bmsy2 <- if (rel.axes) (Bmsy[2] %||% 1) else bscal
  Fmsy2 <- if (rel.axes) (Fmsy[2] %||% 1) else fscal

  clm <- as.matrix(cl)
  if (ncol(clm) == 2 && nrow(clm) > 1) {
    df_ellipse <- data.frame(x = exp(clm[, 1]) / Bmsy2,
                             y = exp(clm[, 2]) / Fmsy2)
  } else {
    df_ellipse <- data.frame(x = exp(clm[1]) / Bmsy2,
                             y = exp(clm[2]) / Fmsy2)
  }

  epsx <- 0.005 * diff(xlim)
  epsy <- 0.005 * diff(ylim)
  clamp_in <- function(v, lo, hi, eps) pmin(pmax(v, lo + eps), hi - eps)

  n_hist <- length(bbb)
  df_traj <- data.frame(
    x = clamp_in(bbb, xminq, xmaxq, epsx),
    y = clamp_in(fff, yminq, ymaxq, epsy),
    ord = seq_len(n_hist)
  )

  df_pred <- NULL
  df_EBseg <- NULL
  if (!(min(rep$inp$dtc) < 1) && !("man" %in% names(rep))) {
    xpred <- Best[rep$inp$indpred, 2] / (if (rel.axes) (Bmsy[2] %||% 1) else bscal)
    ypred <- Fest[rep$inp$indpred, 2] / (if (rel.axes) (Fmsy[2] %||% 1) else fscal)
    if (length(xpred) && length(ypred)) df_pred <- data.frame(x = xpred, y = ypred)
    if (!is.null(df_pred) && nrow(df_pred)) {
      Bll <- utils::tail(df_pred$x, 1)
      Fll <- utils::tail(df_pred$y, 1)
      df_EBseg <- data.frame(x = Bll, xend = EBinf, y = Fll, yend = Fll)
    }
  }

  df_first <- data.frame(x = bbb[1], y = fff[1], lab = round(fbtime[1], 2))
  df_last  <- data.frame(x = utils::tail(bbb,1), y = utils::tail(fff,1), lab = round(utils::tail(fbtime, 1), 2))

  nr <- length(rep$inp$ini$logr %||% numeric())
  df_msy_prev <- df_msy_curr <- NULL
  if (nr > 1) {
    df_msy_prev <- data.frame(x = Bmsy_all[1:(nr-1), 2] / (if (rel.axes) (Bmsy[2] %||% 1) else bscal),
                              y = Fmsy_all[1:(nr-1), 2] / (if (rel.axes) (Fmsy[2] %||% 1) else fscal))
    df_msy_curr <- data.frame(x = Bmsy_all[nr, 2] / (if (rel.axes) (Bmsy[2] %||% 1) else bscal),
                              y = Fmsy_all[nr, 2] / (if (rel.axes) (Fmsy[2] %||% 1) else fscal))
  }

  df_true <- NULL
  if ("true" %in% names(rep$inp)) {
    df_true <- data.frame(x = (rep$inp$true$Bmsy %||% NA_real_) / (if (rel.axes) (Bmsy[2] %||% 1) else bscal),
                          y = (rep$inp$true$Fmsy %||% NA_real_) / (if (rel.axes) (Fmsy[2] %||% 1) else fscal))
  }

  # Canonical scenario order/labels
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

  df_man <- df_man_int <- df_man_jump <- NULL
  leg_man <- NULL
  scs_final <- character(0); lab_final <- character(0); colv <- NULL

  if ("man" %in% names(rep) && length(rep$man)) {
    nman <- length(rep$man)
    leg_man <- names(rep$man)
    if (is.null(leg_man) || !length(leg_man)) leg_man <- paste0("Scenario ", seq_len(nman))

    sc_present <- leg_man
    sc_core    <- intersect(sc_order, sc_present)
    sc_other   <- setdiff(sc_present, sc_core)
    scs_final  <- c(sc_core, sort(sc_other))
    lab_final  <- ifelse(scs_final %in% names(label_map), unname(label_map[scs_final]), scs_final)

    df_list      <- vector("list", nman)
    df_int_list  <- vector("list", nman)
    df_jump_list <- vector("list", nman)

    same_pt <- function(a, b) isTRUE(all.equal(a, b, tolerance = 1e-10))

    for (i in seq_len(nman)) {
      rp <- rep$man[[i]]
      estF <- get.par("logF", rp, exp = TRUE)[, 2]
      estB <- get.par("logB", rp, exp = TRUE)[, 2]
      time <- rp$inp$time
      manint <- rp$inp$maninterval
      lastobs <- rp$inp$timerangeObs[2]

      if (min(rep$inp$dtc) < 1) {
        alb2 <- list(anntime = unique(floor(time)),
                     annvec  = tapply(log(estB), floor(time), mean))
        alf2 <- list(anntime = unique(floor(time)),
                     annvec  = tapply(log(estF), floor(time), mean))
        bx <- exp(unname(alb2$annvec)) / (if (rel.axes) (Bmsy[2] %||% 1) else bscal); ttA <- as.numeric(names(alb2$annvec))
        fy <- exp(unname(alf2$annvec)) / (if (rel.axes) (Fmsy[2] %||% 1) else fscal)
        keep_man <- (ttA >= manint[1])
        keep_int <- (ttA >= lastobs) & (ttA < manint[1])
        pre_mask <- ttA < manint[1]
      } else {
        bx <- estB / (if (rel.axes) (Bmsy[2] %||% 1) else bscal); ttA <- time
        fy <- estF / (if (rel.axes) (Fmsy[2] %||% 1) else fscal)
        keep_man <- (time >= manint[1])
        keep_int <- (time >= lastobs) & (time < manint[1])
        pre_mask <- time < manint[1]
      }

      if (any(keep_man, na.rm = TRUE)) {
        maninds <- which(keep_man)
        df_list[[i]] <- data.frame(x = bx[maninds], y = fy[maninds], scenario = leg_man[i])
      } else df_list[[i]] <- data.frame(x = numeric(0), y = numeric(0), scenario = character(0))

      if (any(keep_int, na.rm = TRUE)) {
        df_int_list[[i]] <- data.frame(x = bx[which(keep_int)], y = fy[which(keep_int)], scenario = leg_man[i])
      } else df_int_list[[i]] <- data.frame(x = numeric(0), y = numeric(0), scenario = character(0))

      x0 <- if (any(pre_mask, na.rm = TRUE)) utils::tail(bx[pre_mask], 1) else NA_real_
      y0 <- if (any(pre_mask, na.rm = TRUE)) utils::tail(fy[pre_mask], 1) else NA_real_
      x1 <- if (any(keep_man, na.rm = TRUE)) bx[which(keep_man)[1]] else NA_real_
      y1 <- if (any(keep_man, na.rm = TRUE)) fy[which(keep_man)[1]] else NA_real_

      if (is.finite(x0) && is.finite(y0) && is.finite(x1) && is.finite(y1) &&
          (abs(x1 - x0) > 0 || abs(y1 - y0) > 0)) {
        df_jump_list[[i]] <- data.frame(x = x0, y = y0, xend = x1, yend = y1, scenario = leg_man[i])
      } else {
        df_jump_list[[i]] <- data.frame(x = numeric(0), y = numeric(0), xend = numeric(0), yend = numeric(0),
                                        scenario = character(0))
      }

      if (is.finite(x0) && is.finite(y0)) {
        if (!nrow(df_int_list[[i]]) ||
            !same_pt(utils::tail(df_int_list[[i]]$x, 1), x0) ||
            !same_pt(utils::tail(df_int_list[[i]]$y, 1), y0)) {
          df_int_list[[i]] <- rbind(df_int_list[[i]], data.frame(x = x0, y = y0, scenario = leg_man[i]))
        }
      }
    }

    df_man       <- do.call(rbind, df_list)
    df_man_int   <- do.call(rbind, df_int_list)
    df_man_jump  <- do.call(rbind, df_jump_list)

    if (!is.null(df_man) && nrow(df_man) && length(scs_final)) {
      df_man$scenario <- factor(df_man$scenario, levels = scs_final)
    }
    if (!is.null(df_man_int) && nrow(df_man_int) && length(scs_final)) {
      df_man_int$scenario <- factor(df_man_int$scenario, levels = scs_final)
    }
    if (!is.null(df_man_jump) && nrow(df_man_jump) && length(scs_final)) {
      df_man_jump$scenario <- factor(df_man_jump$scenario, levels = scs_final)
    }

    colv <- man_cols_local()[seq_along(scs_final)]
    names(colv) <- scs_final
  }

  ## --------------------------- Build ggplot ------------------------------ ##
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
                                   fill = grDevices::adjustcolor("gray90", alpha.f = 0.5),
                                   color = grDevices::adjustcolor("gray90", alpha.f = 0.5), linewidth = 0.8)
  } else {
    p <- p + ggplot2::geom_point(data = df_ellipse, ggplot2::aes(x = x, y = y),
                                 color = "gray40", size = 2)
  }

  p <- p + ggplot2::geom_path(data = df_traj, ggplot2::aes(x = x, y = y),
                              linewidth = 0.4, color = grDevices::adjustcolor("blue", alpha.f = 0.8))

  if (!is.null(df_pred) && nrow(df_pred) > 1) {
    p <- p + ggplot2::geom_path(
      data = df_pred,
      ggplot2::aes(x = x, y = y),
      linetype = "33",
      color = grDevices::adjustcolor("blue", alpha.f = 0.8)
    )
  } else if (!is.null(df_pred) && nrow(df_pred) == 1) {
    p <- p + ggplot2::geom_point(
      data = df_pred,
      ggplot2::aes(x = x, y = y),
      shape = 16, size = 1.8,
      color = grDevices::adjustcolor("blue", alpha.f = 0.8)
    )
  }

  if (!is.null(df_EBseg) && nrow(df_EBseg)) {
    p <- p + ggplot2::geom_segment(data = df_EBseg,
                                   ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
                                   color = "blue", linetype = "33", linewidth = 1.0)
    p <- p + ggplot2::geom_point(data = data.frame(x = df_EBseg$xend, y = df_EBseg$yend),
                                 ggplot2::aes(x = x, y = y),
                                 shape = 23, fill = "gold", color = "black", size = 3, stroke = 0.5)
  }

  if (!is.null(df_man_int) && nrow(df_man_int)) {
    dfmi <- df_man_int
    if (!"scenario" %in% names(dfmi)) dfmi$scenario <- "Scenario"
    dfmi$..n.. <- ave(dfmi$x, dfmi$scenario, FUN = length)

    if (any(dfmi$..n.. > 1)) {
      p <- p + ggplot2::geom_path(
        data = dfmi[dfmi$..n.. > 1, , drop = FALSE],
        ggplot2::aes(x = x, y = y, color = scenario),
        linewidth = 0.7, linetype = "33",
        inherit.aes = FALSE, show.legend = man.legend
      )
    }
    if (any(dfmi$..n.. == 1)) {
      p <- p + ggplot2::geom_point(
        data = dfmi[dfmi$..n.. == 1, , drop = FALSE],
        ggplot2::aes(x = x, y = y, color = scenario),
        size = 1.8, inherit.aes = FALSE, show.legend = man.legend
      )
    }
  }

  if (!is.null(df_man) && nrow(df_man)) {
    dfm <- df_man
    if (!"scenario" %in% names(dfm)) dfm$scenario <- "Scenario"
    dfm$..n.. <- ave(dfm$x, dfm$scenario, FUN = length)

    if (any(dfm$..n.. > 1)) {
      p <- p + ggplot2::geom_path(
        data = dfm[dfm$..n.. > 1, , drop = FALSE],
        ggplot2::aes(x = x, y = y, color = scenario),
        linewidth = 0.7,
        inherit.aes = FALSE, show.legend = man.legend
      )
    }
    if (any(dfm$..n.. == 1)) {
      p <- p + ggplot2::geom_point(
        data = dfm[dfm$..n.. == 1, , drop = FALSE],
        ggplot2::aes(x = x, y = y, color = scenario),
        size = 1.8, inherit.aes = FALSE, show.legend = man.legend
      )
    }
  }

  if (!is.null(df_man_jump) && nrow(df_man_jump)) {
    p <- p + ggplot2::geom_segment(data = df_man_jump,
                                   ggplot2::aes(x = x, y = y, xend = xend, yend = yend, color = scenario),
                                   linewidth = 0.7, inherit.aes = FALSE, show.legend = man.legend)
  }

  if (!is.null(df_msy_prev) && nrow(df_msy_prev)) {
    p <- p + ggplot2::geom_point(data = df_msy_prev, ggplot2::aes(x = x, y = y),
                                 shape = 3, color = "magenta", size = 2)
  }
  if (!is.null(df_msy_curr) && nrow(df_msy_curr)) {
    p <- p + ggplot2::geom_point(data = df_msy_curr, ggplot2::aes(x = x, y = y),
                                 shape = 3, color = "black", size = 2)
  }

  if (!is.null(df_true) && nrow(df_true)) {
    p <- p + ggplot2::geom_point(data = df_true, ggplot2::aes(x = x, y = y),
                                 shape = 25, fill = grDevices::adjustcolor(grDevices::rgb(1, 165/255, 0), alpha.f = 0.7),
                                 color = "black", size = 3)
  }

  p <- p + ggplot2::geom_point(data = df_first, ggplot2::aes(x = x, y = y),
                               shape = 21, fill = "white", color = "black", size = 2.5, stroke = 0.5)
  p <- p + ggplot2::geom_text(data = df_first, ggplot2::aes(x = x, y = y, label = lab),
                              vjust = -0.8, size = 3)
  p <- p + ggplot2::geom_point(data = df_last, ggplot2::aes(x = x, y = y),
                               shape = 22, fill = "white", color = "black", size = 2.5, stroke = 0.5)
  p <- p + ggplot2::geom_text(data = df_last, ggplot2::aes(x = x, y = y, label = lab),
                              vjust = -0.8, size = 3)

  if (!isTRUE(logax) && (0 >= xlim[1]) && (0 <= xlim[2])) {
    p <- p + ggplot2::geom_vline(xintercept = 0, color = "darkred", linetype = 2)
  }

  # Labels & scales
  xlab_final <- if (is.null(xlabel)) xlab_expr else xlabel
  ylab_final <- ylab_expr

  Bmsy_den <- if (is.finite(Bmsy[2]) && (Bmsy[2] > 0)) Bmsy[2] else 1
  Fmsy_den <- if (is.finite(Fmsy[2]) && (Fmsy[2] > 0)) Fmsy[2] else 1

  if (isTRUE(logax)) {
    p <- p + ggplot2::scale_x_log10(
      limits = xlim,
      name   = xlab_final,
      labels = function(x) formatC(x, format = "f", digits = 0),
      expand = c(0, 0),
      sec.axis = if (!rel.axes && isTRUE(ext)) ggplot2::sec_axis(
        trans  = ~ . / Bmsy_den,
        name   = expression(B[t]/B[MSY]),
        breaks = function(lims) { br <- scales::pretty_breaks(n = 6)(lims); sort(unique(c(br, 1))) },
        labels = lab_rel_1
      ) else ggplot2::waiver()
    )

    p <- p + ggplot2::scale_y_log10(
      limits = ylim,
      name   = ylab_final,
      labels = fmt1,
      expand = ggplot2::expansion(mult = c(0.03, 0.03)),
      sec.axis = if (!rel.axes && isTRUE(ext)) ggplot2::sec_axis(
        trans  = ~ . / Fmsy_den,
        name   = expression(F[t]/F[MSY]),
        breaks = function(lims) { br <- scales::pretty_breaks(n = 6)(lims); sort(unique(c(br, 1))) },
        labels = lab_rel_1
      ) else ggplot2::waiver()
    )

  } else {

    p <- p + ggplot2::scale_x_continuous(
      limits = xlim,
      name   = xlab_final,
      labels = function(x) formatC(x, format = "f", digits = 0),
      breaks = function(lims) pretty(lims, n = 6),
      expand = c(0, 0),
      sec.axis = if (!rel.axes && isTRUE(ext)) ggplot2::sec_axis(
        trans  = ~ . / Bmsy_den,
        name   = expression(B[t]/B[MSY]),
        breaks = function(lims) { br <- scales::pretty_breaks(n = 6)(lims); sort(unique(c(br, 1))) },
        labels = lab_rel_1
      ) else ggplot2::waiver()
    )

    p <- p + ggplot2::scale_y_continuous(
      limits = ylim,
      name   = ylab_final,
      labels = fmt1,
      expand = ggplot2::expansion(mult = c(0.03, 0.03)),
      sec.axis = if (!rel.axes && isTRUE(ext)) ggplot2::sec_axis(
        trans  = ~ . / Fmsy_den,
        name   = expression(F[t]/F[MSY]),
        breaks = function(lims) { br <- scales::pretty_breaks(n = 6)(lims); sort(unique(c(br, 1))) },
        labels = lab_rel_1
      ) else ggplot2::waiver()
    )
  }

  if (!is.null(colv) && length(colv) && man.legend) {
    p <- p + ggplot2::scale_color_manual(values = colv,
                                         breaks = names(colv),
                                         labels = if (length(lab_final)) lab_final else names(colv),
                                         name   = "Scenario")
  } else if (!is.null(colv) && length(colv) && !man.legend) {
    p <- p + ggplot2::scale_color_manual(values = colv, guide = "none")
  } else {
    p <- p + ggplot2::scale_color_discrete(guide = "none")
  }

  p <- p + ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.6),
      legend.background = ggplot2::element_rect(fill = grDevices::adjustcolor("white", alpha.f = 0.8), color = NA),
      legend.position = if (man.legend) "top" else "none",
      legend.direction = "horizontal",
      axis.text.x  = ggplot2::element_text(size = 10, face = "bold"),
      axis.text.y  = ggplot2::element_text(size = 10, face = "bold"),
      axis.title.x = ggplot2::element_text(size = 14, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 14, face = "bold"),
      plot.margin = ggplot2::margin(4, 4, 4, 4)
    )

  if (plot.legend && !(min(rep$inp$dtc) < 1)) {
    px <- xlim[2] - 0.02 * diff(xlim)
    py <- ylim[2] - 0.04 * diff(ylim)
    gap <- 0.12 * diff(xlim)
    p <- p + ggplot2::annotate("point", x = px - gap, y = py,
                               shape = 23, size = 2.8, fill = "gold", color = "black", stroke = 0.5)
    p <- p + ggplot2::annotate("text", x = px, y = py,
                               label = "E(B[infinity])", parse = TRUE,
                               hjust = 1, vjust = 0.5, fontface = "bold", size = 5.5)
  }

  if (!is.null(rep$opt$convergence) && rep$opt$convergence != 0) {
    p <- p + ggplot2::annotate("text",
                               x = xlim[1] + 0.01 * diff(xlim),
                               y = ylim[2] + 0.02 * diff(ylim),
                               label = "▲", color = "black", size = 6)
    p <- p + ggplot2::annotate("text",
                               x = xlim[1] + 0.01 * diff(xlim),
                               y = ylim[2] + 0.02 * diff(ylim),
                               label = "!", color = "black", size = 3, vjust = 0.35)
  }

  if (!is.null(stamp) && nzchar(stamp)) {
    p <- p + ggplot2::annotate("text",
                               x = xlim[2] - 0.01 * diff(xlim),
                               y = ylim[1] - 0.04 * diff(ylim),
                               label = stamp, hjust = 1, vjust = 0, size = 3)
  }

  p <- p + ggplot2::guides(color = "none")
  p <- p + ggplot2::theme(legend.position = "none")

  # ----- Auto-print so direct calls render a figure -----
  .ap_ok <- isTRUE(getOption("elu2.autoprint", TRUE)) &&
    (interactive() || isTRUE(getOption("knitr.in.progress", FALSE)))
  if (.ap_ok) {
    try(print(p), silent = TRUE)
  }
  invisible(p)
}


# ===================================================================
# SHARED INTERNALS (EXPORTED/INTERNAL AS NEEDED)
# ===================================================================

#' Reusable colour palette for management scenarios
#'
#' @description
#' Returns a repeatable vector of colours to map across an arbitrary number
#' of management scenarios. Colours recycle to the requested length.
#'
#' @param n Integer. Number of colours required.
#'
#' @return A character vector of length `n` with hex colour strings.
#'
#' @examples
#' man_cols(5)
#'
#' @keywords internal
#' @noRd
man_cols <- function(n) {
  base <- c('darkmagenta','cyan3','darkgreen','coral1','black',
            'magenta','gold','green','cadetblue3','chocolate3',
            'darkolivegreen3','cyan','darkred')
  rep(base, length.out = n)
}

#' Observed indices on the relative B/Bmsy scale
#'
#' @description
#' Builds a tidy data frame of observed index points scaled by \eqn{q} and
#' divided by \eqn{B_{MSY}} so they can be plotted directly against B/Bmsy.
#'
#' @param rep_main A fitted SPiCT report object (`spictcls`).
#' @return A data frame with columns `index`, `time`, `obsrel`, or `NULL`.
#' @keywords internal
#' @noRd
obs_indices_rel_df <- function(rep_main) {
  qest  <- get.par("logq", rep_main, exp = TRUE)
  Bmsy  <- get.par("logBmsy", rep_main, exp = TRUE)
  Bmsy2 <- if (is.null(nrow(Bmsy))) Bmsy[2] else Bmsy[1, 2]
  inp   <- rep_main$inp

  out_list <- list()
  k <- 0L

  if (!is.null(inp$timeI) && length(inp$timeI) >= 1L &&
      !is.null(inp$obsI)  && length(inp$obsI)  >= 1L) {
    k <- k + 1L
    out_list[[k]] <- data.frame(
      index  = "Index1",
      time   = inp$timeI[[1]],
      obsrel = (inp$obsI[[1]] / qest[inp$mapq[1], 2]) / Bmsy2
    )
  }
  if (!is.null(inp$timeI) && length(inp$timeI) >= 2L &&
      !is.null(inp$obsI)  && length(inp$obsI)  >= 2L) {
    k <- k + 1L
    out_list[[k]] <- data.frame(
      index  = "Index2",
      time   = inp$timeI[[2]],
      obsrel = (inp$obsI[[2]] / qest[inp$mapq[2], 2]) / Bmsy2
    )
  }

  if (length(out_list) == 0L) NULL else do.call(rbind, out_list)
}

#' Historical B/Bmsy and F/Fmsy up to the last observed time
#' @keywords internal
#' @noRd
get_base_BB_FF_pre <- function(rep_man, CI = 0.95) {
  stopifnot(inherits(rep_man, "spictcls"))
  t_last_obs <- rep_man$inp$timerangeObs[2]
  tt <- rep_man$inp$time
  keep <- which(tt <= t_last_obs)

  BB <- get.par("logBBmsy", rep_man, exp = TRUE)
  FF <- get.par("logFFmsy", rep_man, exp = TRUE)

  BB <- BB[keep, , drop = FALSE]
  FF <- FF[keep, , drop = FALSE]

  list(
    t_last_obs = t_last_obs,
    BB = data.frame(
      time = tt[keep], lwr = BB[, 1], est = BB[, 2], upr = BB[, 3],
      row.names = NULL
    ),
    FF = data.frame(
      time = tt[keep], lwr = FF[, 1], est = FF[, 2], upr = FF[, 3],
      row.names = NULL
    )
  )
}

#' Prepare data frames for management panels (fit or manage)
#' @keywords internal
#' @noRd
prepare_manage_panel_data <- function(rep_man, CI = 0.95) {
  stopifnot(inherits(rep_man, "spictcls"))

  has_man <- ("man" %in% names(rep_man)) && length(rep_man$man)

  co_df <- if (!is.null(rep_man$inp$timeC) && !is.null(rep_man$inp$obsC)) {
    data.frame(time = rep_man$inp$timeC, catch = rep_man$inp$obsC, catch_type = "Observed")
  } else NULL

  if (has_man) {
    scenarios <- names(rep_man$man)
    ns        <- length(scenarios)

    bb_list <- vector("list", ns)
    ff_list <- vector("list", ns)
    cp_list <- vector("list", ns)
    eval_list <- vector("list", ns)

    for (i in seq_len(ns)) {
      sc <- scenarios[i]; rp <- rep_man$man[[sc]]

      BB <- get.par("logBBmsy", rp, exp = TRUE)
      FF <- get.par("logFFmsy", rp, exp = TRUE)

      bb_list[[i]] <- data.frame(
        time = as.numeric(rownames(BB)), lwr = BB[, 1], est = BB[, 2], upr = BB[, 3], scenario = sc
      )
      ff_list[[i]] <- data.frame(
        time = as.numeric(rownames(FF)), lwr = FF[, 1], est = FF[, 2], upr = FF[, 3], scenario = sc
      )

      CP <- get.par("logCpred", rp, exp = TRUE)
      cp_list[[i]] <- data.frame(
        time = rp$inp$timeCpred, lwr = CP[, 1], catch = CP[, 2], upr = CP[, 3],
        scenario = sc, catch_type = "Predicted"
      )

      ti  <- rp$inp$time; tE  <- rp$inp$maneval; idx <- which(ti == tE)
      if (length(idx) == 1L) {
        eval_list[[i]] <- data.frame(
          scenario = sc, maneval = tE,
          BB_est   = BB[idx, 2], FF_est = FF[idx, 2]
        )
      }
    }

    t0 <- if (!is.null(rep_man$inp$maninterval) && length(rep_man$inp$maninterval) >= 1L)
      rep_man$inp$maninterval[1] else NA_real_

    list(
      bbmsy      = do.call(rbind, bb_list),
      ffmsy      = do.call(rbind, ff_list),
      catch_pred = do.call(rbind, cp_list),
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
  mi <- rep$inp$maninterval
  if (is.null(mi) || length(mi) < 1L || !is.finite(mi[1])) return(list())

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

#' Management vertical lines for Catch panel
#' @keywords internal
#' @noRd
add_management_vlines_catch_good <- function(rep,
                                             color = "grey30",
                                             linetype = "dashed",
                                             linewidth = 0.6,
                                             lineend = "butt") {
  mi <- rep$inp$maninterval
  last_obs_start <- if (!is.null(rep$inp$timeC) && length(rep$inp$timeC)) {
    utils::tail(rep$inp$timeC, 1)
  } else NA_real_

  if (is.null(mi) || length(mi) < 1L || !is.finite(mi[1])) {
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

# ===================================================================
# PANELS: BBMSY / FFMSY / CATCH
# ===================================================================

#' B/Bmsy panel (ggplot2) with optional scenario overlays
#' @export
my_plot_manage_bbmsy_panel <- function(rep_man,
                                       scenario_color = NULL,
                                       show_CIs = TRUE,
                                       CI = 0.95,
                                       show_legend = TRUE) {
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

  dat <- prepare_manage_panel_data(rep_man, CI = CI)

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

  if (!is.null(dat$bbmsy) && nrow(dat$bbmsy) && length(scs_final)) {
    dat$bbmsy$scenario <- factor(dat$bbmsy$scenario, levels = scs_final)
  }

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
    { if (show_CIs && nrow(base_pre$BB)) ggplot2::geom_ribbon(
      data = base_pre$BB,
      ggplot2::aes(x = time, ymin = lwr, ymax = upr),
      inherit.aes = FALSE, fill = spict_blue_ci_fill, color = NA
    ) } +
    { if (show_CIs && nrow(base_pre$BB)) list(
      ggplot2::geom_line(data = base_pre$BB, ggplot2::aes(x = time, y = lwr),
                         inherit.aes = FALSE, color = spict_blue_ci_line, linewidth = 0.6),
      ggplot2::geom_line(data = base_pre$BB, ggplot2::aes(x = time, y = upr),
                         inherit.aes = FALSE, color = spict_blue_ci_line, linewidth = 0.6)
    ) } +
    { if (nrow(base_pre$BB)) ggplot2::geom_line(
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
    } else {
      ggplot2::geom_vline(xintercept = base_pre$t_last_obs,
                          color = "grey75", linetype = "solid", linewidth = 0.4)
    }
    } +
    ggplot2::labs(x = "Year", y = expression(bold(B/B[MSY]))) +
    { if (length(scs_final)) ggplot2::scale_color_manual(
      values = cols, breaks = scs_final, labels = lab_final, name = NULL,
      guide = if (show_legend) "legend" else "none"
    ) else ggplot2::scale_color_discrete(guide = "none") } +
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

  # One-step-ahead for raw fit (no $man)
  if (!length(scs_final)) {
    BB_all <- get.par("logBBmsy", rep_man, exp = TRUE, CI = CI)
    tt <- suppressWarnings(as.numeric(rownames(BB_all)))
    if (!length(tt) || any(!is.finite(tt))) {
      tt_full <- as.numeric(rep_man$inp$time)
      tt <- if (length(tt_full) >= nrow(BB_all)) tt_full[seq_len(nrow(BB_all))] else rep_len(tt_full, nrow(BB_all))
    }
    df_all <- data.frame(time = tt, lwr = BB_all[,1], est = BB_all[,2], upr = BB_all[,3])

    ind_pr <- if (!is.null(rep_man$inp$indpred)) rep_man$inp$indpred else integer(0)
    df_pr  <- if (length(ind_pr)) df_all[ind_pr, , drop = FALSE] else df_all[0, ]

    if (nrow(df_pr)) {
      if (isTRUE(show_CIs)) {
        pB <- pB +
          ggplot2::geom_ribbon(
            data = transform(df_pr, ymin = lwr, ymax = upr),
            ggplot2::aes(x = time, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = spict_blue_ci_fill, color = NA
          ) +
          ggplot2::geom_line(
            data = df_pr, ggplot2::aes(x = time, y = lwr),
            inherit.aes = FALSE, color = spict_blue_ci_line, linewidth = 0.6
          ) +
          ggplot2::geom_line(
            data = df_pr, ggplot2::aes(x = time, y = upr),
            inherit.aes = FALSE, color = spict_blue_ci_line, linewidth = 0.6
          )
      }
      pB <- pB +
        ggplot2::geom_line(
          data = df_pr, ggplot2::aes(x = time, y = est),
          linetype = "dashed", color = spict_blue_mean, linewidth = 0.8
        )
    }
  }

  pB
}

#' F/Fmsy panel (ggplot2) with optional scenario overlays
#' @export
my_plot_manage_ffmsy_panel <- function(rep_man,
                                       scenario_color = NULL,
                                       show_CIs = TRUE,
                                       CI = 0.95,
                                       show_legend = FALSE) {
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

  dat <- prepare_manage_panel_data(rep_man, CI = CI)

  sc_order <- c("currentCatch","currentF","Fmsy","noF",
                "reduceF25","increaseF25","msyHockeyStick","ices")
  sc_present <- if (("man" %in% names(rep_man)) && length(rep_man$man)) names(rep_man$man) else character(0)
  sc_core   <- intersect(sc_order, sc_present)
  sc_other  <- setdiff(sc_present, sc_core)
  scs_final <- c(sc_core, sort(sc_other))

  if (!is.null(dat$ffmsy) && nrow(dat$ffmsy) && length(scs_final)) {
    dat$ffmsy$scenario <- factor(dat$ffmsy$scenario, levels = scs_final)
  }

  if (is.null(scenario_color)) {
    cols <- if (length(scs_final)) { x <- man_cols(length(scs_final)); names(x) <- scs_final; x } else NULL
  } else {
    cols <- scenario_color
    if (length(scs_final)) {
      if (is.null(names(cols))) names(cols) <- scs_final
      cols <- cols[match(scs_final, names(cols))]; names(cols) <- scs_final
    }
  }

  spict_blue_mean    <- "blue"
  spict_blue_ci_line <- grDevices::rgb(0, 0, 1, 0.20)
  spict_blue_ci_fill <- grDevices::rgb(0, 0, 1, 0.10)

  base_pre <- get_base_BB_FF_pre(rep_man, CI = CI)
  has_man  <- length(scs_final) > 0

  pF <- ggplot2::ggplot() +
    { if (show_CIs && nrow(base_pre$FF)) ggplot2::geom_ribbon(
      data = base_pre$FF,
      ggplot2::aes(x = time, ymin = lwr, ymax = upr),
      inherit.aes = FALSE, fill = spict_blue_ci_fill, color = NA
    ) } +
    { if (show_CIs && nrow(base_pre$FF)) list(
      ggplot2::geom_line(data = base_pre$FF, ggplot2::aes(x = time, y = lwr),
                         inherit.aes = FALSE, color = spict_blue_ci_line, linewidth = 0.7),
      ggplot2::geom_line(data = base_pre$FF, ggplot2::aes(x = time, y = upr),
                         inherit.aes = FALSE, color = spict_blue_ci_line, linewidth = 0.7)
    ) } +
    { if (nrow(base_pre$FF)) ggplot2::geom_line(
      data = base_pre$FF, ggplot2::aes(x = time, y = est),
      inherit.aes = FALSE,
      color = spict_blue_mean,
      linewidth = 0.8
    ) } +
    { if (has_man) {
      add_management_vlines_BF_good(rep_man, color = "grey75",
                                    linetype = "solid", linewidth = 0.4, lineend = "butt")
    } else {
      ggplot2::geom_vline(xintercept = base_pre$t_last_obs,
                          color = "grey75", linetype = "solid", linewidth = 0.4)
    }
    } +
    ggplot2::geom_line(
      data = dat$ffmsy,
      ggplot2::aes(x = time, y = est, color = scenario),
      linewidth = 0.8, na.rm = TRUE
    ) +
    ggplot2::geom_hline(yintercept = 1, linetype = "solid") +
    ggplot2::labs(x = "Year", y = expression(bold(F/F[MSY]))) +
    { if (length(scs_final)) ggplot2::scale_color_manual(
      values = cols, guide = if (show_legend) "legend" else "none"
    ) else ggplot2::scale_color_discrete(guide = "none") } +
    ggplot2::scale_fill_manual(values = if (is.null(cols)) NA else cols, guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.02))) +
    ggplot2::geom_line(
      data = base_pre$FF,
      mapping = ggplot2::aes(x = time, y = est, group = 1),
      inherit.aes = FALSE,
      colour = I(spict_blue_mean),
      linewidth = 0.9,
      show.legend = FALSE
    ) +
    theme_minimal_compact2_good_local()

  # Raw fit: one-ahead dotted + ribbon/edges
  if (!has_man) {
    FF_all <- get.par("logFFmsy", rep_man, exp = TRUE, CI = CI)
    tt <- suppressWarnings(as.numeric(rownames(FF_all)))
    if (!length(tt) || any(!is.finite(tt))) {
      tt_full <- as.numeric(rep_man$inp$time)
      tt <- if (length(tt_full) >= nrow(FF_all)) tt_full[seq_len(nrow(FF_all))] else rep_len(tt_full, nrow(FF_all))
    }
    df_all <- data.frame(time = tt, lwr = FF_all[,1], est = FF_all[,2], upr = FF_all[,3])

    ind_pr <- if (!is.null(rep_man$inp$indpred)) rep_man$inp$indpred else integer(0)
    df_pr  <- if (length(ind_pr)) df_all[ind_pr, , drop = FALSE] else df_all[0, ]

    obs_end <- base_pre$t_last_obs

    if (nrow(df_pr)) {
      if (isTRUE(show_CIs)) {
        pF <- pF +
          ggplot2::geom_ribbon(
            data = transform(df_pr, ymin = lwr, ymax = upr),
            ggplot2::aes(x = time, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = spict_blue_ci_fill, color = NA
          ) +
          ggplot2::geom_line(
            data = df_pr, ggplot2::aes(x = time, y = lwr),
            inherit.aes = FALSE, color = spict_blue_ci_line, linewidth = 0.7
          ) +
          ggplot2::geom_line(
            data = df_pr, ggplot2::aes(x = time, y = upr),
            inherit.aes = FALSE, color = spict_blue_ci_line, linewidth = 0.7
          )
      }
      pF <- pF +
        ggplot2::geom_line(
          data = df_pr, ggplot2::aes(x = time, y = est),
          linetype = "dotted", color = spict_blue_mean, linewidth = 0.8
        )

      i0 <- if (any(tt <= obs_end)) max(which(tt <= obs_end)) else NA_integer_
      i1 <- if (any(tt  > obs_end)) min(which(tt  > obs_end)) else NA_integer_
      if (is.finite(i0) && is.finite(i1)) {
        pF <- pF +
          ggplot2::geom_segment(
            ggplot2::aes(x = tt[i0], xend = tt[i1], y = df_all$est[i0], yend = df_all$est[i1]),
            linetype = "dotted", color = spict_blue_mean, linewidth = 0.8, lineend = "round"
          )
      }
    }
  }

  pF
}

#' Catch panel (ggplot2) with MSY band and optional scenario overlays
#' @export
my_plot_manage_catch_panel <- function(rep_man,
                                       scenario_color = NULL,
                                       show_CIs = TRUE,
                                       CI = 0.95,
                                       show_legend = FALSE) {
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

  dat <- prepare_manage_panel_data(rep_man, CI = CI)

  sc_order <- c("currentCatch","currentF","Fmsy","noF",
                "reduceF25","increaseF25","msyHockeyStick","ices")
  sc_present <- if (("man" %in% names(rep_man)) && length(rep_man$man)) names(rep_man$man) else character(0)
  sc_core   <- intersect(sc_order, sc_present)
  sc_other  <- setdiff(sc_present, sc_core)
  scs_final <- c(sc_core, sort(sc_other))
  has_man   <- length(scs_final) > 0

  if (!is.null(dat$catch_pred) && nrow(dat$catch_pred) && length(scs_final)) {
    dat$catch_pred$scenario <- factor(dat$catch_pred$scenario, levels = scs_final)
  }

  if (is.null(scenario_color)) {
    cols <- if (length(scs_final)) { x <- man_cols(length(scs_final)); names(x) <- scs_final; x } else NULL
  } else {
    cols <- scenario_color
    if (length(scs_final)) {
      if (is.null(names(cols))) names(cols) <- scs_final
      cols <- cols[match(scs_final, names(cols))]; names(cols) <- scs_final
    }
  }

  spict_blue_mean    <- "blue"
  spict_blue_ci_line <- grDevices::rgb(0, 0, 1, 0.20)

  t0 <- dat$t0

  CP0 <- try(get.par("logCpred", rep_man, exp = TRUE), silent = TRUE)
  base_C_t0 <- NULL
  if (!inherits(CP0, "try-error") && !is.null(rep_man$inp$timeCpred)) {
    tc <- rep_man$inp$timeCpred
    maskC <- is.finite(tc)
    if (is.finite(t0)) maskC <- maskC & (tc < t0)
    if (any(maskC)) {
      base_C_t0 <- data.frame(
        time = tc[maskC], lwr  = CP0[maskC, 1], est  = CP0[maskC, 2], upr  = CP0[maskC, 3]
      )
    }
  }

  tvgflag <- isTRUE(rep_man$inp$timevaryinggrowth) || isTRUE(rep_man$inp$logmcovflag)
  if (tvgflag) {
    MSYmat <- get.par("logMSYvec", rep_man, exp = TRUE, CI = CI)
    t_msy  <- suppressWarnings(as.numeric(rownames(MSYmat)))
    if (!length(t_msy) || any(!is.finite(t_msy))) t_msy <- rep_man$inp$time
    msy_df <- data.frame(time = t_msy, lwr = MSYmat[, 1], est = MSYmat[, 2], upr = MSYmat[, 3])
  } else {
    MSYone <- get.par("logMSY", rep_man, exp = TRUE, CI = CI) # c(ll, est, ul)
    t_msy  <- rep_man$inp$time
    msy_df <- data.frame(time = t_msy,
                         lwr = rep(MSYone[1], length(t_msy)),
                         est = rep(MSYone[2], length(t_msy)),
                         upr = rep(MSYone[3], length(t_msy)))
  }
  msy_df$lwr[msy_df$lwr < 0] <- 0
  msy_df$est[msy_df$est < 0] <- 0
  msy_df$upr[msy_df$upr < 0] <- 0

  darker_color <- function(col, f = 0.8) {
    m <- grDevices::col2rgb(col) / 255
    grDevices::rgb(m[1]*f, m[2]*f, m[3]*f)
  }
  dot_blue <- darker_color(spict_blue_mean, 0.6)

  obs_end <- if (!is.null(rep_man$inp$timeC) && length(rep_man$inp$timeC))
    utils::tail(rep_man$inp$timeC, 1) else NA_real_

  pC <- ggplot2::ggplot() +
    { if (nrow(msy_df)) ggplot2::geom_ribbon(
      data = msy_df,
      ggplot2::aes(x = time, ymin = lwr, ymax = upr),
      inherit.aes = FALSE, fill = "lightgray"
    ) } +
    { if (nrow(msy_df)) ggplot2::geom_line(
      data = msy_df,
      ggplot2::aes(x = time, y = est),
      inherit.aes = FALSE, color = "black", linewidth = 0.8
    ) } +
    { if (has_man)
      add_management_vlines_catch_good(rep_man, color = "grey75", linetype = "solid",
                                       linewidth = 0.4, lineend = "butt")
      else
        add_management_vlines_catch_good(rep_man, color = "grey75", linetype = "solid",
                                         linewidth = 0.4, lineend = "butt")[1]
    } +
    ggplot2::geom_line(
      data = dat$catch_pred,
      ggplot2::aes(x = time, y = catch, color = scenario),
      linewidth = 0.5, na.rm = TRUE
    ) +
    { if (!is.null(dat$catch_obs) && nrow(dat$catch_obs)) ggplot2::geom_point(
      data = dat$catch_obs,
      ggplot2::aes(x = time, y = catch),
      color = dot_blue, shape = 16, size = 3, na.rm = TRUE
    ) } +
    { if (show_CIs && !is.null(base_C_t0) && nrow(base_C_t0)) list(
      ggplot2::geom_line(
        data = base_C_t0, ggplot2::aes(x = time, y = lwr),
        color = spict_blue_ci_line, linetype = "dashed", linewidth = 0.7, alpha = 0.95, inherit.aes = FALSE
      ),
      ggplot2::geom_line(
        data = base_C_t0, ggplot2::aes(x = time, y = upr),
        color = spict_blue_ci_line, linetype = "dashed", linewidth = 0.7, alpha = 0.95, inherit.aes = FALSE
      )
    ) } +
    { if (!is.null(base_C_t0) && nrow(base_C_t0)) {
      if (has_man || !is.finite(obs_end)) {
        ggplot2::geom_line(
          data = base_C_t0, ggplot2::aes(x = time, y = est),
          color = spict_blue_mean, linewidth = 0.5, inherit.aes = FALSE
        )
      } else {
        list(
          ggplot2::geom_line(
            data = subset(base_C_t0, time <= obs_end),
            ggplot2::aes(x = time, y = est),
            color = spict_blue_mean, linewidth = 0.5, inherit.aes = FALSE,
          ),
          ggplot2::geom_line(
            data = subset(base_C_t0, time >= obs_end),
            ggplot2::aes(x = time, y = est),
            color = spict_blue_mean, linewidth = 0.5, inherit.aes = FALSE,
            linetype = "dotted"
          )
        )
      }
    }
    } +
    ggplot2::labs(x = "Year", y = "Catch (tons)") +
    { if (length(scs_final)) ggplot2::scale_color_manual(values = cols, guide = if (show_legend) "legend" else "none")
      else ggplot2::scale_color_discrete(guide = "none") } +
    ggplot2::scale_fill_manual(values = if (is.null(cols)) NA else cols, guide = "none") +
    theme_minimal_compact2_good_local()

  pC
}

# ===================================================================
# SMALL INLINED HELPERS USED BY BIOMASS/F PANELS
# ===================================================================

#' Minimal compact ggplot2 theme (ELU2)
#' @keywords internal
#' @noRd
.spict_theme_minimal_compact2 <- function(base_size = 10, base_family = "") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12),
      axis.title = ggplot2::element_text(face = "bold", size = 12),
      axis.text  = ggplot2::element_text(size = 10, face = "bold"),
      legend.position = "none",
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(fill = NA, colour = "grey35", linewidth = 2),
      axis.ticks = ggplot2::element_line(linewidth = 0.5, color = "grey35"),
      axis.ticks.length = grid::unit(3, "pt"),
      strip.background = ggplot2::element_rect(fill = "grey35", color = "grey35", linewidth = 0.5),
      strip.text = ggplot2::element_text(face = "bold", size = ggplot2::rel(1)),
      text = ggplot2::element_text(face = "bold", size = 10),
      plot.margin = ggplot2::margin(3, 3, 3, 3)
    )
}

#' End of observation time across components
#' @keywords internal
#' @noRd
.spict_obs_end_overall <- function(model) {
  inp <- model$inp
  tr <- inp$timerangeObs
  if (!is.null(tr) && length(tr) >= 2 && is.finite(tr[2])) return(tr[2])
  if (!is.null(inp$timeC) && length(inp$timeC))       return(utils::tail(inp$timeC, 1))
  if (!is.null(inp$timeI) && length(inp$timeI) && length(inp$timeI[[1]]) > 0) return(utils::tail(inp$timeI[[1]], 1))
  if (!is.null(inp$time)  && length(inp$time))        return(max(inp$time, na.rm = TRUE))
  NA_real_
}

#' Build a ribbon dataframe from lwr/upr vectors
#' @keywords internal
#' @noRd
.spict_make_ribbon_df <- function(time, lwr, upr) {
  ok <- is.finite(time) & is.finite(lwr) & is.finite(upr)
  data.frame(time = time[ok], ymin = lwr[ok], ymax = upr[ok])
}

#' q-scaled index points for first two indices (if present)
#' @keywords internal
#' @noRd
.spict_index_points <- function(model) {
  out <- list()
  inp <- model$inp
  if (!is.null(inp$timeI) && length(inp$timeI) >= 1) {
    qest <- try(get.par("logq", model, exp = TRUE), silent = TRUE)
    if (!inherits(qest, "try-error")) {
      for (i in seq_along(inp$timeI)) {
        qrow <- if (!is.null(inp$mapq) && length(inp$mapq) >= i) inp$mapq[i] else i
        qfac <- if (is.finite(qest[qrow, 2])) qest[qrow, 2] else 1
        out[[i]] <- data.frame(
          time = inp$timeI[[i]],
          obs  = inp$obsI[[i]] / qfac,
          idx  = i
        )
      }
    }
  }
  out
}

#' Prefer rownames(time) else fallback to inp$time
#' @keywords internal
#' @noRd
.spict_time_from_par <- function(model, par_matrix) {
  rn <- suppressWarnings(as.numeric(rownames(par_matrix)))
  if (length(rn) == nrow(par_matrix) && !all(is.na(rn))) return(rn)
  t_full <- as.numeric(model$inp$time)
  if (length(t_full) >= nrow(par_matrix)) return(t_full[seq_len(nrow(par_matrix))])
  return(rep_len(t_full, nrow(par_matrix)))
}

# ===================================================================
# BIOMASS (ABSOLUTE) PANEL
# ===================================================================

#' Absolute biomass panel (managed style)
#' @export
my_plot_manage_biomass_panel <- function(rep_man,
                                         scenario_color = NULL,
                                         show_CIs = TRUE,
                                         CI = 0.95,
                                         show_legend = TRUE) {
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

  Best <- get.par("logB",     rep_man, exp = TRUE, CI = CI)
  BB   <- get.par("logBBmsy", rep_man, exp = TRUE, CI = CI)

  t_B <- .spict_time_from_par(rep_man, Best)
  df_B <- data.frame(time = t_B, lwr = Best[,1], est = Best[,2], upr = Best[,3])

  inp <- rep_man$inp
  ind_in <- inp$indest
  if (length(ind_in) > 0) ind_in <- ind_in[seq_len(length(ind_in) - 1)]
  ind_pr <- if (!length(scs_final)) inp$indpred else integer(0)

  df_B_in <- if (length(ind_in)) df_B[ind_in, , drop = FALSE] else df_B[0, ]
  df_B_pr <- if (length(ind_pr)) df_B[ind_pr, , drop = FALSE] else df_B[0, ]

  repmax   <- if (exists("get.manmax")) get.manmax(rep_man) else rep_man
  Bmsy_all <- get.par("logBmsy", repmax, exp = TRUE, CI = CI)
  Bmsyvec  <- get.msyvec(repmax$inp, Bmsy_all)
  Bmsy     <- if (!is.null(nrow(Bmsy_all))) Bmsy_all[1, ] else Bmsy_all

  df_Bmsy_band <- data.frame(time = repmax$inp$time, ymin = Bmsyvec$ll, ymax = Bmsyvec$ul)
  df_Bmsy_line <- data.frame(time = repmax$inp$time, y = Bmsyvec$msy)

  Bmsy_est <- Bmsy[2]
  df_BB_rib <- .spict_make_ribbon_df(
    time = .spict_time_from_par(rep_man, BB),
    lwr  = BB[, 1] * Bmsy_est,
    upr  = BB[, 3] * Bmsy_est
  )

  df_scen_B <- NULL
  if (length(scs_final)) {
    lst <- vector("list", length(scs_final))
    for (i in seq_along(scs_final)) {
      sc <- scs_final[i]; rp <- rep_man$man[[sc]]
      Bi <- get.par("logB", rp, exp = TRUE, CI = CI)
      ti <- .spict_time_from_par(rp, Bi)
      lst[[i]] <- data.frame(time = ti, est = Bi[,2], scenario = sc)
    }
    df_scen_B <- do.call(rbind, lst)
    if (!is.null(df_scen_B) && nrow(df_scen_B) && length(scs_final)) {
      df_scen_B$scenario <- factor(df_scen_B$scenario, levels = scs_final)
    }
  }

  spict_blue_mean    <- "#0000FF"
  spict_blue_ci_line <- grDevices::rgb(0, 0, 1, 0.20)
  spict_blue_ci_fill <- grDevices::rgb(0, 0, 1, 0.10)

  p <- ggplot2::ggplot()

  if (isTRUE(show_CIs) && nrow(df_Bmsy_band)) {
    p <- p + ggplot2::geom_ribbon(
      data = df_Bmsy_band,
      ggplot2::aes(x = time, ymin = pmax(ymin, 0), ymax = pmax(ymax, 0)),
      fill = "lightgray", colour = NA
    )
  }
  p <- p + ggplot2::geom_line(
    data = df_Bmsy_line, ggplot2::aes(x = time, y = y),
    color = "black", linewidth = 0.7
  )

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
        data = transform(df_BB_rib, y = ymax), ggplot2::aes(x = time, y = y),
        color = spict_blue_ci_line, linewidth = 0.6, linetype = "solid"
      )
  }

  if (length(scs_final)) {
    p <- p + add_management_vlines_BF_good(
      rep_man, color = "grey75",
      linetype = "solid", linewidth = 0.4, lineend = "butt"
    )
  } else {
    obs_end <- .spict_obs_end_overall(rep_man)
    if (is.finite(obs_end)) {
      p <- p + ggplot2::geom_vline(
        xintercept = obs_end, color = "grey75",
        linetype = "solid", linewidth = 0.4
      )
    }
  }

  if (!is.null(df_scen_B) && nrow(df_scen_B)) {
    p <- p + ggplot2::geom_line(
      data = df_scen_B, ggplot2::aes(x = time, y = est, color = scenario),
      linewidth = 0.8, na.rm = TRUE
    )
  }

  if (nrow(df_B_in)) {
    if (isTRUE(show_CIs)) {
      p <- p +
        ggplot2::geom_line(
          data = df_B_in, ggplot2::aes(x = time, y = lwr),
          linetype = "dashed", linewidth = 0.6,
          colour = spict_blue_mean, inherit.aes = FALSE, show.legend = FALSE
        ) +
        ggplot2::geom_line(
          data = df_B_in, ggplot2::aes(x = time, y = upr),
          linetype = "dashed", linewidth = 0.6,
          colour = spict_blue_mean, inherit.aes = FALSE, show.legend = FALSE
        )
    }
    p <- p + ggplot2::geom_line(
      data = df_B_in, ggplot2::aes(x = time, y = est),
      linewidth = 0.9,
      colour = spict_blue_mean, inherit.aes = FALSE, show.legend = FALSE
    )
  }

  if (!length(scs_final) && nrow(df_B_pr)) {
    if (isTRUE(show_CIs)) {
      p <- p +
        ggplot2::geom_line(
          data = df_B_pr, ggplot2::aes(x = time, y = lwr),
          linetype = "dashed", linewidth = 0.6,
          colour = spict_blue_mean, inherit.aes = FALSE, show.legend = FALSE
        ) +
        ggplot2::geom_line(
          data = df_B_pr, ggplot2::aes(x = time, y = upr),
          linetype = "dashed", linewidth = 0.6,
          colour = spict_blue_mean, inherit.aes = FALSE, show.legend = FALSE
        )
    }
    p <- p + ggplot2::geom_line(
      data = df_B_pr, ggplot2::aes(x = time, y = est),
      linetype = "dotted", linewidth = 0.8,
      colour = spict_blue_mean, inherit.aes = FALSE, show.legend = FALSE
    )
  }

  idx <- .spict_index_points(rep_man)
  if (length(idx) >= 1) p <- p + ggplot2::geom_point(
    data = idx[[1]], ggplot2::aes(x = time, y = obs),
    color = "blue", shape = 16, size = 2, inherit.aes = FALSE
  )
  if (length(idx) >= 2) p <- p + ggplot2::geom_point(
    data = idx[[2]], ggplot2::aes(x = time, y = obs),
    shape = 22, color = "black", fill = "green", size = 2, stroke = 0.5, inherit.aes = FALSE
  )

  p +
    ggplot2::labs(title = "Absolute biomass", x = "Year", y = expression(B[t])) +
    { if (length(scs_final)) ggplot2::scale_color_manual(
      values = cols, breaks = scs_final, labels = lab_final, name = NULL,
      guide = if (show_legend) "legend" else "none"
    ) else ggplot2::scale_color_discrete(guide = "none") } +
    ggplot2::scale_fill_manual(values = if (is.null(cols)) NA else cols, guide = "none") +
    ggplot2::scale_y_continuous(
      sec.axis = ggplot2::sec_axis(~ . / Bmsy_est, name = expression(B[t]/B[MSY]))
    ) +
    theme_minimal_compact2_good_local()
}

# ===================================================================
# ABSOLUTE F[t] PANEL
# ===================================================================

#' Absolute F[t] panel (with Fmsy band and optional scenarios)
#' @export
my_plot_manage_f_panel <- function(rep_man,
                                   scenario_color = NULL,
                                   show_CIs = TRUE,
                                   CI = 0.95,
                                   show_legend = FALSE) {
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

  vline_col  <- "grey75"; vline_size <- 0.4
  manflag <- ("man" %in% names(rep_man)) && length(rep_man$man) > 0
  inp <- rep_man$inp

  repmax  <- if (exists("get.manmax")) get.manmax(rep_man) else rep_man
  tvgflag <- isTRUE(repmax$inp$timevaryinggrowth) || isTRUE(repmax$inp$logmcovflag)

  # Absolute Ft
  Fest <- get.par("logFnotS", rep_man, exp = TRUE, CI = CI)
  dfF <- data.frame(
    time = .spict_time_from_par(rep_man, Fest),
    lwr  = Fest[, 1], est = Fest[, 2], upr  = Fest[, 3]
  )

  ind_in <- inp$indest
  if (length(ind_in) > 0) ind_in <- ind_in[seq_len(length(ind_in) - 1)]
  ind_pr <- if (!manflag) inp$indpred else integer(0)

  dfF_in <- if (length(ind_in)) dfF[ind_in, , drop = FALSE] else dfF[0, ]
  dfF_pr <- if (length(ind_pr)) dfF[ind_pr, , drop = FALSE] else dfF[0, ]

  # Fmsy band + secondary axis
  if (tvgflag) {
    Fmsy_all <- get.par("logFmsyvec", repmax, exp = TRUE, CI = CI)
    t_Fmsy   <- .spict_time_from_par(repmax, Fmsy_all)
    Fmsy_band <- data.frame(time = t_Fmsy, ymin = Fmsy_all[,1], ymax = Fmsy_all[,3])
    Fmsy_line <- data.frame(time = t_Fmsy, y = Fmsy_all[,2])
    fmsy_sec  <- ggplot2::waiver()
    fmsy_scale <- approx(x = t_Fmsy, y = Fmsy_all[,2], xout = dfF$time, rule = 2)$y
    rib_pre <- .spict_make_ribbon_df(dfF$time, dfF$lwr, dfF$upr)
  } else {
    Fmsy_one <- get.par("logFmsy", repmax, exp = TRUE, CI = CI)
    t_Fmsy   <- repmax$inp$time
    Fmsy_band <- data.frame(time = t_Fmsy, ymin = rep(Fmsy_one[1], length(t_Fmsy)),
                            ymax = rep(Fmsy_one[3], length(t_Fmsy)))
    Fmsy_line <- data.frame(time = t_Fmsy, y = rep(Fmsy_one[2], length(t_Fmsy)))
    fmsy_sec <- ggplot2::sec_axis(~ . / Fmsy_one[2], name = expression(F[t]/F[MSY]))
    fmsy_scale <- rep(Fmsy_one[2], length(dfF$time))
    rib_pre <- .spict_make_ribbon_df(dfF$time, dfF$lwr, dfF$upr)
  }

  spict_blue_mean    <- "blue"
  spict_blue_ci_line <- grDevices::rgb(0, 0, 1, 0.20)
  spict_blue_ci_fill <- grDevices::rgb(0, 0, 1, 0.10)

  # Scenario colours
  sc_order <- c("currentCatch","currentF","Fmsy","noF",
                "reduceF25","increaseF25","msyHockeyStick","ices")
  sc_present <- if (manflag) names(rep_man$man) else character(0)
  sc_core   <- intersect(sc_order, sc_present)
  sc_other  <- setdiff(sc_present, sc_core)
  scs_final <- c(sc_core, sort(sc_other))
  if (is.null(scenario_color)) {
    cols <- if (length(scs_final)) { x <- man_cols(length(scs_final)); names(x) <- scs_final; x } else NULL
  } else {
    cols <- scenario_color
    if (length(scs_final)) {
      if (is.null(names(cols))) names(cols) <- scs_final
      cols <- cols[match(scs_final, names(cols))]; names(cols) <- scs_final
    }
  }

  # Build plot
  p <- ggplot2::ggplot()

  # Fmsy band + line
  if (nrow(Fmsy_band)) {
    p <- p + ggplot2::geom_ribbon(
      data = Fmsy_band,
      ggplot2::aes(x = time, ymin = pmax(ymin, 0), ymax = pmax(ymax, 0)),
      fill = "lightgray", colour = NA
    )
  }
  if (nrow(Fmsy_line)) {
    p <- p + ggplot2::geom_line(
      data = Fmsy_line, ggplot2::aes(x = time, y = y),
      color = "black", linewidth = 0.7
    )
  }

  # Vlines
  if (manflag) {
    p <- p + add_management_vlines_BF_good(rep_man, color = vline_col,
                                           linetype = "solid", linewidth = vline_size, lineend = "butt")
  } else {
    obs_end <- .spict_obs_end_overall(rep_man)
    if (is.finite(obs_end)) {
      p <- p + ggplot2::geom_vline(xintercept = obs_end, color = vline_col,
                                   linetype = "solid", linewidth = vline_size)
    }
  }

  # In-sample Ft: ribbon edges + mean
  if (isTRUE(show_CIs) && nrow(rib_pre)) {
    p <- p +
      ggplot2::geom_ribbon(
        data = data.frame(time = rib_pre$time, ymin = rib_pre$ymin, ymax = rib_pre$ymax),
        ggplot2::aes(x = time, ymin = ymin, ymax = ymax),
        fill = spict_blue_ci_fill, colour = NA
      ) +
      ggplot2::geom_line(
        data = transform(rib_pre, y = ymin), ggplot2::aes(x = time, y = y),
        color = spict_blue_ci_line, linewidth = 0.6, linetype = "solid"
      ) +
      ggplot2::geom_line(
        data = transform(rib_pre, y = ymax), ggplot2::aes(x = time, y = y),
        color = spict_blue_ci_line, linewidth = 0.6, linetype = "solid"
      )
  }

  if (nrow(dfF_in)) {
    p <- p + ggplot2::geom_line(
      data = dfF_in, ggplot2::aes(x = time, y = est),
      colour = spict_blue_mean, linewidth = 0.9
    )
  }

  # Scenarios (absolute Ft) if available
  if (manflag) {
    lst <- vector("list", length(scs_final))
    for (i in seq_along(scs_final)) {
      sc <- scs_final[i]; rp <- rep_man$man[[sc]]
      Fi <- get.par("logFnotS", rp, exp = TRUE, CI = CI)
      ti <- .spict_time_from_par(rp, Fi)
      lst[[i]] <- data.frame(time = ti, est = Fi[,2], scenario = sc)
    }
    df_scen_F <- do.call(rbind, lst)
    if (!is.null(df_scen_F) && nrow(df_scen_F) && length(scs_final)) {
      df_scen_F$scenario <- factor(df_scen_F$scenario, levels = scs_final)
    }
    p <- p + ggplot2::geom_line(
      data = df_scen_F, ggplot2::aes(x = time, y = est, color = scenario),
      linewidth = 0.8, na.rm = TRUE
    )
  }

  # One-ahead (raw fit only): dotted mean + dashed CI edges
  if (!manflag && nrow(dfF_pr)) {
    if (isTRUE(show_CIs)) {
      p <- p +
        ggplot2::geom_line(
          data = dfF_pr, ggplot2::aes(x = time, y = lwr),
          linetype = "dashed", linewidth = 0.6, colour = spict_blue_mean
        ) +
        ggplot2::geom_line(
          data = dfF_pr, ggplot2::aes(x = time, y = upr),
          linetype = "dashed", linewidth = 0.6, colour = spict_blue_mean
        )
    }
    p <- p + ggplot2::geom_line(
      data = dfF_pr, ggplot2::aes(x = time, y = est),
      linetype = "dotted", linewidth = 0.8, colour = spict_blue_mean
    )
  }

  # Axis + theme
  p +
    ggplot2::labs(x = "Year", y = expression(F[t])) +
    ggplot2::scale_y_continuous(sec.axis = fmsy_sec) +
    { if (manflag) ggplot2::scale_color_manual(values = cols, guide = if (show_legend) "legend" else "none")
      else ggplot2::scale_color_discrete(guide = "none") } +
    theme_minimal_compact2_good_local()
}
