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
my_plot_kobe_kita <- function(rep,
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

  # ----------------------- Local helpers ---------------------- #
  `%||%` <- function(a, b) if (!is.null(a)) a else b

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

  lab_rel_1 <- function(vals) {
    out <- fmt1(vals)
    is_one <- abs(vals - 1) < 1e-9
    out[is_one] <- "1"
    out
  }

  clamp <- function(val, lo, hi) pmin(pmax(val, lo), hi)
  clamp_in <- function(v, lo, hi, eps) pmin(pmax(v, lo + eps), hi - eps)

  # ------------------ Validate & derive quantities --------------- #
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

  if (isTRUE(logax)) {
    if (!is.finite(xlim[1]) || xlim[1] <= 0) xlim[1] <- logminval
    if (!is.finite(ylim[1]) || ylim[1] <= 0) ylim[1] <- logminval
    if (xlim[2] <= xlim[1]) xlim[2] <- xlim[1] * 10
    if (ylim[2] <= ylim[1]) ylim[2] <- ylim[1] * 10
  }

  # ------------------------- Data frames for ggplot ------------------------- #
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
      estF <- get_par("logF", rp, exp = TRUE)[, 2]
      estB <- get_par("logB", rp, exp = TRUE)[, 2]
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
