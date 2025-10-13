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
#'   Useful when the dynamic range is large. Defaults to `FALSE`.
#' @param plot.legend Logical; if `TRUE`, draws a small legend label for
#'   \eqn{E(B_\infty)} in the top-right when available (not for sub-annual `dtc`).
#'   Defaults to `TRUE`.
#' @param man.legend Logical; if `TRUE`, shows a colour legend for management
#'   scenarios when present. Defaults to `TRUE`. (Internally the plot ultimately
#'   hides the legend for a clean Kobe panel; keep `TRUE` to ensure consistent
#'   colour mapping in downstream compositions.)
#' @param ext Logical; if `TRUE` and `rel.axes = FALSE`, adds secondary axes
#'   for the relative scales (\eqn{B_t/B_{MSY}}, \eqn{F_t/F_{MSY}}). Defaults to `TRUE`.
#' @param rel.axes Logical; if `TRUE`, the primary axes are the relative scales
#'   \eqn{B_t/B_{MSY}} and \eqn{F_t/F_{MSY}} (forced when time-varying growth is on).
#'   Defaults to `FALSE`.
#' @param xlim Numeric length-2 vector with x-axis limits. If `NULL`, limits are
#'   computed from data, the RP ellipse, and (if available) \eqn{E(B_\infty)}.
#' @param ylim Numeric length-2 vector with y-axis limits. If `NULL`, limits are
#'   computed from data and the RP ellipse.
#' @param labpos Length-2 numeric giving label anchor for first/last time labels.
#'   (Reserved; kept for compatibility.) Default `c(1, 1)`.
#' @param xlabel Optional x-axis label. If `NULL`, a math expression is used:
#'   `B[t]` (with catch unit appended) or `B[t]/B[MSY]` on relative axes.
#' @param stamp Optional small text (e.g., version string) stamped at the
#'   bottom-right margin of the plot. Default `NULL`.
#' @param verbose Logical; if `TRUE`, enables warnings for parameter extraction
#'   fallbacks in helper routines. Default `TRUE`.
#' @param CI Confidence level in (0, 1) used when constructing parameter
#'   intervals via internal extractors. Default `0.95`.
#'
#' @details
#' Internally, the function:
#' \enumerate{
#'   \item Validates `rep` and enforces relative axes when time-varying growth
#'         or log-`m` covariates are active.
#'   \item Extracts \eqn{B_{MSY}} and \eqn{F_{MSY}} (final or time-varying),
#'         \eqn{B_t} and \eqn{F_t}, and constructs an RP ellipse if the
#'         \code{ellipse} package is available (falls back gracefully otherwise).
#'   \item Optionally aggregates to annual steps when `dtc < 1`.
#'   \item Builds the phase rectangles (green/yellow/red) and overlays:
#'         historical path, one-step prediction, \eqn{E(B_\infty)} segment
#'         (if computable), management scenario paths (with dashed pre-interval
#'         segments and solid post-interval segments), and MSY markers.
#'   \item Styles axes (linear/log10), optionally adds secondary relative axes,
#'         and formats tick labels to emphasise the value 1 on relative scales.
#' }
#'
#' Management scenarios (if present in `rep$man`) are drawn with a consistent,
#' repeatable palette and ordered using a canonical scenario name order when possible.
#'
#' @return
#' Invisibly returns the `ggplot` object after constructing the Kobe panel
#' (so it can be added to patchwork layouts, saved, etc.).
#'
#' @section Dependencies:
#' Requires **ggplot2**. If available, **ellipse** is used to draw the
#' reference point ellipse; otherwise a simple fallback is used.
#'
#' @seealso
#' [spict::manage()], [spict::fit.spict()]
#'
#' @examples
#' \dontrun{
#' # Assume `fit` is a fitted SPiCT object (class 'spictcls'):
#' p <- my_plot_kobe_all_management_scenario(fit)
#' print(p)
#'
#' # With management scenarios (after spict::manage()):
#' # fit$man <- spict::manage(fit, scenarios = c(1,2,3))
#' my_plot_kobe_all_management_scenario(fit, rel.axes = TRUE)
#' }
#'
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
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  ## ----------------------- Local helpers ---------------------- ##
  check_rep <- function(rep, reportmode0 = TRUE) {
    if (!inherits(rep, "spictcls") || !"opt" %in% names(rep)) {
      stop("The argument 'rep' must be a fitted spict object (fit.spict()).")
    } else if (reportmode0 && rep$inp$reportmode != 0) {
      stop("All states must be reported! Set 'inp$reportmode <- 0' and refit.")
    }
    TRUE
  }

  get_par <- function(parname, rep, exp = FALSE, CI = 0.95) {
    if (CI <= 0 || CI >= 1) stop("CI must be in (0,1).")
    z <- qnorm(CI + (1 - CI) / 2)

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
      sdv <- sqrt(diag(rep$cov.fixed))[ind_fix]
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
        sdv <- rep(NA_real_, length(est))
        ll  <- rep(NA_real_, length(est))
        ul  <- rep(NA_real_, length(est))
      } else if ("phases" %in% names(rep$inp) &&
                 parname %in% names(rep$inp$phases) &&
                 rep$inp$phases[[parname]] == -1) {
        est <- rep$inp$parlist[[parname]]
        sdv <- rep(0, length(est))
        ll  <- est; ul <- est
      } else if (!is.na(parname) && identical(parname, "P")) {
        B <- get_par("logB",     rep, exp = TRUE, CI = CI)
        C <- get_par("logCpred", rep, exp = TRUE, CI = CI)
        ic <- rep$inp$ic; nc <- rep$inp$nc
        B0 <- B[ic, 2]; B1 <- B[ic + nc, 2]
        T0 <- rep$inp$time[ic]; T1 <- rep$inp$time[ic + nc]
        est <- (B1 - B0 + C[, 2]) / (T1 - T0)
        sdv <- rep(NA_real_, length(est))
        ll  <- rep(NA_real_, length(est))
        ul  <- rep(NA_real_, length(est))
      } else {
        if (verbose) warning("get_par: could not extract '", parname, "'. Returning NA.")
        est <- NA_real_; sdv <- NA_real_; ll <- NA_real_; ul <- NA_real_
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
      rownames(out) <- rep$inp$time
    }
    out
  }

  add_catchunit <- function(lab, cu) {
    cu <- as.character(cu)
    if (nzchar(cu)) eval(bquote(.(lab[[1]]) * ',' ~ .(cu))) else lab
  }

  make_rp_ellipse <- function(rep) {
    inds <- c(max(which(names(rep$value) == "logBmsy")),
              max(which(names(rep$value) == "logFmsy")))
    sds  <- rep$sd[inds]
    s_sum <- rep$sd[which(names(rep$value) == "logBmsyPluslogFmsy")]
    cova <- (s_sum^2 - sds[1]^2 - sds[2]^2) / 2
    covBF <- matrix(c(sds[1]^2, cova, cova, sds[2]^2), 2, 2, byrow = TRUE)
    parBF <- rep$value[inds]

    if (requireNamespace("ellipse", quietly = TRUE)) {
      ellipse::ellipse(cov2cor(covBF)[1, 2],
                       scale = sqrt(diag(covBF)),
                       centre = parBF, npoints = 300)
    } else {
      matrix(c(parBF[1], parBF[2]), ncol = 2,
             dimnames = list(NULL, c("x", "y")))
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

  man_cols <- function() {
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
  rel_breaks_fun_x <- function(primary_lims) {
    pr <- scales::pretty_breaks(n = 6)(primary_lims)
    sort(unique(c(pr / (Bmsy[2] / 1), 1)))
  }
  rel_breaks_fun_y <- function(primary_lims) {
    pr <- scales::pretty_breaks(n = 6)(primary_lims)
    sort(unique(c(pr / (Fmsy[2] / 1), 1)))
  }

  clamp <- function(val, lo, hi) pmin(pmax(val, lo), hi)

  ## ------------------ Validate & derive quantities --------------- ##
  check_rep(rep)

  tvgflag <- rep$inp$timevaryinggrowth | rep$inp$logmcovflag
  if (tvgflag) rel.axes <- TRUE

  Bmsy_all <- get_par("logBmsy", rep, exp = TRUE, CI = CI)
  Fmsy_all <- get_par("logFmsy", rep, exp = TRUE, CI = CI)
  Bmsy <- tail(Bmsy_all, 1L)
  Fmsy <- tail(Fmsy_all, 1L)

  if (rel.axes) {
    ext <- FALSE
    bscal <- Bmsy[2]
    fscal <- Fmsy[2]
    xlab_expr <- expression(B[t]/B[MSY])
    ylab_expr <- expression(F[t]/F[MSY])
  } else {
    bscal <- 1
    fscal <- 1
    xlab_expr <- expression(B[t])
    xlab_expr <- add_catchunit(xlab_expr, rep$inp$catchunit)
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
    cl <- matrix(c(log(Bmsy[2]), log(Fmsy[2])), ncol = 2)
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

  Fl <- tail(unname(fff), 1)
  Bl <- tail(unname(bbb), 1)
  EBinf <- get_EBinf(rep) / bscal

  if (is.null(xlim)) {
    xlim <- range(c(exp(cl[, 1]), Best[, 2], EBinf) / bscal, na.rm = TRUE)
    if (min(rep$inp$dtc) < 1)
      xlim <- range(c(exp(alb$annvec), exp(cl[, 1]), EBinf) / bscal, na.rm = TRUE)
    xlim[2] <- min(c(xlim[2], 8 * Bmsy[2] / bscal), 2.2 * max(bbb), na.rm = TRUE)
    xlim[2] <- max(c(xlim[2], Bmsy[2] / bscal), na.rm = TRUE)
  }
  if (is.null(ylim)) {
    ylim <- range(c(exp(cl[, 2]), Fest[, 2]) / fscal, na.rm = TRUE)
    if (min(rep$inp$dtc) < 1)
      ylim <- range(c(exp(logFest[, 2]) / fscal, exp(cl[, 2]) / fscal), na.rm = TRUE)
    ylim[2] <- min(c(ylim[2], 8 * Fmsy[2] / fscal), 2.2 * max(fff), na.rm = TRUE)
    ylim[2] <- max(c(ylim[2], Fmsy[2] / fscal), na.rm = TRUE)
    if ("man" %in% names(rep)) ylim <- range(ylim, 0)
  }
  logminval <- 1e-4
  if (isTRUE(logax)) {
    if (xlim[1] < logminval) xlim[1] <- logminval
    if (ylim[1] < logminval) ylim[1] <- logminval
  }

  ## Ensure a little gap above the x-axis on linear scales
  if (!isTRUE(logax)) {
    if (is.finite(ylim[2]) && (ylim[1] >= 0 || !is.finite(ylim[1]))) {
      ylim[1] <- -max(1e-6, 0.02 * (ylim[2] - 0))
    }
  }

  ## ---- Slight padding of hard limits (keeps fills to the frame) ----
  pad_frac <- 0.02
  xpad <- pad_frac * diff(xlim)
  ypad <- pad_frac * diff(ylim)

  xlim <- c(
    max(if (isTRUE(logax)) logminval else 0, xlim[1] - xpad),
    xlim[2] + xpad
  )
  ylim <- c(
    max(if (isTRUE(logax)) logminval else 0, ylim[1] - ypad),
    ylim[2] + ypad
  )

  ## ------------------------- Data frames for ggplot ------------------------- ##
  xminq <- xlim[1]; xmaxq <- xlim[2]
  yminq <- ylim[1]; ymaxq <- ylim[2]

  bx <- if (rel.axes) 1 else Bmsy[2] / bscal
  bx <- clamp(bx, xminq, xmaxq)
  fy <- if (rel.axes) 1 else Fmsy[2] / fscal
  fy <- clamp(fy, yminq, ymaxq)

  df_green_rect   <- data.frame(xmin = bx,   xmax =  Inf, ymin = -Inf, ymax =  fy)
  df_yellowL_rect <- data.frame(xmin = -Inf, xmax =  bx,  ymin = -Inf, ymax =  fy)
  df_yellowT_rect <- data.frame(xmin = bx,   xmax =  Inf, ymin =  fy,  ymax =  Inf)
  df_red_rect     <- data.frame(xmin = -Inf, xmax =  bx,  ymin =  fy,  ymax =  Inf)

  df_ellipse <- if (ncol(as.matrix(cl)) == 2 && nrow(as.matrix(cl)) > 1) {
    data.frame(x = exp(cl[, 1]) / (if (rel.axes) Bmsy[2] else bscal),
               y = exp(cl[, 2]) / (if (rel.axes) Fmsy[2] else fscal))
  } else {
    data.frame(x = exp(cl[1]) / (if (rel.axes) Bmsy[2] else bscal),
               y = exp(cl[2]) / (if (rel.axes) Fmsy[2] else fscal))
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
    xpred <- Best[rep$inp$indpred, 2] / (if (rel.axes) Bmsy[2] else bscal)
    ypred <- Fest[rep$inp$indpred, 2] / (if (rel.axes) Fmsy[2] else fscal)
    df_pred <- data.frame(x = xpred, y = ypred)
    Bll <- tail(xpred, 1)
    Fll <- tail(ypred, 1)
    df_EBseg <- data.frame(x = Bll, xend = EBinf, y = Fll, yend = Fll)
  }

  df_first <- data.frame(x = bbb[1], y = fff[1], lab = round(fbtime[1], 2))
  df_last  <- data.frame(x = tail(bbb,1), y = tail(fff,1), lab = round(tail(fbtime, 1), 2))

  nr <- length(rep$inp$ini$logr)
  df_msy_prev <- df_msy_curr <- NULL
  if (nr > 1) {
    df_msy_prev <- data.frame(x = Bmsy_all[1:(nr-1), 2] / (if (rel.axes) Bmsy[2] else bscal),
                              y = Fmsy_all[1:(nr-1), 2] / (if (rel.axes) Fmsy[2] else fscal))
    df_msy_curr <- data.frame(x = Bmsy_all[nr, 2] / (if (rel.axes) Bmsy[2] else bscal),
                              y = Fmsy_all[nr, 2] / (if (rel.axes) Fmsy[2] else fscal))
  }

  df_true <- NULL
  if ("true" %in% names(rep$inp)) {
    df_true <- data.frame(x = rep$inp$true$Bmsy / (if (rel.axes) Bmsy[2] else bscal),
                          y = rep$inp$true$Fmsy / (if (rel.axes) Fmsy[2] else fscal))
  }

  ## ========== NEW: canonical scenario order, labels, and jump segment ==========
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
  if ("man" %in% names(rep) && length(rep$man)) {
    nman <- length(rep$man)
    leg_man <- names(rep$man)
    if (is.null(leg_man) || !length(leg_man)) leg_man <- paste0("Scenario ", seq_len(nman))

    sc_present <- leg_man
    sc_core    <- intersect(sc_order, sc_present)
    sc_other   <- setdiff(sc_present, sc_core)
    scs_final  <- c(sc_core, sort(sc_other))
    lab_final  <- ifelse(scs_final %in% names(label_map),
                         unname(label_map[scs_final]),
                         scs_final)

    df_list      <- vector("list", nman)
    df_int_list  <- vector("list", nman)
    df_jump_list <- vector("list", nman)

    same_pt <- function(a, b) isTRUE(all.equal(a, b, tolerance = 1e-10))

    i <- 1L
    while (i <= nman) {
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
        bx <- exp(unname(alb2$annvec)) / (if (rel.axes) Bmsy[2] else bscal); ttA <- as.numeric(names(alb2$annvec))
        fy <- exp(unname(alf2$annvec)) / (if (rel.axes) Fmsy[2] else fscal)
        keep_man <- (ttA >= manint[1])
        keep_int <- (ttA >= lastobs) & (ttA < manint[1])
        pre_mask <- ttA < manint[1]
      } else {
        bx <- estB / (if (rel.axes) Bmsy[2] else bscal); ttA <- time
        fy <- estF / (if (rel.axes) Fmsy[2] else fscal)
        keep_man <- (time >= manint[1])
        keep_int <- (time >= lastobs) & (time < manint[1])
        pre_mask <- time < manint[1]
      }

      if (any(keep_man, na.rm = TRUE)) {
        maninds <- which(keep_man)
        df_list[[i]] <- data.frame(
          x = bx[maninds],
          y = fy[maninds],
          scenario = leg_man[i]
        )
      } else {
        df_list[[i]] <- data.frame(x = numeric(0), y = numeric(0), scenario = character(0))
      }

      if (any(keep_int, na.rm = TRUE)) {
        df_int_list[[i]] <- data.frame(
          x = bx[which(keep_int)],
          y = fy[which(keep_int)],
          scenario = leg_man[i]
        )
      } else {
        df_int_list[[i]] <- data.frame(x = numeric(0), y = numeric(0), scenario = character(0))
      }

      x0 <- if (any(pre_mask, na.rm = TRUE)) tail(bx[pre_mask], 1) else NA_real_
      y0 <- if (any(pre_mask, na.rm = TRUE)) tail(fy[pre_mask], 1) else NA_real_
      x1 <- if (any(keep_man, na.rm = TRUE)) bx[which(keep_man)[1]] else NA_real_
      y1 <- if (any(keep_man, na.rm = TRUE)) fy[which(keep_man)[1]] else NA_real_

      if (is.finite(x0) && is.finite(y0) && is.finite(x1) && is.finite(y1) &&
          (abs(x1 - x0) > 0 || abs(y1 - y0) > 0)) {
        df_jump_list[[i]] <- data.frame(
          x = x0, y = y0, xend = x1, yend = y1, scenario = leg_man[i]
        )
      } else {
        df_jump_list[[i]] <- data.frame(
          x = numeric(0), y = numeric(0), xend = numeric(0), yend = numeric(0),
          scenario = character(0)
        )
      }

      if (is.finite(x0) && is.finite(y0)) {
        if (!nrow(df_int_list[[i]]) ||
            !same_pt(tail(df_int_list[[i]]$x, 1), x0) ||
            !same_pt(tail(df_int_list[[i]]$y, 1), y0)) {
          df_int_list[[i]] <- rbind(df_int_list[[i]],
                                    data.frame(x = x0, y = y0, scenario = leg_man[i]))
        }
      }

      i <- i + 1L
    }

    df_man       <- do.call(rbind, df_list)
    df_man_int   <- do.call(rbind, df_int_list)
    df_man_jump  <- do.call(rbind, df_jump_list)

    for (dname in c("df_man","df_man_int","df_man_jump")) {
      dd <- get(dname)
      if (!is.null(dd) && nrow(dd)) {
        dd$scenario <- factor(dd$scenario, levels = scs_final)
        assign(dname, dd)
      }
    }

    colv <- man_cols()[seq_along(scs_final)]
    names(colv) <- scs_final
  } else {
    scs_final <- character(0)
    lab_final <- character(0)
    colv <- NULL
  }
  ## ===================== END new management settings block =====================

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

  # One-step prediction (path-or-point to avoid warnings)
  if (!is.null(df_pred)) {
    if (nrow(df_pred) > 1) {
      p <- p + ggplot2::geom_path(
        data = df_pred,
        ggplot2::aes(x = x, y = y),
        linetype = "33",
        color = grDevices::adjustcolor("blue", alpha.f = 0.8)
      )
    } else {
      p <- p + ggplot2::geom_point(
        data = df_pred,
        ggplot2::aes(x = x, y = y),
        shape = 16, size = 1.8,
        color = grDevices::adjustcolor("blue", alpha.f = 0.8)
      )
    }
  }

  if (!is.null(df_EBseg)) {
    p <- p + ggplot2::geom_segment(data = df_EBseg,
                                   ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
                                   color = "blue", linetype = "33", linewidth = 1.0)
    p <- p + ggplot2::geom_point(data = data.frame(x = df_EBseg$xend, y = df_EBseg$yend),
                                 ggplot2::aes(x = x, y = y),
                                 shape = 23, fill = "gold", color = "black", size = 3, stroke = 0.5)
  }

  # ----- Pre-management dashed path (path-or-point) -----
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

  # ----- Solid management path (path-or-point) -----
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

  # ----- Jump segment -----
  if (!is.null(df_man_jump) && nrow(df_man_jump)) {
    p <- p + ggplot2::geom_segment(data = df_man_jump,
                                   ggplot2::aes(x = x, y = y, xend = xend, yend = yend, color = scenario),
                                   linewidth = 0.7, inherit.aes = FALSE, show.legend = man.legend)
  }

  # MSY markers
  if (!is.null(df_msy_prev) && nrow(df_msy_prev)) {
    p <- p + ggplot2::geom_point(data = df_msy_prev, ggplot2::aes(x = x, y = y),
                                 shape = 3, color = "magenta", size = 2)
  }
  if (!is.null(df_msy_curr) && nrow(df_msy_curr)) {
    p <- p + ggplot2::geom_point(data = df_msy_curr, ggplot2::aes(x = x, y = y),
                                 shape = 3, color = "black", size = 2)
  }

  # True point
  if (!is.null(df_true)) {
    p <- p + ggplot2::geom_point(data = df_true, ggplot2::aes(x = x, y = y),
                                 shape = 25, fill = grDevices::adjustcolor(rgb(1, 165/255, 0), alpha.f = 0.7),
                                 color = "black", size = 3)
  }

  # First/last time labels
  p <- p + ggplot2::geom_point(data = df_first, ggplot2::aes(x = x, y = y),
                               shape = 21, fill = "white", color = "black", size = 2.5, stroke = 0.5)
  p <- p + ggplot2::geom_text(data = df_first, ggplot2::aes(x = x, y = y, label = lab),
                              vjust = -0.8, size = 3)
  p <- p + ggplot2::geom_point(data = df_last, ggplot2::aes(x = x, y = y),
                               shape = 22, fill = "white", color = "black", size = 2.5, stroke = 0.5)
  p <- p + ggplot2::geom_text(data = df_last, ggplot2::aes(x = x, y = y, label = lab),
                              vjust = -0.8, size = 3)

  # Reference vline at 0 only if inside limits
  if (!isTRUE(logax) && (0 >= xlim[1]) && (0 <= xlim[2])) {
    p <- p + ggplot2::geom_vline(xintercept = 0, color = "darkred", linetype = 2)
  }

  ## Small additional padding block (kept from your version)
  pad_frac <- 0.02
  xpad <- pad_frac * diff(xlim)
  ypad <- pad_frac * diff(ylim)
  xlim <- c(max(if (isTRUE(logax)) logminval else 0, xlim[1] - xpad), xlim[2] + xpad)
  ylim <- c(max(if (isTRUE(logax)) logminval else 0, ylim[1] - ypad), ylim[2] + ypad)

  # Labels & scales
  xlab_final <- if (is.null(xlabel)) xlab_expr else xlabel
  ylab_final <- ylab_expr

  if (!rel.axes && isTRUE(ext)) {
    sec_x <- ggplot2::dup_axis(name = expression(B[t]/B[MSY]),
                               labels = function(x) fmt1(x / (Bmsy[2] / 1)))
    sec_y <- ggplot2::dup_axis(name = expression(F[t]/F[MSY]),
                               labels = function(y) fmt1(y / (Fmsy[2] / 1)))
  } else {
    sec_x <- ggplot2::waiver()
    sec_y <- ggplot2::waiver()
  }

  if (isTRUE(logax)) {
    p <- p + ggplot2::scale_x_log10(
      limits = xlim,
      name   = xlab_final,
      labels = function(x) formatC(x, format = "f", digits = 0),
      expand = c(0, 0),
      sec.axis = if (!rel.axes && isTRUE(ext)) ggplot2::sec_axis(
        trans  = ~ . / (Bmsy[2] / 1),
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
        trans  = ~ . / (Fmsy[2] / 1),
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
        trans  = ~ . / (Bmsy[2] / 1),
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
        trans  = ~ . / (Fmsy[2] / 1),
        name   = expression(F[t]/F[MSY]),
        breaks = function(lims) { br <- scales::pretty_breaks(n = 6)(lims); sort(unique(c(br, 1))) },
        labels = lab_rel_1
      ) else ggplot2::waiver()
    )
  }

  # ===== Unified scenario scale (kept, but we'll hide legend below) =====
  if (!is.null(colv) && length(colv) && man.legend) {
    p <- p + ggplot2::scale_color_manual(values = colv,
                                         breaks = names(colv),
                                         labels = if (exists("lab_final") && length(lab_final)) lab_final else names(colv),
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
      axis.text.x  = element_text(size = 10, face = "bold"),
      axis.text.y  = element_text(size = 10, face = "bold"),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      plot.margin = margin(4, 4, 4, 4)
    )

  # Top-right mini-legend for E(B∞)
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

  # Non-convergence warning mark
  if (rep$opt$convergence != 0) {
    p <- p + ggplot2::annotate("text",
                               x = xlim[1] + 0.01 * diff(xlim),
                               y = ylim[2] + 0.02 * diff(ylim),
                               label = "▲", color = "black", size = 6)
    p <- p + ggplot2::annotate("text",
                               x = xlim[1] + 0.01 * diff(xlim),
                               y = ylim[2] + 0.02 * diff(ylim),
                               label = "!", color = "black", size = 3, vjust = 0.35)
  }

  # Version stamp
  if (!is.null(stamp) && nzchar(stamp)) {
    p <- p + ggplot2::annotate("text",
                               x = xlim[2] - 0.01 * diff(xlim),
                               y = ylim[1] - 0.04 * diff(ylim),
                               label = stamp, hjust = 1, vjust = 0, size = 3)
  }

  # --- Force-hide scenario legend for the Kobe plot ---
  p <- p + ggplot2::guides(color = "none")
  p <- p + ggplot2::theme(legend.position = "none")

  invisible(p)
}

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
#' Works on both fitted and managed SPiCT report objects.
#'
#' @param rep_main A fitted SPiCT report object (`spictcls`).
#'
#' @return
#' A data frame with columns `index`, `time`, and `obsrel`, or `NULL` if no
#' indices are present.
#'
#' @details
#' Uses `get.par("logq", ...)` and `get.par("logBmsy", ...)` to scale raw
#' observations in `rep_main$inp$obsI`. Currently supports up to two indices
#' (`Index1`, `Index2`) if available.
#'
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
#'
#' @description
#' Extracts B/Bmsy and F/Fmsy time series for the base fit and truncates them
#' at the last observed time (from `timerangeObs`).
#'
#' @param rep_man A fitted SPiCT report object (`spictcls`), with or without
#'   management scenarios.
#' @param CI Confidence level in (0, 1) to retrieve lower/upper bounds.
#'
#' @return
#' A list with elements:
#' \itemize{
#'   \item `t_last_obs` – numeric, last observed time,
#'   \item `BB` – data frame (`time`, `lwr`, `est`, `upr`),
#'   \item `FF` – data frame (`time`, `lwr`, `est`, `upr`).
#' }
#'
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
#'
#' @description
#' Creates tidy data frames for plotting B/Bmsy, F/Fmsy, and Catch panels.
#' If `rep_man$man` is not present, the function returns empty scenario frames
#' while still providing observed catch and observed indices (relative scale).
#'
#' @param rep_man A fitted SPiCT report object (`spictcls`), possibly containing
#'   a `man` list from [spict::manage()].
#' @param CI Confidence level in (0, 1) used for intervals. Default `0.95`.
#'
#' @return
#' A list with elements:
#' \itemize{
#'   \item `bbmsy`, `ffmsy` – scenario time series with `lwr`, `est`, `upr`,
#'   \item `catch_pred` – predicted catch per scenario (`time`, `lwr`, `catch`, `upr`),
#'   \item `catch_obs` – observed catch (`time`, `catch`),
#'   \item `eval_pts` – optional evaluation points per scenario,
#'   \item `obsI_rel` – observed indices on B/Bmsy scale (or `NULL`),
#'   \item `t0` – left boundary of `maninterval` (or `NA_real_`).
#' }
#'
#' @keywords internal
#' @noRd
prepare_manage_panel_data <- function(rep_man, CI = 0.95) {
  stopifnot(inherits(rep_man, "spictcls"))

  has_man <- ("man" %in% names(rep_man)) && length(rep_man$man)

  # Shared: observed catch series (if present)
  co_df <- if (!is.null(rep_man$inp$timeC) && !is.null(rep_man$inp$obsC)) {
    data.frame(
      time       = rep_man$inp$timeC,
      catch      = rep_man$inp$obsC,
      catch_type = "Observed"
    )
  } else {
    NULL
  }

  if (has_man) {
    scenarios <- names(rep_man$man)
    ns        <- length(scenarios)

    bb_list <- vector("list", ns)
    ff_list <- vector("list", ns)
    cp_list <- vector("list", ns)
    eval_list <- vector("list", ns)

    for (i in seq_len(ns)) {
      sc <- scenarios[i]
      rp <- rep_man$man[[sc]]

      BB <- get.par("logBBmsy", rp, exp = TRUE)
      FF <- get.par("logFFmsy", rp, exp = TRUE)

      bb_list[[i]] <- data.frame(
        time = as.numeric(rownames(BB)),
        lwr  = BB[, 1],
        est  = BB[, 2],
        upr  = BB[, 3],
        scenario = sc
      )
      ff_list[[i]] <- data.frame(
        time = as.numeric(rownames(FF)),
        lwr  = FF[, 1],
        est  = FF[, 2],
        upr  = FF[, 3],
        scenario = sc
      )

      CP <- get.par("logCpred", rp, exp = TRUE)
      cp_list[[i]] <- data.frame(
        time       = rp$inp$timeCpred,
        lwr        = CP[, 1],
        catch      = CP[, 2],
        upr        = CP[, 3],
        scenario   = sc,
        catch_type = "Predicted"
      )

      ti  <- rp$inp$time
      tE  <- rp$inp$maneval
      idx <- which(ti == tE)
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
    # No scenarios available — return empty scenario frames but keep base info
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
#'
#' @description
#' Computes a reasonable time step (dt) from the observed catch time series.
#' Falls back to the minimum of `inp$dtc` and then to `1` if necessary.
#'
#' @param rep A fitted SPiCT report object (`spictcls`).
#'
#' @return A positive numeric time step.
#'
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
#'
#' @description
#' Returns `geom_vline()` layers marking the left (and, if present, right)
#' edges of the management interval `inp$maninterval`.
#'
#' @param rep A fitted SPiCT report object (`spictcls`).
#' @param color Line colour for the vline(s). Default `"grey30"`.
#' @param linetype Line type for the vline(s). Default `"dashed"`.
#' @param linewidth Line width for the vline(s). Default `0.6`.
#' @param lineend Line end style (`"butt"`, `"round"`, `"square"`). Default `"butt"`.
#'
#' @return A list of ggplot2 layers (possibly empty).
#'
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
#'
#' @description
#' Returns one (or two) `geom_vline()` layers for the Catch panel:
#' \itemize{
#'   \item If no `maninterval` is available, draws a single line at the last
#'         observed catch time (if known).
#'   \item If `maninterval` exists, draws one line at the last observed catch
#'         time (or `left - dt`) and one at the left boundary of `maninterval`.
#' }
#'
#' @param rep A fitted SPiCT report object (`spictcls`).
#' @param color Line colour. Default `"grey30"`.
#' @param linetype Line type. Default `"dashed"`.
#' @param linewidth Line width. Default `0.6`.
#' @param lineend Line end style. Default `"butt"`.
#'
#' @return A list of ggplot2 layers (possibly empty).
#'
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
    # No management interval: draw only the last observed catch start if known
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
#' @description
#' Plots the historical B/Bmsy trajectory with optional confidence bands and,
#' when management scenarios are present in `rep_man$man`, overlays each
#' scenario trajectory using a consistent colour scheme and canonical ordering.
#' Observed indices are shown on the B/Bmsy scale when available.
#'
#' @param rep_man A fitted SPiCT report object (`spictcls`), possibly with
#'   management scenarios (i.e., `rep_man$man` from [spict::manage()]).
#' @param scenario_color Optional named vector of colours for scenarios.
#'   If `NULL`, a repeatable palette is used internally.
#' @param show_CIs Logical; if `TRUE`, draws the base-fit 95% CI ribbon and
#'   boundary lines for B/Bmsy. Default `TRUE`.
#' @param CI Confidence level in (0, 1) for parameter intervals. Default `0.95`.
#' @param show_legend Logical; show/hide the scenario legend. Default `TRUE`.
#'
#' @return A `ggplot` object.
#'
#' @details
#' If no scenarios are present, the function still renders the base B/Bmsy
#' (with optional CI ribbon) and, when available, adds a one-step-ahead segment
#' in a dashed style following ELU2 conventions.
#'
#' Vertical lines:
#' \itemize{
#'   \item With scenarios: `add_management_vlines_BF_good()` marks the
#'         management interval.
#'   \item Without scenarios: a single vline at the last observed time is drawn.
#' }
#'
#' @examples
#' \dontrun{
#' # Base fit
#' p <- my_plot_manage_bbmsy_panel(fit)
#' # With scenarios
#' fit$man <- spict::manage(fit, scenarios = c(1,2,3))
#' my_plot_manage_bbmsy_panel(fit, show_legend = TRUE)
#' }
#'
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

  if (!is.null(dat$bbmsy) && nrow(dat$bbmsy)) {
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
    ## vlines: managed = your helper; raw fit = single line at obs end
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

  # ---------- add ELU2-style 1-step-ahead for raw fit (no $man) ----------
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
      # CI band for the prediction segment - same style as the plot's CIs
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
      # One-ahead mean line as dashed (ELU2 style)
      pB <- pB +
        ggplot2::geom_line(
          data = df_pr, ggplot2::aes(x = time, y = est),
          linetype = "dashed", color = spict_blue_mean, linewidth = 0.8
        )
    }
  }
  # -----------------------------------------------------------------------

  pB
}


#' F/Fmsy panel (ggplot2) with optional scenario overlays
#'
#' @description
#' Plots the historical \eqn{F/F_{MSY}} trajectory (up to the last observation)
#' with optional confidence ribbons. When management scenarios are present in
#' `rep_man$man`, overlays each scenario’s \eqn{F/F_{MSY}} path using a stable,
#' repeatable colour mapping (via [man_cols()]) and a canonical scenario order.
#'
#' @param rep_man A fitted SPiCT report object (`spictcls`), possibly containing
#'   a `man` list from \code{spict::manage()}.
#' @param scenario_color Optional named character vector of colours to use for
#'   scenarios. If `NULL`, a repeatable palette from [man_cols()] is used.
#' @param show_CIs Logical; if `TRUE`, draws the base-fit CI ribbon and boundary
#'   lines for \eqn{F/F_{MSY}} (≤ last observation). Default `TRUE`.
#' @param CI Confidence level in (0, 1) used when extracting intervals from the
#'   SPiCT report. Default `0.95`.
#' @param show_legend Logical; if `TRUE`, shows the scenario legend when
#'   scenarios exist. Default `FALSE`.
#'
#' @return A `ggplot` object.
#'
#' @details
#' The panel always shows the base fit up to the last observed time
#' (`rep_man$inp$timerangeObs[2]`). Vertical reference lines are added as:
#' \itemize{
#'   \item With scenarios: [add_management_vlines_BF_good()] marks the management interval.
#'   \item Without scenarios: a single vertical line at the last observed time.
#' }
#'
#' If no scenarios are present, a one-step-ahead prediction segment is added
#' after the last observed time in a dotted style, with a matching CI ribbon
#' (when `show_CIs = TRUE`), bridging smoothly from the in-sample endpoint.
#'
#' Internally, the function uses:
#' \itemize{
#'   \item [prepare_manage_panel_data()] to assemble tidy scenario data,
#'   \item [get_base_BB_FF_pre()] to obtain in-sample \eqn{F/F_{MSY}} and limits,
#'   \item [man_cols()] for a repeatable scenario colour palette.
#' }
#'
#' @examples
#' \dontrun{
#' # Base fit only
#' pF <- my_plot_manage_ffmsy_panel(fit)
#' print(pF)
#'
#' # With management scenarios
#' fit$man <- spict::manage(fit, scenarios = c("currentF","Fmsy","noF"))
#' my_plot_manage_ffmsy_panel(fit, show_legend = TRUE)
#' }
#'
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

  if (!is.null(dat$ffmsy) && nrow(dat$ffmsy)) {
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

  base_pre <- get_base_BB_FF_pre(rep_man, CI = CI)   # ≤ last observation
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
    ## vlines: managed = your helper; raw fit = single vline at obs end
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

  # ---------- RAW FIT ONLY: one-ahead = dotted mean from the vline + ribbon CI ----------
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

    # Find last in-sample point (<= obs end) and first predicted (> obs end)
    obs_end <- base_pre$t_last_obs
    i0 <- if (any(tt <= obs_end)) max(which(tt <= obs_end)) else NA_integer_
    i1 <- if (any(tt  > obs_end)) min(which(tt  > obs_end)) else NA_integer_

    if (nrow(df_pr)) {
      # CI band & solid CI edges over the prediction times
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

      # Dotted mean on prediction times
      pF <- pF +
        ggplot2::geom_line(
          data = df_pr, ggplot2::aes(x = time, y = est),
          linetype = "dotted", color = spict_blue_mean, linewidth = 0.8
        )

      # Bridge from the vline (obs_end) to the first predicted point so the dotted starts at the vline
      if (is.finite(i0) && is.finite(i1)) {
        pF <- pF +
          ggplot2::geom_segment(
            ggplot2::aes(x = tt[i0], xend = tt[i1], y = df_all$est[i0], yend = df_all$est[i1]),
            linetype = "dotted", color = spict_blue_mean, linewidth = 0.8, lineend = "round"
          )
      }
    }
  }
  # --------------------------------------------------------------------------------------

  pF
}


#' Catch panel (ggplot2) with MSY band and optional scenario overlays
#'
#' @description
#' Plots observed and predicted catch time series with an MSY reference band.
#' If management scenarios are present in `rep_man$man`, their predicted catch
#' trajectories are overlaid with a stable colour mapping (via [man_cols()])
#' and canonical ordering. For raw fits (no scenarios), the base predicted
#' catch is shown with a solid mean up to the end of observed catches and a
#' dotted mean thereafter, with optional CI edges prior to the management start.
#'
#' @param rep_man A fitted SPiCT report object (`spictcls`), optionally
#'   containing a `man` list from \code{spict::manage()}.
#' @param scenario_color Optional named character vector of colours for
#'   scenarios. If `NULL`, a repeatable palette from [man_cols()] is used.
#' @param show_CIs Logical; if `TRUE`, draws dashed CI edge lines for the base
#'   predicted catch up to the management-start mask (`t0`). Default `TRUE`.
#' @param CI Confidence level in (0, 1) used when extracting intervals from the
#'   SPiCT object. Default `0.95`.
#' @param show_legend Logical; show the scenario legend when scenarios exist.
#'   Default `FALSE`.
#'
#' @return A `ggplot` object.
#'
#' @details
#' Data are assembled via [prepare_manage_panel_data()], which provides:
#' \itemize{
#'   \item `catch_obs` – observed catch points (blue dots),
#'   \item `catch_pred` – scenario predicted catch series (lines),
#'   \item `t0` – left boundary of `inp$maninterval` used to mask base Cpred.
#' }
#'
#' MSY reference:
#' \itemize{
#'   \item If time-varying growth or log-\eqn{m} covariates are enabled, the
#'         MSY band comes from `logMSYvec` (time-varying).
#'   \item Otherwise, a constant MSY (from `logMSY`) is repeated across time.
#' }
#' Negative MSY values are clamped to zero before plotting.
#'
#' Vertical lines:
#' \itemize{
#'   \item With scenarios: [add_management_vlines_catch_good()] draws a pair:
#'         last observed catch start (or `left - dt`) and the left boundary of
#'         `maninterval`.
#'   \item Without scenarios: the same helper returns at most one vline (the
#'         last observed catch start), which is used here.
#' }
#'
#' Styling:
#' \itemize{
#'   \item Observed catch uses filled blue dots.
#'   \item Base predicted catch mean is solid up to the last observed catch
#'         start and dotted thereafter (raw fit only).
#'   \item Scenario lines sit above the vlines; legend is optional.
#' }
#'
#' @examples
#' \dontrun{
#' # Base fit only
#' pC <- my_plot_manage_catch_panel(fit)
#' print(pC)
#'
#' # With management scenarios
#' fit$man <- spict::manage(fit, scenarios = c("currentCatch","Fmsy","noF"))
#' my_plot_manage_catch_panel(fit, show_legend = TRUE)
#' }
#'
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
  has_man   <- length(scs_final) > 0  # <— added

  if (!is.null(dat$catch_pred) && nrow(dat$catch_pred)) {
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

  # Management start (for masking) — might be NA for raw fit
  t0 <- dat$t0

  # Base (rep) Cpred up to t0 (available for both fit/manage)
  CP0 <- try(get.par("logCpred", rep_man, exp = TRUE), silent = TRUE)
  base_C_t0 <- NULL
  if (!inherits(CP0, "try-error") && !is.null(rep_man$inp$timeCpred)) {
    tc <- rep_man$inp$timeCpred
    maskC <- is.finite(tc)
    if (is.finite(t0)) maskC <- maskC & (tc < t0)
    if (any(maskC)) {
      base_C_t0 <- data.frame(
        time = tc[maskC],
        lwr  = CP0[maskC, 1],
        est  = CP0[maskC, 2],
        upr  = CP0[maskC, 3]
      )
    }
  }

  # MSY series (time-varying vs constant)
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

  # End of observed catch start (vertical line location)
  obs_end <- if (!is.null(rep_man$inp$timeC) && length(rep_man$inp$timeC))
    utils::tail(rep_man$inp$timeC, 1) else NA_real_   # <— added

  pC <- ggplot2::ggplot() +
    # MSY ribbon + mean (CI band is at the very back)
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
    # ---- VLINES MOVED HERE: above CI band, below scenario paths ----
  { if (has_man)
    add_management_vlines_catch_good(rep_man, color = "grey75", linetype = "solid",
                                     linewidth = 0.4, lineend = "butt")
    else
      add_management_vlines_catch_good(rep_man, color = "grey75", linetype = "solid",
                                       linewidth = 0.4, lineend = "butt")[1]
  } +
    # Scenario predicted catch (no-op for raw fit) — sits on top of vlines
    ggplot2::geom_line(
      data = dat$catch_pred,
      ggplot2::aes(x = time, y = catch, color = scenario),
      linewidth = 0.5, na.rm = TRUE
    ) +
    # Observed catch (works for both)
    { if (!is.null(dat$catch_obs) && nrow(dat$catch_obs)) ggplot2::geom_point(
      data = dat$catch_obs,
      ggplot2::aes(x = time, y = catch),
      color = dot_blue, shape = 16, size = 3, na.rm = TRUE
    ) } +
    # Base predicted catch (≤ t0 if present): CI edges unchanged (dashed)
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
    # Mean line: managed = your original solid;
    #            raw fit = solid ≤ obs_end, dotted ≥ obs_end
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


# ---- Inlined helpers (only those used by this panel) ----

#' Minimal compact ggplot2 theme (ELU2)
#' @param base_size Base font size
#' @param base_family Base font family
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
#' @param model Fitted SPiCT object
#' @keywords internal
#' @noRd
.spict_obs_end_overall <- function(model) {
  inp <- model$inp
  tr <- inp$timerangeObs
  if (!is.null(tr) && length(tr) >= 2 && is.finite(tr[2])) return(tr[2])
  if (!is.null(inp$timeC) && length(inp$timeC))       return(tail(inp$timeC, 1))
  if (!is.null(inp$timeI) && length(inp$timeI) && length(inp$timeI[[1]]) > 0) return(tail(inp$timeI[[1]], 1))
  if (!is.null(inp$time)  && length(inp$time))        return(max(inp$time, na.rm = TRUE))
  NA_real_
}

#' Build a ribbon dataframe from lwr/upr vectors
#' @param time Time vector
#' @param lwr Lower bound
#' @param upr Upper bound
#' @keywords internal
#' @noRd
.spict_make_ribbon_df <- function(time, lwr, upr) {
  ok <- is.finite(time) & is.finite(lwr) & is.finite(upr)
  data.frame(time = time[ok], ymin = lwr[ok], ymax = upr[ok])
}

#' q-scaled index points for first two indices (if present)
#' @param model Fitted SPiCT object
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
#' @param model Fitted SPiCT object
#' @param par_matrix Matrix returned by get.par(...)
#' @keywords internal
#' @noRd
.spict_time_from_par <- function(model, par_matrix) {
  rn <- suppressWarnings(as.numeric(rownames(par_matrix)))
  if (length(rn) == nrow(par_matrix) && !all(is.na(rn))) return(rn)
  t_full <- as.numeric(model$inp$time)
  if (length(t_full) >= nrow(par_matrix)) return(t_full[seq_len(nrow(par_matrix))])
  return(rep_len(t_full, nrow(par_matrix)))
}

#' Biomass panel (managed style, with Bmsy band and B/Bmsy-on-B ribbon)
#'
#' @description
#' Draws an absolute-biomass panel consistent with your managed plotting
#' style: a time-varying \eqn{B_{MSY}} band (ribbon + mean line), a
#' \eqn{B/B_{MSY}}-derived uncertainty ribbon mapped back to biomass,
#' optional scenario biomass paths, vertical markers for management intervals
#' or last observation, and in-sample vs. one-ahead styling for the base fit.
#'
#' @details
#' The function keeps visual conventions aligned with your other panels:
#' \itemize{
#'   \item B\eqn{_{MSY}} band from \code{get.msyvec()} using \code{logBmsy}.
#'   \item Base biomass mean in solid blue; CI edges dashed (in-sample).
#'   \item Raw-fit one-ahead segment in dotted blue with dashed CI edges.
#'   \item Scenario lines are drawn beneath the solid blue mean so it remains visible.
#'   \item Vertical lines via \code{add_management_vlines_BF_good()} when scenarios
#'         exist; otherwise a single vline at the overall observation end.
#' }
#'
#' @inheritParams plot_manage_bbmsy_panel
#'
#' @param rep_man A fitted SPiCT report object (`spictcls`), optionally with
#'   management scenarios under \code{$man}.
#' @param scenario_color Optional named vector of colours for scenarios. When
#'   \code{NULL}, colours are provided by \code{man_cols()} in canonical order.
#' @param show_CIs Logical; draw CI ribbons/edges for B\eqn{/}B\eqn{_{MSY}}
#'   (mapped to biomass) and B\eqn{_{MSY}} band. Default \code{TRUE}.
#' @param CI Confidence level in (0, 1) for interval extraction. Default \code{0.95}.
#' @param show_legend Logical; show the scenario legend if scenarios exist.
#'   Default \code{TRUE}.
#'
#' @return A \code{ggplot} object representing the biomass panel.
#'
#' @seealso
#' \code{\link{my_plot_manage_bbmsy_panel}},
#' \code{\link{add_management_vlines_BF_good}},
#' \code{\link{man_cols}}
#'
#' @examples
#' \dontrun{
#' # Base fit (no scenarios)
#' p <- my_plot_manage_biomass_panel(fit)
#' print(p)
#'
#' # With management scenarios
#' fit$man <- spict::manage(fit, scenarios = c("currentCatch","Fmsy","noF"))
#' my_plot_manage_biomass_panel(fit, show_legend = TRUE)
#' }
#'
#' @export
my_plot_manage_biomass_panel <- function(rep_man,
                                         scenario_color = NULL,
                                         show_CIs = TRUE,
                                         CI = 0.95,
                                         show_legend = TRUE) {
  stopifnot(inherits(rep_man, "spictcls"))

  # Same compact theme used in your managed panels
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

  # ---------- scenario ordering & legend text (same as bbmsy panel) ----------
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

  # Colors for scenarios (kept identical to your other panels)
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

  # ---------- SPiCT series ----------
  Best <- get.par("logB",     rep_man, exp = TRUE, CI = CI)
  BB   <- get.par("logBBmsy", rep_man, exp = TRUE, CI = CI)

  t_B <- .spict_time_from_par(rep_man, Best)
  df_B <- data.frame(time = t_B, lwr = Best[,1], est = Best[,2], upr = Best[,3])

  inp <- rep_man$inp
  ind_in <- inp$indest
  if (length(ind_in) > 0) ind_in <- ind_in[seq_len(length(ind_in) - 1)]  # solid ends at vline
  ind_pr <- if (!length(scs_final)) inp$indpred else integer(0)

  df_B_in <- if (length(ind_in)) df_B[ind_in, , drop = FALSE] else df_B[0, ]
  df_B_pr <- if (length(ind_pr)) df_B[ind_pr, , drop = FALSE] else df_B[0, ]

  # Bmsy(t) band & mean line using get.manmax/get.msyvec (as in your biomass panel)
  repmax   <- if (exists("get.manmax")) get.manmax(rep_man) else rep_man
  Bmsy_all <- get.par("logBmsy", repmax, exp = TRUE, CI = CI)
  Bmsyvec  <- get.msyvec(repmax$inp, Bmsy_all)
  Bmsy     <- if (!is.null(nrow(Bmsy_all))) Bmsy_all[1, ] else Bmsy_all

  df_Bmsy_band <- data.frame(time = repmax$inp$time, ymin = Bmsyvec$ll, ymax = Bmsyvec$ul)
  df_Bmsy_line <- data.frame(time = repmax$inp$time, y = Bmsyvec$msy)

  # B/Bmsy-on-B ribbon (use Bmsy_est to map the relative CI to biomass scale)
  Bmsy_est <- Bmsy[2]
  df_BB_rib <- .spict_make_ribbon_df(
    time = .spict_time_from_par(rep_man, BB),
    lwr  = BB[, 1] * Bmsy_est,
    upr  = BB[, 3] * Bmsy_est
  )

  # Scenario biomass paths (managed case)  ---- drawn BEFORE the main mean
  df_scen_B <- NULL
  if (length(scs_final)) {
    lst <- vector("list", length(scs_final))
    for (i in seq_along(scs_final)) {
      sc <- scs_final[i]
      rp <- rep_man$man[[sc]]
      Bi <- get.par("logB", rp, exp = TRUE, CI = CI)
      ti <- .spict_time_from_par(rp, Bi)
      lst[[i]] <- data.frame(time = ti, est = Bi[,2], scenario = sc)
    }
    df_scen_B <- do.call(rbind, lst)
    df_scen_B$scenario <- factor(df_scen_B$scenario, levels = scs_final)
  }

  # ---------- colours and lines identical to your managed panels ----------
  spict_blue_mean    <- "#0000FF"  # force pure blue so it's identical across prints
  spict_blue_ci_line <- grDevices::rgb(0, 0, 1, 0.20)
  spict_blue_ci_fill <- grDevices::rgb(0, 0, 1, 0.10)

  # ---------- Build plot ----------
  p <- ggplot2::ggplot()

  # Bmsy band & mean
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

  # B/Bmsy-on-B ribbon with solid edges (your ELU2 style)
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

  # --- VLINES HERE (after ribbons, before scenario lines) ---
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

  # Scenario biomass lines (managed) -- drawn BEFORE the main mean so blue stays visible
  if (!is.null(df_scen_B) && nrow(df_scen_B)) {
    p <- p + ggplot2::geom_line(
      data = df_scen_B, ggplot2::aes(x = time, y = est, color = scenario),
      linewidth = 0.8, na.rm = TRUE
    )
  }

  # In-sample Bt: solid mean (blue) OVER the scenarios; dashed CI (blue)
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

  # One-ahead (raw fit only): mean dotted; CI dashed; blue and independent of scales
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

  # Optional q-scaled index points (kept)
  idx <- .spict_index_points(rep_man)
  if (length(idx) >= 1) p <- p + ggplot2::geom_point(
    data = idx[[1]], ggplot2::aes(x = time, y = obs),
    color = "blue", shape = 16, size = 2, inherit.aes = FALSE
  )
  if (length(idx) >= 2) p <- p + ggplot2::geom_point(
    data = idx[[2]], ggplot2::aes(x = time, y = obs),
    shape = 22, color = "black", fill = "green", size = 2, stroke = 0.5, inherit.aes = FALSE
  )

  # Labels, legend, sec axis, theme (consistent with your managed panels)
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


#' Absolute F[t] panel for managed/fit runs (with Fmsy band and optional scenarios)
#'
#' @description
#' Plots absolute fishing mortality \eqn{F_t} through time with a
#' \eqn{F_{MSY}} reference band (ribbon + mean line). If management scenarios
#' are present in `rep_man$man`, their absolute-\eqn{F_t} trajectories are
#' overlaid in a canonical scenario order with stable colours (via `man_cols()`).
#' For raw fits (no scenarios), a one-step-ahead dotted segment is shown after
#' the last observation, with matching CI treatment.
#'
#' @param rep_man A fitted SPiCT object (`spictcls`). If `rep_man$man` exists,
#'   scenario paths are drawn; otherwise a one-ahead dotted segment is shown.
#' @param scenario_color Optional named vector of colours for scenarios. When
#'   `NULL`, colours are generated by `man_cols()` following the canonical order.
#' @param show_CIs Logical; draw ribbons and CI edges for in-sample segments
#'   (and for the one-ahead part in raw fits). Default `TRUE`.
#' @param CI Confidence level for interval extraction in `(0, 1)`. Default `0.95`.
#' @param show_legend Logical; whether to display the scenario legend (only
#'   relevant when scenarios exist). Default `FALSE`.
#'
#' @return A `ggplot` object.
#'
#' @details
#' \itemize{
#'   \item \strong{Fmsy band and line}: If time-varying growth or \code{logmcovflag}
#'         is enabled, the band comes from `logFmsyvec`; otherwise a constant
#'         band is derived from `logFmsy` via `get.msyvec()`.
#'   \item \strong{(F/Fmsy)×Fmsy ribbon}: The \eqn{F/F_{MSY}} uncertainty is mapped
#'         onto absolute \eqn{F_t} using the (time-varying or constant) \eqn{F_{MSY}}.
#'         Pre-management (or pre-one-ahead) segments are shown as a ribbon with
#'         solid CI edges; for raw fits, a post-observation one-ahead ribbon is
#'         also drawn when available.
#'   \item \strong{Vertical lines}: With scenarios, management interval lines are
#'         added via `add_management_vlines_BF_good()`. Otherwise a single vline
#'         is placed at the overall observation end (`.spict_obs_end_overall()`).
#'   \item \strong{In-sample vs one-ahead styling}: In-sample \eqn{F_t} is solid
#'         blue with dashed CI edges; the one-ahead mean is dotted blue with a
#'         small bridging segment so the dotted line starts from the vline.
#'   \item \strong{Sub-annual \eqn{F_t}}: If `inp$dtc < 1`, a semi-transparent
#'         blue line shows sub-annual \eqn{F_t} up to the in-sample boundary.
#'   \item A secondary y-axis displays \eqn{F_t/F_{MSY}} when \eqn{F_{MSY}} is
#'         constant; it is omitted for time-varying \eqn{F_{MSY}}.
#' }
#'
#' @seealso
#' `my_plot_manage_ffmsy_panel()`, `add_management_vlines_BF_good()`, `man_cols()`
#'
#' @examples
#' \dontrun{
#' # Base fit only (no scenarios)
#' p <- my_plot_manage_f_panel(fit)
#' print(p)
#'
#' # With management scenarios and legend
#' fit$man <- spict::manage(fit, scenarios = c("currentF","Fmsy","noF"))
#' my_plot_manage_f_panel(fit, show_legend = TRUE)
#' }
#'
#' @export
my_plot_manage_f_panel <- function(rep_man,
                                   scenario_color = NULL,
                                   show_CIs = TRUE,
                                   CI = 0.95,
                                   show_legend = FALSE) {
  stopifnot(inherits(rep_man, "spictcls"))

  # ---- Local theme (managed panel style) ----
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

  # ---- Core series: Ft (absolute) ----
  Fest <- get.par("logFnotS", rep_man, exp = TRUE, CI = CI)
  dfF <- data.frame(
    time = .spict_time_from_par(rep_man, Fest),
    lwr  = Fest[, 1], est = Fest[, 2], upr = Fest[, 3]
  )

  # In-sample vs one-ahead split (raw fit only)
  ind_in <- inp$indest
  if (length(ind_in) > 0) ind_in <- ind_in[seq_len(length(ind_in) - 1)]
  ind_pr <- if (!manflag) inp$indpred else integer(0)

  dfF_in <- if (length(ind_in)) dfF[ind_in, , drop = FALSE] else dfF[0, ]
  dfF_pr <- if (length(ind_pr)) dfF[ind_pr, , drop = FALSE] else dfF[0, ]

  # ---- Fmsy band + line + secondary axis setup ----
  if (tvgflag) {
    Fmsy_all <- get.par("logFmsyvec", repmax, exp = TRUE, CI = CI)
    Fmsyvec  <- as.data.frame(Fmsy_all); Fmsyvec$msy <- Fmsyvec$est
    Fmsy_band <- data.frame(time = repmax$inp$time, ymin = Fmsyvec$ll, ymax = Fmsyvec$ul)
    Fmsy_line <- data.frame(time = repmax$inp$time, y = Fmsyvec$msy)
    fmsy_mult <- matrix(rep(Fmsyvec$msy, each = 3), ncol = 3, byrow = TRUE)
    sec_axis_obj <- ggplot2::waiver()
  } else {
    Fmsy_all <- get.par("logFmsy", repmax, exp = TRUE, CI = CI)
    if (any(is.na(Fmsy_all))) Fmsy_all <- get.par("logFmsyd", repmax, exp = TRUE, CI = CI)
    xm <- as.matrix(Fmsy_all)
    col_est <- if (!is.null(colnames(xm))) which(colnames(xm) %in% c("est","Est","EST"))[1] else NA_integer_
    if (is.na(col_est)) col_est <- 2
    Fmsy_est <- as.numeric(xm[1, col_est])

    Fmsyvec  <- get.msyvec(repmax$inp, Fmsy_all)
    Fmsy_band <- data.frame(time = repmax$inp$time, ymin = Fmsyvec$ll, ymax = Fmsyvec$ul)
    Fmsy_line <- data.frame(time = repmax$inp$time, y = Fmsyvec$msy)
    fmsy_mult <- Fmsy_est
    sec_axis_obj <- ggplot2::sec_axis(~ . / Fmsy_est, name = expression(F[t]/F[MSY]))
  }

  # ---- (F/Fmsy) × Fmsy ribbon mapped to absolute F ----
  FFrel_all <- get.par("logFFmsynotS", rep_man, exp = TRUE, CI = CI)[, 1:3]
  tt_all    <- .spict_time_from_par(rep_man, FFrel_all)

  last_obs_time <- .spict_obs_end_overall(rep_man)
  ix_last <- which(repmax$inp$time == last_obs_time); if (!length(ix_last)) ix_last <- length(repmax$inp$time)
  indxmax <- max(ix_last - if (manflag) 1 else 0, 1)

  time_pre <- repmax$inp$time[seq_len(indxmax)]
  if (is.matrix(fmsy_mult)) FFabs_all <- FFrel_all * fmsy_mult else FFabs_all <- FFrel_all * as.numeric(fmsy_mult)
  ribFF_pre  <- .spict_make_ribbon_df(time_pre, FFabs_all[seq_len(indxmax), 1], FFabs_all[seq_len(indxmax), 3])

  ribFF_post <- data.frame(time = numeric(0), ymin = numeric(0), ymax = numeric(0))
  if (!manflag && length(ind_pr)) {
    ribFF_post <- .spict_make_ribbon_df(tt_all[ind_pr], FFabs_all[ind_pr, 1], FFabs_all[ind_pr, 3])
  }

  # ---- Scenario Ft paths (absolute) ----
  sc_order <- c("currentCatch","currentF","Fmsy","noF",
                "reduceF25","increaseF25","msyHockeyStick","ices")
  sc_present <- if (manflag) names(rep_man$man) else character(0)
  sc_core   <- intersect(sc_order, sc_present)
  sc_other  <- setdiff(sc_present, sc_core)
  scs_final <- c(sc_core, sort(sc_other))

  df_scen_F <- NULL
  if (length(scs_final)) {
    lst <- vector("list", length(scs_final))
    for (i in seq_along(scs_final)) {
      sc <- scs_final[i]
      rp <- rep_man$man[[sc]]
      Fi <- get.par("logFnotS", rp, exp = TRUE, CI = CI)
      ti <- .spict_time_from_par(rp, Fi)
      lst[[i]] <- data.frame(time = ti, est = Fi[,2], scenario = sc)
    }
    df_scen_F <- do.call(rbind, lst)
    df_scen_F$scenario <- factor(df_scen_F$scenario, levels = scs_final)
  }

  # ---- Colors ----
  if (is.null(scenario_color)) {
    cols <- if (length(scs_final)) { x <- man_cols(length(scs_final)); names(x) <- scs_final; x } else NULL
  } else {
    cols <- scenario_color
    if (length(scs_final)) {
      if (is.null(names(cols))) names(cols) <- scs_final
      cols <- cols[match(scs_final, names(cols))]; names(cols) <- scs_final
    }
  }
  spict_blue_mean    <- "#0000FF"
  spict_blue_ci_line <- grDevices::rgb(0, 0, 1, 0.20)
  spict_blue_ci_fill <- grDevices::rgb(0, 0, 1, 0.10)

  # ---- Build plot ----
  p <- ggplot2::ggplot()

  # Fmsy band & line (background)
  if (isTRUE(show_CIs) && nrow(Fmsy_band)) {
    p <- p + ggplot2::geom_ribbon(
      data = Fmsy_band, ggplot2::aes(x = time, ymin = ymin, ymax = ymax),
      fill = "lightgray", colour = NA
    )
  }
  p <- p + ggplot2::geom_line(
    data = Fmsy_line, ggplot2::aes(x = time, y = y),
    color = "black", linewidth = 0.7
  )

  # (F/Fmsy)×Fmsy ribbons (pre + post)
  if (isTRUE(show_CIs) && nrow(ribFF_pre)) {
    p <- p +
      ggplot2::geom_ribbon(
        data = ribFF_pre, ggplot2::aes(x = time, ymin = ymin, ymax = ymax),
        fill = spict_blue_ci_fill, colour = NA
      ) +
      ggplot2::geom_line(
        data = transform(ribFF_pre, y = ymin), ggplot2::aes(x = time, y = y),
        color = spict_blue_ci_line, linewidth = 0.6
      ) +
      ggplot2::geom_line(
        data = transform(ribFF_pre, y = ymax), ggplot2::aes(x = time, y = y),
        color = spict_blue_ci_line, linewidth = 0.6
      )
  }
  if (!manflag && isTRUE(show_CIs) && nrow(ribFF_post)) {
    p <- p +
      ggplot2::geom_ribbon(
        data = ribFF_post, ggplot2::aes(x = time, ymin = ymin, ymax = ymax),
        fill = spict_blue_ci_fill, colour = NA
      ) +
      ggplot2::geom_line(
        data = transform(ribFF_post, y = ymin), ggplot2::aes(x = time, y = y),
        color = spict_blue_ci_line, linewidth = 0.6
      ) +
      ggplot2::geom_line(
        data = transform(ribFF_post, y = ymax), ggplot2::aes(x = time, y = y),
        color = spict_blue_ci_line, linewidth = 0.6
      )
  }

  # --- VLINES HERE (after ribbons, before scenario lines) ---
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

  # Scenario Ft paths (draw AFTER vlines so scenarios sit on top)
  if (!is.null(df_scen_F) && nrow(df_scen_F)) {
    p <- p + ggplot2::geom_line(
      data = df_scen_F, ggplot2::aes(x = time, y = est, color = scenario),
      linewidth = 0.8, na.rm = TRUE
    )
  }

  # Sub-annual Ft (semi-transparent blue), if present (pre window)
  if (min(inp$dtc, na.rm = TRUE) < 1 && indxmax >= 1) {
    sF <- get.par("logFs", rep_man, exp = TRUE, CI = CI)
    p <- p + ggplot2::geom_line(
      data = data.frame(time = repmax$inp$time[seq_len(indxmax)], est = sF[seq_len(indxmax), 2]),
      ggplot2::aes(x = time, y = est),
      color = grDevices::rgb(0,0,1,0.4), linewidth = 0.6
    )
  }

  # In-sample Ft mean/CI (blue, on top)
  if (nrow(dfF_in)) {
    if (isTRUE(show_CIs)) {
      p <- p +
        ggplot2::geom_line(
          data = dfF_in, ggplot2::aes(x = time, y = lwr),
          linetype = "dashed", linewidth = 0.6,
          colour = "#0000FF", inherit.aes = FALSE, show.legend = FALSE
        ) +
        ggplot2::geom_line(
          data = dfF_in, ggplot2::aes(x = time, y = upr),
          linetype = "dashed", linewidth = 0.6,
          colour = "#0000FF", inherit.aes = FALSE, show.legend = FALSE
        )
    }
    p <- p + ggplot2::geom_line(
      data = dfF_in, ggplot2::aes(x = time, y = est),
      linewidth = 0.9, colour = "#0000FF",
      inherit.aes = FALSE, show.legend = FALSE
    )
  }

  # RAW FIT ONLY: one-ahead dotted mean (blue) + bridge
  if (!manflag && nrow(dfF_pr)) {
    p <- p + ggplot2::geom_line(
      data = dfF_pr, ggplot2::aes(x = time, y = est),
      linetype = "dotted", linewidth = 0.8,
      colour = "#0000FF", inherit.aes = FALSE, show.legend = FALSE
    )
    obs_end <- .spict_obs_end_overall(rep_man)
    ttF <- dfF$time
    i0 <- if (any(ttF <= obs_end)) max(which(ttF <= obs_end)) else NA_integer_
    i1 <- if (any(ttF  > obs_end)) min(which(ttF  > obs_end)) else NA_integer_
    if (is.finite(i0) && is.finite(i1)) {
      p <- p + ggplot2::geom_segment(
        ggplot2::aes(x = ttF[i0], xend = ttF[i1], y = dfF$est[i0], yend = dfF$est[i1]),
        linetype = "dotted", color = "#0000FF", linewidth = 0.8, lineend = "round"
      )
    }
  }

  # Labels, scales, theme
  p +
    ggplot2::labs(title = "Absolute fishing mortality", x = "Year", y = expression(F[t])) +
    { if (length(scs_final)) ggplot2::scale_color_manual(
      values = cols, guide = if (show_legend) "legend" else "none"
    ) else ggplot2::scale_color_discrete(guide = "none") } +
    ggplot2::scale_fill_manual(values = if (is.null(cols)) NA else cols, guide = "none") +
    ggplot2::scale_y_continuous(sec.axis = sec_axis_obj,
                                expand = ggplot2::expansion(mult = c(0.02, 0.02))) +
    theme_minimal_compact2_good_local()
}


