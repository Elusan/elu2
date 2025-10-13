
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
#' @importFrom stats qnorm cov2cor
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

  # ---- tiny safe-field helpers ----
  .sf  <- function(x, nm, default = NULL) if (!is.null(x) && nm %in% names(x)) x[[nm]] else default
  .sfb <- function(x, nm) isTRUE(.sf(x, nm, FALSE))
  .sfi <- function(x, nm, default = integer(0)) .sf(x, nm, default)
  .sfn <- function(x, nm, default = numeric(0)) .sf(x, nm, default)

  # ---- stats wrappers (explicitly from stats) ----
  .qnorm   <- function(...) stats::qnorm(...)
  .cov2cor <- function(...) stats::cov2cor(...)

  # ---- sec_axis() compatibility (ggplot2 >= 3.5 uses `transform=`) ----
  .sec_axis_compat <- function(.transform, ...) {
    fa <- formals(ggplot2::sec_axis)
    if ("transform" %in% names(fa)) {
      ggplot2::sec_axis(transform = .transform, ...)
    } else {
      ggplot2::sec_axis(trans = .transform, ...)
    }
  }

  get_par <- function(parname, rep, exp = FALSE, CI = 0.95) {
    if (CI <= 0 || CI >= 1) stop("CI must be in (0,1).")
    z <- .qnorm(CI + (1 - CI) / 2)

    ind_ran <- which(names(rep$par.random) == parname)
    ind_fix <- which(names(rep$par.fixed)  == parname)
    ind_sdr <- which(names(rep$value)      == parname)
    ind_opt <- which(names(rep$opt$par)    == parname)

    est <- ll <- ul <- sdv <- NULL

    if (length(ind_ran)) {
      est <- rep$par.random[ind_ran]
      if (!is.null(rep$diag.cov.random)) {
        sdv <- sqrt(rep$diag.cov.random[ind_ran])
      } else if (!is.null(rep$cov.random)) {
        sdv <- sqrt(diag(rep$cov.random))[ind_ran]
      } else {
        sdv <- rep(NA_real_, length(ind_ran))
      }
      ll  <- est - z * sdv; ul <- est + z * sdv
    }
    if (length(ind_fix)) {
      est <- rep$par.fixed[ind_fix]
      sd_all <- diag(rep$cov.fixed)
      sdv <- if (length(sd_all)) sqrt(sd_all)[ind_fix] else rep(NA_real_, length(ind_fix))
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
    if (parname %in% c("logB", "logF", "logBBmsy", "logFFmsy", "logFs")) {
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
      ellipse::ellipse(.cov2cor(covBF)[1, 2],
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
    K     <- get_par("logK",   rep, exp = TRUE)[1, 2]
    n     <- get_par("logn",   rep, exp = TRUE)[1, 2]
    sdb2  <- get_par("logsdb", rep, exp = TRUE)[1, 2]^2
    Fmsy  <- tail(get_par("logFmsy", rep, exp = TRUE), 1)[, 2]
    logFs <- get_par("logFs",  rep,        CI = CI)
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

  clamp <- function(val, lo, hi) pmin(pmax(val, lo), hi)

  ## ------------------ Validate & derive quantities --------------- ##
  check_rep(rep)

  tvgflag <- .sfb(rep$inp, "timevaryinggrowth") | .sfb(rep$inp, "logmcovflag")
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

  indest  <- .sfi(rep$inp, "indest",  seq_along(rep$inp$time))
  indpred <- .sfi(rep$inp, "indpred", integer(0))

  if (min(rep$inp$dtc) < 1) {
    alb <- annual_avg(rep$inp$time, logBest[, 2])
    alf <- annual_avg(rep$inp$time, logFest[, 2])
    bbb <- exp(alb$annvec) / bscal
    fff <- exp(alf$annvec) / fscal
    fbtime <- alb$anntime
  } else {
    bbb <- Best[indest, 2] / bscal
    fff <- Fest[indest, 2] / fscal
    fbtime <- rep$inp$time[indest]
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
  if (!(min(rep$inp$dtc) < 1) && !("man" %in% names(rep)) && length(indpred)) {
    xpred <- Best[indpred, 2] / (if (rel.axes) Bmsy[2] else bscal)
    ypred <- Fest[indpred, 2] / (if (rel.axes) Fmsy[2] else fscal)
    df_pred <- data.frame(x = xpred, y = ypred)
    Bll <- tail(xpred, 1)
    Fll <- tail(ypred, 1)
    df_EBseg <- data.frame(x = Bll, xend = EBinf, y = Fll, yend = Fll)
  }

  df_first <- data.frame(x = bbb[1], y = fff[1], lab = round(fbtime[1], 2))
  df_last  <- data.frame(x = tail(bbb,1), y = tail(fff,1), lab = round(tail(fbtime, 1), 2))

  logr_ini <- .sfn(rep$inp$ini, "logr", 1)
  nr <- length(logr_ini)
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

  ## ========== Scenario order/labels and jump segment ==========
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
      estF <- get_par("logFs", rp, exp = TRUE)[, 2]   # <-- FIXED (logFs)
      estB <- get_par("logB",  rp, exp = TRUE)[, 2]
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

  ## Small additional padding block
  pad_frac <- 0.02
  xpad <- pad_frac * diff(xlim)
  ypad <- pad_frac * diff(ylim)
  xlim <- c(max(if (isTRUE(logax)) logminval else 0, xlim[1] - xpad), xlim[2] + xpad)
  ylim <- c(max(if (isTRUE(logax)) logminval else 0, ylim[1] - ypad), ylim[2] + ypad)

  # Labels & scales
  xlab_final <- if (is.null(xlabel)) xlab_expr else xlabel
  ylab_final <- ylab_expr

  # Base 'pretty' breaks with forced 1 in labels (no 'scales' dependency)
  pretty_rel_breaks <- function(lims) {
    br <- pretty(lims, n = 6)
    sort(unique(c(br, 1)))
  }

  if (isTRUE(logax)) {
    p <- p + ggplot2::scale_x_log10(
      limits = xlim,
      name   = xlab_final,
      labels = function(x) formatC(x, format = "f", digits = 0),
      expand = c(0, 0),
      sec.axis = if (!rel.axes && isTRUE(ext)) .sec_axis_compat(
        ~ . / (Bmsy[2] / 1),
        name   = expression(B[t]/B[MSY]),
        breaks = pretty_rel_breaks,
        labels = lab_rel_1
      ) else ggplot2::waiver()
    )

    p <- p + ggplot2::scale_y_log10(
      limits = ylim,
      name   = ylab_final,
      labels = fmt1,
      expand = ggplot2::expansion(mult = c(0.03, 0.03)),
      sec.axis = if (!rel.axes && isTRUE(ext)) .sec_axis_compat(
        ~ . / (Fmsy[2] / 1),
        name   = expression(F[t]/F[MSY]),
        breaks = pretty_rel_breaks,
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
      sec.axis = if (!rel.axes && isTRUE(ext)) .sec_axis_compat(
        ~ . / (Bmsy[2] / 1),
        name   = expression(B[t]/B[MSY]),
        breaks = pretty_rel_breaks,
        labels = lab_rel_1
      ) else ggplot2::waiver()
    )

    p <- p + ggplot2::scale_y_continuous(
      limits = ylim,
      name   = ylab_final,
      labels = fmt1,
      expand = ggplot2::expansion(mult = c(0.03, 0.03)),
      sec.axis = if (!rel.axes && isTRUE(ext)) .sec_axis_compat(
        ~ . / (Fmsy[2] / 1),
        name   = expression(F[t]/F[MSY]),
        breaks = pretty_rel_breaks,
        labels = lab_rel_1
      ) else ggplot2::waiver()
    )
  }

  # ===== Unified scenario scale (hidden legend for clean panel) =====
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
      plot.margin  = ggplot2::margin(4, 4, 4, 4)
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

  p
}
