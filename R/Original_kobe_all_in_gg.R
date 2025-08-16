#' All-in-one Kobe plot (ggplot2)
#'
#' @description
#' Produces a Kobe plot of fishing mortality vs. biomass using a fitted
#' **SPiCT** result, rendered with **ggplot2**. The plot includes colored
#' quadrants, an MSY reference-point “banana” (ellipse from
#' \eqn{\log B_{\mathrm{MSY}}}–\eqn{\log F_{\mathrm{MSY}}} covariance), the
#' biomass–F trajectory, optional management-scenario overlays, first/last-year
#' markers, and optional secondary axes in relative units.
#'
#' @details
#' The function expects a fitted object of class `spictcls` (from
#' \code{elu2::fit.elu2()}) with `inp$reportmode == 0` so that all states are
#' reported. When `rel.axes = TRUE`, the plot uses relative axes
#' \eqn{B_t/B_{\mathrm{MSY}}} vs. \eqn{F_t/F_{\mathrm{MSY}}} (and the extra
#' top/right secondary axes are suppressed). If either `inp$timevaryinggrowth`
#' or `inp$logmcovflag` is `TRUE`, relative axes are enforced.
#'
#' For sub-annual time steps (`min(inp$dtc) < 1`), biomass and fishing
#' mortality series are aggregated to annual means for plotting. The MSY
#' “banana” is computed from the estimated covariance of
#' \eqn{\log B_{\mathrm{MSY}}} and \eqn{\log F_{\mathrm{MSY}}}. If the
#' **ellipse** package is available, a proper ellipse is used; otherwise, a
#' fallback point at (\eqn{B_{\mathrm{MSY}}}, \eqn{F_{\mathrm{MSY}}}) is drawn.
#'
#' If management scenarios are present in `rep$man`, their trajectories (and
#' the intermediate period up to the management start) are overlaid and can be
#' shown in a legend. Otherwise, a one-step prediction segment and the
#' \eqn{\mathrm{E}(B_{\infty})} marker are drawn.
#'
#' @param rep A fitted SPiCT report object of class `spictcls`
#'   (result of \code{elu2::fit.elu2()}). Must have `inp$reportmode == 0`.
#' @param logax Logical; if `TRUE`, both axes are log-scaled (log10). Default `FALSE`.
#' @param plot.legend Logical; show the small top-right legend for
#'   \eqn{\mathrm{E}(B_{\infty})} (and “True” if present). Default `TRUE`.
#' @param man.legend Logical; when management scenarios exist in `rep$man`, show
#'   their legend. Default `TRUE`.
#' @param ext Logical; add secondary axes with relative units on the top/right
#'   when `rel.axes = FALSE`. Ignored when `rel.axes = TRUE`. Default `TRUE`.
#' @param rel.axes Logical; plot axes as relative units
#'   \eqn{B_t/B_{\mathrm{MSY}}} vs. \eqn{F_t/F_{\mathrm{MSY}}}. Forces `ext = FALSE`.
#'   Default `FALSE`.
#' @param xlim,ylim Numeric length-2 vectors giving x/y limits. If `NULL`
#'   (default), limits are computed from the ellipse, time series, and reference
#'   points with a small padding to avoid border overlap.
#' @param labpos Integer length-2; relative positions for the first/last-year
#'   labels (used internally when computing offsets). Default `c(1, 1)`.
#' @param xlabel Optional x-axis label (overrides the default). Default `NULL`.
#' @param stamp Optional character string stamped at the lower-right margin
#'   (e.g., a version tag). Default `NULL`.
#' @param verbose Logical; print notes/warnings (e.g., if a parameter cannot be
#'   extracted). Default `TRUE`.
#' @param CI Numeric in `(0,1)`; confidence level used when extracting
#'   parameters and constructing intervals (e.g., `0.95`). Default `0.95`.
#'
#' @return
#' (Invisibly) returns the `ggplot` object, and prints it as a side effect.
#'
#' @section Secondary axes:
#' When `rel.axes = FALSE` and `ext = TRUE`, duplicated axes are added on the
#' top and right that convert absolute units to relative units using the
#' estimated \eqn{B_{\mathrm{MSY}}} and \eqn{F_{\mathrm{MSY}}}. Labels format `1`
#' exactly as `"1"` for easy identification of the MSY thresholds.
#'
#' @section Dependencies:
#' Requires **ggplot2**. If available, **ellipse** is used for the MSY
#' covariance ellipse. Uses \code{stats::qnorm} and \code{scales::pretty_breaks}
#' (referenced with explicit namespaces).
#'
#' @seealso
#' \code{elu2::fit.elu2()} for model fitting; \code{spict::plotspict()} for base plots.
#'
#' @examples
#' \dontrun{
#' # (Pseudo-code) Fit a SPiCT model, then:
#' # rep <- elu2::fit.elu2(inp)
#'
#' # Default: absolute axes
#' kobe_all_in_one_gg(rep)
#'
#' # Relative axes (B/Bmsy vs F/Fmsy) with log scaling
#' kobe_all_in_one_gg(rep, rel.axes = TRUE, logax = TRUE)
#'
#' # Show management overlays but hide their legend, add a stamp
#' kobe_all_in_one_gg(rep, man.legend = FALSE, stamp = "Analysis v1.0")
#' }
#'
#' @author
#' Elhadji Ndiaye & Contributors
#'
#' @encoding UTF-8
#'
#' @importFrom stats qnorm
#' @import ggplot2
#' @export
Original_kobe_all_in_gg <- function(rep,
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
    # Use your fmt1() (1 decimal) but rewrite exactly-1 to "1"
    out <- fmt1(vals)
    is_one <- abs(vals - 1) < 1e-9
    out[is_one] <- "1"
    out
  }
  rel_breaks_fun_x <- function(primary_lims) {
    # Convert pretty primary ticks to relative units, add 1, ensure unique & sorted
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
  # Full-panel bounds
  xminq <- xlim[1]; xmaxq <- xlim[2]
  yminq <- ylim[1]; ymaxq <- ylim[2]

  # Split lines, clamped to panel to avoid "removed rows" warnings
  bx <- if (rel.axes) 1 else Bmsy[2] / bscal
  bx <- clamp(bx, xminq, xmaxq)
  fy <- if (rel.axes) 1 else Fmsy[2] / fscal
  fy <- clamp(fy, yminq, ymaxq)

  # Quadrant rectangles (fill everything, no gaps)
  # Quadrant rectangles — tile panel to borders using -Inf / Inf
  df_green_rect   <- data.frame(xmin = bx,   xmax =  Inf, ymin = -Inf, ymax =  fy)  # lower-right
  df_yellowL_rect <- data.frame(xmin = -Inf, xmax =  bx,  ymin = -Inf, ymax =  fy)  # lower-left
  df_yellowT_rect <- data.frame(xmin = bx,   xmax =  Inf, ymin =  fy,  ymax =  Inf) # upper-right
  df_red_rect     <- data.frame(xmin = -Inf, xmax =  bx,  ymin =  fy,  ymax =  Inf) # upper-left


  df_ellipse <- if (ncol(as.matrix(cl)) == 2 && nrow(as.matrix(cl)) > 1) {
    data.frame(x = exp(cl[, 1]) / (if (rel.axes) Bmsy[2] else bscal),
               y = exp(cl[, 2]) / (if (rel.axes) Fmsy[2] else fscal))
  } else {
    data.frame(x = exp(cl[1]) / (if (rel.axes) Bmsy[2] else bscal),
               y = exp(cl[2]) / (if (rel.axes) Fmsy[2] else fscal))
  }

  # Clamp trajectory just inside the frame so it never touches borders
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

  df_man <- df_man_int <- NULL
  leg_man <- NULL
  if ("man" %in% names(rep) && length(rep$man)) {
    nman <- length(rep$man)
    leg_man <- names(rep$man)
    if (is.null(leg_man) || !length(leg_man)) leg_man <- paste0("Scenario ", seq_len(nman))
    df_list <- vector("list", nman)
    df_int_list <- vector("list", nman)
    i <- 1L
    while (i <= nman) {
      rp <- rep$man[[i]]
      estF <- get_par("logF", rp, exp = TRUE)[, 2]
      estB <- get_par("logB", rp, exp = TRUE)[, 2]
      time <- rp$inp$time
      manint <- rp$inp$maninterval
      indmanstart <- which(time >= manint[1])
      if (length(indmanstart)) {
        maninds <- indmanstart[1]:tail(indmanstart, 1)
        df_list[[i]] <- data.frame(x = estB[maninds] / (if (rel.axes) Bmsy[2] else bscal),
                                   y = estF[maninds] / (if (rel.axes) Fmsy[2] else fscal),
                                   scenario = leg_man[i])
      } else {
        df_list[[i]] <- data.frame(x = numeric(0), y = numeric(0), scenario = character(0))
      }
      lastobs <- rp$inp$timerangeObs[2]
      ind <- which(time >= lastobs & time < manint[1])
      if (length(ind)) {
        df_int_list[[i]] <- data.frame(x = estB[ind] / (if (rel.axes) Bmsy[2] else bscal),
                                       y = estF[ind] / (if (rel.axes) Fmsy[2] else fscal),
                                       scenario = leg_man[i])
      } else {
        df_int_list[[i]] <- data.frame(x = numeric(0), y = numeric(0), scenario = character(0))
      }
      i <- i + 1L
    }
    df_man <- do.call(rbind, df_list)
    df_man_int <- do.call(rbind, df_int_list)
  }

  ## --------------------------- Build ggplot ------------------------------ ##
  p <- ggplot2::ggplot()

  # Quadrants (fill entire panel; no gaps)
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

  # Ellipse (banana CI) — transparent gray
  if (nrow(df_ellipse) > 1) {
    p <- p + ggplot2::geom_polygon(data = df_ellipse, ggplot2::aes(x = x, y = y),
                                   fill = grDevices::adjustcolor("gray80", alpha.f = 0.7),
                                   color = grDevices::adjustcolor("gray80", alpha.f = 0.7), linewidth = 0.5)
  } else {
    p <- p + ggplot2::geom_point(data = df_ellipse, ggplot2::aes(x = x, y = y),
                                 color = "gray40", size = 2)
  }

  # Historical trajectory (thinner and safely inside the frame)
  p <- p + ggplot2::geom_path(data = df_traj, ggplot2::aes(x = x, y = y),
                              linewidth = 0.4, color = grDevices::adjustcolor("blue", alpha.f = 0.8))

  # One-step prediction
  if (!is.null(df_pred)) {
    p <- p + ggplot2::geom_path(data = df_pred, ggplot2::aes(x = x, y = y),
                                linetype = "33", color = grDevices::adjustcolor("blue", alpha.f = 0.8))
  }
  if (!is.null(df_EBseg)) {
    p <- p + ggplot2::geom_segment(data = df_EBseg,
                                   ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
                                   color = "blue", linetype = "33", linewidth = 1.0)
    p <- p + ggplot2::geom_point(data = data.frame(x = df_EBseg$xend, y = df_EBseg$yend),
                                 ggplot2::aes(x = x, y = y),
                                 shape = 23, fill = "gold", color = "black", size = 3, stroke = 0.7)
  }

  # Management lines
  if (!is.null(df_man) && nrow(df_man)) {
    p <- p + ggplot2::geom_path(data = df_man, ggplot2::aes(x = x, y = y, color = scenario),
                                linewidth = 0.7, inherit.aes = FALSE, show.legend = man.legend)
  }
  if (!is.null(df_man_int) && nrow(df_man_int)) {
    p <- p + ggplot2::geom_path(data = df_man_int, ggplot2::aes(x = x, y = y, color = scenario),
                                linewidth = 0.7, linetype = "33", inherit.aes = FALSE, show.legend = man.legend)
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
                               shape = 21, fill = "white", color = "black", size = 3, stroke = 0.5)
  p <- p + ggplot2::geom_text(data = df_first, ggplot2::aes(x = x, y = y, label = lab),
                              vjust = -0.8, size = 3)
  p <- p + ggplot2::geom_point(data = df_last, ggplot2::aes(x = x, y = y),
                               shape = 22, fill = "white", color = "black", size = 3, stroke = 0.5)
  p <- p + ggplot2::geom_text(data = df_last, ggplot2::aes(x = x, y = y, label = lab),
                              vjust = -0.8, size = 3)

  # Reference vline at 0 only if inside limits
  if (!isTRUE(logax) && (0 >= xlim[1]) && (0 <= xlim[2])) {
    p <- p + ggplot2::geom_vline(xintercept = 0, color = "darkred", linetype = 2)
  }

  ## Label function for secondary axes: force exact 1 → "1" (no decimal)
  lab_rel_1 <- function(vals) {
    out <- fmt1(vals)                 # your existing formatter (1 decimal)
    out[abs(vals - 1) < 1e-9] <- "1"  # but print 1 as "1"
    out
  }

  ## (Optional but recommended) Small padding of limits to avoid overlap with frame
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
        trans  = ~ . / (Bmsy[2] / 1),                           # absolute → relative
        name   = expression(B[t]/B[MSY]),
        breaks = function(lims) {
          br <- scales::pretty_breaks(n = 6)(lims)
          sort(unique(c(br, 1)))
        },
        labels = lab_rel_1
      ) else ggplot2::waiver()
    )

    p <- p + ggplot2::scale_y_log10(
      limits = ylim,
      name   = ylab_final,
      labels = fmt1,
      expand = c(0, 0),
      sec.axis = if (!rel.axes && isTRUE(ext)) ggplot2::sec_axis(
        trans  = ~ . / (Fmsy[2] / 1),
        name   = expression(F[t]/F[MSY]),
        breaks = function(lims) {
          br <- scales::pretty_breaks(n = 6)(lims)
          sort(unique(c(br, 1)))
        },
        labels = lab_rel_1
      ) else ggplot2::waiver()
    )

  } else {

    p <- p + ggplot2::scale_x_continuous(
      limits = xlim,
      name   = xlab_final,
      labels = function(x) formatC(x, format = "f", digits = 0),
      expand = c(0, 0),
      sec.axis = if (!rel.axes && isTRUE(ext)) ggplot2::sec_axis(
        trans  = ~ . / (Bmsy[2] / 1),
        name   = expression(B[t]/B[MSY]),
        breaks = function(lims) {
          br <- scales::pretty_breaks(n = 6)(lims)
          sort(unique(c(br, 1)))
        },
        labels = lab_rel_1
      ) else ggplot2::waiver()
    )

    p <- p + ggplot2::scale_y_continuous(
      limits = ylim,
      name   = ylab_final,
      labels = fmt1,
      expand = c(0, 0),
      sec.axis = if (!rel.axes && isTRUE(ext)) ggplot2::sec_axis(
        trans  = ~ . / (Fmsy[2] / 1),
        name   = expression(F[t]/F[MSY]),
        breaks = function(lims) {
          br <- scales::pretty_breaks(n = 6)(lims)
          sort(unique(c(br, 1)))
        },
        labels = lab_rel_1
      ) else ggplot2::waiver()
    )
  }




  if (!is.null(df_man) && nrow(df_man) && man.legend) {
    uniq_sc <- unique(df_man$scenario)
    n_sc <- length(uniq_sc)
    colv <- man_cols()[seq_len(n_sc)]
    names(colv) <- uniq_sc
    p <- p + ggplot2::scale_color_manual(values = colv, name = "Scenario")
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

  # Top-right mini-legend for E(B∞): smaller diamond, more gap, math label bigger
  if (plot.legend && !(min(rep$inp$dtc) < 1)) {
    px <- xlim[2] - 0.02 * diff(xlim)
    py <- ylim[2] - 0.04 * diff(ylim)
    gap <- 0.15 * diff(xlim)   # more horizontal spacing
    p <- p + ggplot2::annotate("point", x = px - gap, y = py,
                               shape = 23, size = 2.5, fill = "gold", color = "black", stroke = 0.6)
    p <- p + ggplot2::annotate("text", x = px, y = py,
                               label = "E(B[infinity])", parse = TRUE,
                               hjust = 1, vjust = 0.5, fontface = "bold", size = 4)
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

  print(p)
  invisible(p)
}
