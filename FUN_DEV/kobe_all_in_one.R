#' All-in-one Kobe plot (base graphics)
#'
#' @description
#' Draws a Kobe plot of fishing mortality versus biomass using a fitted
#' **SPiCT** result. The plot shows colored management quadrants, the
#' MSY reference-point “banana” (from the covariance of \eqn{\log B_{\mathrm{MSY}}}
#' and \eqn{\log F_{\mathrm{MSY}}}), the biomass–fishing mortality trajectory,
#' optional management-scenario overlays, and stamped first/last years.
#'
#' @details
#' The function expects a fitted object of class `spictcls` (from
#' `fit.spict()`) with `inp$reportmode == 0` so that all states are reported.
#' When \code{rel.axes = TRUE}, axes are drawn in relative units
#' \eqn{B_t/B_{\mathrm{MSY}}} and \eqn{F_t/F_{\mathrm{MSY}}} (and the extra top/right
#' relative axes are suppressed). If \code{inp$timevaryinggrowth} or
#' \code{inp$logmcovflag} is \code{TRUE}, relative axes are enforced. When time steps
#' are sub-annual (\code{min(inp$dtc) < 1}), biomass and fishing mortality time
#' series are aggregated to annual means for plotting.
#'
#' The “banana” is produced from the estimated covariance between
#' \eqn{\log B_{\mathrm{MSY}}} and \eqn{\log F_{\mathrm{MSY}}}. If the
#' \pkg{ellipse} package is available, a proper ellipse is drawn; otherwise a
#' fallback marker at (\eqn{B_{\mathrm{MSY}}}, \eqn{F_{\mathrm{MSY}}}) is used.
#'
#' If management scenarios are present in \code{rep$man}, corresponding
#' trajectories are overlaid and optionally labeled. Otherwise, a one-step
#' prediction segment and the equilibrium biomass marker \eqn{\mathrm{E}(B_{\infty})}
#' are shown.
#'
#' @param rep A fitted SPiCT report object of class \code{spictcls}
#'   (the result of \code{fit.spict()}). Must have \code{inp$reportmode == 0}.
#' @param logax Logical; if \code{TRUE}, both axes are log-scaled (equivalent to
#'   base-graphics \code{log = "xy"}). Default \code{FALSE}.
#' @param plot.legend Logical; show the legend for MSY markers,
#'   \eqn{\mathrm{E}(B_{\infty})}, and/or "True" point if present. Default \code{TRUE}.
#' @param man.legend Logical; when management scenarios exist in \code{rep$man},
#'   show their legend. Default \code{TRUE}.
#' @param ext Logical; add relative axes (\eqn{B_t/B_{\mathrm{MSY}}} on top,
#'   \eqn{F_t/F_{\mathrm{MSY}}} on right) when \code{rel.axes = FALSE}. Ignored when
#'   \code{rel.axes = TRUE}. Default \code{TRUE}.
#' @param rel.axes Logical; plot axes in relative units
#'   \eqn{B_t/B_{\mathrm{MSY}}} vs \eqn{F_t/F_{\mathrm{MSY}}}. Forces \code{ext = FALSE}.
#'   Default \code{FALSE}.
#' @param xlim,ylim Numeric length-2 vectors giving x/y limits. If \code{NULL}
#'   (default), limits are computed from the ellipse, time series, and
#'   reference points with safety expansion.
#' @param labpos Integer length-2 vector passed to \code{text(pos = ...)} to position
#'   the first and last year labels. Default \code{c(1, 1)}.
#' @param xlabel Optional x-axis label (overrides the default). Default \code{NULL}.
#' @param stamp Optional character string to stamp in the bottom-right margin.
#'   If \code{NULL}, no extra stamp is drawn (the internal version stamp helper can
#'   be used by passing a string). Default \code{NULL}.
#' @param verbose Logical; print notes (e.g., when a management period is
#'   shorter than one year). Default \code{TRUE}.
#' @param CI Numeric in \code{(0,1)}; confidence level used when extracting
#'   parameters and constructing intervals (e.g., \code{0.95}). Default \code{0.95}.
#'
#' @return Invisibly returns \code{NULL}. The function is called for its side
#'   effect of drawing a plot on the active device (base graphics).
#'
#' @section Quadrants and markers:
#' The background is divided into the usual Kobe quadrants by
#' \eqn{B_{\mathrm{MSY}}} and \eqn{F_{\mathrm{MSY}}}:
#' green (safe: \eqn{B \ge B_{\mathrm{MSY}}} and \eqn{F \le F_{\mathrm{MSY}}}),
#' yellow (transition), and red (unsafe). Current and previous MSY estimates
#' (for multiple initializations) are marked with \code{pch = 3}; a “True” point
#' (if supplied in \code{rep$inp$true}) is drawn with \code{pch = 25}.
#'
#' @section Dependencies:
#' Uses base R graphics. If available, \pkg{ellipse} is used to draw the MSY
#' covariance ellipse. Utilities from \pkg{utils} are used for version stamping.
#'
#' @seealso
#' \code{\link[spict]{fit.spict}}, \code{\link[spict]{plotspict}} for the original
#' plotting functions and model fitting.
#'
#' @examples
#' \dontrun{
#' # Fit a SPiCT model (pseudo-code)
#' # inp  <- make.spict.input(obsC = ..., timeC = ..., obsI = ..., timeI = ...)
#' # rep  <- fit.spict(inp)
#'
#' # Default Kobe plot (absolute axes)
#' kobe_all_in_one(rep)
#'
#' # Relative axes (B/Bmsy vs F/Fmsy) with log scaling
#' kobe_all_in_one(rep, rel.axes = TRUE, logax = TRUE)
#'
#' # Suppress scenario legend, provide a custom stamp
#' kobe_all_in_one(rep, man.legend = FALSE, stamp = "MyAnalysis v1.0")
#' }
#'
#' @author
#' Elhadji Ndiaye & Contributors
#'
#' @encoding UTF-8
#'
#' @importFrom stats qnorm
#' @importFrom graphics plot points lines polygon axis legend mtext text box par
#' @importFrom grDevices rgb col2rgb
#' @importFrom utils packageDescription
#'
#' @export
kobe_all_in_one <- function(rep,
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

  ## ------------------------- Local helpers (not exported) ------------------ ##
  check_rep <- function(rep, reportmode0 = TRUE) {
    if (!inherits(rep, "spictcls") || !"opt" %in% names(rep)) {
      stop("The argument 'rep' must be a fitted spict object (fit.spict()).")
    } else if (reportmode0 && rep$inp$reportmode != 0) {
      stop("All states must be reported! Set inp$reportmode <- 0 and refit.")
    }
    TRUE
  }

  get_par <- function(parname, rep, exp = FALSE, CI = 0.95) {
    if (CI <= 0 || CI >= 1) stop("CI must be in (0,1).")
    z <- qnorm(CI + (1 - CI) / 2)

    # Try to find the parameter in different slots
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

    # Nothing found yet?
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
        # Special case: annual production P
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
        warning("get_par: could not extract '", parname, "'. Returning NA.")
        est <- NA_real_; sdv <- NA_real_; ll <- NA_real_; ul <- NA_real_
      }
    }

    # Ensure all vectors have same length
    n <- length(est)
    if (length(sdv) == 0) sdv <- rep(NA_real_, n)
    if (length(ll)  == 0) ll  <- rep(NA_real_, n)
    if (length(ul)  == 0) ul  <- rep(NA_real_, n)

    if (isTRUE(exp)) {
      # log-normal transform with proper CI handling
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


  list_quantities <- function(rep) {
    sort(unique(c(names(rep$value), names(rep$par.fixed), names(rep$par.random))))
  }

  add_catchunit <- function(lab, cu) {
    cu <- as.character(cu)
    if (nzchar(cu)) eval(bquote(.(lab[[1]]) * ',' ~ .(cu))) else lab
  }

  make_rp_ellipse <- function(rep) {
    # Compute covariance of (logBmsy, logFmsy) from stored sd and sd of sum
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
      # Fallback: single point at (logBmsy, logFmsy)
      matrix(c(parBF[1], parBF[2]), ncol = 2,
             dimnames = list(NULL, c("x", "y")))
    }
  }

  calc_EBinf <- function(K, n, Fl, Fmsy, sdb2) {
    # Pella-Tomlinson equilibrium with process variance correction
    base <- 1 - (n - 1) / n * (Fl / Fmsy)
    base <- max(0, base)
    corr <- 1 - n / 2 / (1 - (1 - n * Fmsy + (n - 1) * Fl))
    EB  <- K * (base)^(1 / (n - 1)) * (1 - corr * sdb2)
    max(0, EB)
  }

  annual_avg <- function(intime, vec, type = "mean") {
    fun <- match.fun(type)
    anntime   <- unique(floor(intime))
    floortime <- floor(intime)
    nstepvec  <- vapply(anntime, function(a) sum(a == floortime), numeric(1))
    # Keep only years with full complement of within-year steps
    anntime <- anntime[which(nstepvec == max(nstepvec))]
    annvec  <- vapply(anntime, function(a) fun(vec[which(a == floortime)]), numeric(1))
    list(anntime = anntime, annvec = annvec)
  }

  get_EBinf <- function(rep) {
    K     <- get_par("logK",   rep, exp = TRUE)[2]
    n     <- get_par("logn",   rep, exp = TRUE)[2]
    sdb2  <- get_par("logsdb", rep, exp = TRUE)[2]^2
    Fmsy  <- tail(get_par("logFmsy", rep, exp = TRUE), 1)[2]
    logFs <- get_par("logFs", rep, CI = CI)
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

  add_man_lines <- function(rep, par, par2 = NULL, index.shift = 0,
                            plot.legend = TRUE, verbose = TRUE, ...) {
    if (!("man" %in% names(rep))) return(invisible(NULL))
    nman <- length(rep$man)
    errflag <- rep(FALSE, nman)

    for (i in seq_len(nman)) {
      rp <- rep$man[[i]]
      ip <- rp$inp
      manint <- ip$maninterval
      mandiff <- diff(manint)
      est <- get_par(par, rp, exp = TRUE)[, 2]
      dtcp <- ip$dtcp

      if (identical(par, "logCpred")) {
        time <- rp$inp$timeCpred
        lastobs <- ip$lastCatchObs
        indlastobs <- which(time == lastobs)
      } else {
        time <- rp$inp$time
        lastobs <- ip$timerangeObs[2]
        indlastobs <- which(time == lastobs)
      }

      indmanstart <- which(time >= manint[1])

      if (identical(par, "logCpred")) {
        if (mandiff < 1) {
          errflag[i] <- TRUE
          mantime <- NULL; manc <- NULL
        } else if (mandiff > 1) {
          mantime <- seq(manint[1], manint[2], 1)
          mantime <- mantime[-length(mantime)]
          manc <- rep(est[indmanstart] / mandiff, length(mantime))
        } else {
          mantime <- manint[1]
          manc <- est[indmanstart]
        }
        if (any(dtcp[indmanstart] < 1)) {
          alo <- annual_avg(time[indmanstart], est[indmanstart] / dtcp[indmanstart], mean)
          mantime <- alo$anntime
          manc <- alo$annvec
        }
        ind <- which(time < manint[1])
        preind <- tail(which(dtcp[ind] == 1), 1)
        if (any(dtcp[ind] < 1)) {
          alo <- annual_avg(time[ind], est[ind] / dtcp[ind], mean)
          time <- alo$anntime
          est <- alo$annvec
          ind <- which(time < manint[1])
          dtcp <- diff(unique(c(time, manint[1])))
          preind <- which.max(time[which(dtcp >= 1)])
        }
        premantime <- time[preind]
        premanc <- est[preind]
        x <- c(premantime, mantime); y <- c(premanc, manc)
      } else {
        maninds <- (indmanstart[1] - index.shift):tail(indmanstart, 1)
        x <- if (is.null(par2)) time[maninds] else get_par(par2, rp, exp = TRUE)[maninds, 2]
        y <- est[maninds]
      }
      lines(x, y, col = man_cols()[i], lwd = 1.5, ...)

      # Intermediate period
      intdiff <- manint[1] - lastobs
      if (intdiff > 0 && (!identical(par, "logCpred") || (identical(par, "logCpred") && intdiff %% 1 == 0))) {
        ind <- which(time >= lastobs & time < manint[1])
        if (identical(par, "logCpred")) {
          if (intdiff %% 1 == 0) {
            intinds <- (ind[1] - index.shift):preind
          } else intinds <- integer(0)
        } else {
          intinds <- (ind[1] - index.shift):maninds[1]
        }
        if (length(intinds)) {
          x <- if (is.null(par2)) time[intinds] else get_par(par2, rp, exp = TRUE)[intinds, 2]
          y <- est[intinds]
          lines(x, y, col = man_cols()[i], lwd = 1.5, lty = 3)
        }
      }
    }

    if (plot.legend) {
      leg <- names(rep$man)
      if (is.null(leg) || !length(leg)) leg <- paste0("Scenario ", seq_len(nman))
      legend("topleft", legend = leg, lty = 1, col = man_cols()[seq_len(nman)],
             bg = "transparent", cex = 0.8)
    }

    if (any(errflag) && verbose) {
      bad <- paste0(names(rep$man)[errflag], collapse = ", ")
      message("Note: management period of scenario(s) [", bad,
              "] < 1 year; annual catch display may be inaccurate.")
    }
    invisible(NULL)
  }

  get_version <- function(pkg = "spict") {
    pd <- utils::packageDescription(pkg)
    v <- paste0(pd$Package, "_v", pd$Version)
    if (is.null(pd$GithubRef)) v else paste0(v, "@", substr(pd$GithubSHA1, 1, 6))
  }

  warning_stamp <- function() {
    opar <- par(yaxt = "s", xaxt = "s", xpd = NA); on.exit(par(opar))
    usr <- par("usr")
    xcoord <- usr[1]
    ycoord <- usr[4] + 0.035 / diff(par()$fig[3:4]) * diff(usr[3:4])
    if (par("xlog")) xcoord <- 10^(xcoord)
    if (par("ylog")) ycoord <- 10^(ycoord)
    points(xcoord, ycoord, pch = 24, bg = "yellow", col = "black", cex = 2, lwd = 1.5)
    text(xcoord, ycoord, "!", cex = 0.8)
  }

  txt_stamp <- function(string = get_version(), cex = 0.5) {
    if (is.null(string) || !nzchar(string)) return(invisible(NULL))
    if (mean(par()$mfrow) > 1) return(invisible(NULL))
    opar <- par(new = TRUE, plt = c(0, 1, 0, 1), mfrow = c(1, 1), xpd = FALSE)
    on.exit(par(opar))
    plot(1, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
    plt <- par("plt"); usr <- par("usr")
    xcoord <- usr[2] + (usr[2] - usr[1]) / (plt[2] - plt[1]) * (1 - plt[2]) - 0.4 * strwidth("m")
    ycoord <- usr[3] - diff(usr[3:4]) / diff(plt[3:4]) * (plt[3]) + 0.4 * strheight("m")
    if (par("xlog")) xcoord <- 10^(xcoord)
    if (par("ylog")) ycoord <- 10^(ycoord)
    text(xcoord, ycoord, string, adj = 1, cex = cex)
    invisible(NULL)
  }
  ## ------------------------------------------------------------------------ ##

  ## Validate input
  check_rep(rep)

  ## Margins safeguard for single-plot devices
  if (!"sderr" %in% names(rep)) {
    mar <- c(5.1, 4.3, 4.1, 4.1)
    if (dev.cur() == 1L || (dev.cur() == 2L && all(par()$mfrow == c(1, 1)))) {
      opar <- par(mar = mar)
      on.exit(par(opar), add = TRUE)
    }
  }

  # Rel-axes forced if time-varying growth / mcovflag
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
    xlab <- expression(B[t]/B[MSY])
    ylab <- expression(F[t]/F[MSY])
  } else {
    bscal <- 1
    fscal <- 1
    xlab <- expression(B[t]); xlab <- add_catchunit(xlab, rep$inp$catchunit)
    ylab <- expression(F[t])
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

  # RP ellipse (or fallback)
  cl <- try(make_rp_ellipse(rep), silent = TRUE)
  if (inherits(cl, "try-error")) {
    cl <- matrix(c(log(Bmsy[2]), log(Fmsy[2])), ncol = 2)
  }

  # Build time series in annual steps if needed
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

  # Axis limits
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
    xlim[1] <- max(xlim[1], logminval)
    ylim[1] <- max(ylim[1], logminval)
  }

  # Plot
  if (!is.null(xlabel)) xlab <- xlabel
  logarg <- if (logax) "xy" else ""

  plot(Bmsy[2] / bscal, Fmsy[2] / fscal, type = "n",
       xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, log = logarg)

  # External relative axes (top/right)
  if (isTRUE(ext)) {
    if (logax) {
      expx <- pretty(log10(xlim / (Bmsy[2] / bscal)))
      expy <- pretty(log10(ylim / (Fmsy[2] / fscal)))
      labx <- 10^expx; laby <- 10^expy
    } else {
      labx <- pretty(xlim / (Bmsy[2] / bscal))
      laby <- pretty(ylim / (Fmsy[2] / fscal))
    }
    atx <- labx * (Bmsy[2] / bscal)
    aty <- laby * (Fmsy[2] / fscal)
    axis(3, labels = labx, at = atx)
    mtext(expression(B[t]/B[MSY]), side = 3, line = 2, cex = par("cex"))
    axis(4, labels = laby, at = aty)
    mtext(expression(F[t]/F[MSY]), side = 4, line = 2.5, cex = par("cex"))
  }

  # Kobe quadrants
  ymin <- if (logax) logminval * 1e-2 else -10
  xmin <- if (logax) logminval * 1e-2 else xlim[1] - xlim[2]
  # Green (B>=Bmsy, F<=Fmsy)
  polygon(c(Bmsy[2]/bscal, Bmsy[2]/bscal, xlim[2]*2, xlim[2]*2),
          c(Fmsy[2], ymin, ymin, Fmsy[2])/fscal, col = rgb(0.5, 0.8, 0.4, 1), border = NA)
  yel <- rgb(1, 0.925, 0.55, 1)
  # Yellow zones
  polygon(c(Bmsy[2]/bscal, Bmsy[2]/bscal, xmin, xmin),
          c(Fmsy[2], ymin, ymin, Fmsy[2])/fscal, col = yel, border = NA)
  polygon(c(Bmsy[2]/bscal, Bmsy[2]/bscal, xlim[2]*2, xlim[2]*2),
          c(Fmsy[2], (ylim[2]+1)*2, (ylim[2]+1)*2, Fmsy[2])/fscal, col = yel, border = NA)
  # Red (B<Bmsy, F>Fmsy)
  polygon(c(Bmsy[2]/bscal, Bmsy[2]/bscal, xmin, xmin),
          c(Fmsy[2], (ylim[2]+1)*2, (ylim[2]+1)*2, Fmsy[2])/fscal, col = rgb(1, 0.188, 0.188, 1), border = NA)
  abline(v = 0, col = "darkred", lty = 2)

  # RP ellipse (banana shape)
  cicol  <- "lightgray"; cicol2 <- "gray"
  cicolu <- do.call(rgb, as.list(c(as.numeric(col2rgb(cicol))/255, alpha = 0.7)))
  cicol2u<- do.call(rgb, as.list(c(as.numeric(col2rgb(cicol2))/255, alpha = 0.7)))
  polygon(exp(cl[, 1]) / bscal, exp(cl[, 2]) / fscal, col = cicolu, border = cicol2u)

  # True point if present
  if ("true" %in% names(rep$inp)) {
    points(rep$inp$true$Bmsy / bscal, rep$inp$true$Fmsy / fscal, pch = 25, bg = rgb(1, 165/255, 0, alpha = 0.7))
  }

  # Trajectory and management overlays
  maincol <- rgb(0, 0, 1, 0.8)
  lines(bbb, fff, col = maincol, lwd = 1.5)

  if (!("man" %in% names(rep))) {
    # one-step prediction path + EBinf marker
    if (!(min(rep$inp$dtc) < 1)) {
      lines(Best[rep$inp$indpred, 2] / bscal, Fest[rep$inp$indpred, 2] / fscal, col = maincol, lty = 3)
      Bll <- tail(Best[rep$inp$indpred, 2] / bscal, 1)
      Fll <- tail(Fest[rep$inp$indpred, 2] / fscal, 1)
      lines(c(Bll, EBinf), rep(Fll, 2), lwd = 1.5, lty = 3, col = "blue")
      points(EBinf, Fll, pch = 23, bg = "gold")
    }
  } else {
    add_man_lines(rep, "logF", par2 = "logB", index.shift = 1, plot.legend = man.legend, verbose = verbose)
  }

  # Previous vs current MSY markers (if multiple initial values)
  nr <- length(rep$inp$ini$logr)
  if (nr > 1) {
    points(Bmsy_all[1:(nr-1), 2] / bscal, Fmsy_all[1:(nr-1), 2] / fscal, pch = 3, col = "magenta")
    points(Bmsy_all[nr, 2] / bscal,      Fmsy_all[nr, 2] / fscal,      pch = 3, col = "black")
  }

  # Legends
  if (isTRUE(plot.legend)) {
    if (nr > 1) {
      legend("topright", c("Current MSY", "Previous MSY"), pch = 3,
             col = c("black", "magenta"), bg = "transparent")
    } else if ("true" %in% names(rep$inp)) {
      if (min(rep$inp$dtc) < 1) {
        legend("topright", "True", pch = 25, pt.bg = rgb(1,165/255,0,0.7), bg = "white")
      } else {
        legend("topright", c(expression('E(B'[infinity]*')'), "True"),
               pch = c(23, 25), pt.bg = c("gold", rgb(1,165/255,0,0.7)), bg = "transparent")
      }
    } else if (!(min(rep$inp$dtc) < 1)) {
      legend("topright", expression('E(B'[infinity]*')'), pch = 23, pt.bg = "gold", bg = "transparent")
    }
  }

  # Stamp first/last years
  points(bbb[1], fff[1], pch = 21, bg = "white")
  text(bbb[1], fff[1], round(fbtime[1], 2), pos = labpos[1], cex = 0.75, xpd = TRUE, offset = 0.5)
  points(Bl, Fl, pch = 22, bg = "white")
  text(Bl, Fl, round(tail(fbtime, 1), 2), pos = labpos[2], cex = 0.75, xpd = TRUE, offset = 0.5)

  box(lwd = 1.5)
  if (rep$opt$convergence != 0) warning_stamp()
  if (!is.null(stamp)) txt_stamp(stamp)

  invisible(NULL)
}
