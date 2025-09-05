#' Truth overlay color (semi-transparent orange)
#'
#' Returns a semi-transparent orange color used to draw "true" trajectories
#' and reference levels in simulation/OSA overlays.
#'
#' @return A length-1 character string with RGBA hex color.
#' @keywords internal
#' @noRd
true.col <- function() grDevices::rgb(1, 165/255, 0, alpha = 0.7)  # orange, semi-transp.

#' Add "truth" overlays to a panel ggplot
#'
#' Conditionally overlays simulated **true** trajectories / levels onto a
#' panel ggplot, based on fields present in \code{rep$inp$true}.
#' Supported panels: \code{"biomass"}, \code{"f"}, \code{"bbmsy"},
#' \code{"ffmsy"}, \code{"catch"}.
#'
#' - \strong{biomass}: draws true \code{B(t)} line and \code{Bmsy} (as a line or series).
#' - \strong{f}: draws true \code{F(t)} line and \code{Fmsy} (as a line or series).
#' - \strong{bbmsy}: draws true \code{B(t)/Bmsy(t)}.
#' - \strong{ffmsy}: draws true \code{F(t)/Fmsy(t)}.
#' - \strong{catch}: draws horizontal \code{MSY} line.
#'
#' If no \code{rep$inp$true} is present, the plot is returned unchanged.
#'
#' @param p A ggplot object built by a panel function.
#' @param rep A fitted \code{spictcls} object; may contain \code{rep$inp$true}.
#' @param panel Single string, one of \code{"biomass"}, \code{"f"},
#'   \code{"bbmsy"}, \code{"ffmsy"}, \code{"catch"}.
#'
#' @return The input ggplot \code{p} with overlays added (or unchanged).
#' @examples
#' \dontrun{
#' p <- plot_elu2_panel_biomass(rep)
#' p <- .elu2_overlay_truth(p, rep, panel = "biomass")
#' }
#' @keywords internal
#' @noRd
.elu2_overlay_truth <- function(p, rep, panel) {
  if (!("true" %in% names(rep$inp)) || is.null(rep$inp$true)) return(p)
  tr <- rep$inp$true
  has <- function(nm) !is.null(tr[[nm]]) && length(tr[[nm]]) > 0

  if (panel == "biomass") {
    if (has("time") && has("B")) {
      p <- p + ggplot2::geom_line(
        data = data.frame(time = tr$time, y = tr$B),
        ggplot2::aes(x = time, y = y),
        inherit.aes = FALSE, color = true.col(), linewidth = 0.8
      )
    }
    if (has("Bmsy")) {
      if (length(tr$Bmsy) == 1L) {
        p <- p + ggplot2::geom_hline(yintercept = tr$Bmsy, color = true.col(), linewidth = 0.8) +
          ggplot2::geom_hline(yintercept = tr$Bmsy, color = "black", linetype = "dotted")
      } else if (has("time")) {
        p <- p + ggplot2::geom_line(
          data = data.frame(time = tr$time, y = tr$Bmsy),
          ggplot2::aes(x = time, y = y),
          inherit.aes = FALSE, color = true.col(), linewidth = 0.8
        )
      }
    }
  }

  if (panel == "f") {
    if (has("time") && has("Fs")) {
      p <- p + ggplot2::geom_line(
        data = data.frame(time = tr$time, y = tr$Fs),
        ggplot2::aes(x = time, y = y),
        inherit.aes = FALSE, color = true.col(), linewidth = 0.8
      )
    }
    if (has("Fmsy")) {
      if (length(tr$Fmsy) == 1L) {
        p <- p + ggplot2::geom_hline(yintercept = tr$Fmsy, color = true.col(), linewidth = 0.8) +
          ggplot2::geom_hline(yintercept = tr$Fmsy, color = "black", linetype = "dotted")
      } else if (has("time")) {
        p <- p + ggplot2::geom_line(
          data = data.frame(time = tr$time, y = tr$Fmsy),
          ggplot2::aes(x = time, y = y),
          inherit.aes = FALSE, color = true.col(), linewidth = 0.8
        )
      }
    }
  }

  if (panel == "bbmsy") {
    if (has("time") && has("B") && has("Bmsy")) {
      p <- p + ggplot2::geom_line(
        data = data.frame(time = tr$time, y = tr$B / tr$Bmsy),
        ggplot2::aes(x = time, y = y),
        inherit.aes = FALSE, color = true.col(), linewidth = 0.8
      )
    }
  }

  if (panel == "ffmsy") {
    if (has("time") && has("Fs") && has("Fmsy")) {
      p <- p + ggplot2::geom_line(
        data = data.frame(time = tr$time, y = tr$Fs / tr$Fmsy),
        ggplot2::aes(x = time, y = y),
        inherit.aes = FALSE, color = true.col(), linewidth = 0.8
      )
    }
  }

  if (panel == "catch") {
    if (has("MSY")) {
      p <- p + ggplot2::geom_hline(yintercept = tr$MSY, color = true.col(), linewidth = 0.8) +
        ggplot2::geom_hline(yintercept = tr$MSY, color = "black", linetype = "dotted")
    }
  }

  p
}

#' Kobe panel (ggplot2): OSA on, MSY ellipse, scenarios, and left-aligned mini-legend
#'
#' Draws a phase plot of biomass vs. fishing mortality (absolute or relative to
#' MSY), including the MSY ellipse, historical and predicted trajectories,
#' optional management scenarios, and a compact, left-aligned mini-legend that
#' shows \eqn{E(B_\infty)} and the \strong{True} MSY marker when available.
#'
#' The function supports:
#' \itemize{
#'   \item Absolute axes (\eqn{B_t}, \eqn{F_t}) with optional secondary relative axes.
#'   \item Relative axes (\eqn{B_t/B_{MSY}}, \eqn{F_t/F_{MSY}}) via \code{rel.axes = TRUE}.
#'   \item Annual averaging when sub-annual time steps are present.
#'   \item Optional scenario overlays from \code{model$man} (solid in-management,
#'         dashed pre-management).
#'   \item A gold diamond for \eqn{E(B_\infty)} and an orange triangle (not in legend)
#'         for the true MSY point if \code{model$inp$true} exists.
#' }
#'
#' @param model A fitted \code{spictcls} object (\code{fit.spict()}) with
#'   \code{model$inp$reportmode == 0}.
#' @param logax Logical; use log10 scales for both axes. Default \code{FALSE}.
#' @param plot.legend Logical; draw the compact mini-legend
#'   (\eqn{E(B_\infty)} and \strong{True}). Default \code{TRUE}.
#' @param man.legend Logical; show scenario legend when \code{model$man} exists.
#'   Default \code{TRUE}.
#' @param ext Logical; add secondary relative axes when plotting absolute axes.
#'   Ignored when \code{rel.axes = TRUE}. Default \code{TRUE}.
#' @param rel.axes Logical; plot relative axes \eqn{B_t/B_{MSY}} vs \eqn{F_t/F_{MSY}}.
#'   Forces \code{ext = FALSE}. Default \code{FALSE}.
#' @param xlim,ylim Optional numeric length-2 limits for x and y. Computed if \code{NULL}.
#' @param labpos Numeric length-2; reserved for future label placement (currently unused).
#' @param xlabel Optional x-axis label override (expression or character).
#' @param stamp Optional character stamp (currently unused).
#' @param verbose Logical; show extraction fallbacks warnings. Default \code{TRUE}.
#' @param CI Confidence level for intervals (used in internal extractors). Default \code{0.95}.
#' @param print_it Logical; if \code{TRUE}, also prints the ggplot. Default \code{FALSE}.
#'
#' @return A \code{ggplot} object representing the Kobe panel.
#'
#' @section Dependencies:
#' Uses \pkg{ggplot2} geoms/scales and base \pkg{stats}/\pkg{grDevices}.
#' If \pkg{ellipse} is available, the MSY ellipse is drawn; otherwise a point
#' at \eqn{(B_{MSY}, F_{MSY})} is shown.
#'
#' @examples
#' \dontrun{
#' p <- plot_elu2_panel_kobe(fit, rel.axes = TRUE)
#' print(p)
#' }
#'
#' @import ggplot2
#' @importFrom grDevices adjustcolor rgb
#' @importFrom stats qnorm cov2cor
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

  # ----------------------- file-local helpers ---------------------- #
  .elu2_check_rep <- function(rep, reportmode0 = TRUE) {
    if (!inherits(rep, "spictcls") || !"opt" %in% names(rep)) {
      stop("The argument 'model' must be a fitted spict object (fit.spict()).")
    } else if (reportmode0 && rep$inp$reportmode != 0) {
      stop("All states must be reported! Set 'inp$reportmode <- 0' and refit.")
    }
    TRUE
  }
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
  .elu2_add_catchunit <- function(lab, cu) {
    cu <- as.character(cu)
    if (nzchar(cu)) eval(bquote(.(lab[[1]]) * ',' ~ .(cu))) else lab
  }
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
  .elu2_annual_avg <- function(intime, vec, type = "mean") {
    fun <- match.fun(type)
    anntime   <- unique(floor(intime))
    floortime <- floor(intime)
    nstepvec  <- vapply(anntime, function(a) sum(a == floortime), numeric(1))
    anntime   <- anntime[which(nstepvec == max(nstepvec))]
    annvec    <- vapply(anntime, function(a) fun(vec[which(a == floortime)]), numeric(1))
    list(anntime = anntime, annvec = annvec)
  }
  .elu2_calc_EBinf <- function(K, n, Fl, Fmsy, sdb2) {
    base <- 1 - (n - 1)/n * (Fl/Fmsy)
    base <- max(0, base)
    corr <- 1 - n/2 / (1 - (1 - n*Fmsy + (n - 1)*Fl))
    EB   <- K * (base)^(1/(n - 1)) * (1 - corr*sdb2)
    max(0, EB)
  }
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
  .elu2_man_cols <- function() {
    colvec <- c('darkmagenta','cyan3','darkgreen','coral1','black','magenta',
                'gold','green','cadetblue3','chocolate3','darkolivegreen3',
                'cyan','darkred')
    rep(colvec, 3)
  }
  .elu2_fmt1 <- function(x) { x <- as.numeric(x); ifelse(is.finite(x),
                                                         formatC(x, format = "f", digits = 1),
                                                         NA_character_) }
  .elu2_lab_rel_1 <- function(vals) {
    out <- .elu2_fmt1(vals); out[abs(vals - 1) < 1e-9] <- "1"; out
  }
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

  pad_frac <- 0.02; xpad <- pad_frac * diff(xlim); ypad <- pad_frac * diff(ylim)
  xlim <- c(max(if (isTRUE(logax)) logminval else 0, xlim[1] - xpad), xlim[2] + xpad)
  ylim <- c(max(if (isTRUE(logax)) logminval else 0, ylim[1] - ypad), ylim[2] + ypad)

  # ------------------------- Data frames for ggplot ------------------------- #
  xminq <- xlim[1]; xmaxq <- xlim[2]; yminq <- ylim[1]; ymaxq <- ylim[2]
  bx <- if (rel.axes) 1 else Bmsy[2]; bx <- .elu2_clamp(bx, xminq, xmaxq)
  fy <- if (rel.axes) 1 else Fmsy[2]; fy <- .elu2_clamp(fy, yminq, ymaxq)

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

  # True MSY point (plotted in the phase space; no ggplot legend)
  if (!is.null(df_true) && nrow(df_true)) {
    p <- p + ggplot2::geom_point(
      data = df_true,
      ggplot2::aes(x = x, y = y),
      shape = 25, fill = true.col(),
      color = "black", size = 3, stroke = 0.5, show.legend = FALSE
    )
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

  # Axes
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

  # Scenario legend (color)
  if (!is.null(df_man) && nrow(df_man) && man.legend) {
    uniq_sc <- unique(df_man$scenario)
    colv <- .elu2_man_cols()[seq_along(uniq_sc)]; names(colv) <- uniq_sc
    p <- p + ggplot2::scale_color_manual(values = colv, name = "Scenario")
  } else {
    p <- p + ggplot2::scale_color_discrete(guide = "none")
  }

  # Theme
  show_scen_legend <- (!is.null(df_man) && nrow(df_man) && man.legend)
  legend_pos <- if (show_scen_legend) "top" else "none"

  p <- p + ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "gray25", fill = NA, linewidth = 2),
      legend.background = ggplot2::element_rect(fill = grDevices::adjustcolor("white", alpha.f = 0.8), color = NA),
      legend.position = legend_pos,
      legend.direction = "horizontal",
      axis.text.x  = ggplot2::element_text(size = 10, face = "bold", color = "gray25"),
      axis.text.y  = ggplot2::element_text(size = 10, face = "bold", color = "gray25"),
      axis.title.x = ggplot2::element_text(size = 14, face = "bold", color = "gray25"),
      axis.title.y = ggplot2::element_text(size = 14, face = "bold", color = "gray25"),
      plot.margin = ggplot2::margin(4, 4, 4, 4)
    )

  # === Mini-legend for E(B∞) + "True" (left-aligned text; tighter spacing) ===
  if (plot.legend) {
    dx <- diff(xlim); dy <- diff(ylim)

    # a little closer to the right edge than before
    px_right <- xlim[2] - 0.015 * dx

    # tighter horizontal spacing (marker from right edge + text from marker)
    gapx   <- 0.1 * dx
    txtpad <- 0.02 * dx

    # tighter vertical spacing between the two rows
    vgap <- 0.06 * dy

    x_marker <- px_right - gapx
    x_text   <- x_marker + txtpad  # LEFT-ALIGNED label anchor for both rows

    py_top <- ylim[2] - 0.04 * dy

    # E(B∞) mini-legend (only when seasonal off), LEFT-aligned text
    eb_inf_shown <- !(min(rep$inp$dtc) < 1)
    if (eb_inf_shown) {
      p <- p +
        ggplot2::annotate("point",
                          x = x_marker, y = py_top,
                          shape = 23, size = 2.4, fill = "gold", color = "black", stroke = 0.6
        ) +
        ggplot2::annotate("text",
                          x = x_text, y = py_top,
                          label = "bold(E(B[infinity]))", parse = TRUE,
                          hjust = 0, vjust = 0.5, size = 4
        )
    }

    # True MSY mini-legend entry (if available), text LEFT-aligned and bold
    if ("true" %in% names(model$inp) && !is.null(model$inp$true)) {
      y_true <- if (eb_inf_shown) py_top - vgap else py_top
      p <- p +
        ggplot2::annotate("point",
                          x = x_marker, y = y_true,
                          shape = 25, size = 2.4, fill = true.col(), color = "black", stroke = 0.6
        ) +
        ggplot2::annotate("text",
                          x = x_text, y = y_true,
                          label = "True",
                          hjust = 0, vjust = 0.5, size = 4, fontface = "bold"
        )
    }
  }

  if (isTRUE(print_it)) print(p)
  p
}
