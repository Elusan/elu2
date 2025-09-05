#' Plot SPiCT Scenario (one model) — 6 panels with Kobe (OSA on, Original-style mini-legend)
#'
#' Produces a 3×2 patchwork for a single fitted SPiCT model: Biomass, B/Bmsy, Catch,
#' F, F/Fmsy, and a **Kobe** panel. Confidence intervals are drawn as **dotted** lines
#' and, where relevant, **filled ribbons** (B/Bmsy-on-B scale, B/Bmsy, F/Fmsy, MSY band,
#' Absolute F band) when `show_CIs = TRUE`. The biomass panel overlays observed indices
#' scaled by \eqn{\hat q}. Thin **solid grey** vertical lines mark the end of observed data
#' (SPiCT convention for Catch).
#'
#' The Kobe panel is built by calling an internal, non-printing copy of
#' `Original_kobe_all_in_gg()` using the **same** `model` (with `$man` stripped),
#' so the OSA segment + E(B∞) marker are shown and the mini-legend matches
#' `Original_kobe_all_in_gg()`.
#'
#' @param model A fitted SPiCT object (class `spictcls`).
#' @param extract_catch_data Optional function returning a data.frame with columns
#'   `time, catch, lwr, upr, catch_type` (Observed/Predicted). If `NULL`, a built-in
#'   extractor is used.
#' @param line_color Single color for the series (default `#1b9e77`).
#' @param return_patchwork If `TRUE` (default), return a patchwork object; if `FALSE`,
#'   return a named list of ggplot objects.
#' @param show_CIs If `TRUE` (default), draw dotted lower/upper CI lines and filled ribbons
#'   (where applicable). If `FALSE`, CIs and ribbons are suppressed.
#'
#' @return A patchwork object (default) or `list(biomass, bbmsy, catch, f, ffmsy, kobe)`.
#' @export
plot_spict_scenarios_by_model_NEW4_one_MODEL <- function(model,
                                                         extract_catch_data = NULL,
                                                         line_color = "#1b9e77",
                                                         return_patchwork = TRUE,
                                                         show_CIs = TRUE) {
  stopifnot(inherits(model, "spictcls"))
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("Package 'patchwork' is required.")
  library(ggplot2)
  library(patchwork)
  library(grid)

  ## ─────────────────────────── Kobe (internal copy, non-printing) ───────────────────────────
  Original_kobe_all_in_gg_internal <- function(rep,
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
                                               CI = 0.95,
                                               print_it = FALSE) {
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
      z <- stats::qnorm(CI + (1 - CI) / 2)
      ind_ran <- which(names(rep$par.random) == parname)
      ind_fix <- which(names(rep$par.fixed)  == parname)
      ind_sdr <- which(names(rep$value)      == parname)
      ind_opt <- which(names(rep$opt$par)    == parname)
      est <- ll <- ul <- sdv <- NULL
      if (length(ind_ran)) { est <- rep$par.random[ind_ran]; sdv <- sqrt(rep$diag.cov.random[ind_ran]); ll <- est - z*sdv; ul <- est + z*sdv }
      if (length(ind_fix)) { est <- rep$par.fixed[ind_fix];  sdv <- sqrt(diag(rep$cov.fixed))[ind_fix]; ll <- est - z*sdv; ul <- est + z*sdv }
      if (length(ind_sdr)) { est <- rep$value[ind_sdr];       sdv <- rep$sd[ind_sdr];                   ll <- est - z*sdv; ul <- est + z*sdv }
      if (length(est) == 0) {
        if (length(ind_opt)) {
          est <- rep$opt$par[ind_opt]; sdv <- ll <- ul <- rep(NA_real_, length(est))
        } else if ("phases" %in% names(rep$inp) && parname %in% names(rep$inp$phases) && rep$inp$phases[[parname]] == -1) {
          est <- rep$inp$parlist[[parname]]; sdv <- rep(0, length(est)); ll <- est; ul <- est
        } else if (!is.na(parname) && identical(parname, "P")) {
          B <- get_par("logB", rep, exp = TRUE, CI = CI); C <- get_par("logCpred", rep, exp = TRUE, CI = CI)
          ic <- rep$inp$ic; nc <- rep$inp$nc; B0 <- B[ic, 2]; B1 <- B[ic + nc, 2]; T0 <- rep$inp$time[ic]; T1 <- rep$inp$time[ic + nc]
          est <- (B1 - B0 + C[, 2]) / (T1 - T0); sdv <- ll <- ul <- rep(NA_real_, length(est))
        } else {
          if (verbose) warning("get_par: could not extract '", parname, "'. Returning NA.")
          est <- sdv <- ll <- ul <- NA_real_
        }
      }
      n <- length(est); if (length(sdv) == 0) sdv <- rep(NA_real_, n); if (length(ll) == 0) ll <- rep(NA_real_, n); if (length(ul) == 0) ul <- rep(NA_real_, n)
      if (isTRUE(exp)) { cv <- ifelse(is.finite(sdv), sqrt(exp(sdv^2) - 1), NA_real_); ll <- exp(ll); ul <- exp(ul); est <- exp(est); ul[is.infinite(ul)] <- exp(705) }
      else             { cv <- sdv / est }
      out <- cbind(ll, est, ul, sdv, cv)
      if (parname %in% c("logB","logF","logBBmsy","logFFmsy")) rownames(out) <- rep$inp$time
      out
    }
    add_catchunit <- function(lab, cu) { cu <- as.character(cu); if (nzchar(cu)) eval(bquote(.(lab[[1]]) * ',' ~ .(cu))) else lab }
    make_rp_ellipse <- function(rep) {
      inds <- c(max(which(names(rep$value) == "logBmsy")), max(which(names(rep$value) == "logFmsy")))
      sds  <- rep$sd[inds]; s_sum <- rep$sd[which(names(rep$value) == "logBmsyPluslogFmsy")]
      cova <- (s_sum^2 - sds[1]^2 - sds[2]^2) / 2
      covBF <- matrix(c(sds[1]^2, cova, cova, sds[2]^2), 2, 2, byrow = TRUE); parBF <- rep$value[inds]
      if (requireNamespace("ellipse", quietly = TRUE)) ellipse::ellipse(stats::cov2cor(covBF)[1, 2], scale = sqrt(diag(covBF)), centre = parBF, npoints = 300)
      else matrix(c(parBF[1], parBF[2]), ncol = 2, dimnames = list(NULL, c("x","y")))
    }
    annual_avg <- function(intime, vec, type = "mean") {
      fun <- match.fun(type); anntime <- unique(floor(intime)); floortime <- floor(intime)
      nstepvec <- vapply(anntime, function(a) sum(a == floortime), numeric(1))
      anntime  <- anntime[which(nstepvec == max(nstepvec))]
      annvec   <- vapply(anntime, function(a) fun(vec[which(a == floortime)]), numeric(1))
      list(anntime = anntime, annvec = annvec)
    }
    calc_EBinf <- function(K, n, Fl, Fmsy, sdb2) {
      base <- 1 - (n - 1)/n * (Fl/Fmsy); base <- max(0, base)
      corr <- 1 - n/2 / (1 - (1 - n*Fmsy + (n - 1)*Fl))
      EB  <- K * (base)^(1/(n - 1)) * (1 - corr*sdb2); max(0, EB)
    }
    get_EBinf <- function(rep) {
      K <- get_par("logK", rep, exp = TRUE)[2]; n <- get_par("logn", rep, exp = TRUE)[2]; sdb2 <- get_par("logsdb", rep, exp = TRUE)[2]^2
      Fmsy <- tail(get_par("logFmsy", rep, exp = TRUE), 1)[2]; logFs <- get_par("logFs", rep, CI = CI)
      if (min(rep$inp$dtc) < 1) { alf <- annual_avg(rep$inp$time, logFs[, 2]); fff <- exp(alf$annvec) } else { fff <- exp(logFs[, 2]) }
      Fl <- tail(unname(fff), 1); calc_EBinf(K, n, Fl, Fmsy, sdb2)
    }
    man_cols <- function() {
      colvec <- c('darkmagenta','cyan3','darkgreen','coral1','black','magenta','gold','green','cadetblue3','chocolate3','darkolivegreen3','cyan','darkred')
      rep(colvec, 3)
    }
    fmt1 <- function(x) { x <- as.numeric(x); ifelse(is.finite(x), formatC(x, format = "f", digits = 1), NA_character_) }
    lab_rel_1 <- function(vals) { out <- fmt1(vals); out[abs(vals - 1) < 1e-9] <- "1"; out }
    clamp <- function(val, lo, hi) pmin(pmax(val, lo), hi)

    ## Validate & derive
    check_rep(rep)
    tvgflag <- rep$inp$timevaryinggrowth | rep$inp$logmcovflag
    if (tvgflag) rel.axes <- TRUE

    Bmsy_all <- get_par("logBmsy", rep, exp = TRUE, CI = CI)
    Fmsy_all <- get_par("logFmsy", rep, exp = TRUE, CI = CI)
    Bmsy <- tail(Bmsy_all, 1L); Fmsy <- tail(Fmsy_all, 1L)

    if (rel.axes) { ext <- FALSE; bscal <- Bmsy[2]; fscal <- Fmsy[2]; xlab_expr <- expression(B[t]/B[MSY]); ylab_expr <- expression(F[t]/F[MSY]) }
    else          { bscal <- 1; fscal <- 1; xlab_expr <- add_catchunit(expression(B[t]), rep$inp$catchunit); ylab_expr <- expression(F[t]) }

    Best <- get_par("logB", rep, exp = TRUE, CI = CI)
    logBest <- get_par("logB", rep, CI = CI)
    if (tvgflag) { Fest <- get_par("logFFmsy", rep, exp = TRUE, CI = CI); fscal <- 1; Fmsy <- c(1,1) }
    else         { Fest <- get_par("logFs",    rep, exp = TRUE, CI = CI) }
    logFest <- get_par("logFs", rep, CI = CI)

    cl <- try(make_rp_ellipse(rep), silent = TRUE)
    if (inherits(cl,"try-error")) cl <- matrix(c(log(Bmsy[2]), log(Fmsy[2])), ncol = 2)

    if (min(rep$inp$dtc) < 1) {
      alb <- annual_avg(rep$inp$time, logBest[,2]); alf <- annual_avg(rep$inp$time, logFest[,2])
      bbb <- exp(alb$annvec) / bscal; fff <- exp(alf$annvec) / fscal; fbtime <- alb$anntime
    } else {
      bbb <- Best[rep$inp$indest, 2] / bscal; fff <- Fest[rep$inp$indest, 2] / fscal; fbtime <- rep$inp$time[rep$inp$indest]
    }

    Fl <- tail(unname(fff), 1); Bl <- tail(unname(bbb), 1); EBinf <- get_EBinf(rep) / bscal

    if (is.null(xlim)) {
      xlim <- range(c(exp(cl[,1]), Best[,2], EBinf) / bscal, na.rm = TRUE)
      if (min(rep$inp$dtc) < 1) xlim <- range(c(exp(alb$annvec), exp(cl[,1]), EBinf) / bscal, na.rm = TRUE)
      xlim[2] <- min(c(xlim[2], 8*Bmsy[2]/bscal), 2.2*max(bbb), na.rm = TRUE); xlim[2] <- max(c(xlim[2], Bmsy[2]/bscal), na.rm = TRUE)
    }
    if (is.null(ylim)) {
      ylim <- range(c(exp(cl[,2]), Fest[,2]) / fscal, na.rm = TRUE)
      if (min(rep$inp$dtc) < 1) ylim <- range(c(exp(logFest[,2]) / fscal, exp(cl[,2]) / fscal), na.rm = TRUE)
      ylim[2] <- min(c(ylim[2], 8*Fmsy[2]/fscal), 2.2*max(fff), na.rm = TRUE); ylim[2] <- max(c(ylim[2], Fmsy[2]/fscal), na.rm = TRUE)
      if ("man" %in% names(rep)) ylim <- range(ylim, 0)
    }
    logminval <- 1e-4
    if (isTRUE(logax)) { if (xlim[1] < logminval) xlim[1] <- logminval; if (ylim[1] < logminval) ylim[1] <- logminval }

    pad_frac <- 0.02; xpad <- pad_frac * diff(xlim); ypad <- pad_frac * diff(ylim)
    xlim <- c(max(if (isTRUE(logax)) logminval else 0, xlim[1] - xpad), xlim[2] + xpad)
    ylim <- c(max(if (isTRUE(logax)) logminval else 0, ylim[1] - ypad), ylim[2] + ypad)

    xminq <- xlim[1]; xmaxq <- xlim[2]; yminq <- ylim[1]; ymaxq <- ylim[2]
    bx <- if (rel.axes) 1 else Bmsy[2] / (if (rel.axes) Bmsy[2] else 1); bx <- pmin(pmax(bx, xminq), xmaxq)
    fy <- if (rel.axes) 1 else Fmsy[2] / (if (rel.axes) Fmsy[2] else 1); fy <- pmin(pmax(fy, yminq), ymaxq)

    df_green_rect   <- data.frame(xmin = bx,   xmax =  Inf, ymin = -Inf, ymax =  fy)
    df_yellowL_rect <- data.frame(xmin = -Inf, xmax =  bx,  ymin = -Inf, ymax =  fy)
    df_yellowT_rect <- data.frame(xmin = bx,   xmax =  Inf, ymin =  fy,  ymax =  Inf)
    df_red_rect     <- data.frame(xmin = -Inf, xmax =  bx,  ymin =  fy,  ymax =  Inf)

    df_ellipse <- if (ncol(as.matrix(cl)) == 2 && nrow(as.matrix(cl)) > 1) {
      data.frame(x = exp(cl[,1]) / (if (rel.axes) Bmsy[2] else 1), y = exp(cl[,2]) / (if (rel.axes) Fmsy[2] else 1))
    } else {
      data.frame(x = exp(cl[1]) / (if (rel.axes) Bmsy[2] else 1), y = exp(cl[2]) / (if (rel.axes) Fmsy[2] else 1))
    }

    epsx <- 0.005*diff(xlim); epsy <- 0.005*diff(ylim)
    clamp_in <- function(v, lo, hi, eps) pmin(pmax(v, lo + eps), hi - eps)

    n_hist <- length(bbb)
    df_traj <- data.frame(x = clamp_in(bbb, xminq, xmaxq, epsx), y = clamp_in(fff, yminq, ymaxq, epsy), ord = seq_len(n_hist))

    df_pred <- NULL; df_EBseg <- NULL
    if (!(min(rep$inp$dtc) < 1) && !("man" %in% names(rep))) {
      xpred <- get_par("logB", rep, exp = TRUE)[rep$inp$indpred, 2] / (if (rel.axes) Bmsy[2] else 1)
      ypred <- get_par("logFs", rep, exp = TRUE)[rep$inp$indpred, 2] / (if (rel.axes) Fmsy[2] else 1)
      df_pred <- data.frame(x = xpred, y = ypred)
      Bll <- tail(xpred, 1); Fll <- tail(ypred, 1); EBinf <- get_EBinf(rep) / (if (rel.axes) Bmsy[2] else 1)
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
      nman <- length(rep$man); leg_man <- names(rep$man); if (is.null(leg_man) || !length(leg_man)) leg_man <- paste0("Scenario ", seq_len(nman))
      df_list <- vector("list", nman); df_int_list <- vector("list", nman)
      for (i in seq_len(nman)) {
        rp <- rep$man[[i]]; estF <- get_par("logF", rp, exp = TRUE)[,2]; estB <- get_par("logB", rp, exp = TRUE)[,2]
        time <- rp$inp$time; manint <- rp$inp$maninterval; indmanstart <- which(time >= manint[1])
        if (length(indmanstart)) { maninds <- indmanstart[1]:tail(indmanstart, 1)
        df_list[[i]] <- data.frame(x = estB[maninds] / (if (rel.axes) Bmsy[2] else 1),
                                   y = estF[maninds] / (if (rel.axes) Fmsy[2] else 1),
                                   scenario = leg_man[i]) }
        else df_list[[i]] <- data.frame(x = numeric(0), y = numeric(0), scenario = character(0))
        lastobs <- rp$inp$timerangeObs[2]; ind <- which(time >= lastobs & time < manint[1])
        if (length(ind)) {
          df_int_list[[i]] <- data.frame(x = estB[ind] / (if (rel.axes) Bmsy[2] else 1),
                                         y = estF[ind] / (if (rel.axes) Fmsy[2] else 1),
                                         scenario = leg_man[i])
        } else df_int_list[[i]] <- data.frame(x = numeric(0), y = numeric(0), scenario = character(0))
      }
      df_man <- do.call(rbind, df_list); df_man_int <- do.call(rbind, df_int_list)
    }

    p <- ggplot()
    p <- p + geom_rect(data = df_green_rect,   aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                       fill = grDevices::adjustcolor("darkseagreen3", alpha.f = 1), color = NA)
    p <- p + geom_rect(data = df_yellowL_rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                       fill = grDevices::adjustcolor("khaki1", alpha.f = 1), color = NA)
    p <- p + geom_rect(data = df_yellowT_rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                       fill = grDevices::adjustcolor("khaki1", alpha.f = 1), color = NA)
    p <- p + geom_rect(data = df_red_rect,     aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                       fill = grDevices::adjustcolor("indianred1", alpha.f = 1), color = NA)

    if (nrow(df_ellipse) > 1) {
      p <- p + geom_polygon(data = df_ellipse, aes(x = x, y = y),
                            fill = grDevices::adjustcolor("gray80", alpha.f = 0.8),
                            color = grDevices::adjustcolor("gray80", alpha.f = 0.8), linewidth = 0.3)
    } else {
      p <- p + geom_point(data = df_ellipse, aes(x = x, y = y),
                          color = "gray40", size = 2)
    }

    p <- p + geom_path(data = df_traj, aes(x = x, y = y),
                       linewidth = 0.4, color = grDevices::adjustcolor("blue", alpha.f = 0.8))

    if (!is.null(df_pred)) {
      p <- p + geom_path(data = df_pred, aes(x = x, y = y),
                         linetype = "11", color = grDevices::adjustcolor("blue", alpha.f = 0.8))
    }
    if (!is.null(df_EBseg)) {
      p <- p + geom_segment(data = df_EBseg,
                            aes(x = x, xend = xend, y = y, yend = yend),
                            color = "blue", linetype = "11", linewidth = 0.7) +
        geom_point(data = data.frame(x = df_EBseg$xend, y = df_EBseg$yend),
                   aes(x = x, y = y), shape = 23, fill = "gold", color = "black", size = 2, stroke = 0.5)
    }

    if (!is.null(df_man) && nrow(df_man)) {
      p <- p + geom_path(data = df_man, aes(x = x, y = y, color = scenario),
                         linewidth = 0.7, inherit.aes = FALSE, show.legend = man.legend)
    }
    if (!is.null(df_man_int) && nrow(df_man_int)) {
      p <- p + geom_path(data = df_man_int, aes(x = x, y = y, color = scenario),
                         linewidth = 0.7, linetype = "33", inherit.aes = FALSE, show.legend = man.legend)
    }

    if (!is.null(df_msy_prev) && nrow(df_msy_prev)) p <- p + geom_point(data = df_msy_prev, aes(x = x, y = y), shape = 3, color = "magenta", size = 2)
    if (!is.null(df_msy_curr) && nrow(df_msy_curr)) p <- p + geom_point(data = df_msy_curr, aes(x = x, y = y), shape = 3, color = "black", size = 2)

    if (!is.null(df_true)) {
      p <- p + geom_point(data = df_true, aes(x = x, y = y),
                          shape = 25, fill = grDevices::adjustcolor(rgb(1, 165/255, 0), alpha.f = 0.7),
                          color = "black", size = 3)
    }

    p <- p + geom_point(data = df_first, aes(x = x, y = y),
                        shape = 21, fill = "white", color = "black", size = 2, stroke = 0.5) +
      geom_text(data = df_first, aes(x = x, y = y, label = lab), vjust = -0.8, size = 3) +
      geom_point(data = df_last, aes(x = x, y = y),
                 shape = 22, fill = "white", color = "black", size = 2, stroke = 0.5) +
      geom_text(data = df_last, aes(x = x, y = y, label = lab), vjust = -0.8, size = 3)

    if (!isTRUE(logax) && (0 >= xlim[1]) && (0 <= xlim[2])) p <- p + geom_vline(xintercept = 0, color = "darkred", linetype = 2)

    xlab_final <- if (is.null(xlabel)) xlab_expr else xlabel; ylab_final <- ylab_expr

    if (isTRUE(logax)) {
      p <- p + scale_x_log10(limits = xlim, name = xlab_final,
                             labels = function(x) formatC(x, format = "f", digits = 0), expand = c(0, 0)) +
        scale_y_log10(limits = ylim, name = ylab_final,
                      labels = fmt1, expand = c(0, 0))
    } else {
      p <- p + scale_x_continuous(limits = xlim, name = xlab_final,
                                  labels = function(x) formatC(x, format = "f", digits = 0), expand = c(0, 0)) +
        scale_y_continuous(limits = ylim, name = ylab_final,
                           labels = fmt1, expand = c(0, 0))
    }

    if (!is.null(df_man) && nrow(df_man) && man.legend) {
      uniq_sc <- unique(df_man$scenario)
      colv <- c('darkmagenta','cyan3','darkgreen','coral1','black','magenta','gold','green','cadetblue3','chocolate3','darkolivegreen3','cyan','darkred')
      colv <- rep(colv, 3)[seq_along(uniq_sc)]; names(colv) <- uniq_sc
      p <- p + scale_color_manual(values = colv, name = "Scenario")
    } else p <- p + scale_color_discrete(guide = "none")

    p <- p + theme_bw() +
      theme(panel.grid = element_blank(),
            panel.border = element_rect(color = "gray25", fill = NA, linewidth = 0.6),
            legend.background = element_rect(fill = grDevices::adjustcolor("white", alpha.f = 0.8), color = NA),
            legend.position = if (man.legend) "top" else "none",
            legend.direction = "horizontal",
            axis.text.x  = element_text(size = 10, face = "bold", color = "gray25"),
            axis.text.y  = element_text(size = 10, face = "bold", color = "gray25"),
            axis.title.x = element_text(size = 14, face = "bold", color = "gray25"),
            axis.title.y = element_text(size = 14, face = "bold", color = "gray25"),
            plot.margin = margin(4, 4, 4, 4))

    if (plot.legend && !(min(rep$inp$dtc) < 1)) {
      px <- xlim[2] - 0.02 * diff(xlim); py <- ylim[2] - 0.04 * diff(ylim); gap <- 0.12 * diff(xlim)
      p <- p + annotate("point", x = px - gap, y = py, shape = 23, size = 2.5, fill = "gold", color = "black", stroke = 0.6) +
        annotate("text", x = px, y = py, label = "bold(E(B[infinity]))", parse = TRUE, hjust = 1, vjust = 0.5, size = 4)
    }

    if (rep$opt$convergence != 0) {
      p <- p + annotate("text", x = xlim[1] + 0.01 * diff(xlim), y = ylim[2] + 0.02 * diff(ylim),
                        label = "▲", color = "black", size = 6) +
        annotate("text", x = xlim[1] + 0.01 * diff(xlim), y = ylim[2] + 0.02 * diff(ylim),
                 label = "!", color = "black", size = 3, vjust = 0.35)
    }

    if (!is.null(stamp) && nzchar(stamp)) {
      p <- p + annotate("text", x = xlim[2] - 0.01 * diff(xlim), y = ylim[1] - 0.04 * diff(ylim),
                        label = stamp, hjust = 1, vjust = 0, size = 3)
    }

    if (isTRUE(print_it)) print(p)
    return(p)
  }
  ## ────────────────────────────────────────────────────────────────────────────

  ## Shared styles & helpers for time-series panels
  vline_col  <- "grey50"
  vline_size <- 0.2
  theme_minimal_compact2 <- function(base_size = 10, base_family = "") {
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

  ## Observation end helpers
  get_obs_end_generic <- function(m) {
    tr <- m$inp$timerangeObs
    if (!is.null(tr) && length(tr) >= 2 && is.finite(tr[2])) return(tr[2])
    if (!is.null(m$inp$timeC) && length(m$inp$timeC))       return(tail(m$inp$timeC, 1))
    if (!is.null(m$inp$timeI) && length(m$inp$timeI) && length(m$inp$timeI[[1]]) > 0) return(tail(m$inp$timeI[[1]], 1))
    if (!is.null(m$inp$time)  && length(m$inp$time))        return(max(m$inp$time, na.rm = TRUE))
    NA_real_
  }
  obs_end_overall <- get_obs_end_generic(model)

  get_catch_obs_end <- function(m) {
    inp <- m$inp
    if (is.null(inp$timeC) || !length(inp$timeC)) return(NA_real_)
    if (!is.null(inp$dtc) && length(inp$dtc) && min(inp$dtc, na.rm = TRUE) < 1) {
      ix <- which((inp$timeC %% 1) == 0)
      if (length(ix)) return(tail(inp$timeC[ix], 1))
      return(tail(inp$timeC, 1))
    } else return(tail(inp$timeC, 1))
  }
  catch_obs_end <- get_catch_obs_end(model)

  ## Convenience extractors
  ## --- replace your current get_series() with this one ---
  get_series <- function(parname) {
    par <- get.par(parname, model, exp = TRUE)

    # turn into matrix safely
    pm <- as.matrix(par)

    # robust time: prefer rownames if they parse cleanly, otherwise fall back to inp$time
    t_row <- suppressWarnings(as.numeric(rownames(pm)))
    if (length(t_row) != nrow(pm) || all(is.na(t_row))) {
      # fall back to model$inp$time, truncated/expanded to nrow(pm) if needed
      t_full <- as.numeric(model$inp$time)
      if (length(t_full) >= nrow(pm)) {
        t_row <- t_full[seq_len(nrow(pm))]
      } else {
        # in odd cases, recycle last time to match rows
        t_row <- rep_len(t_full, nrow(pm))
      }
    }

    # standardize column names
    colnames(pm) <- c("lwr","est","upr","sd","cv")[seq_len(ncol(pm))]

    data.frame(
      time = t_row,
      lwr  = pm[, "lwr"],
      est  = pm[, "est"],
      upr  = pm[, "upr"],
      sd   = if ("sd" %in% colnames(pm))  pm[, "sd"] else NA_real_,
      cv   = if ("cv" %in% colnames(pm))  pm[, "cv"] else NA_real_
    )
  }


  ## Indices for in-sample & prediction portions (match SPiCT convention)
  CI <- 0.95
  manflag <- ("man" %in% names(model))
  inp <- model$inp
  ind_in <- inp$indest
  if (length(ind_in) > 0) ind_in <- ind_in[seq_len(length(ind_in) - 1)]
  ind_pr <- if (!manflag) inp$indpred else integer(0)

  make_ribbon_df <- function(time, lwr, upr) {
    ok <- is.finite(lwr) & is.finite(upr) & is.finite(time)
    data.frame(time = time[ok], ymin = lwr[ok], ymax = upr[ok])
  }

  ## q-scaled index points (first two indices if available)
  get_index_points <- function() {
    out <- list()
    if (!is.null(inp$timeI) && length(inp$timeI) >= 1) {
      qest <- try(get.par("logq", model, exp = TRUE), silent = TRUE)
      if (!inherits(qest, "try-error")) {
        for (i in seq_along(inp$timeI)) {
          out[[i]] <- data.frame(
            time = inp$timeI[[i]],
            obs  = inp$obsI[[i]] / qest[inp$mapq[i], 2],
            idx  = i
          )
        }
      }
    }
    out
  }

  ## ───────────────────────── Biomass panel (Bmsy band + B/Bmsy-on-B ribbon) ─────────────────
  make_biomass_plot <- function() {
    Best <- get.par("logB",     model, exp = TRUE, CI = CI)
    BB   <- get.par("logBBmsy", model, exp = TRUE, CI = CI)

    repmax  <- if (manflag) get.manmax(model) else model
    Bmsy_all <- get.par("logBmsy", repmax, exp = TRUE, CI = CI)  # keep matrix form first
    Bmsyvec  <- get.msyvec(repmax$inp, Bmsy_all)
    Bmsy     <- if (!is.null(nrow(Bmsy_all))) Bmsy_all[1, ] else Bmsy_all

    df_B <- data.frame(time = as.numeric(rownames(Best)),
                       lwr  = Best[, 1], est = Best[, 2], upr = Best[, 3])
    df_B_in <- if (length(ind_in)) df_B[ind_in, , drop = FALSE] else df_B[0, ]
    df_B_pr <- if (length(ind_pr)) df_B[ind_pr, , drop = FALSE] else df_B[0, ]

    Bmsy_est <- Bmsy[2]
    df_BB_rib <- make_ribbon_df(
      time = as.numeric(rownames(BB)),
      lwr  = BB[, 1] * Bmsy_est,
      upr  = BB[, 3] * Bmsy_est
    )

    df_Bmsy_band <- data.frame(time = repmax$inp$time, ymin = Bmsyvec$ll, ymax = Bmsyvec$ul)
    df_Bmsy_line <- data.frame(time = repmax$inp$time, y = Bmsyvec$msy)

    p <- ggplot()
    if (isTRUE(show_CIs)) {
      p <- p + geom_ribbon(data = df_Bmsy_band, aes(x = time, ymin = ymin, ymax = ymax), fill = "grey80")
      p <- p + geom_ribbon(data = df_BB_rib,    aes(x = time, ymin = ymin, ymax = ymax), fill = rgb(0, 0, 1, 0.10)) +
        geom_line(data = transform(df_BB_rib, y = ymin), aes(x = time, y = y), color = rgb(0, 0, 1, 0.20), linewidth = 0.6) +
        geom_line(data = transform(df_BB_rib, y = ymax), aes(x = time, y = y), color = rgb(0, 0, 1, 0.20), linewidth = 0.6)
    }
    p <- p + geom_line(data = df_Bmsy_line, aes(x = time, y = y), color = "black", linewidth = 0.7)

    if (isTRUE(show_CIs) && nrow(df_B_in)) {
      p <- p + geom_line(data = df_B_in, aes(x = time, y = lwr), linetype = "dotted", linewidth = 0.6, color = line_color) +
        geom_line(data = df_B_in, aes(x = time, y = upr), linetype = "dotted", linewidth = 0.6, color = line_color)
    }
    if (nrow(df_B_in)) p <- p + geom_line(data = df_B_in, aes(x = time, y = est), linewidth = 0.8, color = line_color)
    if (!manflag && nrow(df_B_pr)) {
      p <- p + geom_line(data = df_B_pr, aes(x = time, y = est), linetype = "dashed", linewidth = 0.8, color = line_color)
      if (isTRUE(show_CIs)) {
        p <- p + geom_line(data = df_B_pr, aes(x = time, y = lwr), linetype = "dotted", linewidth = 0.6, color = line_color) +
          geom_line(data = df_B_pr, aes(x = time, y = upr), linetype = "dotted", linewidth = 0.6, color = line_color)
      }
    }
    if (is.finite(obs_end_overall)) p <- p + geom_vline(xintercept = obs_end_overall, color = vline_col, linewidth = vline_size)

    p <- p + labs(x = "Year", y = "Biomass") + theme_minimal_compact2()
    idx <- get_index_points()
    if (length(idx) >= 1) p <- p + geom_point(data = idx[[1]], aes(x = time, y = obs), color = "blue", shape = 16, size = 2, inherit.aes = FALSE)
    if (length(idx) >= 2) p <- p + geom_point(data = idx[[2]], aes(x = time, y = obs), shape = 22, color = "black", fill = "green", size = 2, stroke = 0.5, inherit.aes = FALSE)

    p <- p + scale_y_continuous(sec.axis = sec_axis(~ . / Bmsy_est, name = expression(B[t]/B[MSY])))
    p
  }

  ## ───────────────────────────── B/Bmsy panel (ribbon + outlines) ───────────────────────────
  make_bbmsy_plot <- function() {
    BB <- get.par("logBBmsy", model, exp = TRUE, CI = CI)
    df <- data.frame(time = as.numeric(rownames(BB)), lwr = BB[,1], est = BB[,2], upr = BB[,3])
    df_in <- if (length(ind_in)) df[ind_in, , drop = FALSE] else df[0, ]
    df_pr <- if (!manflag && length(ind_pr)) df[ind_pr, , drop = FALSE] else df[0, ]
    rib <- make_ribbon_df(df$time, df$lwr, df$upr)

    p <- ggplot()
    if (isTRUE(show_CIs)) {
      p <- p + geom_ribbon(data = rib, aes(x = time, ymin = ymin, ymax = ymax), fill = rgb(0,0,1,0.10)) +
        geom_line(data = transform(rib, y = ymin), aes(x = time, y = y), color = rgb(0,0,1,0.20), linewidth = 0.6) +
        geom_line(data = transform(rib, y = ymax), aes(x = time, y = y), color = rgb(0,0,1,0.20), linewidth = 0.6)
    }
    if (nrow(df_in)) p <- p + geom_line(data = df_in, aes(x = time, y = est), color = line_color, linewidth = 0.8)
    if (nrow(df_pr)) p <- p + geom_line(data = df_pr, aes(x = time, y = est), linetype = "dashed", color = line_color, linewidth = 0.8)
    if (is.finite(obs_end_overall)) p <- p + geom_vline(xintercept = obs_end_overall, color = vline_col, linewidth = vline_size)
    p <- p + geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
      labs(x = "Year", y = expression(bold(B/B[MSY]))) + theme_minimal_compact2()
    p
  }

  ## ─────────────────────────── F/Fmsy panel (ribbon + outlines) ─────────────────────────────
  make_ffmsy_plot <- function() {
    FF <- get.par("logFFmsy", model, exp = TRUE, CI = CI)
    df <- data.frame(time = as.numeric(rownames(FF)), lwr = FF[,1], est = FF[,2], upr = FF[,3])
    df_in <- if (length(ind_in)) df[ind_in, , drop = FALSE] else df[0, ]
    df_pr <- if (!manflag && length(ind_pr)) df[ind_pr, , drop = FALSE] else df[0, ]
    rib <- make_ribbon_df(df$time, df$lwr, df$upr)

    p <- ggplot()
    if (isTRUE(show_CIs)) {
      p <- p + geom_ribbon(data = rib, aes(x = time, ymin = ymin, ymax = ymax), fill = rgb(0,0,1,0.10)) +
        geom_line(data = transform(rib, y = ymin), aes(x = time, y = y), color = rgb(0,0,1,0.20), linewidth = 0.6) +
        geom_line(data = transform(rib, y = ymax), aes(x = time, y = y), color = rgb(0,0,1,0.20), linewidth = 0.6)
    }
    if (nrow(df_in)) p <- p + geom_line(data = df_in, aes(x = time, y = est), color = line_color, linewidth = 0.8)
    if (nrow(df_pr)) p <- p + geom_line(data = df_pr, aes(x = time, y = est), linetype = "dashed", color = line_color, linewidth = 0.8)
    if (is.finite(obs_end_overall)) p <- p + geom_vline(xintercept = obs_end_overall, color = vline_col, linewidth = vline_size)
    p <- p + geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
      labs(x = "Year", y = expression(bold(F/F[MSY]))) + theme_minimal_compact2()

    if (min(inp$dtc, na.rm = TRUE) < 1 && nrow(df_in)) {
      FFs <- get.par("logFFmsy", model, exp = TRUE, CI = CI)
      p <- p + geom_line(data = data.frame(time = as.numeric(rownames(FFs))[ind_in], est = FFs[ind_in, 2]),
                         aes(x = time, y = est), color = rgb(0,0,1,0.4), linewidth = 0.6)
    }
    p
  }

  ## ─────────────────────────── Absolute F panel (Ft + CI + bands) ───────────────────────────
  ## --- replace your current make_f_absolute_plot() with this one ---
  make_f_absolute_plot <- function() {
    repmax  <- if (manflag) get.manmax(model) else model
    tvgflag <- isTRUE(repmax$inp$timevaryinggrowth) || isTRUE(repmax$inp$logmcovflag)

    ## Absolute Ft median + CI (direct)
    Fest <- get.par("logFnotS", model, exp = TRUE, CI = CI)
    pmF  <- as.matrix(Fest)

    # time for Ft: do NOT rely on rownames
    t_full <- as.numeric(model$inp$time)
    if (length(t_full) >= nrow(pmF)) {
      timeF <- t_full[seq_len(nrow(pmF))]
    } else {
      timeF <- rep_len(t_full, nrow(pmF))
    }

    dfF <- data.frame(
      time = timeF,
      lwr  = pmF[, 1],
      est  = pmF[, 2],
      upr  = pmF[, 3]
    )
    dfF_in <- if (length(ind_in)) dfF[ind_in, , drop = FALSE] else dfF[0, ]
    dfF_pr <- if (!manflag && length(ind_pr)) dfF[ind_pr, , drop = FALSE] else dfF[0, ]

    ## Fmsy band (+ line)
    if (tvgflag) {
      Fmsy_all <- get.par("logFmsyvec", repmax, exp = TRUE, CI = CI)
      Fmsyvec  <- as.data.frame(Fmsy_all); Fmsyvec$msy <- Fmsyvec$est
      Fmsy_band <- data.frame(time = repmax$inp$time, ymin = Fmsyvec$ll, ymax = Fmsyvec$ul)
      Fmsy_line <- data.frame(time = repmax$inp$time, y = Fmsyvec$msy)
      fmsy_mult <- matrix(rep(Fmsyvec$msy, each = 3), ncol = 3, byrow = TRUE)
      sec_axis_obj <- waiver()
    } else {
      Fmsy_all <- get.par("logFmsy", repmax, exp = TRUE, CI = CI)
      if (any(is.na(Fmsy_all))) Fmsy_all <- get.par("logFmsyd", repmax, exp = TRUE, CI = CI)

      # make it a 1×5 matrix and take “est” robustly
      xm <- as.matrix(Fmsy_all)
      col_est <- if (!is.null(colnames(xm))) {
        which(colnames(xm) %in% c("est","Est","EST"))[1]
      } else NA_integer_
      if (is.na(col_est)) col_est <- 2
      Fmsy_est <- as.numeric(xm[1, col_est])

      Fmsyvec  <- get.msyvec(repmax$inp, Fmsy_all)
      Fmsy_band <- data.frame(time = repmax$inp$time, ymin = Fmsyvec$ll, ymax = Fmsyvec$ul)
      Fmsy_line <- data.frame(time = repmax$inp$time, y = Fmsyvec$msy)
      fmsy_mult <- Fmsy_est
      sec_axis_obj <- sec_axis(~ . / Fmsy_est, name = expression(F[t]/F[MSY]))
    }

    ## Absolute-F ribbon from (FFmsy_notS × Fmsy)
    FFrel <- get.par("logFFmsynotS", model, exp = TRUE, CI = CI)[, 1:3]
    # align time for ribbon up to last observed
    last_obs_time <- if (manflag) repmax$inp$timerange[2] else max(repmax$inp$time, na.rm = TRUE)
    ix_last <- which(repmax$inp$time == last_obs_time)
    if (!length(ix_last)) ix_last <- length(repmax$inp$time)
    indxmax <- max(ix_last - if (manflag) 1 else 0, 1)

    timef <- repmax$inp$time[seq_len(indxmax)]
    if (is.matrix(fmsy_mult)) {
      FFabs <- FFrel * fmsy_mult
    } else {
      FFabs <- FFrel * as.numeric(fmsy_mult)
    }
    clf <- FFabs[seq_len(indxmax), 1]
    cuf <- FFabs[seq_len(indxmax), 3]
    ribFF <- {
      ok <- is.finite(timef) & is.finite(clf) & is.finite(cuf)
      data.frame(time = timef[ok], ymin = clf[ok], ymax = cuf[ok])
    }

    p <- ggplot()

    # Fmsy grey band + mid line
    if (isTRUE(show_CIs)) {
      p <- p + geom_ribbon(data = Fmsy_band, aes(x = time, ymin = ymin, ymax = ymax), fill = "grey80")
    }
    p <- p + geom_line(data = Fmsy_line, aes(x = time, y = y), color = "black", linewidth = 0.7)

    # Absolute-F (FF×Fmsy) ribbon + faint edges
    if (isTRUE(show_CIs) && nrow(ribFF)) {
      p <- p + geom_ribbon(data = ribFF, aes(x = time, ymin = ymin, ymax = ymax), fill = rgb(0,0,1,0.10)) +
        geom_line(data = transform(ribFF, y = ymin), aes(x = time, y = y), color = rgb(0,0,1,0.20), linewidth = 0.6) +
        geom_line(data = transform(ribFF, y = ymax), aes(x = time, y = y), color = rgb(0,0,1,0.20), linewidth = 0.6)
    }

    # Sub-annual Ft
    if (min(inp$dtc, na.rm = TRUE) < 1) {
      sF <- get.par("logFs", model, exp = TRUE, CI = CI)
      p <- p + geom_line(
        data = data.frame(time = repmax$inp$time[seq_len(indxmax)], est = sF[seq_len(indxmax), 2]),
        aes(x = time, y = est),
        color = rgb(0,0,1,0.4), linewidth = 0.6
      )
    }

    # In-sample Ft median + dotted CI
    if (isTRUE(show_CIs) && nrow(dfF_in)) {
      p <- p + geom_line(data = dfF_in, aes(x = time, y = lwr), linetype = "dotted", linewidth = 0.6, color = line_color) +
        geom_line(data = dfF_in, aes(x = time, y = upr), linetype = "dotted", linewidth = 0.6, color = line_color)
    }
    if (nrow(dfF_in)) {
      p <- p + geom_line(data = dfF_in, aes(x = time, y = est), color = line_color, linewidth = 0.8)
    }

    # Prediction (no scenarios)
    if (!manflag && nrow(dfF_pr)) {
      if (isTRUE(show_CIs)) {
        p <- p + geom_line(data = dfF_pr, aes(x = time, y = lwr), linetype = "dotted", linewidth = 0.6, color = line_color) +
          geom_line(data = dfF_pr, aes(x = time, y = upr), linetype = "dotted", linewidth = 0.6, color = line_color)
      }
      p <- p + geom_line(data = dfF_pr, aes(x = time, y = est), linetype = "dashed", color = line_color, linewidth = 0.8)
    }

    # Observed F-proxy
    if (!is.null(inp$timeE) && length(inp$timeE) && !is.null(inp$obsE) && !is.null(inp$dte)) {
      qf <- try(get.par("logqf", model, exp = TRUE, CI = CI), silent = TRUE)
      if (!inherits(qf, "try-error") && is.finite(qf[2])) {
        df_obsF <- data.frame(time = inp$timeE, val = inp$obsE / inp$dte * qf[2])
        p <- p + geom_point(data = df_obsF, aes(x = time, y = val), color = "black", size = 1.8)
      }
    }

    if (is.finite(obs_end_overall)) p <- p + geom_vline(xintercept = obs_end_overall, color = vline_col, linewidth = vline_size)
    p <- p + labs(x = "Year", y = "Fishing mortality") + theme_minimal_compact2()

    # Right-hand axis for F/Fmsy only if Fmsy is constant
    p <- p + scale_y_continuous(sec.axis = sec_axis_obj)
    p
  }

  ## ─────────────────────────── Build panels ───────────────────────────
  plots <- list(
    biomass = make_biomass_plot(),
    bbmsy   = make_bbmsy_plot(),
    catch   = NULL,  # set below after MSY band prep
    f       = make_f_absolute_plot(),
    ffmsy   = make_ffmsy_plot()
  )

  ## ─────────────────────────── Catch panel with MSY band ───────────────────────
  if (is.null(extract_catch_data)) {
    CP <- get.par("logCpred", model, exp = TRUE, CI = CI)
    predicted <- data.frame(time = model$inp$timeCpred, lwr = CP[,1], catch = CP[,2], upr = CP[,3])
    observed  <- data.frame(time = model$inp$timeC,    catch = model$inp$obsC)
  } else {
    cd <- extract_catch_data(model, scenario_name = "Model")
    predicted <- subset(cd, catch_type == "Predicted")
    observed  <- subset(cd, catch_type == "Observed")
  }

  repmax <- if (manflag) get.manmax(model) else model
  tvgflag <- repmax$inp$timevaryinggrowth || repmax$inp$logmcovflag
  if (tvgflag) {
    MSY <- get.par("logMSYvec", repmax, exp = TRUE, CI = CI)
    MSYvec <- transform(as.data.frame(MSY), msy = est)
    MSY_band <- data.frame(time = repmax$inp$time, ymin = MSYvec$ll, ymax = MSYvec$ul)
    MSY_line <- data.frame(time = repmax$inp$time, y = MSYvec$msy)
  } else {
    MSY <- get.par("logMSY", repmax, exp = TRUE, CI = CI)
    MSYvec <- get.msyvec(repmax$inp, MSY)
    MSY_band <- data.frame(time = repmax$inp$time, ymin = MSYvec$ll, ymax = MSYvec$ul)
    MSY_line <- data.frame(time = repmax$inp$time, y = MSYvec$msy)
  }

  plots$catch <- ggplot() +
    { if (isTRUE(show_CIs)) geom_ribbon(data = MSY_band, aes(x = time, ymin = ymin, ymax = ymax), fill = "grey80") } +
    { if (isTRUE(show_CIs) && nrow(predicted)) geom_line(data = predicted, aes(x = time, y = lwr), linetype = "dotted", linewidth = 0.6, color = line_color) } +
    { if (isTRUE(show_CIs) && nrow(predicted)) geom_line(data = predicted, aes(x = time, y = upr), linetype = "dotted", linewidth = 0.6, color = line_color) } +
    { if (nrow(predicted)) geom_line(data = predicted, aes(x = time, y = catch), linewidth = 0.8, color = line_color) } +
    { if (nrow(observed))  geom_point(data = observed,  aes(x = time, y = catch), color = "black", size = 1.3) } +
    { if (is.finite(catch_obs_end)) geom_vline(xintercept = catch_obs_end, color = vline_col, linewidth = vline_size) } +
    geom_line(data = MSY_line, aes(x = time, y = y)) +
    labs(x = "Year", y = "Catch") +
    theme_minimal_compact2()

  ## ─────────────────────────── Kobe panel ───────────────────────────
  rep_kobe <- model
  if ("man" %in% names(rep_kobe)) rep_kobe$man <- NULL
  plots$kobe <- Original_kobe_all_in_gg_internal(
    rep         = rep_kobe,
    logax       = FALSE,
    plot.legend = TRUE,
    man.legend  = FALSE,
    ext         = TRUE,
    rel.axes    = FALSE,
    CI          = CI,
    print_it    = FALSE
  )

  ## ─────────────────────────── Layout / return ───────────────────────────
  if (isTRUE(return_patchwork)) {
    (plots$biomass | plots$bbmsy | plots$catch) /
      patchwork::plot_spacer() /
      (plots$f | plots$ffmsy | plots$kobe) +
      patchwork::plot_layout(heights = c(1, 0.05, 1)) &
      theme(plot.margin = margin(4, 4, 4, 4))
  } else {vtoo
    plots
  }
}
