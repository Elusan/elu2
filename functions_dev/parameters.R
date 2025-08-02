#' Calculate gamma from the Pella-Tomlinson exponent
#'
#' This function computes the scaling factor gamma used in the Pella-Tomlinson surplus production function.
#'
#' @param n Numeric; exponent of the Pella-Tomlinson production function.
#'
#' @return A numeric value of gamma.
#'
#' @examples
#' calc.gamma(2)
#'
#' @export
calc.gamma <- function(n) {
  n^(n / (n - 1)) / (n - 1)
}

#' Extract parameter estimates and uncertainty from a SPiCT model result
#'
#' Extracts the estimate, standard deviation, confidence intervals, and coefficient of variation
#' for a specified parameter from a fitted SPiCT result object.
#'
#' @param parname Character; name of the parameter to extract (e.g. `"logBmsy"`).
#' @param rep_obj A fitted SPiCT model object (i.e., the output of `fit.spict()`).
#' @param exp Logical; whether to return values on the exponentiated scale (default: `FALSE`).
#' @param CI Numeric; confidence level for intervals (default: `0.95`).
#' @param random Deprecated. Ignored.
#' @param fixed Deprecated. Ignored.
#'
#' @return A matrix with columns: lower bound, estimate, upper bound, standard deviation, and coefficient of variation.
#' If `parname` is time-varying (e.g., `"logB"`, `"logF"`), row names will be the time vector.
#'
#' @examples
#' \dontrun{
#' get.par("logBmsy", rep_obj, exp = TRUE)
#' }
#'
#' @export
get.par <- function(parname, rep_obj, exp = FALSE, random = FALSE, fixed = FALSE, CI = 0.95) {
  if (CI > 1 || CI < 0) stop("CI has to be between 0 and 1!")
  zscore <- qnorm(CI + (1 - CI) / 2)

  if (!"sderr" %in% names(rep_obj)) {
    est <- sd <- ll <- ul <- NA
    indran  <- which(names(rep_obj$par.random) == parname)
    indfix  <- which(names(rep_obj$par.fixed) == parname)
    indsdr  <- which(names(rep_obj$value) == parname)
    indopt  <- which(names(rep_obj$opt$par) == parname)

    if (length(indran) > 0) {
      est <- rep_obj$par.random[indran]
      sd  <- sqrt(rep_obj$diag.cov.random[indran])
      ll <- est - zscore * sd
      ul <- est + zscore * sd
    }
    if (length(indfix) > 0) {
      est <- rep_obj$par.fixed[indfix]
      sd  <- sqrt(diag(rep_obj$cov.fixed))[indfix]
      ll <- est - zscore * sd
      ul <- est + zscore * sd
    }
    if (length(indsdr) > 0) {
      est <- rep_obj$value[indsdr]
      sd  <- rep_obj$sd[indsdr]
      ll <- est - zscore * sd
      ul <- est + zscore * sd
    }

    if (length(est) == 0) {
      if (length(indopt) > 0) {
        est <- rep_obj$opt$par[indopt]
      } else if ("phases" %in% names(rep_obj$inp) &&
                 parname %in% names(rep_obj$inp$phases) &&
                 rep_obj$inp$phases[[parname]] == -1) {
        est <- rep_obj$inp$parlist[[parname]]
        ll <- ul <- est
      } else if (!is.na(parname) && parname == "P") {
        B <- get.par("logB", rep_obj, exp = TRUE, CI = CI)
        C <- get.par("logCpred", rep_obj, exp = TRUE, CI = CI)
        ic <- rep_obj$inp$ic
        nc <- rep_obj$inp$nc
        B0 <- B[ic, 2]
        B1 <- B[ic + nc, 2]
        T0 <- rep_obj$inp$time[ic]
        T1 <- rep_obj$inp$time[ic + nc]
        est <- (B1 - B0 + C[, 2]) / (T1 - T0)
      } else {
        warning("get.par WARNING: could not extract ", parname)
      }
    }

    if (exp) {
      cv <- sqrt(exp(sd^2) - 1)
      ll <- exp(ll)
      ul <- exp(ul)
      ul[ul == Inf] <- exp(705)  # Replace Inf with a large number
      est <- exp(est)
    } else {
      cv <- sd / est
    }

    out <- cbind(ll, est, ul, sd, cv)
    if (parname %in% c("logB", "logF", "logBBmsy", "logFFmsy")) {
      rownames(out) <- rep_obj$inp$time
    }
    return(out)
  }
}
