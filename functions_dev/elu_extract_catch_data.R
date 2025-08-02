#' Extract Observed and Predicted Catch Time Series from SPiCT Model Fit
#'
#' Returns a tidy data frame containing both observed and model-predicted catches (and their confidence intervals) for a given scenario from a SPiCT model fit. Handles both annual and seasonal data.
#'
#' @param fit A fitted SPiCT model object (typically the output of \code{fit.spict()} or a compatible function). Must include input list \code{inp} and prediction results.
#' @param scenario_name Character. Name of the scenario/model for labeling in the output (default: "Base").
#' @param CI Numeric. Width of confidence interval to extract for predictions (default: 0.95).
#'
#' @details
#' This function extracts both observed and predicted catch time series from a SPiCT model fit. For seasonal data (where \code{dtc < 1}), the function annualizes both observed and predicted catches and their time points using the \code{annual()} function, then calculates mean or median values across each year. For standard annual models, catches are used as-is.
#'
#' Output columns include time, catch (tons), lower and upper CI (for predicted only), type ("Observed" or "Predicted"), and scenario label.
#'
#' @return
#' A data frame with columns:
#'   \describe{
#'     \item{time}{Time of observation or prediction (year or fractional year).}
#'     \item{catch}{Catch value (tons, annualized if appropriate).}
#'     \item{lwr}{Lower confidence interval (only for predicted catches; NA for observed).}
#'     \item{upr}{Upper confidence interval (only for predicted catches; NA for observed).}
#'     \item{catch_type}{Type: "Observed" or "Predicted".}
#'     \item{scenario}{Name of the scenario/model for labeling.}
#'   }
#'
#' @export
elu_extract_catch_data <- function(fit, scenario_name = "Base", CI = 0.95) {
  # Defensive checks
  if (!("inp" %in% names(fit))) stop("Input must be a SPiCT fit object.")

  inp <- fit$inp
  # Predicted catch: rows = times, cols = [lwr, est, upr]
  Cp <- get.par("logCpred", fit, exp = TRUE, CI = CI)
  Cp[Cp < 0] <- 0

  # Handle annual/seasonal data if dtc < 1
  if (min(inp$dtc) < 1) {
    # Use annualized time/values for observed catch
    alo <- annual(inp$timeC, inp$obsC / inp$dtc)
    time_obs <- alo$anntime
    obsC <- alo$annvec

    # Annualized predictions (use median/est column for prediction, lwr/upr for CI)
    al1 <- annual(inp$timeCpred, Cp[, 1] / inp$dtcp)
    al2 <- annual(inp$timeCpred, Cp[, 2] / inp$dtcp)
    al3 <- annual(inp$timeCpred, Cp[, 3] / inp$dtcp)
    inds <- which(!is.na(al2$annvec))
    time_pred <- al2$anntime[inds]
    catch_pred <- al2$annvec[inds]
    lwr <- al1$annvec[inds]
    upr <- al3$annvec[inds]
  } else {
    # Standard case (annual)
    time_obs <- inp$timeC
    obsC <- inp$obsC / inp$dtc
    time_pred <- inp$timeCpred
    catch_pred <- Cp[, 2] / inp$dtcp
    lwr <- Cp[, 1] / inp$dtcp
    upr <- Cp[, 3] / inp$dtcp
  }

  # Observed
  obs_df <- data.frame(
    time = time_obs,
    catch = obsC,
    lwr = NA,
    upr = NA,
    catch_type = "Observed",
    scenario = scenario_name
  )

  # Predicted
  pred_df <- data.frame(
    time = time_pred,
    catch = catch_pred,
    lwr = lwr,
    upr = upr,
    catch_type = "Predicted",
    scenario = scenario_name
  )

  # Combine and sort
  out <- rbind(obs_df, pred_df)
  out <- out[order(out$time, out$catch_type), ]
  rownames(out) <- NULL
  return(out)
}
