# ===========================
# 3/6 — Catch panel
# ===========================

#' Catch panel (Observed vs Predicted) with MSY band
#'
#' Plots predicted catches (median with dotted CI bounds) as a line,
#' observed catches as points, and overlays the \eqn{MSY} band/line. A vertical
#' line marks the end of observed catches. Optionally accepts a custom extractor.
#'
#' @param model A fitted SPiCT object (\code{spictcls}).
#' @param extract_catch_data Optional function \code{function(model, scenario_name)}
#'   returning a data frame with columns \code{time}, \code{catch}, \code{lwr},
#'   \code{upr}, and \code{catch_type} in \{ "Observed","Predicted" \}.
#' @param line_color Colour for the predicted series.
#' @param show_CIs Logical; draw dotted CI edges for predictions and MSY band.
#' @param CI Confidence level for intervals (default \code{0.95}).
#'
#' @return A \code{ggplot} object.
#' @family ELU2-panels
#' @export
#' @examples
#' \dontrun{
#' m1 <- all_models$S1$S1P.SDM
#' plot_elu2_panel_catch(m1)
#' }
plot_elu2_panel_catch <- function(model,
                                  extract_catch_data = NULL,
                                  line_color = "blue",
                                  show_CIs   = TRUE,
                                  CI         = 0.95) {

  stopifnot(inherits(model, "spictcls"))
  vline_col  <- "grey50"; vline_size <- 0.2
  manflag <- ("man" %in% names(model))
  inp <- model$inp

  if (is.null(extract_catch_data)) {
    CP <- get.par("logCpred", model, exp = TRUE, CI = CI)
    predicted <- data.frame(time = inp$timeCpred, lwr = CP[,1], catch = CP[,2], upr = CP[,3])
    observed  <- data.frame(time = inp$timeC,    catch = inp$obsC)
  } else {
    cd <- extract_catch_data(model, scenario_name = "Model")
    predicted <- subset(cd, catch_type == "Predicted")
    observed  <- subset(cd, catch_type == "Observed")
  }

  repmax <- if (manflag) get.manmax(model) else model
  tvgflag <- isTRUE(repmax$inp$timevaryinggrowth) || isTRUE(repmax$inp$logmcovflag)
  if (tvgflag) {
    MSY <- get.par("logMSYvec", repmax, exp = TRUE, CI = CI)
    MSYvec   <- transform(as.data.frame(MSY), msy = est)
    MSY_band <- data.frame(time = repmax$inp$time, ymin = MSYvec$ll, ymax = MSYvec$ul)
    MSY_line <- data.frame(time = repmax$inp$time, y = MSYvec$msy)
  } else {
    MSY <- get.par("logMSY", repmax, exp = TRUE, CI = CI)
    MSYvec   <- get.msyvec(repmax$inp, MSY)
    MSY_band <- data.frame(time = repmax$inp$time, ymin = MSYvec$ll, ymax = MSYvec$ul)
    MSY_line <- data.frame(time = repmax$inp$time, y = MSYvec$msy)
  }

  p <- ggplot2::ggplot()
  if (isTRUE(show_CIs)) {
    p <- p + ggplot2::geom_ribbon(data = MSY_band, ggplot2::aes(x = time, ymin = ymin, ymax = ymax), fill = "grey80")
  }
  if (isTRUE(show_CIs) && nrow(predicted)) {
    p <- p +
      ggplot2::geom_line(data = predicted, ggplot2::aes(x = time, y = lwr), linetype = "dotted", linewidth = 0.6, color = line_color) +
      ggplot2::geom_line(data = predicted, ggplot2::aes(x = time, y = upr), linetype = "dotted", linewidth = 0.6, color = line_color)
  }
  if (nrow(predicted)) p <- p + ggplot2::geom_line(data = predicted, ggplot2::aes(x = time, y = catch), linewidth = 0.8, color = line_color)
  if (nrow(observed))  p <- p + ggplot2::geom_point(data = observed,  ggplot2::aes(x = time, y = catch), color = "black", size = 1.3)

  catch_obs_end <- .spict_catch_obs_end(model)
  if (is.finite(catch_obs_end)) p <- p + ggplot2::geom_vline(xintercept = catch_obs_end, color = vline_col, linewidth = vline_size)

  p <- p + ggplot2::geom_line(data = MSY_line, ggplot2::aes(x = time, y = y), linetype = "solid") +
    ggplot2::labs(title = "Catch", x = "Year", y = "Catch") +
    .spict_theme_minimal_compact2()

  p

}  # (body unchanged)
