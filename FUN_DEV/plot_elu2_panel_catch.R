#' Catch panel (Observed vs Predicted) with MSY band
#'
#' @param model A fitted SPiCT object (`spictcls`).
#' @param extract_catch_data Optional function(model, scenario_name) -> data.frame(time, catch, lwr, upr, catch_type).
#' @param line_color Line color for predicted series.
#' @param show_CIs Logical; draw ribbons & dotted CI edges.
#' @param CI Confidence level for intervals (default 0.95).
#' @return A ggplot object.
#' @export
#' @import ggplot2
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

  p <- ggplot()
  if (isTRUE(show_CIs)) {
    p <- p + geom_ribbon(data = MSY_band, aes(x = time, ymin = ymin, ymax = ymax), fill = "grey80")
  }
  if (isTRUE(show_CIs) && nrow(predicted)) {
    p <- p +
      geom_line(data = predicted, aes(x = time, y = lwr), linetype = "dotted", linewidth = 0.6, color = line_color) +
      geom_line(data = predicted, aes(x = time, y = upr), linetype = "dotted", linewidth = 0.6, color = line_color)
  }
  if (nrow(predicted)) p <- p + geom_line(data = predicted, aes(x = time, y = catch), linewidth = 0.8, color = line_color)
  if (nrow(observed))  p <- p + geom_point(data = observed,  aes(x = time, y = catch), color = "black", size = 1.3)

  catch_obs_end <- .spict_catch_obs_end(model)
  if (is.finite(catch_obs_end)) p <- p + geom_vline(xintercept = catch_obs_end, color = vline_col, linewidth = vline_size)

  p <- p + geom_line(data = MSY_line, aes(x = time, y = y), linetype = "solid") +
    labs(title = "Catch", x = "Year", y = "Catch") +
    .spict_theme_minimal_compact2()

  p
}

# ---- Inlined helpers (only those used by this panel) ----

#' Minimal compact ggplot2 theme (ELU2)
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

#' End of observation time for catches (handles sub-annual dtc)
#' @keywords internal
#' @noRd
.spict_catch_obs_end <- function(model) {
  inp <- model$inp
  if (is.null(inp$timeC) || !length(inp$timeC)) return(NA_real_)
  if (!is.null(inp$dtc) && length(inp$dtc) && min(inp$dtc, na.rm = TRUE) < 1) {
    ix <- which((inp$timeC %% 1) == 0)
    if (length(ix)) return(tail(inp$timeC[ix], 1))
    return(tail(inp$timeC, 1))
  } else {
    return(tail(inp$timeC, 1))
  }
}
