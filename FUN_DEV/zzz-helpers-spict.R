# ===========================
# Shared helpers (internal)
# ===========================

#' Compact panel theme for ELU2 plots (internal)
#'
#' Minimal, bold, boxed theme used across ELU2 panels. Relies on
#' \pkg{ggplot2} geoms/labels and \pkg{grid} units. Intended for internal use.
#'
#' @param base_size Base font size.
#' @param base_family Base font family.
#'
#' @return A \code{ggplot2} theme object.
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

  }  # (body unchanged)

#' Last observation time across available series (internal)
#'
#' Returns the most reasonable "end of observations" time using, in order:
#' \code{inp$timerangeObs[2]}, last \code{timeC}, last first-index \code{timeI[[1]]},
#' or the maximum of \code{inp$time}. Falls back to \code{NA_real_}.
#'
#' @param model A fitted SPiCT result (\code{spictcls}).
#' @return Numeric scalar (time) or \code{NA_real_}.
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

  }  # (body unchanged)

#' Last observation time for Catch (internal)
#'
#' Returns the final observed catch time, respecting seasonal models (\code{dtc < 1})
#' by preferring the last whole-year time if available.
#'
#' @param model A fitted SPiCT result (\code{spictcls}).
#' @return Numeric scalar (time) or \code{NA_real_}.
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

  }  # (body unchanged)

#' Build a ribbon dataframe from time and CI bounds (internal)
#'
#' Filters to finite rows and returns a data frame with
#' \code{time}, \code{ymin}, \code{ymax}.
#'
#' @param time Numeric vector of times.
#' @param lwr Numeric vector of lower bounds.
#' @param upr Numeric vector of upper bounds.
#' @return Data frame with columns \code{time}, \code{ymin}, \code{ymax}.
#' @keywords internal
#' @noRd
.spict_make_ribbon_df <- function(time, lwr, upr) {
  ok <- is.finite(time) & is.finite(lwr) & is.finite(upr)
  data.frame(time = time[ok], ymin = lwr[ok], ymax = upr[ok])
  }  # (body unchanged)

#' Extract q-scaled index points for overlay (internal)
#'
#' Builds a list of data frames (one per index series) with \code{time}, \code{obs},
#' and \code{idx}, dividing observed indices by \eqn{\hat q} for comparability.
#' Respects \code{inp$mapq} if present.
#'
#' @param model A fitted SPiCT result (\code{spictcls}).
#' @return List of data frames (possibly empty).
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
  }  # (body unchanged)

#' Derive a "time" vector matching a parameter matrix (internal)
#'
#' Uses numeric row names if present; otherwise falls back to \code{inp$time}
#' (trimming or recycling as needed).
#'
#' @param model A fitted SPiCT result (\code{spictcls}).
#' @param par_matrix Parameter matrix returned by \code{get.par()}.
#' @return Numeric vector of times aligned to \code{par_matrix} rows.
#' @keywords internal
#' @noRd
.spict_time_from_par <- function(model, par_matrix) {
  # Prefer clean numeric rownames; fallback to inp$time (trim/recycle).
  rn <- suppressWarnings(as.numeric(rownames(par_matrix)))
  if (length(rn) == nrow(par_matrix) && !all(is.na(rn))) return(rn)
  t_full <- as.numeric(model$inp$time)
  if (length(t_full) >= nrow(par_matrix)) return(t_full[seq_len(nrow(par_matrix))])
  return(rep_len(t_full, nrow(par_matrix)))
}  # (body unchanged)
