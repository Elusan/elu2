#' ELU2 3×3 summary grid (ggplot2)
#'
#' @description
#' Builds a compact 3×3 grid of ELU2/SPiCT panels (ggplot2) and prints it.
#' The grid contains: Biomass (Bt), Fishing mortality (Ft), Catch, B/Bmsy,
#' F/Fmsy, Kobe phase plot, Production curve, Time to Bmsy, and a single
#' prior–posterior density panel (first active prior).
#'
#' When simulation truth is available in \code{rep$inp$true}, the function
#' overlays the appropriate “true” trajectories/levels on relevant panels via
#' an internal helper \code{.elu2_overlay_truth()} (not exported).
#'
#' @details
#' \strong{Panels included}
#' \itemize{
#'   \item \code{plot_elu2_panel_biomass()}
#'   \item \code{plot_elu2_panel_f()}
#'   \item \code{plot_elu2_panel_catch()}
#'   \item \code{plot_elu2_panel_bbmsy()}
#'   \item \code{plot_elu2_panel_ffmsy()}
#'   \item \code{plot_elu2_panel_kobe()}
#'   \item \code{plot_elu2_panel_production()}
#'   \item \code{plot_elu2_panel_time_to_bmsy()}
#'   \item \code{plot_elu2_panel_per_model(do.plot = 1)} (falls back to a blank placeholder if unavailable)
#' }
#'
#' \strong{Dependencies}
#' Uses \pkg{patchwork} for the layout (checked at runtime with
#' \code{requireNamespace("patchwork", quietly = TRUE)}). You should list
#' \pkg{patchwork} in \strong{Suggests} of your DESCRIPTION, or import it if you
#' prefer a hard dependency.
#'
#' @param x A fitted SPiCT object (class \code{spictcls}); typically returned by \code{fit.spict()}.
#' @param stamp Character stamp for the figure (kept for API compatibility; not used).
#' @param verbose Logical; if \code{TRUE}, enables informational messages in sub-panels that support it.
#' @param CI Numeric in \code{(0, 1)}; confidence level passed through to the panel helpers. Default \code{0.95}.
#' @param ... Dots for future expansion (currently unused).
#'
#' @return A \pkg{patchwork} object representing the 3×3 grid, returned \emph{invisibly}
#' after being printed.
#'
#' @seealso
#' \code{\link{plot_elu2_panel_biomass}},
#' \code{\link{plot_elu2_panel_f}},
#' \code{\link{plot_elu2_panel_catch}},
#' \code{\link{plot_elu2_panel_bbmsy}},
#' \code{\link{plot_elu2_panel_ffmsy}},
#' \code{\link{plot_elu2_panel_kobe}},
#' \code{\link{plot_elu2_panel_production}},
#' \code{\link{plot_elu2_panel_time_to_bmsy}},
#' \code{\link{plot_elu2_panel_per_model}}
#'
#' @examples
#' \dontrun{
#' # Assuming `fit` is a fitted spictcls object:
#' g <- plot_spictcls_elu2_gg(fit, CI = 0.95)
#' # `g` is returned invisibly; print explicitly if needed:
#' print(g)
#' }
#'
#' @export
plot_spictcls_elu2_gg <- function(x, stamp = get.version(), verbose = TRUE, CI = 0.95, ...) {
  rep <- x
  if (!inherits(rep, "spictcls")) stop("x must be a fitted SPiCT object (class 'spictcls').")
  if (!isTRUE(rep$inp$reportall)) message("inp$reportall = FALSE: showing compact layout.")

  p_biomass <- plot_elu2_panel_biomass(rep, CI = CI)
  p_f       <- plot_elu2_panel_f(rep,       CI = CI)
  p_catch   <- plot_elu2_panel_catch(rep,   CI = CI)
  p_bbmsy   <- plot_elu2_panel_bbmsy(rep,   CI = CI)
  p_ffmsy   <- plot_elu2_panel_ffmsy(rep,   CI = CI)
  p_kobe    <- plot_elu2_panel_kobe(rep,    CI = CI)
  p_prod    <- plot_elu2_panel_production(rep, CI = CI)
  p_time    <- plot_elu2_panel_time_to_bmsy(rep, CI = CI)
  p_priors  <- try(plot_elu2_panel_per_model(rep, do.plot = 1, CI = CI), silent = TRUE)
  if (inherits(p_priors, "try-error") || is.null(p_priors)) {
    p_priors <- ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = "Prior–posterior")
  }

  p_biomass <- .elu2_overlay_truth(p_biomass, rep, "biomass")
  p_f       <- .elu2_overlay_truth(p_f,       rep, "f")
  p_catch   <- .elu2_overlay_truth(p_catch,   rep, "catch")
  p_bbmsy   <- .elu2_overlay_truth(p_bbmsy,   rep, "bbmsy")
  p_ffmsy   <- .elu2_overlay_truth(p_ffmsy,   rep, "ffmsy")

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required for the 3×3 grid. Please install it.")
  }
  g <- patchwork::wrap_plots(
    p_biomass, p_f, p_catch,
    p_bbmsy,   p_ffmsy, p_kobe,
    p_prod,    p_time,  p_priors,
    ncol = 3
  )

  print(g)
  invisible(g)
}
