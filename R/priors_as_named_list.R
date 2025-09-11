#' @title Priors as Named ggplot List (identical visuals to \code{priors.elu10})
#'
#' @description
#' Builds one prior–posterior density panel per **active** prior in a fitted ELU/SPiCT
#' object and returns them as a **named list of ggplot objects** (rather than a
#' single patchwork). Visual settings are aligned with \code{priors.elu10}:
#'
#' \itemize{
#'   \item No horizontal padding on the x-axis
#'         (\code{expand = expansion(mult = 0, add = 0)}), so curves “kiss” the frame.
#'   \item Exact x-limits per panel (\code{limits = range(df_plot$x, na.rm = TRUE)}).
#'   \item Gamma prior fix: x is \emph{not} exponentiated for gamma-scale priors/posteriors;
#'         for log-scale priors, x is exponentiated to data scale.
#'   \item Scale-aware annotation padding (top-left), using
#'         \code{pad_x = 0.0126} of the x-range and \code{pad_y = 0.035} of the y-range.
#'   \item Automatic choice of log vs. linear x-scale (log if span ratio > 10).
#'   \item Optional vertical dashed line at the posterior mean (on the data scale)
#'         when available.
#' }
#'
#' This adapter is handy when you want to arrange panels yourself (e.g., with
#' \pkg{patchwork}) or save each parameter plot separately. List element names are
#' parameter names; for multi-row priors (e.g., \code{logsdi1}, \code{logsdi2}) a numeric
#' suffix is appended.
#'
#' @param rep A fitted ELU/SPiCT model object (e.g., from \code{fit.elu2()}), with
#'   \code{$inp$} components \code{priors}, \code{priorsuseflags}, etc.
#' @param model_id Optional character label drawn inside each panel (top-left).
#'   Defaults to \code{rep$name} if present, otherwise \code{"Model"}.
#' @param do.plot Optional integer; if provided, stop after producing this many panels
#'   (useful for previews).
#' @param stamp Optional string printed via \code{message()} (e.g., version tag).
#'   Default: \code{get.version()}.
#' @param CI Numeric in (0,1). Credible interval width passed to \code{get.par()}.
#'   Default \code{0.95}.
#' @param show_suffix_in_title Logical; if \code{TRUE}, append a numeric suffix to the
#'   panel title for multi-row priors (e.g., \code{q1}, \code{q2}). Default \code{TRUE}.
#'
#' @return A **named list** of \pkg{ggplot2} objects. If no active priors, returns an empty list.
#'
#' @details
#' Relies on package-internal helpers such as \code{get.par()}, \code{fd()}, and
#' \code{add.catchunit()}. For gamma-scale priors (names containing \code{"gamma"}),
#' the x-values are kept on their natural scale; for log-scale priors, x is exponentiated.
#'
#' @examples
#' \dontrun{
#' gp_list <- priors_as_named_list(fit, model_id = "S1F.SDM")
#' # Print one panel:
#' gp_list[["logr"]]
#' # Arrange a subset with patchwork:
#' library(patchwork)
#' wrap_plots(gp_list[c("logsdi1","logsdi2","logr")])
#' }
#'
#' @seealso \code{\link{priors.elu10}}, \code{\link{theme_minimal_compact}}
#'
#' @importFrom dplyr bind_rows
#' @importFrom ggplot2 ggplot aes geom_line geom_area labs annotate theme element_text
#' @importFrom ggplot2 element_rect element_blank scale_x_log10 scale_x_continuous
#' @importFrom ggplot2 scale_color_manual scale_fill_manual guides guide_legend expansion margin
#' @importFrom grid unit
#' @importFrom stats dnorm dgamma qgamma
#' @export
# --- Adapter: identical visuals to priors.elu10(), but returns a NAMED list of ggplots ---
priors_as_named_list <- function(rep, model_id = NULL, do.plot = NULL,
                                 stamp = get.version(), CI = 0.95,
                                 show_suffix_in_title = TRUE) {
  inp <- rep$inp
  useflags <- inp$priorsuseflags
  inds <- which(useflags == 1)
  plots <- list()
  n_made <- 0L

  if (is.null(model_id)) model_id <- if (!is.null(rep$name)) rep$name else "Model"
  if (length(inds) == 0L) return(plots)

  for (j in inds) {
    priorvec <- inp$priors[[j]]
    nm <- names(inp$priors)[j]     # e.g. "logsdi", "logq", "logsdc", "logr"
    isGamma <- FALSE

    nmpl <- sub('log', '', nm)
    nmpl <- sub('gamma', '', nmpl)
    par <- get.par(nm, rep, exp = FALSE, CI = CI)

    if (grepl('gamma', nm)) {
      isGamma <- TRUE
      if (nm == "logngamma") par <- get.par("logn", rep, exp = FALSE, CI = CI)
    }

    repriors <- c('logB','logF','logBBmsy','logFFmsy')
    if (nm %in% repriors) {
      par <- par[priorvec[5], , drop = FALSE]
      nmpl <- paste0(nmpl, fd(priorvec[4]))
      if (nm == 'logB') nmpl <- add.catchunit(nmpl, inp$catchunit)
    }

    nrows <- if (is.null(dim(par))) 1L else nrow(par)
    for (rr in seq_len(nrows)) {
      # Key with suffix when multi-row (e.g., "logsdi1", "logq2")
      current_key <- if (nrows > 1L) paste0(nm, rr) else nm

      prvec <- if (is.list(priorvec)) priorvec[[rr]] else priorvec
      mu <- ifelse(is.na(par[rr, 4]), prvec[1], par[rr, 2])
      sd <- ifelse(is.na(par[rr, 4]), prvec[2], par[rr, 4])

      # support on working scale
      if (isGamma && is.na(par[rr, 4])) {
        xmin <- 1e-12; xmax <- qgamma(0.99, shape = mu, rate = sd)
      } else {
        xmin <- mu - 3*sd; xmax <- mu + 3*sd
      }

      xpr <- xpo <- seq(xmin, xmax, length.out = 200)
      priorvals <- if (!isGamma) dnorm(xpr, prvec[1], prvec[2]) else dgamma(xpr, prvec[1], prvec[2])

      # --- gamma prior fix: do NOT exponentiate gamma-scale x; exponentiate log-scale x
      x_prior <- if (isGamma) xpr else exp(xpr)
      df_prior <- data.frame(x = x_prior, density = priorvals, type = "Prior")

      # posterior
      df_post <- NULL
      if (!is.na(par[rr, 4])) {
        if (isGamma) xpo <- seq(mu - 3*sd, mu + 3*sd, length.out = 200)
        posteriorvals <- dnorm(xpo, par[rr, 2], par[rr, 4])
        x_post <- if (isGamma) xpo else exp(xpo)
        df_post <- data.frame(x = x_post, density = posteriorvals, type = "Posterior")
      }

      df_plot <- dplyr::bind_rows(df_prior, df_post)
      if (!"Prior" %in% df_plot$type)     df_plot <- rbind(df_plot, data.frame(x=NA, density=NA, type="Prior"))
      if (!"Posterior" %in% df_plot$type) df_plot <- rbind(df_plot, data.frame(x=NA, density=NA, type="Posterior"))

      # choose x scale
      x_range_ratio <- max(df_plot$x, na.rm = TRUE) / max(min(df_plot$x, na.rm = TRUE), 1e-10)
      use_log_scale <- x_range_ratio > 10

      # --- exact x-limits + no horizontal padding (curves kiss borders)
      x_limits <- range(df_plot$x, na.rm = TRUE)
      no_exp   <- ggplot2::expansion(mult = 0, add = 0)

      # tiny insets for annotation (avoid touching panel), scale-aware
      pad_x <- 0.0126
      pad_y <- 0.035
      if (use_log_scale) {
        logL <- log10(x_limits[1]); logR <- log10(x_limits[2])
        x_inset <- 10^(logL + pad_x * (logR - logL))
      } else {
        x_inset <- x_limits[1] + pad_x * diff(x_limits)
      }
      y_limits <- range(df_plot$density, na.rm = TRUE)
      y_inset  <- y_limits[2] - pad_y * diff(y_limits)

      # Title: add suffix for clarity when multi-row (q1/q2, sdi1/sdi2)
      panel_title <- if (show_suffix_in_title && nrows > 1L) {
        paste0(sub('log','', nm), rr)
      } else {
        sub('log','', nm)
      }

      p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = x, y = density, color = type)) +
        ggplot2::geom_line(linewidth = 0.7, linetype = "solid", na.rm = TRUE, show.legend = TRUE) +
        ggplot2::geom_area(ggplot2::aes(fill = type), alpha = 0.14, position = "identity", show.legend = FALSE) +
        ggplot2::labs(title = panel_title, x = NULL, y = "Density", color = NULL) +
        scale_color_manual(
          values = c("Prior" = "black", "Posterior" = "red"),
          breaks = c("Prior", "Posterior"),
          labels = c("Priors", "Post.")
        ) +
        scale_fill_manual(  # legend hidden; kept for consistency
          values = c("Prior" = "black", "Posterior" = "red"),
          breaks = c("Prior", "Posterior"),
          labels = c("Priors", "Post.")
        ) +
        theme_minimal_compact() +
        ggplot2::theme(
          legend.position = c(0.98, 0.98),
          legend.justification = c("right", "top"),
          legend.background = ggplot2::element_rect(fill = "white", color = NA, linewidth = 0),
          legend.box.background = ggplot2::element_rect(color = "grey60", linewidth = 0),
          legend.key.width  = grid::unit(0.5, "lines"),
          legend.key.height = grid::unit(0.2, "lines"),
          legend.text = ggplot2::element_text(face = "bold", size = 8),
          legend.title = ggplot2::element_blank(),
          legend.spacing.y = grid::unit(0, "pt"),
          plot.margin = ggplot2::margin(8, 8, 8, 8)
        ) +
        ggplot2::guides(
          color = ggplot2::guide_legend(override.aes = list(linewidth = 2.5), nrow = 2, byrow = TRUE),
          fill  = "none"
        ) +
        # inset annotation (avoids overlap with panel border)
        ggplot2::annotate("text",
                          x = x_inset, y = y_inset,
                          label = model_id, hjust = 0, vjust = 1,
                          fontface = "bold", size = 3, color = "grey20"
        )

      # apply x scale with exact limits + no expansion
      if (use_log_scale) {
        p <- p + ggplot2::scale_x_log10(limits = x_limits, expand = no_exp)
      } else {
        p <- p + ggplot2::scale_x_continuous(limits = x_limits, expand = no_exp)
      }
      # keep a little vertical headroom (default y expansion).
      # If you ever want zero vertical padding too, add:
      # p <- p + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0))

      # dashed line at posterior mean (on data scale)
      if (is.na(par[rr, 4]) && !is.na(par[rr, 2])) {
        vline_x <- if (isGamma) par[rr, 2] else exp(par[rr, 2])
        p <- p + ggplot2::geom_vline(xintercept = vline_x, linetype = "dashed", color = "red", linewidth = 1)
      }

      plots[[current_key]] <- p
      n_made <- n_made + 1L
      if (!is.null(do.plot) && n_made >= do.plot) return(plots)
    }
  }

  if (!is.null(stamp)) message("elu custom spict")
  plots
}
