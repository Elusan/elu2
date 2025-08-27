#' @title Minimal Compact Theme for ggplot2
#' @description
#' A lightweight theme tuned for small-panel, publication-style figures.
#' It removes grid lines, draws a visible panel border, and uses compact,
#' bold text suitable for dense layouts (e.g., prior–posterior grids).
#'
#' @param base_size Numeric. Base font size. Default: `8`.
#' @param base_family Character. Base font family. Default: `""` (use default).
#'
#' @return A \code{ggplot2} theme object to be added with \code{+}.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) +
#'   geom_point() +
#'   theme_minimal_compact()
#' }
#' @export
library(ggplot2)
library(grid)       # unit()
library(patchwork)  # wrap_plots()

theme_minimal_compact <- function(base_size = 8, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.title = element_text(face = "bold", size = 8),
      axis.text = element_text(size = 8, face = "bold"),
      #legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 10, face = "bold"),
      legend.spacing.y = unit(0, "pt"),
      legend.spacing.x = unit(0, "pt"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, colour = "grey45", linewidth = 1.5),
      axis.ticks = element_line(linewidth = 0.5, color = "grey45"),
      axis.ticks.length = unit(3, "pt"),
      strip.background = element_rect(fill = "grey45", color = "grey45", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = rel(1)),
      text = element_text(face = "bold", size = 10),
      plot.margin = margin(2, 2, 2, 2)
    )
}

# --- helpers -------------------------------------------------------

#' @title Smart numeric labels
#' @description
#' Formats numbers for axis labels: integers are shown with one decimal
#' (e.g., \code{2.0}); non-integers with two decimals (e.g., \code{2.75}).
#' Returns \code{NA} for \code{NA}/non-finite inputs.
#'
#' @param x Numeric vector.
#'
#' @return Character vector of formatted labels.
#'
#' @examples
#' smart_label(c(1, 2.5, 3, NA, Inf))
#' #> "1.0" "2.50" "3.0" NA NA
#' @noRd
#' @keywords internal
smart_label <- function(x) {
  sapply(x, function(xx) {
    if (is.na(xx) || !is.finite(xx)) return(NA_character_)
    xx_round <- round(xx, 2)
    if (abs(xx_round - round(xx_round)) < .Machine$double.eps^0.5) {
      sprintf("%.1f", xx_round)
    } else sprintf("%.2f", xx_round)
  })
}

#' @title Scientific notation only when needed
#' @description
#' For each value in \code{x}, if the number of meaningful decimal digits
#' (after trimming trailing zeros) exceeds \code{threshold}, format in
#' scientific notation (e.g., \code{"3e-04"}). Otherwise, use
#' \code{smart_label()}.
#'
#' @param x Numeric vector.
#' @param threshold Integer. Max allowed decimal digits before switching
#'   to scientific notation. Default: \code{3}.
#'
#' @return Character vector of formatted labels.
#' @noRd
#' @keywords internal
sci_if_many_decimals <- function(x, threshold = 3) {
  sapply(x, function(xx) {
    if (is.na(xx) || !is.finite(xx)) return(NA_character_)
    parts    <- strsplit(sprintf("%.10f", xx), ".", fixed = TRUE)[[1]]
    dec_part <- sub("0+$", "", parts[2])
    if (nchar(dec_part) > threshold) sprintf("%.0e", xx) else smart_label(xx)
  })
}

# --- main ----------------------------------------------------------

#' @title Prior–posterior density panels for ELU/SPiCT priors
#' @description
#' Draws density curves comparing prior and posterior distributions for all
#' active priors in a fitted ELU/SPiCT model object (\code{rep}). Panels
#' auto-select linear vs log x-scale, remove horizontal expansion so curves
#' touch panel borders, and inset the panel annotation (e.g., model ID) slightly
#' from the top-left border to avoid overlap.
#'
#' @param rep A fitted model object (e.g., result of \code{fit.elu2()}),
#'   containing \code{$inp$} fields such as \code{priors}, \code{priorsuseflags},
#'   and the usual posterior summary available via \code{get.par()}.
#' @param model_id Optional character label shown inside each panel (top-left).
#'   Defaults to \code{rep$name} if available, otherwise \code{"Model"}.
#' @param do.plot Optional integer. If provided, stop after plotting this many
#'   priors (useful for quick previews).
#' @param stamp Optional provenance/version string (printed via \code{message()}).
#'   Default: \code{get.version()}.
#' @param CI Numeric in (0,1). Posterior CI width forwarded to \code{get.par()}.
#'   Default: \code{0.95}.
#'
#' @details
#' This function relies on package-internal helpers and objects:
#' \itemize{
#' \item \code{get.par()} to retrieve posterior summaries on the working scale,
#' \item \code{fd()} and \code{add.catchunit()} for display-name decoration,
#' \item \code{rep$inp$priors}, \code{rep$inp$priorsuseflags} following SPiCT conventions.
#' }
#' For gamma-type priors (e.g., \code{"...gamma"}), x is kept on the natural scale.
#' For log-scale priors, x-values are exponentiated so the panel shows the parameter
#' on its interpretable scale.
#'
#' The annotation is in-set by a small fraction of the x/y ranges
#' (\code{pad_x = 0.0126}, \code{pad_y = 0.035}) and is consistent on log or linear scales.
#'
#' @return A \pkg{patchwork} object combining the individual panels. If no priors
#'   are active, returns \code{NULL}.
#'
#' @examples
#' \dontrun{
#' p <- priors.elu10(fit, model_id = "S1F.SDM")
#' p
#' }
#'
#' @importFrom ggplot2 ggplot geom_line geom_area scale_x_log10 scale_x_continuous
#' @importFrom ggplot2 labs annotate theme element_text element_rect element_blank
#' @importFrom ggplot2 guide_legend unit margin
#' @importFrom patchwork wrap_plots
#' @export
priors.elu10 <- function(rep, model_id = NULL, do.plot = NULL, stamp = get.version(), CI = 0.95) {
  inp <- rep$inp
  useflags <- inp$priorsuseflags
  inds <- which(useflags == 1)
  plots <- list()
  counter <- 0

  if (is.null(model_id)) model_id <- if (!is.null(rep$name)) rep$name else "Model"

  if (length(inds) > 0) {
    for (j in inds) {
      priorvec <- inp$priors[[j]]
      nm <- names(inp$priors)[j]
      isGamma <- FALSE
      nmpl <- sub('log', '', nm)
      nmpl <- sub('gamma', '', nmpl)
      par <- get.par(nm, rep, exp = FALSE, CI = CI)

      if (grepl('gamma', nm)) {
        isGamma <- TRUE
        if (nm == "logngamma") par <- get.par("logn", rep, exp = FALSE, CI = CI)
      }

      repriors <- c('logB', 'logF', 'logBBmsy', 'logFFmsy')
      if (nm %in% repriors) {
        par <- par[priorvec[5], , drop = FALSE]
        nmpl <- paste0(nmpl, fd(priorvec[4]))
        if (nm == 'logB') nmpl <- add.catchunit(nmpl, inp$catchunit)
      }

      for (rr in seq_len(nrow(par))) {
        current_nmpl <- if (nrow(par) > 1) paste0(nmpl, rr) else nmpl
        prvec <- if (is.list(priorvec)) priorvec[[rr]] else priorvec

        mu <- ifelse(is.na(par[rr, 4]), prvec[1], par[rr, 2])
        sd <- ifelse(is.na(par[rr, 4]), prvec[2], par[rr, 4])

        # density support (working scale)
        if (isGamma && is.na(par[rr, 4])) {
          xmin <- 1e-12
          xmax <- qgamma(0.99, shape = mu, rate = sd)
        } else {
          xmin <- mu - 3 * sd
          xmax <- mu + 3 * sd
        }
        xpr <- xpo <- seq(xmin, xmax, length.out = 200)

        # prior curve
        priorvals <- if (!isGamma) dnorm(xpr, prvec[1], prvec[2]) else dgamma(xpr, prvec[1], prvec[2])
        x_prior   <- if (isGamma) xpr else exp(xpr)
        df_prior <- data.frame(x = x_prior, density = priorvals, type = "Prior")

        # posterior curve
        df_post <- NULL
        if (!is.na(par[rr, 4])) {
          if (isGamma) xpo <- seq(mu - 3 * sd, mu + 3 * sd, length.out = 200)
          posteriorvals <- dnorm(xpo, par[rr, 2], par[rr, 4])
          x_post        <- if (isGamma) xpo else exp(xpo)
          df_post <- data.frame(x = x_post, density = posteriorvals, type = "Posterior")
        }

        df_plot <- dplyr::bind_rows(df_prior, df_post)
        if (!"Prior" %in% df_plot$type)     df_plot <- rbind(df_plot, data.frame(x = NA, density = NA, type = "Prior"))
        if (!"Posterior" %in% df_plot$type) df_plot <- rbind(df_plot, data.frame(x = NA, density = NA, type = "Posterior"))

        # choose x scale
        x_range_ratio <- max(df_plot$x, na.rm = TRUE) / max(min(df_plot$x, na.rm = TRUE), 1e-10)
        use_log_scale <- x_range_ratio > 10

        # exact limits (no left/right padding for curves)
        x_limits <- range(df_plot$x, na.rm = TRUE)
        no_exp   <- expansion(mult = 0, add = 0)

        # tiny insets for annotation (avoid touching panel)
        pad_x <- 0.0126
        pad_y <- 0.035
        if (use_log_scale) {
          # move a constant FRACTION of the log-range
          logL <- log10(x_limits[1]); logR <- log10(x_limits[2])
          x_inset <- 10^(logL + pad_x * (logR - logL))
        } else {
          x_inset <- x_limits[1] + pad_x * diff(x_limits)
        }
        y_limits <- range(df_plot$density, na.rm = TRUE)
        y_inset  <- y_limits[2] - pad_y * diff(y_limits)

        p <- ggplot(df_plot, aes(x = x, y = density, color = type)) +
          geom_line(linewidth = 0.6, linetype = "solid", na.rm = TRUE, show.legend = TRUE) +
          geom_area(aes(fill = type), alpha = 0.14, position = "identity", show.legend = FALSE) +
          labs(title = current_nmpl, x = NULL, y = "Density", color = NULL) +
          scale_color_manual(values = c("Prior" = "black", "Posterior" = "red")) +
          scale_fill_manual(values = c("Prior" = "black", "Posterior" = "red")) +
          theme_minimal_compact() +
          theme(
            legend.position = c(0.98, 0.98),
            legend.justification = c("right", "top"),
            legend.background = element_rect(fill = "white", color = NA, linewidth = 0),
            legend.box.background = element_rect(color = "grey60", linewidth = 0),
            legend.key.width = unit(0.2, "lines"),
            legend.key.height = unit(0, "lines"),
            legend.text = element_text(face = "bold", size = 6),
            legend.title = element_blank(),
            legend.spacing.y = unit(0, "pt"),
            plot.margin = margin(8, 8, 8, 8)
          ) +
          guides(
            color = guide_legend(override.aes = list(linewidth = 2.5), nrow = 2, byrow = TRUE),
            fill = "none"
          ) +
          # annotation placed a hair inside the frame
          annotate("text",
                   x = x_inset, y = y_inset,
                   label = model_id,
                   hjust = 0, vjust = 1, fontface = "bold", size = 3, color = "grey15"
          )

        # apply x scale with no horizontal expansion
        if (use_log_scale) {
          p <- p + scale_x_log10(limits = x_limits, expand = no_exp)
        } else {
          p <- p + scale_x_continuous(limits = x_limits, expand = no_exp)
        }
        # optional: a whisper of top headroom
        # p <- p + scale_y_continuous(expand = expansion(mult = c(0, 0.02)))

        # optional: dashed line at posterior mean
        if (is.na(par[rr, 4]) && !is.na(par[rr, 2])) {
          vline_x <- if (isGamma) par[rr, 2] else exp(par[rr, 2])
          p <- p + geom_vline(xintercept = vline_x, linetype = "dashed", color = "red", linewidth = 1)
        }

        plots[[length(plots) + 1]] <- p
        counter <- counter + 1
        if (!is.null(do.plot) && counter >= do.plot) break
      }
    }
  }

  if (length(plots) > 0) {
    final_plot <- wrap_plots(plots)
    if (!is.null(stamp)) message("elu custom spict")
    return(final_plot)
  }
}
