library(ggplot2)
library(grid)  # for unit()
theme_minimal_compact <- function(base_size = 8, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.title = element_text(face = "bold", size = 8),
      axis.text = element_text(size = 8, face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 3, face = "bold"),
      legend.key.size = unit(0.1, "lines"),
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


#' Format Numeric Values for Smart Axis Labels
#'
#' Provides rounded numeric labels for plotting, showing one decimal place for integer values (e.g. 2.0) and two for non-integers (e.g. 2.75). Handles NA and Inf gracefully.
#'
#' @param x A numeric vector of values to format.
#'
#' @return A character vector of formatted labels.
#' @examples
#' smart_label(c(1, 2.5, 3, NA, Inf))
#' # [1] "1.0" "2.50" "3.0" NA NA
#' @export
smart_label <- function(x) {
  sapply(x, function(xx) {
    if (is.na(xx) || !is.finite(xx)) return(NA_character_) # Safely handle NA and Inf
    xx_round <- round(xx, 2)
    # If value is an integer after rounding, show with .0
    if (abs(xx_round - round(xx_round)) < .Machine$double.eps^0.5) {
      return(sprintf("%.1f", xx_round)) # e.g., 1.0, 2.0, etc.
    } else {
      return(sprintf("%.2f", xx_round)) # e.g., 2.50, 2.75, etc.
    }
  })
}

#' Plot Priors and Posteriors for ELU Model Parameters
#'
#' Generates density plots comparing prior and posterior distributions for all active priors in an ELU/SPiCT model fit, with intelligent x-axis scaling and plot annotations.
#'
#' @param rep A fitted model object (typically output from \code{fit.elu2()}).
#' @param model_id Optional label for the model, shown on the plot. Defaults to \code{rep$name} if available.
#' @param do.plot Integer; maximum number of priors to plot. Default: \code{NULL} (all).
#' @param stamp Optional text annotation for provenance/version info (used internally).
#' @param CI Confidence interval width for the posterior. Default: 0.95.
#'
#' @details
#' For each parameter with an active prior, this function computes the prior and posterior densities,
#' chooses an appropriate x-axis scale (logarithmic or linear), and plots the results. Plots are
#' combined into a single layout using \code{patchwork::wrap_plots}. If \code{do.plot} is specified,
#' plotting stops after that many priors.
#'
#' @return A \code{patchwork} plot object combining the individual prior-posterior plots.
#'
#' @importFrom dplyr bind_rows
#' @importFrom ggplot2 ggplot geom_line geom_area scale_color_manual scale_fill_manual labs theme_minimal theme annotate
#' @importFrom patchwork wrap_plots
#' @importFrom scales log_breaks pretty_breaks
#' @examples
#' \dontrun{
#' fit <- fit.elu2(inp)
#' priors.elu6(fit, model_id = "Scenario 1")
#' }
#' @export
priors.elu66 <- function(rep, model_id = NULL, do.plot = NULL, stamp = get.version(), CI = 0.95) {
  inp <- rep$inp
  useflags <- inp$priorsuseflags
  inds <- which(useflags == 1)
  ninds <- length(inds)
  plots <- list()
  counter <- 0

  if (is.null(model_id)) {
    model_id <- rep$name %||% "Model"
  }

  if (ninds > 0) {
    for (i in seq_len(ninds)) {
      j <- inds[i]
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
        if (nm == 'logB') {
          nmpl <- add.catchunit(nmpl, inp$catchunit)
        }
      }

      for (rr in seq_len(nrow(par))) {
        current_nmpl <- if (nrow(par) > 1) paste0(nmpl, rr) else nmpl
        prvec <- if (is.list(priorvec)) priorvec[[rr]] else priorvec

        mu <- ifelse(is.na(par[rr, 4]), prvec[1], par[rr, 2])
        sd <- ifelse(is.na(par[rr, 4]), prvec[2], par[rr, 4])

        if (isGamma && is.na(par[rr, 4])) {
          xmin <- 1e-12
          xmax <- qgamma(0.99, shape = mu, rate = sd)
        } else {
          xmin <- mu - 3 * sd
          xmax <- mu + 3 * sd
        }

        xpr <- xpo <- seq(xmin, xmax, length.out = 200)
        priorvals <- if (!isGamma) dnorm(xpr, prvec[1], prvec[2]) else dgamma(xpr, prvec[1], prvec[2])

        df_prior <- data.frame(
          x = exp(xpr),
          density = priorvals,
          type = "Prior"
        )

        df_post <- NULL
        if (!is.na(par[rr, 4])) {
          if (isGamma) xpo <- seq(mu - 3 * sd, mu + 3 * sd, length.out = 200)
          posteriorvals <- dnorm(xpo, par[rr, 2], par[rr, 4])
          df_post <- data.frame(
            x = exp(xpo),
            density = posteriorvals,
            type = "Posterior"
          )
        }

        df_plot <- dplyr::bind_rows(df_prior, df_post)

        if (!"Prior" %in% df_plot$type) df_plot <- rbind(df_plot, data.frame(x = NA, density = NA, type = "Prior"))
        if (!"Posterior" %in% df_plot$type) df_plot <- rbind(df_plot, data.frame(x = NA, density = NA, type = "Posterior"))

        x_range_ratio <- max(df_plot$x, na.rm = TRUE) / max(min(df_plot$x, na.rm = TRUE), 1e-10)
        use_log_scale <- x_range_ratio > 10

        y_top <- max(df_plot$density, na.rm = TRUE)
        x_left <- min(df_plot$x, na.rm = TRUE)

        p <- ggplot(df_plot, aes(x = x, y = density, color = type)) +
          geom_line(linewidth = 0.5, linetype = "solid", na.rm = TRUE) +
          geom_area(aes(fill = type), alpha = 0.14, position = "identity", show.legend = FALSE) +
          labs(title = current_nmpl, x = NULL, y = "Density", color = NULL) +
          scale_color_manual(values = c("Prior" = "black", "Posterior" = "red")) +
          scale_fill_manual(values = c("Prior" = "black", "Posterior" = "red")) +
          theme_minimal_compact() +
          theme(
            legend.position = c(0.98, 0.98),
            legend.justification = c("right", "top"),
            legend.background = element_rect(fill = "white", color = NA),
            legend.box.background = element_rect(color = "grey60"),
            legend.key.width = unit(0.2, "lines"),
            legend.key.height = unit(0.05, "lines"),
            legend.text = element_text(face = "bold", size = 2),
            plot.margin = margin(8, 8, 8, 8)
          ) +
          guides(color = guide_legend(override.aes = list(linewidth = 2.2), nrow = 2)) +
          annotate("text", x = x_left, y = y_top * 1.06, label = model_id,
                   hjust = 0, vjust = 1, fontface = "bold", size = 4, color = "grey20")

        if (use_log_scale) {
          p <- p + scale_x_log10()
        } else {
          p <- p + scale_x_continuous()
        }

        if (is.na(par[rr, 4]) && !is.na(par[rr, 2])) {
          p <- p + geom_vline(xintercept = exp(par[rr, 2]), linetype = "dashed", color = "red", linewidth = 1)
        }

        # Use the parameter name (nm) as key — matches what elu_prior_posterior_grid2() expects
        plots[[nm]] <- p

        counter <- counter + 1
        if (!is.null(do.plot) && counter >= do.plot) break
      }
    }
  }

  if (length(plots) > 0) {
    if (!is.null(stamp)) message("elu custom spict")
    return(plots)  # ✅ named list of ggplot objects
  } else {
    return(list())  # empty list if no priors used
  }
}
