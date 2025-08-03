#' Minimal Compact Theme for ggplot2
#'
#' A clean and compact ggplot2 theme optimized for density plots used in ELU prior–posterior comparisons.
#'
#' @param base_size Base font size. Default is 8.
#' @param base_family Font family. Default is "" (use system default).
#'
#' @return A `ggplot2::theme` object.
#'
#' @importFrom ggplot2 theme_minimal element_text element_blank element_line
#' @importFrom ggplot2 margin element_rect theme
#' @importFrom grid unit
#' @export
theme_minimal_compact <- function(base_size = 8, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.title = element_text(face = "bold", size = 8),
      axis.text = element_text(size = 8, face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 10, face = "bold"),
      legend.key.size = unit(0.6, "lines"),
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

#' Plot Priors and Posteriors for ELU/SPiCT Model Parameters
#'
#' This function generates density plots comparing prior and posterior distributions for all
#' active priors in a fitted ELU/SPiCT model. It supports log- or linear-scaled x-axes,
#' highlights posterior means, and allows for faceted visualization using `patchwork`.
#'
#' @param rep A fitted model object returned by `fit.elu2()` or equivalent. Must contain a valid `inp` list and estimated parameters.
#' @param model_id Optional character string specifying the model label shown in plot annotations. If `NULL`, uses `rep$name` or defaults to `"Model"`.
#' @param do.plot Integer. Optional limit on the number of prior–posterior panels to display (useful for debugging or limiting output). Default is `NULL` (plot all).
#' @param stamp Optional text string used for internal version tracking. If not `NULL`, prints a message when the function runs.
#' @param CI Numeric (default = 0.95). Confidence interval width used when computing posterior uncertainty.
#' @param return_list Logical (default = `FALSE`). If `TRUE`, returns a named list of individual ggplot objects for each prior; if `FALSE`, returns a `patchwork` object.
#'
#' @details
#' The function iterates through all priors actively used in the model (`inp$priorsuseflags == 1`),
#' computes their prior and posterior densities (on the original scale), and plots them using
#' consistent styling and annotation. Posterior means (if available) are marked with a dashed line.
#' Axis scale (log or linear) is chosen heuristically based on the range of each parameter.
#'
#' This function supports both visual inspection (as a patchwork grid) and programmatic extraction
#' (via `return_list = TRUE`) for use in other plotting layouts like `elu_prior_posterior_grid()`.
#'
#' @return A `patchwork` object (default) or a named list of `ggplot` objects if `return_list = TRUE`.
#' Returns `NULL` if no active priors are found.
#'
#' @importFrom dplyr bind_rows
#' @importFrom ggplot2 ggplot geom_line geom_area scale_color_manual scale_fill_manual labs annotate
#' @importFrom patchwork wrap_plots
#' @importFrom scales log_breaks pretty_breaks
#' @examples
#' \dontrun{
#' fit <- fit.elu2(inp1)
#' priors.elu6(fit, model_id = "S1P")
#' # To return list of ggplot objects:
#' plots <- priors.elu6(fit, model_id = "S1P", return_list = TRUE)
#' }
#' @export
priors.elu6 <- function(rep, model_id = NULL, do.plot = NULL, stamp = get.version(), CI = 0.95, return_list = FALSE) {
  inp <- rep$inp
  useflags <- inp$priorsuseflags
  inds <- which(useflags == 1)
  ninds <- length(inds)
  plots <- list()
  plot_titles <- character(0)  # Store names like "logr", "q1", etc.
  counter <- 0

  # Get model name if not provided
  if (is.null(model_id)) {
    model_id <- if (!is.null(rep$name)) rep$name else "Model"
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
        if (nm == 'logB') nmpl <- add.catchunit(nmpl, inp$catchunit)
      }

      for (rr in seq_len(nrow(par))) {
        current_nmpl <- if (nrow(par) > 1) paste0(nmpl, rr) else nmpl
        plot_titles <- c(plot_titles, current_nmpl)

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
        df_prior <- data.frame(x = exp(xpr), density = priorvals, type = "Prior")

        df_post <- NULL
        if (!is.na(par[rr, 4])) {
          if (isGamma) xpo <- seq(mu - 3 * sd, mu + 3 * sd, length.out = 200)
          posteriorvals <- dnorm(xpo, par[rr, 2], par[rr, 4])
          df_post <- data.frame(x = exp(xpo), density = posteriorvals, type = "Posterior")
        }

        df_plot <- dplyr::bind_rows(df_prior, df_post)
        if (!"Prior" %in% df_plot$type) df_plot <- rbind(df_plot, data.frame(x = NA, density = NA, type = "Prior"))
        if (!"Posterior" %in% df_plot$type) df_plot <- rbind(df_plot, data.frame(x = NA, density = NA, type = "Posterior"))

        x_range_ratio <- max(df_plot$x, na.rm = TRUE) / max(min(df_plot$x, na.rm = TRUE), 1e-10)
        use_log_scale <- x_range_ratio > 10

        y_top <- max(df_plot$density, na.rm = TRUE)
        x_left <- min(df_plot$x, na.rm = TRUE)

        p <- ggplot(df_plot, aes(x = x, y = density, color = type)) +
          geom_line(linewidth = 0.8, linetype = "solid", na.rm = TRUE) +
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
            legend.key.width = unit(0.5, "lines"),
            legend.key.height = unit(0.18, "lines"),
            legend.text = element_text(face = "bold", size = 7),
            legend.title = element_blank(),
            plot.margin = margin(8, 8, 8, 8)
          ) +
          annotate("text", x = x_left, y = y_top * 1.06,
                   label = model_id, hjust = 0, vjust = 1,
                   fontface = "bold", size = 4, color = "grey20")

        if (use_log_scale) {
          p <- p + scale_x_log10(
            breaks = scales::log_breaks(n = 8),
            labels = smart_label
          ) + coord_cartesian(xlim = c(min(df_plot$x, na.rm = TRUE), max(df_plot$x, na.rm = TRUE)))
        } else {
          p <- p + scale_x_continuous(
            breaks = scales::pretty_breaks(n = 8),
            labels = smart_label
          ) + coord_cartesian(xlim = c(min(df_plot$x, na.rm = TRUE), max(df_plot$x, na.rm = TRUE)))
        }

        if (is.na(par[rr, 4]) && !is.na(par[rr, 2])) {
          p <- p + geom_vline(xintercept = exp(par[rr, 2]), linetype = "dashed", color = "red", linewidth = 1)
        }

        plots[[length(plots) + 1]] <- p
        counter <- counter + 1
        if (!is.null(do.plot) && counter >= do.plot) break
      }
    }
  }

  if (length(plots) > 0) {
    if (!is.null(stamp)) message("elu custom spict")
    if (return_list) {
      names(plots) <- plot_titles
      return(plots)
    } else {
      return(patchwork::wrap_plots(plots))
    }
  }
}
