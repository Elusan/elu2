#' Prior–posterior density (one model) — named ggplot or patchwork grid
#'
#' @description
#' Builds prior–posterior density panels for each **active** prior in a fitted
#' SPiCT/ELU model. If `do.plot = 1`, returns a single ggplot panel (first active
#' prior) immediately; otherwise returns a patchwork grid of all active priors.
#'
#' @param rep A fitted SPiCT report object (`spictcls`).
#' @param model_id Optional model label embedded in messages; default tries `rep$name`.
#' @param do.plot Optional integer cap on number of panels. If `1`, returns the
#'   single panel immediately. Default `NULL` (all active priors).
#' @param stamp Optional message (emitted when returning a grid).
#' @param CI Confidence level (0,1) for parameter extraction. Default `0.95`.
#'
#' @return If `do.plot = 1`, a single `ggplot`. Otherwise, a `patchwork`
#'   object containing all panels (invisible).
#'
#' @examples
#' \dontrun{
#' rep <- fit.spict(inp)
#' p1  <- plot_elu2_panel_priors.elu6_per_model(rep, do.plot = 1)
#' pA  <- plot_elu2_panel_priors.elu6_per_model(rep)
#' }
#'
#' @import ggplot2
#' @importFrom spict get.par get.version
#' @importFrom grid unit
#' @importFrom stats dnorm dgamma qgamma
#' @export
plot_elu2_panel_per_model <- function(rep, model_id = NULL,
                                                  do.plot = NULL,
                                                  stamp = get.version(),
                                                  CI = 0.95) {
  inp <- rep$inp
  useflags <- inp$priorsuseflags
  inds <- which(useflags == 1)
  plots <- list()
  counter <- 0

  if (is.null(model_id)) model_id <- if (!is.null(rep$name)) rep$name else "Model"

  # ------------- file-local helpers (plain comments) -------------
  # Compact theme for density panels
  theme_minimal_compact <- function(base_size = 8, base_family = "") {
    ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12),
        axis.title = ggplot2::element_text(face = "bold", size = 8),
        axis.text  = ggplot2::element_text(size = 8, face = "bold"),
        legend.position = "bottom",
        legend.title = ggplot2::element_blank(),
        legend.text  = ggplot2::element_text(size = 10, face = "bold"),
        legend.spacing.y = grid::unit(0, "pt"),
        legend.spacing.x = grid::unit(0, "pt"),
        legend.margin    = ggplot2::margin(0, 0, 0, 0),
        legend.box.margin= ggplot2::margin(0, 0, 0, 0),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(fill = NA, colour = "grey35", linewidth = 2),
        axis.ticks = ggplot2::element_line(linewidth = 0.5, color = "grey45"),
        axis.ticks.length = grid::unit(3, "pt"),
        strip.background = ggplot2::element_rect(fill = "grey45", color = "grey45", linewidth = 0.5),
        strip.text = ggplot2::element_text(face = "bold", size = ggplot2::rel(1)),
        text = ggplot2::element_text(face = "bold", size = 10),
        plot.margin = ggplot2::margin(2, 2, 2, 2)
      )
  }

  # Optional index suffix formatter
  fd <- function(x) {
    if (is.null(x) || length(x) == 0 || is.na(x) || x == 0) "" else as.character(x)
  }

  # Append catch unit to label
  add.catchunit <- function(lab, cu) {
    cu <- as.character(cu)
    if (nzchar(cu)) paste0(lab, ", ", cu) else lab
  }

  if (length(inds) > 0) {
    for (i in seq_along(inds)) {
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
        prvec <- if (is.list(priorvec)) priorvec[[rr]] else priorvec

        mu <- ifelse(is.na(par[rr, 4]), prvec[1], par[rr, 2])
        sd <- ifelse(is.na(par[rr, 4]), prvec[2], par[rr, 4])

        if (isGamma && is.na(par[rr, 4])) {
          xmin <- 1e-12
          xmax <- stats::qgamma(0.99, shape = mu, rate = sd)
        } else {
          xmin <- mu - 3 * sd
          xmax <- mu + 3 * sd
        }

        xpr <- xpo <- seq(xmin, xmax, length.out = 200)
        priorvals <- if (!isGamma) {
          stats::dnorm(xpr, prvec[1], prvec[2])
        } else {
          stats::dgamma(xpr, prvec[1], prvec[2])
        }

        df_prior <- data.frame(
          x = exp(xpr),
          density = priorvals,
          type = "Prior"
        )

        df_post <- NULL
        if (!is.na(par[rr, 4])) {
          if (isGamma) xpo <- seq(mu - 3 * sd, mu + 3 * sd, length.out = 200)
          posteriorvals <- stats::dnorm(xpo, par[rr, 2], par[rr, 4])
          df_post <- data.frame(
            x = exp(xpo),
            density = posteriorvals,
            type = "Posterior"
          )
        }

        # bind rows without requiring dplyr
        parts <- Filter(Negate(is.null), list(df_prior, df_post))
        df_plot <- do.call(rbind, parts)

        # ensure both legend levels exist (even if NA rows)
        if (!"Prior" %in% df_plot$type)
          df_plot <- rbind(df_plot, data.frame(x = NA, density = NA, type = "Prior"))
        if (!"Posterior" %in% df_plot$type)
          df_plot <- rbind(df_plot, data.frame(x = NA, density = NA, type = "Posterior"))

        x_range_ratio <- max(df_plot$x, na.rm = TRUE) / max(min(df_plot$x, na.rm = TRUE), 1e-10)
        use_log_scale <- x_range_ratio > 10

        x_limits <- range(df_plot$x, na.rm = TRUE)
        no_exp   <- ggplot2::expansion(mult = 0, add = 0)

        p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = x, y = density, color = type)) +
          ggplot2::geom_line(linewidth = 0.8, linetype = "solid", na.rm = TRUE, show.legend = TRUE) +
          ggplot2::geom_area(ggplot2::aes(fill = type), alpha = 0.14, position = "identity", show.legend = FALSE) +
          ggplot2::labs(title = current_nmpl, x = NULL, y = "Density", color = NULL) +
          ggplot2::scale_color_manual(values = c("Prior" = "black", "Posterior" = "red")) +
          ggplot2::scale_fill_manual(values = c("Prior" = "black", "Posterior" = "red")) +
          theme_minimal_compact() +
          ggplot2::theme(
            legend.position = c(0.98, 0.98),
            legend.justification = c("right", "top"),
            legend.background = ggplot2::element_rect(fill = "white", color = NA, linewidth = 0),
            legend.box.background = ggplot2::element_rect(color = "white", linewidth = 0),
            legend.key.width = grid::unit(0.2, "lines"),
            legend.key.height = grid::unit(0, "lines"),
            legend.text = ggplot2::element_text(face = "bold", size = 6),
            legend.title = ggplot2::element_blank(),
            legend.spacing.y = grid::unit(0, "pt"),
            plot.margin = ggplot2::margin(8, 8, 8, 8)
          ) +
          ggplot2::guides(
            color = ggplot2::guide_legend(override.aes = list(linewidth = 2.5), nrow = 2, byrow = TRUE),
            fill = "none"
          )

        if (use_log_scale) {
          p <- p + ggplot2::scale_x_log10(limits = x_limits, expand = no_exp)
        } else {
          p <- p + ggplot2::scale_x_continuous(limits = x_limits, expand = no_exp)
        }

        if (is.na(par[rr, 4]) && !is.na(par[rr, 2])) {
          p <- p + ggplot2::geom_vline(xintercept = exp(par[rr, 2]),
                                       linetype = "dashed", color = "red", linewidth = 1)
        }

        plots[[length(plots) + 1]] <- p
        counter <- counter + 1

        if (!is.null(do.plot) && counter >= do.plot) {
          if (do.plot == 1) return(p)
          break
        }
      }
      if (!is.null(do.plot) && counter >= do.plot) break
    }
  }

  if (length(plots) > 0) {
    if (!is.null(do.plot) && do.plot == 1) return(plots[[1]])
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      stop("Package 'patchwork' is required to combine multiple prior plots. Install it or set do.plot = 1.")
    }
    final_plot <- patchwork::wrap_plots(plots)
    if (!is.null(stamp)) message("elu custom spict")
    return(final_plot)
  }

  invisible(NULL)
}
