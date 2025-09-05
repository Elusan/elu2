# ===========================
# Priors — single model
# ===========================

#' Prior vs posterior density (single model; single or multi-panel)
#'
#' Builds prior–posterior density panels for all **active** priors in a fitted
#' SPiCT model. If \code{do.plot = 1}, returns the **first** active prior panel
#' directly (useful for embedding a single density). Otherwise returns a
#' patchwork of panels. Gamma priors are handled; x-scale is log when dynamic
#' range is wide.
#'
#' @param rep A fitted SPiCT result (\code{spictcls}).
#' @param model_id Optional character used in labeling (not required).
#' @param do.plot Optional integer cap on number of priors to draw; if \code{1},
#'   returns a single \code{ggplot} panel immediately.
#' @param stamp Optional message emitted when returning a patchwork (unchanged).
#' @param CI Confidence level for parameter extraction (\code{0.95} by default).
#'
#' @return Either a \code{ggplot} (single panel) or a \code{patchwork} object.
#' @family ELU2-priors
#' @export
#' @examples
#' \dontrun{
#' # First prior only:
#' priors.elu6_for_one_model(all_models$S1$S1P.SDM, do.plot = 1)
#' # All active priors:
#' priors.elu6_for_one_model(all_models$S1$S1P.SDM)
#' }
priors.elu6_for_one_model <- function(rep, model_id = NULL, do.plot = NULL, stamp = get.version(), CI = 0.95) {


  inp <- rep$inp
  useflags <- inp$priorsuseflags
  inds <- which(useflags == 1)
  ninds <- length(inds)
  plots <- list()
  counter <- 0

  if (is.null(model_id)) model_id <- if (!is.null(rep$name)) rep$name else "Model"

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
        priorvals <- if (!isGamma) {
          dnorm(xpr, prvec[1], prvec[2])
        } else {
          dgamma(xpr, prvec[1], prvec[2])
        }

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

        if (!"Prior" %in% df_plot$type)     df_plot <- rbind(df_plot, data.frame(x = NA, density = NA, type = "Prior"))
        if (!"Posterior" %in% df_plot$type) df_plot <- rbind(df_plot, data.frame(x = NA, density = NA, type = "Posterior"))

        x_range_ratio <- max(df_plot$x, na.rm = TRUE) / max(min(df_plot$x, na.rm = TRUE), 1e-10)
        use_log_scale <- x_range_ratio > 10

        x_limits <- range(df_plot$x, na.rm = TRUE)
        no_exp   <- ggplot2::expansion(mult = 0, add = 0)

        y_top  <- max(df_plot$density, na.rm = TRUE)
        x_left <- min(df_plot$x, na.rm = TRUE)

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
            legend.box.background = ggplot2::element_rect(color = "grey60", linewidth = 0),
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

        # --- enforce do.plot cap, with early-return for single-panel use ----
        if (!is.null(do.plot) && counter >= do.plot) {
          if (do.plot == 1) return(p)   # return the single panel immediately
          break                         # otherwise stop inner loop
        }
      } # end rr loop

      # also break outer loop if we've reached the cap
      if (!is.null(do.plot) && counter >= do.plot) break
    } # end i loop
  }

  if (length(plots) > 0) {
    if (!is.null(do.plot) && do.plot == 1) return(plots[[1]])
    final_plot <- patchwork::wrap_plots(plots)
    if (!is.null(stamp)) message("elu custom spict")
    return(final_plot)
  }

}  # (body unchanged)
