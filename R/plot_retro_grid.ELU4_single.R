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

#' Plot individual retrospective grids for each SPiCT model in 2x2 layout
#'
#' For each model in the list, generates and optionally saves a 2x2 retrospective grid:
#' Top-left: B[t], Bottom-left: B[t]/B[MSY], Top-right: F[t], Bottom-right: F[t]/F[MSY]
#'
#' @param models A named list of SPiCT model fits (each with $retro)
#' @param width, height, dpi Figure size and resolution
#' @param output_dir If given, saves plots to this folder
#'
#' @return Invisibly returns a list of patchwork plots
#' @export
plot_retro_grid.ELU4_single <- function(
    models,
    width = 8,
    height = 5,
    dpi = 400,
    output_dir = NULL
) {
  require(patchwork)
  require(ggtext)
  require(dplyr)
  require(grid)
  require(ggplot2)

  # Fixed peel color palette (as in plot_retro_grid.G4)
  cols <- function() {
    cs <- c("#000000", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
    c(cs, adjustcolor(cs[-1], 0.5), adjustcolor(cs[-1], 0.2))
  }
  peel_colors_all <- cols()

  plot_list <- list()

  for (model_name in names(models)) {
    model <- models[[model_name]]
    message("Plotting retrospective for: ", model_name)

    if (is.null(model$retro)) {
      warning("Skipping model '", model_name, "' because it lacks a $retro object.")
      next
    }

    runs <- model$retro
    names(runs)[1] <- "All"
    peel_names <- names(runs)
    n_peels <- length(peel_names)

    # Assign peel colors based on order in peel_names
    peel_colors_used <- peel_colors_all[seq_len(n_peels)]
    names(peel_colors_used) <- peel_names

    # Compute Mohn's rho
    rho <- tryCatch(mohns_rho(model, what = c("FFmsy", "BBmsy")) %>% round(3), error = function(e) rep(NA, 2))
    lab_BB <- paste0("Mohn*\"'s\"~rho[B/B[MSY]]==", rho["BBmsy"])
    lab_FF <- paste0("Mohn*\"'s\"~rho[F/F[MSY]]==", rho["FFmsy"])

    # Extract time series with uncertainty
    extract_df <- function(fit, peel, par_name) {
      idx  <- fit$inp$indest
      pars <- get.par(par_name, fit, exp = TRUE, CI = 0.95)[idx, 1:3]
      tibble(
        peel     = peel,
        time     = fit$inp$time[idx],
        estimate = pars[, 2],
        lower    = pars[, 1],
        upper    = pars[, 3]
      )
    }

    df_B    <- purrr::imap_dfr(runs, ~ extract_df(.x, .y, "logB"))
    df_F    <- purrr::imap_dfr(runs, ~ extract_df(.x, .y, "logFnotS"))
    df_Bmsy <- purrr::imap_dfr(runs, ~ extract_df(.x, .y, "logBBmsy"))
    df_Fmsy <- purrr::imap_dfr(runs, ~ extract_df(.x, .y, "logFFmsynotS"))

    # Standardize peel ordering: All → numeric peels (desc) → others
    all_present <- "All" %in% peel_names
    other_peels <- setdiff(peel_names, "All")
    peel_nums   <- suppressWarnings(as.integer(other_peels))
    valid_peels <- other_peels[!is.na(peel_nums)]
    invalid_peels <- other_peels[is.na(peel_nums)]
    ordered_peels <- c(if (all_present) "All", valid_peels[order(as.integer(valid_peels), decreasing = TRUE)], invalid_peels)

    # Apply factor levels
    df_B$peel    <- factor(df_B$peel,    levels = ordered_peels)
    df_F$peel    <- factor(df_F$peel,    levels = ordered_peels)
    df_Bmsy$peel <- factor(df_Bmsy$peel, levels = ordered_peels)
    df_Fmsy$peel <- factor(df_Fmsy$peel, levels = ordered_peels)

    # Build each panel
    make_panel <- function(df, ylab, rho_text = NULL) {
      ggplot(df, aes(time, estimate, group = peel)) +
        geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#B0B0B0", alpha = 0.25) +
        geom_line(aes(color = peel), linewidth = 0.7) +
        labs(x = NULL, y = ylab) +
        theme_minimal_compact() +
        scale_color_manual(values = peel_colors_used, guide = "none") +
        theme(axis.title.y = element_text(face = "bold", size = 11)) +
        {if (!is.null(rho_text)) annotate("text", x = Inf, y = Inf, label = rho_text,
                                          parse = TRUE, hjust = 1.1, vjust = 1.5, size = 3.5)}
    }

    # Create the 2x2 panel grid
    p_B    <- make_panel(df_B,    expression(bold(B[t])))
    p_F    <- make_panel(df_F,    expression(bold(F[t])))
    p_Bmsy <- make_panel(df_Bmsy, expression(bold(B[t]/B[MSY])), lab_BB)
    p_Fmsy <- make_panel(df_Fmsy, expression(bold(F[t]/F[MSY])), lab_FF)

    p_grid <- (p_B + p_F) / (p_Bmsy + p_Fmsy)

    # Final annotated plot
    p_final <- p_grid +
      plot_annotation(
        title = model_name,
        subtitle = make_peel_subtitle(peel_names, peel_colors_used),
        theme = theme(
          plot.title = element_markdown(hjust = 0.5, size = 17, face = "bold"),
          plot.subtitle = element_markdown(hjust = 0.5, size = 15, face = "bold", margin = margin(b = 15))
        )
      )

    # Save if output_dir is set
    if (!is.null(output_dir)) {
      if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
      file <- file.path(output_dir, paste0("retrospective_", model_name, ".png"))
      ggsave(file, plot = p_final, width = width, height = height, dpi = dpi)
      message(" → Saved to: ", file)
    }

    plot_list[[model_name]] <- p_final
  }

  invisible(plot_list)
}
