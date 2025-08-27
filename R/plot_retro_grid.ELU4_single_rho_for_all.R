#' Plot retrospective grids (2×2) with Mohn's rho for B, F, B/Bmsy, F/Fmsy
#'
#' @description
#' For each fitted **SPiCT** model in a named list, generate a 2×2 retrospective
#' grid showing:
#' - **Top-left:** Biomass *B[t]* (with 95% CI ribbons),
#' - **Top-right:** Fishing mortality *F[t]*,
#' - **Bottom-left:** Relative biomass *B[t]/B[MSY]*,
#' - **Bottom-right:** Relative fishing mortality *F[t]/F[MSY]*.
#'
#' Each panel overlays the **baseline** and **peel** runs from `retro()` and
#' annotates **Mohn’s ρ** for its corresponding quantity (B, F, B/Bmsy, F/Fmsy)
#' when available (annotation is omitted if the value is `NA`). Peel ordering is
#' standardized as `"All"` (baseline) followed by numeric peels in **descending**
#' order (e.g., `-5, -4, …, -1`), then any non-numeric peel names.
#'
#' @details
#' This function expects each element of `models` to be a fitted **SPiCT** object
#' (class `spictcls`) that already contains a `retro` component produced by
#' `retro()`. Mohn’s rho is computed via `mohns_rho()` for the four quantities
#' and displayed in each respective panel. Quantity names are resolved robustly:
#' for example, *F* falls back from `"logFnotS"` to `"logF"` if needed, and
#' *F/Fmsy* from `"logFFmsynotS"` to `"logFFmsy"`.
#'
#' If `output_dir` is provided, each model’s figure is saved as a PNG named
#' `"retrospective_<model_name>.png"`. Otherwise, the function returns the plots
#' invisibly for further composition/export.
#'
#' The function relies on a user-provided `theme_minimal_compact()` and
#' `make_peel_subtitle()` (not defined here). Make sure these helpers exist in
#' your package namespace.
#'
#' @param models Named list of fitted **SPiCT** models (class `spictcls`), each
#'   with a `$retro` element as produced by `retro()`. Names are used for plot
#'   titles and output filenames.
#' @param width,height Numeric. Figure width and height in inches for exported
#'   PNGs. Defaults: `width = 8`, `height = 5`.
#' @param dpi Numeric. Resolution for exported PNGs. Default: `400`.
#' @param output_dir Character or `NULL`. If non-`NULL`, directory where PNGs
#'   are written (created recursively if it does not exist). If `NULL`, plots are
#'   not written to disk and are only returned.
#'
#' @return
#' (Invisibly) a named list of **patchwork** plot objects (one per model). If
#' `output_dir` is provided, PNG files are also written to disk.
#'
#' @section Panels:
#' - **B[t]** (biomass) with CI ribbon and Mohn’s ρ\eqn{(B)}.
#' - **F[t]** (fishing mortality) with CI ribbon and Mohn’s ρ\eqn{(F)}.
#' - **B[t]/B[MSY]** with CI ribbon and Mohn’s ρ\eqn{(B/B_{MSY})}.
#' - **F[t]/F[MSY]** with CI ribbon and Mohn’s ρ\eqn{(F/F_{MSY})}.
#'
#' @seealso
#' [retro()], [mohns_rho()], [get.par()], [list.quantities()],
#' **patchwork**::[patchwork::plot_annotation()], **ggtext**::[ggtext::element_markdown()]
#'
#' @examples
#' \dontrun{
#' library(spict)
#' data(pol)
#' rep <- fit.spict(pol$albacore)
#' rep <- retro(rep, nretroyear = 5)
#'
#' # Suppose you have a named list of models with $retro present:
#' models <- list(S1P.SDM = rep)
#'
#' # Just build the plots (not saved):
#' out <- plot_retro_grid.ELU4_single_rho_for_all(models)
#'
#' # Or save each model’s figure to FIG/retro:
#' dir.create("FIG/retro", recursive = TRUE, showWarnings = FALSE)
#' plot_retro_grid.ELU4_single_rho_for_all(models, output_dir = "FIG/retro")
#' }
#'
#' @export
#' @encoding UTF-8
#'
#' @import ggplot2
#' @importFrom patchwork plot_annotation
#' @importFrom ggtext element_markdown
#' @importFrom purrr imap_dfr
#' @importFrom grDevices adjustcolor
#' @importFrom tibble tibble
#' @importFrom stats setNames
plot_retro_grid.ELU4_single_rho_for_all <- function(
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
  require(purrr)

  # Fixed peel color palette (as in plot_retro_grid.G4)
  cols <- function() {
    cs <- c("#000000", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
            "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
    c(cs, adjustcolor(cs[-1], 0.5), adjustcolor(cs[-1], 0.2))
  }
  peel_colors_all <- cols()

  # helper: choose first quantity name that exists in this model's report list
  pick_q <- function(fit, candidates) {
    qs <- list.quantities(fit)
    for (nm in candidates) {
      if (nm %in% qs) return(nm)
    }
    # fall back to first; extract_df will error meaningfully if truly absent
    candidates[[1]]
  }

  # helper: NA-safe label builder (returns NULL if rho is NA so we skip annotate)
  make_rho_label <- function(txt) {
    if (is.na(txt)) return(NULL)
    txt
  }

  # Data frame extractor for one parameter
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

  # Panel builder with optional rho annotation
  make_panel <- function(df, ylab, rho_text = NULL, peel_colors_used = NULL) {
    p <- ggplot(df, aes(time, estimate, group = peel)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#B0B0B0", alpha = 0.25) +
      geom_line(aes(color = peel), linewidth = 0.7) +
      labs(x = NULL, y = ylab) +
      theme_minimal_compact() +
      scale_color_manual(values = peel_colors_used, guide = "none") +
      theme(axis.title.y = element_text(face = "bold", size = 11))
    if (!is.null(rho_text)) {
      p <- p + annotate("text", x = Inf, y = Inf, label = rho_text,
                        parse = TRUE, hjust = 1.1, vjust = 1.5, size = 3.5)
    }
    p
  }

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

    # peel colors
    peel_colors_used <- peel_colors_all[seq_len(n_peels)]
    names(peel_colors_used) <- peel_names

    # --- Mohn's rho (NA-safe) ---
    rho <- tryCatch(
      {
        out <- mohns_rho(model, what = c("FFmsy", "BBmsy", "B", "F"))
        out <- round(out, 3)
        # ensure full set and name order
        out[c("FFmsy", "BBmsy", "B", "F")]
      },
      error = function(e) setNames(rep(NA_real_, 4), c("FFmsy", "BBmsy", "B", "F"))
    )

    lab_BB <- make_rho_label(paste0("Mohn*\"'s\"~rho[B/B[MSY]]==", rho["BBmsy"]))
    lab_FF <- make_rho_label(paste0("Mohn*\"'s\"~rho[F/F[MSY]]==",  rho["FFmsy"]))
    lab_B  <- make_rho_label(paste0("Mohn*\"'s\"~rho[B]==",         rho["B"]))
    lab_F  <- make_rho_label(paste0("Mohn*\"'s\"~rho[F]==",         rho["F"]))

    # quantities (with gentle fallback names)
    q_B    <- "logB"
    q_F    <- pick_q(model, c("logFnotS", "logF"))
    q_BB   <- "logBBmsy"
    q_FF   <- pick_q(model, c("logFFmsynotS", "logFFmsy"))

    # Extract time series with uncertainty
    df_B    <- purrr::imap_dfr(runs, ~ extract_df(.x, .y, q_B))
    df_F    <- purrr::imap_dfr(runs, ~ extract_df(.x, .y, q_F))
    df_Bmsy <- purrr::imap_dfr(runs, ~ extract_df(.x, .y, q_BB))
    df_Fmsy <- purrr::imap_dfr(runs, ~ extract_df(.x, .y, q_FF))

    # Standardize peel ordering: All → numeric peels (desc) → others
    all_present   <- "All" %in% peel_names
    other_peels   <- setdiff(peel_names, "All")
    peel_nums     <- suppressWarnings(as.integer(other_peels))
    valid_peels   <- other_peels[!is.na(peel_nums)]
    invalid_peels <- other_peels[is.na(peel_nums)]
    ordered_peels <- c(if (all_present) "All",
                       valid_peels[order(as.integer(valid_peels), decreasing = TRUE)],
                       invalid_peels)

    # Apply factor levels
    df_B$peel    <- factor(df_B$peel,    levels = ordered_peels)
    df_F$peel    <- factor(df_F$peel,    levels = ordered_peels)
    df_Bmsy$peel <- factor(df_Bmsy$peel, levels = ordered_peels)
    df_Fmsy$peel <- factor(df_Fmsy$peel, levels = ordered_peels)

    # Panels (now with rho also on B and F)
    p_B    <- make_panel(df_B,    expression(bold(B[t])),           rho_text = lab_B,  peel_colors_used = peel_colors_used)
    p_F    <- make_panel(df_F,    expression(bold(F[t])),           rho_text = lab_F,  peel_colors_used = peel_colors_used)
    p_Bmsy <- make_panel(df_Bmsy, expression(bold(B[t]/B[MSY])),    rho_text = lab_BB, peel_colors_used = peel_colors_used)
    p_Fmsy <- make_panel(df_Fmsy, expression(bold(F[t]/F[MSY])),    rho_text = lab_FF, peel_colors_used = peel_colors_used)

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
