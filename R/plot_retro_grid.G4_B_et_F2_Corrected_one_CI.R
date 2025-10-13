#' Retrospective diagnostics grid (B, F, B/BMSY, F/FMSY) for multiple SPiCT models
#'
#' @description
#' Creates a multi-column grid (one column per model) with four panels per model:
#' biomass \eqn{B_t}, fishing mortality \eqn{F_t}, relative biomass \eqn{B_t/B_{MSY}},
#' and relative fishing mortality \eqn{F_t/F_{MSY}}. Each panel overlays the baseline
#' run (`"All"`) and all retrospective peels; a 95% CI ribbon is drawn **only** for the
#' `"All"` run, while all runs (including `"All"`) are shown as colored lines. When
#' available, Mohn’s \eqn{\rho} is annotated in each panel for the corresponding quantity.
#'
#' @details
#' - **Inputs:** `models` must be a **named** list of fitted SPiCT objects with a non-empty
#'   `$retro` component (the output of `retro()`), where list names become the column headers.
#' - **Peel handling:** Peel names are normalized so the first element is `"All"`. Peels are
#'   ordered deterministically: `"All"`
#'   \if{html}{\out{→}}\if{latex}{\eqn{\rightarrow}{->}}
#'   numeric peels in **descending** order
#'   \if{html}{\out{→}}\if{latex}{\eqn{\rightarrow}{->}}
#'   any non-numeric peels.
#' - **Colors:** `peel_colors` may be `NULL`, a character vector, or a function returning
#'   a vector of colors given `n`. If `palette = "JABBA"`, a fixed fallback palette is used
#'   unless `jabba_colors` exists. If **MetBrewer** is installed and `palette` matches a valid
#'   MetBrewer palette, it is used; otherwise a hue palette from **scales** is used.
#' - **Quantities extracted:** `get.par(..., exp = TRUE, CI = 0.95)` is called for
#'   `logB`, `logFnotS`, `logBBmsy`, and `logFFmsynotS` using indices in `inp$indest`.
#' - **Theme & layout:** Uses a compact minimal theme (bold axis text, panel border, no grid).
#'   Title/legend spacing is controlled via margins on `plot.title` and `plot.subtitle`
#'   (space **below** the title and **above/below** the legend). If **ggtext** is available,
#'   title/subtitle are rendered with `element_markdown()`.
#' - **Legend (Peels):** The legend content is provided by an external helper
#'   `make_peel_subtitle(peel_names, peel_colors)` (must exist in the search path) and is
#'   passed to `plot_annotation(subtitle = ...)`.
#'
#' @param models A **named** list of fitted SPiCT models with `$retro` already computed.
#'   Names are used as column headers (e.g., `"S1P"`, `"S1F"`, `"S1S"`).
#' @param scenario_name Character. Title displayed above the grid. Default: `"Scenario 1"`.
#' @param peel_colors `NULL`, a character vector of colors, or a function that returns a
#'   vector of colors given `n` peels (e.g., `function(n) scales::hue_pal()(n)`).
#' @param palette Character. `"JABBA"` for a fixed fallback palette (unless `jabba_colors`
#'   is present). If **MetBrewer** is installed, any valid MetBrewer palette name is accepted.
#'   Otherwise a hue palette is used.
#' @param width,height,dpi Numeric. Device size and resolution **only used when**
#'   `output_dir` is non-`NULL`. Defaults: `width = 20`, `height = 8`, `dpi = 400`.
#' @param output_dir `NULL` or directory path. If non-`NULL`, the plot is saved as a PNG
#'   in this directory (created recursively if needed).
#' @param filename Optional file name for the PNG. If `NULL`, a name is derived from
#'   `scenario_name` (non-alphanumerics replaced with `_`).
#'
#' @return A **patchwork** plot object representing the combined grid. If `output_dir`
#'   is provided, a PNG file is also written to disk.
#'
#' @section Dependencies:
#' Requires **ggplot2** and **patchwork**. Optionally uses **ggtext** (for markdown
#' titles/subtitles), **MetBrewer** (palettes), **scales** (fallback palette),
#' and **spict** for `get.par()` if a user-defined `get.par` is not found.
#'
#' @seealso spict::retro(), `make_peel_subtitle()`, `theme_minimal_compact()`,
#'   and optionally `mohns_rho()`.
#'
#' @examples
#' \dontrun{
#' # Suppose 'mods' is a named list of SPiCT fits with $retro present:
#' p <- plot_retro_grid.G4_B_et_F2_Corrected_one_CI(
#'   models = mods,
#'   scenario_name = "Scenario 1",
#'   palette = "JABBA"
#' )
#' print(p)
#'
#' # Save to file:
#' plot_retro_grid.G4_B_et_F2_Corrected_one_CI(
#'   models = mods,
#'   scenario_name = "Scenario 1",
#'   output_dir = "FIG/retro_grids"
#' )
#' }
#'
#' @export
#' @encoding UTF-8
#' @import ggplot2
#' @importFrom patchwork plot_annotation
#' @importFrom ggtext element_markdown
plot_retro_grid.G4_B_et_F2_Corrected_one_CI <- function(
    models,
    scenario_name = "Scenario 1",
    peel_colors = NULL,
    palette = "JABBA",
    width = 20, height = 8, dpi = 400,
    output_dir = NULL, filename = NULL
) {
  # ---- Basic validation ------------------------------------------------------
  if (is.null(models) || !length(models)) {
    stop("`models` must be a non-empty named list of fitted SPiCT objects with `$retro`.", call. = FALSE)
  }
  if (is.null(names(models)) || any(!nzchar(names(models)))) {
    stop("`models` must be a *named* list (names become column headers).", call. = FALSE)
  }

  .has_pkg <- function(p) isTRUE(requireNamespace(p, quietly = TRUE))
  if (!.has_pkg("ggplot2")) stop("Package 'ggplot2' is required.", call. = FALSE)
  if (!.has_pkg("patchwork")) stop("Package 'patchwork' is required.", call. = FALSE)
  has_ggtext    <- .has_pkg("ggtext")
  has_MetBrewer <- .has_pkg("MetBrewer")
  has_scales    <- .has_pkg("scales")

  .exists_fun <- function(name) exists(name, mode = "function", inherits = TRUE)
  .get_fun    <- function(name) get(name, inherits = TRUE)

  .apply_theme <- local({
    if (.exists_fun("theme_minimal_compact")) {
      .get_fun("theme_minimal_compact")
    } else {
      function(base_size = 8, base_family = "") {
        ggplot2::theme_minimal(base_size = base_size, base_family = "") +
          ggplot2::theme(
            plot.title       = ggplot2::element_text(hjust = 0.5, face = "bold", size = 10),
            axis.title       = ggplot2::element_text(face = "bold", size = 8),
            axis.text        = ggplot2::element_text(size = 7, face = "bold"),
            legend.position  = "bottom",
            legend.title     = ggplot2::element_blank(),
            legend.text      = ggplot2::element_text(size = 8, face = "bold"),
            panel.grid       = ggplot2::element_blank(),
            panel.border     = ggplot2::element_rect(fill = NA, colour = "grey45", linewidth = 0.5),
            strip.background = ggplot2::element_rect(fill = "grey90", color = "grey80", linewidth = 0.5),
            strip.text       = ggplot2::element_text(face = "bold", size = ggplot2::rel(1)),
            plot.margin      = ggplot2::margin(4, 4, 4, 4)
          )
      }
    }
  })

  .get_par <- function(...) {
    if (.exists_fun("get.par")) .get_fun("get.par")(...)
    else if (.has_pkg("spict")) spict::get.par(...)
    else stop("Cannot find `get.par()` (neither your version nor spict::get.par).", call. = FALSE)
  }

  .mohns_rho <- function(...) {
    if (.exists_fun("mohns_rho")) .get_fun("mohns_rho")(...)
    else c(FFmsy = NA_real_, BBmsy = NA_real_, B = NA_real_, F = NA_real_)
  }

  .JABBA_FIXED <- c("#000000","#1B9E77","#D95F02","#7570B3","#E7298A",
                    "#66A61E","#E6AB02","#A6761D","#666666","#1F78B4")

  .coerce_colors <- function(x, n, nm = NULL) {
    .as_chr_len <- function(v, n) { v <- as.character(v); if (!length(v)) v <- rep("", n); rep(v, length.out = n) }
    if (is.null(x)) return(NULL)
    if (is.function(x)) {
      res <- tryCatch({
        fml <- names(formals(x))
        if ("n" %in% fml) x(n = n) else if (length(fml) == 0) x() else x(n)
      }, error = function(e) NULL)
      if (is.null(res)) res <- tryCatch(x(n), error = function(e) NULL)
      if (is.null(res)) res <- .JABBA_FIXED
      res <- .as_chr_len(res, n)
    } else if (is.character(x)) {
      res <- .as_chr_len(x, n)
    } else {
      stop("`peel_colors` must be NULL, a character vector, or a function returning a character vector.", call. = FALSE)
    }
    if (!is.null(nm)) names(res) <- nm
    res
  }

  mycols <- if (is.null(peel_colors)) { if (.exists_fun("cols")) .get_fun("cols") else NULL } else peel_colors

  make_model_panels <- function(model_fit, model_nm, peel_colors_input) {
    if (is.null(model_fit) || is.null(model_fit$retro) || !length(model_fit$retro)) {
      stop("Model '", model_nm, "' lacks a non-empty `$retro`. Run retro() first.", call. = FALSE)
    }
    runs <- model_fit$retro
    if (length(runs) >= 1L) names(runs)[1] <- "All"

    peel_names_raw <- names(runs)
    if (is.null(peel_names_raw) || any(!nzchar(peel_names_raw))) {
      peel_names_raw <- paste0("Peel", seq_along(runs))
      names(runs) <- peel_names_raw
      peel_names_raw[1] <- "All"; names(runs)[1] <- "All"
    }
    n_peels <- length(runs)

    peel_colors_final <- NULL
    if (!is.null(peel_colors_input)) {
      peel_colors_final <- .coerce_colors(peel_colors_input, n_peels, peel_names_raw)
    } else if (toupper(palette) == "JABBA" && exists("jabba_colors", inherits = TRUE)) {
      peel_colors_final <- .coerce_colors(get("jabba_colors", inherits = TRUE), n_peels, peel_names_raw)
    } else if (toupper(palette) == "JABBA") {
      peel_colors_final <- .coerce_colors(.JABBA_FIXED, n_peels, peel_names_raw)
    } else if (!is.null(mycols)) {
      peel_colors_final <- .coerce_colors(mycols, n_peels, peel_names_raw)
    } else if (!is.null(palette) && has_MetBrewer && palette %in% MetBrewer::list_met_palettes()) {
      peel_colors_final <- MetBrewer::met.brewer(palette, n_peels); names(peel_colors_final) <- peel_names_raw
    } else {
      if (!has_scales) stop("Package 'scales' is required for fallback palette.", call. = FALSE)
      peel_colors_final <- scales::hue_pal()(n_peels); names(peel_colors_final) <- peel_names_raw
    }

    extract_df <- function(fit, peel, par_name, exponentiate = TRUE) {
      idx <- fit$inp$indest
      if (is.null(idx) || !length(idx)) stop("`inp$indest` missing/empty in one of the retro fits (peel ", peel, ").", call. = FALSE)
      mat <- tryCatch(.get_par(par_name, fit, exp = exponentiate, CI = 0.95), error = function(e) NULL)
      if (is.null(mat)) stop("get.par('", par_name, "') failed for peel ", peel, ".", call. = FALSE)
      if (!is.matrix(mat) || ncol(mat) < 3) stop("Unexpected structure from get.par('", par_name, "').", call. = FALSE)
      if (length(idx) > nrow(mat)) idx <- idx[seq_len(nrow(mat))]
      data.frame(peel = peel, time = fit$inp$time[idx],
                 lower = mat[idx, 1], estimate = mat[idx, 2], upper = mat[idx, 3],
                 stringsAsFactors = FALSE)
    }

    df_B  <- do.call(rbind, lapply(seq_along(runs), function(i) extract_df(runs[[i]], names(runs)[i], "logB",          TRUE)))
    df_F  <- do.call(rbind, lapply(seq_along(runs), function(i) extract_df(runs[[i]], names(runs)[i], "logFnotS",     TRUE)))
    df_BB <- do.call(rbind, lapply(seq_along(runs), function(i) extract_df(runs[[i]], names(runs)[i], "logBBmsy",     TRUE)))
    df_FF <- do.call(rbind, lapply(seq_along(runs), function(i) extract_df(runs[[i]], names(runs)[i], "logFFmsynotS", TRUE)))

    peel_levels <- {
      all_present <- "All" %in% peel_names_raw
      others <- setdiff(peel_names_raw, "All")
      nums   <- suppressWarnings(as.integer(others))
      valid  <- others[!is.na(nums)]
      inval  <- others[is.na(nums)]
      c(if (all_present) "All", valid[order(as.integer(valid), decreasing = TRUE)], inval)
    }

    df_B$peel  <- factor(df_B$peel,  levels = peel_levels)
    df_F$peel  <- factor(df_F$peel,  levels = peel_levels)
    df_BB$peel <- factor(df_BB$peel, levels = peel_levels)
    df_FF$peel <- factor(df_FF$peel, levels = peel_levels)

    if (!all(peel_levels %in% names(peel_colors_final))) {
      peel_colors_final <- .coerce_colors(peel_colors_final, length(peel_levels), peel_levels)
    } else {
      peel_colors_final <- peel_colors_final[peel_levels]
    }

    rho_vals <- tryCatch(
      .mohns_rho(model_fit, what = c("FFmsy","BBmsy","B","F")),
      error = function(e) c(FFmsy = NA_real_, BBmsy = NA_real_, B = NA_real_, F = NA_real_)
    )
    rho_vals <- suppressWarnings(round(rho_vals, 3))
    lab_BB <- paste0("Mohn*\"'s\"~rho[B/B[MSY]]==", if (is.finite(rho_vals["BBmsy"])) rho_vals["BBmsy"] else "NA")
    lab_FF <- paste0("Mohn*\"'s\"~rho[F/F[MSY]]==", if (is.finite(rho_vals["FFmsy"])) rho_vals["FFmsy"] else "NA")
    lab_B  <- paste0("Mohn*\"'s\"~rho[B]==",        if (is.finite(rho_vals["B"]))      rho_vals["B"]      else "NA")
    lab_F  <- paste0("Mohn*\"'s\"~rho[F]==",        if (is.finite(rho_vals["F"]))      rho_vals["F"]      else "NA")

    ci_gray <- "gray30"

    make_panel <- function(df, ylab, rho_text = NULL) {
      df_all   <- df[df$peel == "All", , drop = FALSE]
      p <- ggplot2::ggplot() +
        { if (nrow(df_all)) ggplot2::geom_ribbon(data = df_all,
                                                 ggplot2::aes(x = .data$time, ymin = .data$lower, ymax = .data$upper),
                                                 fill = ci_gray, alpha = 0.25, inherit.aes = FALSE) else NULL } +
        ggplot2::geom_line(data = df, ggplot2::aes(x = .data$time, y = .data$estimate, color = .data$peel),
                           linewidth = 0.6) +
        ggplot2::labs(x = NULL, y = ylab) +
        .apply_theme() +
        ggplot2::scale_color_manual(values = peel_colors_final, guide = "none") +
        ggplot2::theme(axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 2))) +
        { if (!is.null(rho_text)) ggplot2::annotate("text", x = Inf, y = Inf, label = rho_text,
                                                    parse = TRUE, hjust = 1.05, vjust = 1.5, size = 3) else NULL }
      p
    }

    list(
      B     = make_panel(df_B,  expression(bold(B[t])),        rho_text = lab_B),
      F     = make_panel(df_F,  expression(bold(F[t])),        rho_text = lab_F),
      BBmsy = make_panel(df_BB, expression(bold(B[t]/B[MSY])), rho_text = lab_BB),
      FFmsy = make_panel(df_FF, expression(bold(F[t]/F[MSY])), rho_text = lab_FF),
      peel_names   = peel_names_raw,
      peel_colors  = peel_colors_final
    )
  }

  model_order <- names(models)
  cols_list <- vector("list", length(model_order)); names(cols_list) <- model_order
  subtitle_peel_names <- NULL; subtitle_peel_cols <- NULL

  for (j in seq_along(model_order)) {
    nm <- model_order[j]; fit <- models[[nm]]
    panels <- make_model_panels(fit, nm, mycols)
    if (is.null(subtitle_peel_names)) { subtitle_peel_names <- panels$peel_names; subtitle_peel_cols <- panels$peel_colors }
    header_plot <- ggplot2::ggplot() + ggplot2::theme_void() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = nm,
                        size = 5, fontface = "bold", hjust = 0.5, vjust = 0.5) +
      ggplot2::theme(plot.margin = ggplot2::margin(2.5, 0, 0, 0, "pt"))
    cols_list[[j]] <- patchwork::wrap_plots(header_plot, panels$B, panels$F, panels$BBmsy, panels$FFmsy,
                                            ncol = 1, heights = c(0.16, 1, 1, 1, 1))
  }

  grid_plot <- cols_list[[1]]
  if (length(cols_list) >= 2L) for (k in 2:length(cols_list)) grid_plot <- grid_plot | cols_list[[k]]

  if (!.exists_fun("make_peel_subtitle")) {
    stop("`make_peel_subtitle()` must exist (as in your ELU4 function) to build the legend.", call. = FALSE)
  }
  subtitle_html <- get("make_peel_subtitle", inherits = TRUE)(subtitle_peel_names, subtitle_peel_cols)

  # --- ONLY CHANGE BELOW: add vertical spacing via margins on title/subtitle ---
  final <- grid_plot +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.title = if (has_ggtext)
        ggtext::element_markdown(hjust = 0.5, size = 17, face = "bold",
                                 margin = ggplot2::margin(b = 10)) else
                                   ggplot2::element_text(hjust = 0.5, size = 17, face = "bold",
                                                         margin = ggplot2::margin(b = 10)),
      plot.subtitle = if (has_ggtext)
        ggtext::element_markdown(hjust = 0.5, size = 15, face = "bold",
                                 margin = ggplot2::margin(t = 8, b = 12)) else
                                   ggplot2::element_text(hjust = 0.5, size = 15, face = "bold",
                                                         margin = ggplot2::margin(t = 8, b = 12))
    ) &
    patchwork::plot_annotation(title = scenario_name, subtitle = subtitle_html)

  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    if (is.null(filename)) {
      scen_lab <- gsub("[^A-Za-z0-9]+", "_", scenario_name)
      filename <- paste0("retrospective_", scen_lab, ".png")
    }
    outpath <- file.path(output_dir, filename)
    ggplot2::ggsave(filename = outpath, plot = final, width = width, height = height, dpi = dpi)
    message("Saved to: ", outpath)
  }

  final
}
