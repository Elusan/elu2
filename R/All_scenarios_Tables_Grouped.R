#' Export Combined SDM+GLM SPiCT Scenario Tables with Grouped Separation
#'
#' @description
#' Automatically summarizes all models, creates wide summary tables, and exports
#' a fully formatted LaTeX table for each scenario, with left-aligned columns and
#' group separation (shown in LaTeX as \verb{\addlinespace} and \verb{\hline})
#' after desired parameter groups.
#'
#' @param model_list Named list of SPiCT model fits for all scenarios (names like \code{S1F.SDM}).
#' @param output_dir Directory for \code{.tex} files (must exist).
#' @param scenarios Character vector of scenario names; default auto-detected from \code{names(model_list)}.
#' @param param_groups Optional list of character vectors defining parameter groupings (for separation).
#' @param verbose Logical; print progress messages?
#' @param caption Optional LaTeX caption override (character). If \code{NULL}, a default is used.
#'
#' @return Invisibly returns a named list of \pkg{kableExtra} table objects (one per scenario).
#' @export
All_scenarios_Tables_Grouped <- function(
    model_list,
    output_dir = ".",
    scenarios = NULL,
    param_groups = NULL,
    verbose = TRUE,
    caption = NULL
) {
  if (!requireNamespace("kableExtra", quietly = TRUE)) {
    stop("Package 'kableExtra' is required for table export.")
  }
  if (!dir.exists(output_dir)) stop("Output directory does not exist!")
  if (is.null(names(model_list))) stop("`model_list` must be a *named* list (e.g., 'S1F.SDM').")

  # ---- default groups/order (kept as in your version) ----
  if (is.null(param_groups)) {
    group1 <- c("$r$", "$m$", "$K$",
                "$\\sigma_B$", "$\\sigma_F$",
                "$\\sigma_{I_{1}}$", "$\\sigma_{I_{2}}$", "$\\sigma_C$")
    group2 <- c("$B_{\\mathrm{MSY}}$", "$F_{\\mathrm{MSY}}$", "$\\mathrm{MSY}$")
    group3 <- c("$B_{2022}$", "$F_{2022}$", "$B_{2022}/B_{\\mathrm{MSY}}$", "$F_{2022}/F_{\\mathrm{MSY}}$")
    param_order  <- c(group1, group2, group3)
    param_groups <- list(group1, group2, group3)
  } else {
    param_order <- unlist(param_groups)
  }

  # ---- auto-detect scenarios ----
  if (is.null(scenarios) || length(scenarios) == 0) {
    scen_raw <- gsub("^((S\\d+)).*$", "\\1", names(model_list))
    scen_keep <- scen_raw[grepl("^S\\d+$", scen_raw)]
    scenarios <- unique(scen_keep)
    scenarios <- scenarios[order(as.integer(sub("^S(\\d+)$", "\\1", scenarios)))]
    if (length(scenarios) == 0) stop("Could not auto-detect scenarios from names(model_list).")
  }

  # ---- helpers ----
  pivot_to_wide <- function(df) {
    cols_needed <- c("Parameters",
                     "Fox_estimate", "Fox_low", "Fox_up",
                     "Schaefer_estimate", "Schaefer_low", "Schaefer_up",
                     "Pella_estimate", "Pella_low", "Pella_up")
    df_clean <- df[, cols_needed]
    keep <- !apply(df_clean, 1, function(row) all(is.na(row[-1])))
    df_clean <- df_clean[keep, , drop = FALSE]
    unique_params <- unique(as.character(df_clean$Parameters))
    wide <- data.frame(Parameters = unique_params, stringsAsFactors = FALSE)
    for (param in unique_params) {
      rows <- df_clean[df_clean$Parameters == param, , drop = FALSE]
      for (mod in c("Fox", "Schaefer", "Pella")) {
        est <- rows[[paste0(mod, "_estimate")]]
        low <- rows[[paste0(mod, "_low")]]
        up  <- rows[[paste0(mod, "_up")]]
        idx <- which(!is.na(est))
        if (length(idx) == 0) {
          est_val <- low_val <- up_val <- NA
        } else {
          est_val <- est[idx[1]]; low_val <- low[idx[1]]; up_val <- up[idx[1]]
        }
        wide[[paste0(mod, "_estimate")]][wide$Parameters == param] <- est_val
        wide[[paste0(mod, "_low")]][wide$Parameters == param]      <- low_val
        wide[[paste0(mod, "_up")]][wide$Parameters == param]       <- up_val
      }
    }
    wide
  }

  format_cell <- function(est, low, up) {
    if (is.na(est)) {
      est_str <- ""; ci_str <- ""
    } else {
      est_str <- sprintf("\\textbf{%.3f}", as.numeric(est))
      ci_str  <- if (is.na(low) || is.na(up)) "" else
        sprintf("\\scriptsize(%.1f--%.1f)", as.numeric(low), as.numeric(up))
    }
    paste0("\\begin{tabular}[t]{@{}l@{}}", est_str,
           if (nzchar(ci_str)) " \\\\ " else "", ci_str, "\\end{tabular}")
  }

  table_list <- list()

  for (scen in scenarios) {
    if (verbose) message("Exporting combined SDM+GLM summary for scenario: ", scen)

    get_mod <- function(model, type) {
      nm <- paste0(scen, model, ".", type)
      if (!nm %in% names(model_list)) stop("Model not found: ", nm, " in scenario ", scen)
      model_list[[nm]]
    }
    models_SDM <- list(Fox = get_mod("F", "SDM"),
                       Schaefer = get_mod("S", "SDM"),
                       Pella = get_mod("P", "SDM"))
    models_GLM <- list(Fox = get_mod("F", "GLM"),
                       Schaefer = get_mod("S", "GLM"),
                       Pella = get_mod("P", "GLM"))

    df_SDM_long <- make_spict_summary_df1_nopipe(models_SDM, names(models_SDM))
    df_GLM_long <- make_spict_summary_df1_nopipe(models_GLM, names(models_GLM))
    wide_SDM <- pivot_to_wide(df_SDM_long)
    wide_GLM <- pivot_to_wide(df_GLM_long)

    idx_SDM <- match(param_order, wide_SDM$Parameters)
    idx_GLM <- match(param_order, wide_GLM$Parameters)
    if (any(is.na(idx_SDM)) || any(is.na(idx_GLM))) {
      stop("Some parameters are missing in SDM or GLM summary for scenario ", scen, ".")
    }
    wide_SDM <- wide_SDM[idx_SDM, , drop = FALSE]
    wide_GLM <- wide_GLM[idx_GLM, , drop = FALSE]

    Fox_SDM      <- mapply(format_cell, wide_SDM$Fox_estimate,      wide_SDM$Fox_low,      wide_SDM$Fox_up)
    Fox_GLM      <- mapply(format_cell, wide_GLM$Fox_estimate,      wide_GLM$Fox_low,      wide_GLM$Fox_up)
    Schaefer_SDM <- mapply(format_cell, wide_SDM$Schaefer_estimate, wide_SDM$Schaefer_low, wide_SDM$Schaefer_up)
    Schaefer_GLM <- mapply(format_cell, wide_GLM$Schaefer_estimate, wide_GLM$Schaefer_low, wide_GLM$Schaefer_up)
    Pella_SDM    <- mapply(format_cell, wide_SDM$Pella_estimate,    wide_SDM$Pella_low,    wide_SDM$Pella_up)
    Pella_GLM    <- mapply(format_cell, wide_GLM$Pella_estimate,    wide_GLM$Pella_low,    wide_GLM$Pella_up)

    ParametersOut <- kableExtra::cell_spec(param_order, format = "latex", bold = TRUE, escape = FALSE)
    df_out <- data.frame(
      Scenario      = rep(scen, length(param_order)),
      Parameters    = ParametersOut,
      Fox_SDM       = Fox_SDM,
      Fox_GLM       = Fox_GLM,
      Schaefer_SDM  = Schaefer_SDM,
      Schaefer_GLM  = Schaefer_GLM,
      Pella_SDM     = Pella_SDM,
      Pella_GLM     = Pella_GLM,
      stringsAsFactors = FALSE
    )

    # flags for spacing/rules
    plain_parameters <- param_order
    df_out$group <- NA_integer_
    for (i in seq_along(param_groups)) {
      df_out$group[plain_parameters %in% param_groups[[i]]] <- i
    }
    df_out$add_linespace <- plain_parameters %in% c(tail(param_groups[[1]], 1),
                                                    tail(param_groups[[2]], 1))
    df_out$add_hline <- plain_parameters %in% c(head(param_groups[[2]], 1),
                                                head(param_groups[[3]], 1))
    df_out$add_hline_row <- plain_parameters %in% c("$\\sigma_C$", "$\\mathrm{MSY}$")

    this_caption <- if (!is.null(caption)) caption else
      paste0("\\captionsetup{width=\\textwidth}Parameter estimates and 95\\% confidence intervals for Fox, Schaefer, and Pella-Tomlinson models (SDM and GLM) for scenario ", scen, ".")

    out_file <- file.path(output_dir, paste0("table_scenario", scen, "_SDM_GLM.tex"))

    # IMPORTANT: escape backslashes in R string so LaTeX sees a single backslash
    kbl_out <- kableExtra::kbl(
      df_out[, 1:8],
      format   = "latex",
      escape   = FALSE,
      linesep  = "\\\\addlinespace",
      booktabs = TRUE,
      align    = "llllllll",
      col.names = c(
        "\\textbf{Scenario}", "\\textbf{Parameter}",
        "SDM", "GLM", "SDM", "GLM", "SDM", "GLM"
      ),
      caption  = this_caption
    )
    kbl_out <- kableExtra::add_header_above(
      kbl_out,
      header = c(" " = 2, "Fox" = 2, "Schaefer" = 2, "Pella-Tomlinson" = 2),
      escape = FALSE
    )

    for (i in seq_len(nrow(df_out))) {
      extra_code <- ""
      if (isTRUE(df_out$add_hline_row[i])) extra_code <- paste0(extra_code, "\\hline\n")
      if (isTRUE(df_out$add_linespace[i])) extra_code <- paste0(extra_code, "\\addlinespace\n")
      kbl_out <- kableExtra::row_spec(kbl_out, row = i, extra_latex_after = extra_code)
    }

    kbl_out <- kableExtra::column_spec(kbl_out, 3:8, width = "3cm", latex_column_spec = "p{3cm}")
    kbl_out <- kableExtra::kable_styling(
      kbl_out, latex_options = "hold_position", full_width = FALSE, font_size = 8, position = "center"
    )
    kbl_out <- kableExtra::footnote(
      kbl_out,
      general = "‘Fixed’ indicates that the parameter was not estimated but assigned a fixed value.",
      general_title = "Note: ",
      footnote_as_chunk = TRUE,
      escape = FALSE
    )

    writeLines(as.character(kbl_out), out_file)
    if (verbose) message("Written: ", out_file)
    table_list[[scen]] <- kbl_out
  }

  invisible(table_list)
}
