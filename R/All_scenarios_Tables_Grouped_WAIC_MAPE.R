#' Export Combined SDM+GLM SPiCT Scenario Tables with Grouped Separation
#'
#' Automatically summarizes all models, creates wide summary tables, and exports
#' a fully formatted LaTeX table for each scenario, with left-aligned columns and
#' group separation (\addlinespace and \hline) after desired parameter groups.
#' Adds WAIC and pooled MAPE (%) rows per scenario.
#'
#' @param model_list Named flat list of SPiCT fits with names like "S1F.SDM".
#' @param output_dir Directory for .tex files (must exist).
#' @param scenarios Character vector of scenario names (default auto-detected).
#' @param param_groups Optional parameter grouping (same behavior as before).
#' @param verbose logical
#' @param caption Optional caption override.
#' @param include_waic logical; append WAIC row (default TRUE).
#' @param include_mape logical; append pooled MAPE row (default TRUE).
#' @param waic_digits integer; decimals for WAIC (default 1).
#' @param mape_digits integer; decimals for MAPE (default 1).
#' @param cand_pred_names optional; forwarded to mape_by_scenario() / mape_for_fit_elu2() to help find predictions.
#' @return Invisibly, a named list of kableExtra objects.
#' @export
All_scenarios_Tables_Grouped_WAIC_MAPE <- function(
    model_list,
    output_dir = ".",
    scenarios = NULL,  # auto-detect by default
    param_groups = NULL,
    verbose = TRUE,
    caption = NULL,
    include_waic = TRUE,
    include_mape = TRUE,
    waic_digits = 1,
    mape_digits = 1,
    cand_pred_names = NULL
) {
  if (!dir.exists(output_dir)) stop("Output directory does not exist!")

  # -------- helpers -----------------------------------------------------------
  family_from_code <- function(code) {
    if (identical(code, "F")) "Fox" else if (identical(code, "S")) "Schaefer" else if (identical(code, "P")) "Pella" else NA_character_
  }
  parse_name <- function(nm) {
    base <- sub("\\..*$", "", nm)               # "S1P"
    scen <- sub("^((S\\d+)).*$", "\\1", base)   # "S1"
    type <- sub("^.*\\.", "", nm)               # "SDM"/"GLM"
    code <- sub("^S\\d+", "", base)             # "P" / "S" / "F"
    fam  <- family_from_code(code)
    list(scenario = scen, type = type, family = fam, code = code)
  }
  extract_waic_single <- function(fit) {
    if (!is.null(fit$opt$convergence) && is.finite(fit$opt$convergence) && fit$opt$convergence != 0) return(NA_real_)
    v <- NA_real_
    if (!is.null(fit$report) && !is.null(fit$report$WAIC_plug_total)) {
      v <- fit$report$WAIC_plug_total
    } else if (is.function(fit$report)) {
      rr <- fit$report()
      if (!is.null(rr$WAIC_plug_total)) v <- rr$WAIC_plug_total
    }
    if (!is.finite(v)) NA_real_ else as.numeric(v)
  }
  # compact one-line cell (bold only; no CI line)
  format_single <- function(x, digits = 1) {
    if (!is.finite(x)) return("\\begin{tabular}[t]{@{}l@{}}\\end{tabular}")
    val <- formatC(as.numeric(x), format = "f", digits = digits)
    paste0("\\begin{tabular}[t]{@{}l@{}}\\textbf{", val, "}\\end{tabular}")
  }
  # original 2-line (estimate + CI) cell for parameters
  format_cell <- function(est, low, up) {
    if (is.na(est)) {
      est_str <- ""; ci_str <- ""
    } else {
      est_str <- sprintf("\\textbf{%.3f}", as.numeric(est))
      if (is.na(low) || is.na(up)) {
        ci_str <- ""
      } else {
        ci_str <- sprintf("\\scriptsize(%.1f--%.1f)", as.numeric(low), as.numeric(up))
      }
    }
    paste0("\\begin{tabular}[t]{@{}l@{}}", est_str, if (nzchar(ci_str)) " \\\\ " else "", ci_str, "\\end{tabular}")
  }

  # -------- default parameter groups/order -----------------------------------
  if (is.null(param_groups)) {
    group1 <- c("$r$", "$m$", "$K$",
                "$\\sigma_B$", "$\\sigma_F$",
                "$\\sigma_{I_{1}}$", "$\\sigma_{I_{2}}$", "$\\sigma_C$")
    group2 <- c("$B_{\\mathrm{MSY}}$", "$F_{\\mathrm{MSY}}$", "$\\mathrm{MSY}$")
    group3 <- c("$B_{2022}$", "$F_{2022}$", "$B_{2022}/B_{\\mathrm{MSY}}$", "$F_{2022}/F_{\\mathrm{MSY}}$")
    param_order <- c(group1, group2, group3)
    param_groups <- list(group1, group2, group3)
  } else {
    param_order <- unlist(param_groups)
  }

  # -------- auto-detect scenarios (flat list) --------------------------------
  if (is.null(scenarios) || length(scenarios) == 0) {
    nm <- names(model_list)
    if (is.null(nm)) stop("`model_list` must be a named list like 'S1F.SDM'.")
    scen_raw <- gsub("^((S\\d+)).*$", "\\1", sub("\\..*$", "", nm))
    scen_keep <- scen_raw[grepl("^S\\d+$", scen_raw)]
    scenarios <- unique(scen_keep)
    scenarios <- scenarios[order(as.integer(sub("^S(\\d+)$", "\\1", scenarios)))]
    if (length(scenarios) == 0) stop("Could not auto-detect scenarios from names(model_list).")
  }

  # -------- precompute WAIC + MAPE across all fits ---------------------------
  # WAIC from flat list
  waic_df <- {
    rows <- vector("list", length(model_list)); k <- 0L
    for (nm in names(model_list)) {
      info <- parse_name(nm)
      fit  <- model_list[[nm]]
      k <- k + 1L
      rows[[k]] <- data.frame(
        scenario = info$scenario,
        type     = info$type,
        family   = info$family,
        waic     = extract_waic_single(fit),
        stringsAsFactors = FALSE
      )
    }
    do.call(rbind, rows)
  }

  # Pooled ("AllIndices") MAPE from your helper (works on flat lists too)
  mape_df <- mape_by_scenario(model_list, cand_pred_names = cand_pred_names)

  # --- Normalize keys & map family names (prevents blank MAPE due to mismatched keys) ---
  mape_df$scenario <- sub("^((S\\d+)).*$", "\\1", trimws(mape_df$scenario))
  mape_df$type     <- toupper(trimws(mape_df$type))
  mape_df$model    <- trimws(mape_df$model)
  mape_df$index    <- trimws(mape_df$index)

  # keep only the pooled row
  mape_df <- mape_df[mape_df$index == "AllIndices", , drop = FALSE]

  # map model code -> family
  model_to_family <- function(z) {
    if (z == "F") "Fox" else if (z == "S") "Schaefer" else if (z == "P") "Pella" else NA_character_
  }
  mape_df$family <- vapply(mape_df$model, model_to_family, FUN.VALUE = character(1))

  # -------- inner helpers: scenario lookups + robust MAPE fallback -----------
  # return named numeric vectors (Fox, Schaefer, Pella) for SDM/GLM
  waic_for <- function(scen, typ) {
    sub <- waic_df[waic_df$scenario == scen & waic_df$type == typ, , drop = FALSE]
    out <- setNames(rep(NA_real_, 3L), c("Fox", "Schaefer", "Pella"))
    if (nrow(sub)) {
      for (fam in intersect(names(out), sub$family)) {
        out[fam] <- sub$waic[sub$family == fam][1]
      }
    }
    out
  }

  # compute pooled MAPE directly from a fit if needed
  safe_pooled_mape <- function(fit, scen_lab) {
    out <- NA_real_
    tryCatch({
      df <- mape_for_fit_elu2(
        fit,
        model_name = NA_character_,
        scenario   = scen_lab,
        type       = NA_character_,
        cand_pred_names = cand_pred_names
      )
      v <- df$MAPE[df$index == "AllIndices"]
      if (length(v)) out <- as.numeric(v[1])
    }, error = function(e) { out <<- NA_real_ })
    if (is.finite(out)) out else NA_real_
  }

  # robust getter: try precomputed table first, else compute directly
  get_mape_vec <- function(scen, typ, models_list) {
    sub <- mape_df[mape_df$scenario == scen & mape_df$type == typ, , drop = FALSE]
    vec <- setNames(rep(NA_real_, 3L), c("Fox", "Schaefer", "Pella"))
    if (nrow(sub)) {
      fams_here <- intersect(names(vec), sub$family)
      for (fam in fams_here) {
        vec[fam] <- sub$MAPE[sub$family == fam][1]
      }
    }
    # fallback per family if still NA
    if (!is.finite(vec["Fox"]) && "Fox" %in% names(models_list)) {
      vec["Fox"] <- safe_pooled_mape(models_list$Fox, scen)
    }
    if (!is.finite(vec["Schaefer"]) && "Schaefer" %in% names(models_list)) {
      vec["Schaefer"] <- safe_pooled_mape(models_list$Schaefer, scen)
    }
    if (!is.finite(vec["Pella"]) && "Pella" %in% names(models_list)) {
      vec["Pella"] <- safe_pooled_mape(models_list$Pella, scen)
    }
    vec
  }

  # -------- core summarization (unchanged for parameter rows) ----------------
  pivot_to_wide <- function(df) {
    cols_needed <- c("Parameters",
                     "Fox_estimate", "Fox_low", "Fox_up",
                     "Schaefer_estimate", "Schaefer_low", "Schaefer_up",
                     "Pella_estimate", "Pella_low", "Pella_up")
    df_clean <- df[, cols_needed]
    keep <- !apply(df_clean, 1, function(row) all(is.na(row[-1])))
    df_clean <- df_clean[keep, ]
    unique_params <- unique(as.character(df_clean$Parameters))
    wide <- data.frame(Parameters = unique_params, stringsAsFactors = FALSE)
    for (param in unique_params) {
      rows <- df_clean[df_clean$Parameters == param, ]
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

  table_list <- list()

  for (scen in scenarios) {
    if (verbose) message("Exporting combined SDM+GLM summary for scenario: ", scen)

    # ---- Gather models for the scenario
    get_mod <- function(model, type) {
      nm <- paste0(scen, model, ".", type)
      if (!nm %in% names(model_list)) stop("Model not found: ", nm, " in scenario ", scen)
      model_list[[nm]]
    }
    models_SDM <- list(
      Fox      = get_mod("F", "SDM"),
      Schaefer = get_mod("S", "SDM"),
      Pella    = get_mod("P", "SDM")
    )
    models_GLM <- list(
      Fox      = get_mod("F", "GLM"),
      Schaefer = get_mod("S", "GLM"),
      Pella    = get_mod("P", "GLM")
    )

    # ---- Summaries to wide
    df_SDM_long <- make_spict_summary_df1_nopipe(models_SDM, names(models_SDM))
    df_GLM_long <- make_spict_summary_df1_nopipe(models_GLM, names(models_GLM))
    wide_SDM <- pivot_to_wide(df_SDM_long)
    wide_GLM <- pivot_to_wide(df_GLM_long)

    # ---- Order by parameter group
    idx_SDM <- match(param_order, wide_SDM$Parameters)
    idx_GLM <- match(param_order, wide_GLM$Parameters)
    if (any(is.na(idx_SDM)) || any(is.na(idx_GLM))) {
      stop("Some parameters missing in SDM or GLM summary for scenario ", scen, ".")
    }
    wide_SDM <- wide_SDM[idx_SDM, , drop = FALSE]
    wide_GLM <- wide_GLM[idx_GLM, , drop = FALSE]

    # ---- Build parameter block (cells with CI)
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

    # ---- group flags for parameter rows (original behavior)
    plain_parameters <- param_order
    df_out$group <- NA_integer_
    for (i in seq_along(param_groups)) {
      df_out$group[plain_parameters %in% param_groups[[i]]] <- i
    }
    df_out$add_linespace <- plain_parameters %in% c(tail(param_groups[[1]], 1), tail(param_groups[[2]], 1))
    df_out$add_hline     <- plain_parameters %in% c(head(param_groups[[2]], 1), head(param_groups[[3]], 1))
    df_out$add_hline_row <- plain_parameters %in% c("$\\sigma_C$", "$\\mathrm{MSY}$")

    # ---- Append WAIC + MAPE rows (no CI)
    nb_param_rows <- nrow(df_out)

    # lookup values
    waic_sdm  <- waic_for(scen, "SDM")
    waic_glm  <- waic_for(scen, "GLM")
    mape_sdm  <- get_mape_vec(scen, "SDM", models_SDM)
    mape_glm  <- get_mape_vec(scen, "GLM", models_GLM)

    new_rows <- list()
    if (isTRUE(include_waic)) {
      row_WAIC <- data.frame(
        Scenario     = scen,
        Parameters   = kableExtra::cell_spec("WAIC", format = "latex", bold = TRUE, escape = FALSE),
        Fox_SDM      = format_single(waic_sdm["Fox"],       digits = waic_digits),
        Fox_GLM      = format_single(waic_glm["Fox"],       digits = waic_digits),
        Schaefer_SDM = format_single(waic_sdm["Schaefer"],  digits = waic_digits),
        Schaefer_GLM = format_single(waic_glm["Schaefer"],  digits = waic_digits),
        Pella_SDM    = format_single(waic_sdm["Pella"],     digits = waic_digits),
        Pella_GLM    = format_single(waic_glm["Pella"],     digits = waic_digits),
        stringsAsFactors = FALSE
      )
      new_rows <- c(new_rows, list(row_WAIC))
    }
    if (isTRUE(include_mape)) {
      row_MAPE <- data.frame(
        Scenario     = scen,
        Parameters   = kableExtra::cell_spec("MAPE (\\%)", format = "latex", bold = TRUE, escape = FALSE),
        Fox_SDM      = format_single(mape_sdm["Fox"],       digits = mape_digits),
        Fox_GLM      = format_single(mape_glm["Fox"],       digits = mape_digits),
        Schaefer_SDM = format_single(mape_sdm["Schaefer"],  digits = mape_digits),
        Schaefer_GLM = format_single(mape_glm["Schaefer"],  digits = mape_digits),
        Pella_SDM    = format_single(mape_sdm["Pella"],     digits = mape_digits),
        Pella_GLM    = format_single(mape_glm["Pella"],     digits = mape_digits),
        stringsAsFactors = FALSE
      )
      new_rows <- c(new_rows, list(row_MAPE))
    }

    if (length(new_rows)) {
      # add a linespace BEFORE the WAIC/MAPE block by tagging the last parameter row
      df_out$add_linespace[nb_param_rows] <- TRUE

      # append rows and extend flag columns with FALSE
      add_block <- do.call(rbind, new_rows)
      add_block$group <- NA_integer_
      add_block$add_linespace <- FALSE
      add_block$add_hline     <- FALSE
      add_block$add_hline_row <- FALSE
      df_out <- rbind(df_out, add_block)
    }

    # ---- caption & LaTeX output
    this_caption <- if (!is.null(caption)) caption else
      paste0("\\captionsetup{width=\\textwidth}Parameter estimates and 95\\% confidence intervals for Fox, Schaefer, and Pella-Tomlinson models (SDM and GLM) for scenario ", scen, ".")

    out_file <- file.path(output_dir, paste0("table_scenario", scen, "_SDM_GLM.tex"))
    kbl_out <- kableExtra::kbl(
      df_out[, 1:8],
      format = "latex",
      escape = FALSE,
      linesep = "\\hline",
      booktabs = TRUE,
      align = "llllllll",
      col.names = c(
        "\\textbf{Scenario}", "\\textbf{Parameter}",
        "SDM", "GLM", "SDM", "GLM", "SDM", "GLM"
      ) #,
      #caption = this_caption
    )
    kbl_out <- kableExtra::add_header_above(
      kbl_out,
      header = c(" " = 2, "Fox" = 2, "Schaefer" = 2, "Pella-Tomlinson" = 2),
      escape = FALSE
    )

    # apply per-row extras
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
      #general = "‘Fixed’ indicates that the parameter was not estimated but assigned a fixed value.",
      #general_title = "Note: ",
      footnote_as_chunk = TRUE,
      escape = FALSE
    )

    writeLines(as.character(kbl_out), out_file)
    if (verbose) message("Written: ", out_file)
    table_list[[scen]] <- kbl_out
  }

  invisible(table_list)
}
