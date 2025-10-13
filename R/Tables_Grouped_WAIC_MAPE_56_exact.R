#' Export Combined SDM+GLM SPiCT Scenario Tables with Grouped Separation (+ WAIC & MAPE) — EXACT STYLE
#'
#' Produces LaTeX tables that exactly match your specified style:
#' row banding with \code{\\rowcolors}, bold \emph{Fox/Schaefer/Pella-Tomlinson} headers with
#' custom \code{\\cmidrule(l\{3pt\}r\{3pt\})}, a rule after every row, white section banners,
#' and compact two-line cells (bold estimate on line 1; \code{\\scriptsize} CI on line 2).
#' Appends WAIC and pooled MAPE rows per scenario. Writes one \code{table_scenario<SX>_SDM_GLM.tex}
#' per scenario and invisibly returns the LaTeX as a named list.
#'
#' @param model_list Named flat list of SPiCT fits with names like \code{"S1F.SDM"} and \code{"S1F.GLM"}.
#'   Names must follow the pattern \code{<scenario><family>.<type>} where \code{scenario="S\\d+"},
#'   \code{family in \{"F","S","P"\}} for Fox, Schaefer, Pella, and \code{type in \{"SDM","GLM"\}}.
#' @param output_dir Directory where \code{.tex} files are written. Must exist.
#' @param scenarios Optional character vector of scenario IDs (e.g., \code{c("S1","S2")}). If \code{NULL},
#'   scenarios are auto-detected from \code{names(model_list)}.
#' @param param_groups Optional list of three character vectors defining the row order and sectioning:
#'   group 1 (process & observation SDs), group 2 (stochastic reference points),
#'   group 3 (stock status). If \code{NULL}, sensible defaults are used.
#' @param verbose Logical; print progress messages while writing files. Default \code{TRUE}.
#' @param caption Optional caption string (kept for API symmetry; not inserted in the EXACT style output).
#' @param include_waic Logical; if \code{TRUE} append WAIC row. Default \code{TRUE}.
#' @param include_mape Logical; if \code{TRUE} append pooled MAPE row. Default \code{TRUE}.
#' @param waic_digits Integer; decimals for WAIC cells. Default \code{1}.
#' @param mape_digits Integer; decimals for MAPE cells. Default \code{1}.
#' @param cand_pred_names Optional character vector of candidate report element names to obtain
#'   index predictions for pooled MAPE when \code{fit$osar$Ihat} is unavailable.
#'
#' @details
#' This function writes raw LaTeX (no \pkg{kableExtra}) to guarantee byte-for-byte styling.
#' Ensure your LaTeX preamble includes \code{\\usepackage[table]{xcolor}} and \code{\\usepackage{booktabs}}.
#'
#' @return Invisibly, a named list of character scalars containing the LaTeX for each scenario.
#'   Files are written as \code{file.path(output_dir, paste0("table_scenario", scen, "_SDM_GLM.tex"))}.
#'
#' @examples
#' \dontrun{
#' out <- All_scenarios_Tables_Grouped_WAIC_MAPE_56_exact(model_list, output_dir = "TABLES")
#' }
#'
#' @export
Tables_Grouped_WAIC_MAPE_56_exact <- function(
    model_list,
    output_dir = ".",
    scenarios = NULL,
    param_groups = NULL,
    verbose = TRUE,
    caption = NULL,                # kept for API symmetry (unused in exact style)
    include_waic = TRUE,
    include_mape = TRUE,
    waic_digits = 1,
    mape_digits = 1,
    cand_pred_names = NULL
) {
  if (!dir.exists(output_dir)) stop("Output directory does not exist!")
  if (is.null(names(model_list))) stop("`model_list` must be a *named* list (e.g., 'S1F.SDM').")

  `%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

  # ---------- helpers copied/adapted from your v55 ----------
  family_from_code <- function(code) {
    if (is.null(code) || is.na(code)) return(NA_character_)
    z <- toupper(trimws(as.character(code)))
    if (nchar(z) > 1L) z <- substr(z, 1L, 1L)
    if (z == "F") "Fox" else if (z == "S") "Schaefer" else if (z == "P") "Pella" else NA_character_
  }
  parse_name <- function(nm) {
    base <- sub("\\..*$", "", nm)               # "S1F"
    scen <- sub("^((S\\d+)).*$", "\\1", base)   # "S1"
    type <- sub("^.*\\.", "", nm)               # "SDM"/"GLM"
    code <- substr(base, nchar(base), nchar(base))  # "F"/"S"/"P"
    list(scenario = scen, type = toupper(type), family = family_from_code(code))
  }
  extract_waic_single <- function(fit) {
    if (!is.null(fit$opt$convergence) && is.finite(fit$opt$convergence) && fit$opt$convergence != 0)
      return(NA_real_)
    v <- NA_real_
    if (!is.null(fit$report) && !is.null(fit$report$WAIC_plug_total)) {
      v <- fit$report$WAIC_plug_total
    } else if (is.function(fit$report)) {
      rr <- fit$report()
      if (!is.null(rr$WAIC_plug_total)) v <- rr$WAIC_plug_total
    }
    v <- suppressWarnings(as.numeric(v[1]))
    if (!is.finite(v)) NA_real_ else v
  }

  # two-line cell (estimate + CI) for parameters
  format_cell <- function(est, low, up) {
    if (is.na(est)) return("\\begin{tabular}[t]{@{}l@{}}\\end{tabular}")
    est_str <- sprintf("\\textbf{%.3f}", as.numeric(est))
    ci_str  <- if (is.na(low) || is.na(up)) "" else
      sprintf("\\scriptsize(%.1f--%.1f)", as.numeric(low), as.numeric(up))
    paste0("\\begin{tabular}[t]{@{}l@{}}", est_str,
           if (nzchar(ci_str)) " \\\\ " else "", ci_str, "\\end{tabular}")
  }
  # single-line bold numeric cell (for WAIC & MAPE)
  format_single <- function(x, digits = 1) {
    ax <- suppressWarnings(as.numeric(x))
    if (!is.finite(ax)) return("\\begin{tabular}[t]{@{}l@{}}\\end{tabular}")
    val <- formatC(ax, format = "f", digits = digits)
    paste0("\\begin{tabular}[t]{@{}l@{}}\\textbf{", val, "}\\end{tabular}")
  }

  # robust pooled MAPE (simple fallback)
  pooled_mape_for_fit_simple <- function(fit, cand_pred_names = NULL) {
    if (is.null(fit) || is.null(fit$inp) || is.null(fit$inp$obsI)) return(NA_real_)
    obs_all <- unlist(fit$inp$obsI, use.names = FALSE)
    if (!length(obs_all)) return(NA_real_)
    if (any(!is.finite(obs_all)) || any(obs_all == 0)) return(NA_real_)
    nI    <- length(fit$inp$obsI)
    nobsI <- if (!is.null(fit$inp$nobsI)) as.integer(fit$inp$nobsI) else vapply(fit$inp$obsI, length, 1L)
    totI  <- sum(nobsI)
    iq    <- NULL
    if ("iqin" %in% names(fit$inp)) {
      iq <- fit$inp$iqin
      if (is.list(iq)) iq <- unlist(iq, use.names = FALSE)
      iq <- as.integer(iq)
      if (length(iq) && min(iq, na.rm = TRUE) == 0L) iq <- iq + 1L
      if (!length(iq)) iq <- NULL
    }
    if (is.null(cand_pred_names)) {
      cand_pred_names <- c("Ihat_obs","Ipred_obs","Ip_obs","logIhat_obs","logIpred_obs","logIp_obs",
                           "Ihat","Ipred","Ip","logIhat","logIpred","logIp")
    }
    if (!is.null(fit$osar) && !is.null(fit$osar$Ihat)) {
      v <- fit$osar$Ihat; if (is.matrix(v) || is.data.frame(v)) v <- as.numeric(v)
      if (is.numeric(v)) {
        if (length(v) == totI) return(mean(abs(v - obs_all) / abs(obs_all)) * 100)
        if (!is.null(iq) && max(iq, na.rm = TRUE) <= length(v))
          return(mean(abs(v[iq] - obs_all) / abs(obs_all)) * 100)
      }
    }
    rep <- fit$report; if (is.null(rep)) return(NA_real_)
    nm <- names(rep); if (is.null(nm)) nm <- character()
    use_vec <- function(vec, is_log = FALSE) {
      if (is.matrix(vec) || is.data.frame(vec)) vec <- as.numeric(vec)
      if (!is.numeric(vec)) return(NULL)
      if (length(vec) == totI) return(if (is_log) exp(as.numeric(vec)) else as.numeric(vec))
      if (!is.null(iq) && max(iq, na.rm = TRUE) <= length(vec)) {
        vv <- vec[iq]; return(if (is_log) exp(as.numeric(vv)) else as.numeric(vv))
      }
      if (length(vec) == nI) {
        if (is_log) vec <- exp(as.numeric(vec)) else vec <- as.numeric(vec)
        out <- numeric(totI); pos <- 1L
        for (i in seq_len(nI)) {
          cnt <- nobsI[i]; if (cnt > 0L) { out[pos:(pos+cnt-1L)] <- vec[i]; pos <- pos + cnt }
        }
        return(out)
      }
      NULL
    }
    for (cn in cand_pred_names) {
      if (cn %in% nm) {
        pred_all <- use_vec(rep[[cn]], is_log = startsWith(cn, "log"))
        if (!is.null(pred_all) && length(pred_all) == totI && all(is.finite(pred_all)))
          return(mean(abs(pred_all - obs_all) / abs(obs_all)) * 100)
      }
    }
    likeI <- nm[grepl("^(log)?I(hat|pred|p)?(_obs)?$", nm, perl = TRUE)]
    likeI <- setdiff(likeI, nm[grepl("^(C|E)p", nm)])
    for (cn in likeI) {
      pred_all <- use_vec(rep[[cn]], is_log = startsWith(cn, "log"))
      if (!is.null(pred_all) && length(pred_all) == totI && all(is.finite(pred_all)))
        return(mean(abs(pred_all - obs_all) / abs(obs_all)) * 100)
    }
    NA_real_
  }

  # default groups and order if not provided
  if (is.null(param_groups)) {
    group1 <- c("$r$", "$m$", "$K$", "$\\sigma_B$", "$\\sigma_F$", "$\\sigma_{I_{1}}$", "$\\sigma_{I_{2}}$", "$\\sigma_C$")
    group2 <- c("$B_{\\mathrm{MSY}}$", "$F_{\\mathrm{MSY}}$", "$\\mathrm{MSY}$")
    group3 <- c("$B_{2022}$", "$F_{2022}$", "$B_{2022}/B_{\\mathrm{MSY}}$", "$F_{2022}/F_{\\mathrm{MSY}}$")
    param_order  <- c(group1, group2, group3)
    param_groups <- list(group1, group2, group3)
  } else {
    param_order <- unlist(param_groups)
  }

  # scenarios autodetect
  if (is.null(scenarios) || !length(scenarios)) {
    scen_raw  <- gsub("^((S\\d+)).*$", "\\1", sub("\\..*$", "", names(model_list)))
    scen_keep <- scen_raw[grepl("^S\\d+$", scen_raw)]
    scenarios <- unique(scen_keep)
    scenarios <- scenarios[order(as.integer(sub("^S(\\d+)$", "\\1", scenarios)))]
    if (!length(scenarios)) stop("Could not auto-detect scenarios from names(model_list).")
  }

  # WAIC table over flat list
  waic_df <- {
    rows <- vector("list", length(model_list)); k <- 0L
    for (nm in names(model_list)) {
      info <- parse_name(nm); fit <- model_list[[nm]]; k <- k + 1L
      rows[[k]] <- data.frame(
        scenario = info$scenario, type = info$type, family = info$family,
        waic = extract_waic_single(fit), stringsAsFactors = FALSE
      )
    }
    do.call(rbind, rows)
  }
  waic_for <- function(scen, typ) {
    sub <- waic_df[waic_df$scenario == scen & waic_df$type == typ, , drop = FALSE]
    out <- setNames(rep(NA_real_, 3L), c("Fox", "Schaefer", "Pella"))
    if (nrow(sub)) for (fam in intersect(names(out), sub$family)) out[fam] <- sub$waic[sub$family == fam][1]
    out
  }
  safe_pooled_mape <- function(fit) {
    out <- tryCatch(pooled_mape_for_fit_simple(fit, cand_pred_names = cand_pred_names), error = function(e) NA_real_)
    if (is.finite(out)) out else NA_real_
  }

  # collect outputs
  out_text <- list()

  for (scen in scenarios) {
    if (verbose) message("Exporting combined SDM+GLM summary for scenario: ", scen)

    get_mod <- function(model, type) {
      nm <- paste0(scen, model, ".", type)
      if (!nm %in% names(model_list)) stop("Model not found: ", nm, " in scenario ", scen)
      model_list[[nm]]
    }
    models_SDM <- list(Fox = get_mod("F", "SDM"), Schaefer = get_mod("S", "SDM"), Pella = get_mod("P", "SDM"))
    models_GLM <- list(Fox = get_mod("F", "GLM"), Schaefer = get_mod("S", "GLM"), Pella = get_mod("P", "GLM"))

    df_SDM_long <- make_spict_summary_df1_nopipe(models_SDM, names(models_SDM))
    df_GLM_long <- make_spict_summary_df1_nopipe(models_GLM, names(models_GLM))

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
    wide_SDM <- pivot_to_wide(df_SDM_long)
    wide_GLM <- pivot_to_wide(df_GLM_long)

    # enforce order
    idx_SDM <- match(param_order, wide_SDM$Parameters)
    idx_GLM <- match(param_order, wide_GLM$Parameters)
    if (any(is.na(idx_SDM)) || any(is.na(idx_GLM))) {
      stop("Some parameters are missing in SDM or GLM summary for scenario ", scen, ".")
    }
    wide_SDM <- wide_SDM[idx_SDM, , drop = FALSE]
    wide_GLM <- wide_GLM[idx_GLM, , drop = FALSE]

    # assemble rows for each parameter
    build_row <- function(param_label, row_idx) {
      fox_sdm      <- format_cell(wide_SDM$Fox_estimate[row_idx],      wide_SDM$Fox_low[row_idx],      wide_SDM$Fox_up[row_idx])
      fox_glm      <- format_cell(wide_GLM$Fox_estimate[row_idx],      wide_GLM$Fox_low[row_idx],      wide_GLM$Fox_up[row_idx])
      sch_sdm      <- format_cell(wide_SDM$Schaefer_estimate[row_idx], wide_SDM$Schaefer_low[row_idx], wide_SDM$Schaefer_up[row_idx])
      sch_glm      <- format_cell(wide_GLM$Schaefer_estimate[row_idx], wide_GLM$Schaefer_low[row_idx], wide_GLM$Schaefer_up[row_idx])
      pel_sdm      <- format_cell(wide_SDM$Pella_estimate[row_idx],    wide_SDM$Pella_low[row_idx],    wide_SDM$Pella_up[row_idx])
      pel_glm      <- format_cell(wide_GLM$Pella_estimate[row_idx],    wide_GLM$Pella_low[row_idx],    wide_GLM$Pella_up[row_idx])
      paste0(
        scen, " & \\textbf{", param_label, "} & ",
        fox_sdm, " & ", fox_glm, " & ",
        sch_sdm, " & ", sch_glm, " & ",
        pel_sdm, " & ", pel_glm, "\\\\\n\\hline\n"
      )
    }

    # section banners
    banner_row <- function(title) {
      paste0("\\addlinespace\n\\rowcolor{White} \\multicolumn{8}{c}{\\large{", title, "}}\\\\\n\\addlinespace\n\\hline\n")
    }

    # build the full body in the exact grouping layout
    n_param <- length(param_order)
    body_txt <- character(0)

    # Group 1
    for (i in seq_along(param_groups[[1]])) {
      body_txt <- c(body_txt, build_row(param_groups[[1]][i], which(param_order == param_groups[[1]][i])))
    }
    # Banner: Stochastic reference points
    body_txt <- c(body_txt, banner_row("Stochastic reference points"))
    # Group 2
    for (i in seq_along(param_groups[[2]])) {
      body_txt <- c(body_txt, build_row(param_groups[[2]][i], which(param_order == param_groups[[2]][i])))
    }
    # Banner: Stock status with 95% CI
    body_txt <- c(body_txt, banner_row("Stock status with 95\\% CI"))
    # Group 3
    for (i in seq_along(param_groups[[3]])) {
      body_txt <- c(body_txt, build_row(param_groups[[3]][i], which(param_order == param_groups[[3]][i])))
    }

    # WAIC + MAPE block
    waic_sdm <- waic_for(scen, "SDM"); waic_glm <- waic_for(scen, "GLM")
    mape_sdm <- c(Fox = safe_pooled_mape(models_SDM$Fox),
                  Schaefer = safe_pooled_mape(models_SDM$Schaefer),
                  Pella = safe_pooled_mape(models_SDM$Pella))
    mape_glm <- c(Fox = safe_pooled_mape(models_GLM$Fox),
                  Schaefer = safe_pooled_mape(models_GLM$Schaefer),
                  Pella = safe_pooled_mape(models_GLM$Pella))

    waic_row <- paste0(
      scen, " & \\textbf{WAIC} & ",
      format_single(waic_sdm["Fox"], waic_digits), " & ",
      format_single(waic_glm["Fox"], waic_digits), " & ",
      format_single(waic_sdm["Schaefer"], waic_digits), " & ",
      format_single(waic_glm["Schaefer"], waic_digits), " & ",
      format_single(waic_sdm["Pella"], waic_digits), " & ",
      format_single(waic_glm["Pella"], waic_digits),
      "\\\\\n\\hline\n"
    )
    mape_row <- paste0(
      scen, " & \\textbf{MAPE (\\%)} & ",
      format_single(mape_sdm["Fox"], mape_digits), " & ",
      format_single(mape_glm["Fox"], mape_digits), " & ",
      format_single(mape_sdm["Schaefer"], mape_digits), " & ",
      format_single(mape_glm["Schaefer"], mape_digits), " & ",
      format_single(mape_sdm["Pella"], mape_digits), " & ",
      format_single(mape_glm["Pella"], mape_digits),
      "\\\\\n\\hline\n"
    )

    body_txt <- c(body_txt,
                  banner_row("WAIC and MAPE"),
                  if (isTRUE(include_waic)) waic_row else NULL,
                  if (isTRUE(include_mape)) mape_row else NULL)

    # --------- stitch full EXACT LaTeX -----------
    header_tex <- paste0(
      "\\rowcolors{1}{blue!20}{blue!10}\n",
      "\\renewcommand{\\arraystretch}{1.15}\n",
      "\\setlength{\\extrarowheight}{1pt}\n",
      "\\setlength{\\aboverulesep}{0.6ex}\n",
      "\\setlength{\\belowrulesep}{0.6ex}\n",
      "\\begin{table}[!h]\n",
      "\\centering\\begingroup\\fontsize{8}{10}\\selectfont\n\n",
      "\\begin{tabular}[t]{llp{3cm}p{3cm}p{3cm}p{3cm}p{3cm}p{3cm}}\n",
      "\\hline\n",
      "\\multicolumn{2}{c}{ } & \\multicolumn{2}{c}{\\bf \\large{Fox}} & \\multicolumn{2}{c}{\\bf \\large{Schaefer}} & \\multicolumn{2}{c}{\\bf \\large{Pella-Tomlinson}} \\\\\n",
      "\\cmidrule(l{3pt}r{3pt}){3-4} \\cmidrule(l{3pt}r{3pt}){5-6} \\cmidrule(l{3pt}r{3pt}){7-8}\n",
      "\\textbf{Scenario} & \\textbf{Parameter} & SDM & GLM & SDM & GLM & SDM & GLM\\\\\n",
      "\\hline\n"
    )

    footer_tex <- paste0(
      "\\end{tabular}\n",
      "\\endgroup{}\n",
      "\\end{table}\n"
    )

    latex_full <- paste0(header_tex, paste0(body_txt, collapse = ""), footer_tex)

    # write file
    out_file <- file.path(output_dir, paste0("table_scenario", scen, "_SDM_GLM.tex"))
    writeLines(latex_full, out_file)
    if (verbose) message("Written: ", out_file)

    out_text[[scen]] <- latex_full
  }

  invisible(out_text)
}
