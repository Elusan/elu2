#' Export one LaTeX table comparing parameter values across ALL scenarios
#'
#' @param model_list Named list with elements like "S1F.SDM", "S1F.GLM", "S1S.SDM", ...
#' @param output_dir Directory for the .tex file
#' @param scenarios Character vector (auto-detected from names(model_list) if NULL)
#' @param block_def Named list: names = section headers (e.g., "alpha1"),
#'        values = character vector of parameter labels (exactly as they appear
#'        in your summary `Parameters` column, e.g., "$r$", "$\\alpha_1$").
#'        Each parameter in a block becomes a sub-section with all scenarios listed.
#' @param caption Optional custom caption
#' @param verbose logical
#' @return (invisible) the kableExtra object
#' @export
export_spict_compare_across_scenarios_SDM_GLM_007 <- function(
    model_list,
    output_dir = ".",
    scenarios = NULL,  # CHANGED: auto-detect by default
    block_def = NULL,
    caption = NULL,
    verbose = TRUE
) {
  if (!dir.exists(output_dir)) stop("Output directory does not exist!")

  ## ---- Reuse helpers from your existing function ----
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
          est_val <- est[idx[1]]
          low_val <- low[idx[1]]
          up_val  <- up[idx[1]]
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
      est_str <- ""
      ci_str  <- ""
    } else {
      est_str <- sprintf("\\textbf{%.3f}", as.numeric(est))
      ci_str  <- sprintf("\\scriptsize(%.1f--%.1f)", as.numeric(low), as.numeric(up))
    }
    paste0("\\begin{tabular}[t]{@{}l@{}}", est_str, " \\\\ ", ci_str, "\\end{tabular}")
  }

  ## ---- Default blocks if user does not pass a custom one ----
  ## (You can replace by your alpha1/alpha2/beta/r set as needed.)
  if (is.null(block_def)) {
    block_def <- list(
      "Process & obs. variances" = c("$\\sigma_B$", "$\\sigma_F$", "$\\sigma_{I_{1}}$", "$\\sigma_{I_{2}}$", "$\\sigma_C$"),
      "MSY quantities"           = c("$B_{\\mathrm{MSY}}$", "$F_{\\mathrm{MSY}}$", "$\\mathrm{MSY}$"),
      "Status (2022)"            = c("$B_{2022}$", "$F_{2022}$", "$B_{2022}/B_{\\mathrm{MSY}}$", "$F_{2022}/F_{\\mathrm{MSY}}$"),
      "Dynamics"                  = c("$r$", "$m$", "$K$")
    )
  }

  ## CHANGED: Auto-detect scenarios from names(model_list) if not provided
  if (is.null(scenarios) || length(scenarios) == 0) {
    nm <- names(model_list)
    if (is.null(nm)) stop("`model_list` must be a named list with elements like 'S1F.SDM'.")
    scen_raw  <- gsub("^((S\\d+)).*$", "\\1", nm)
    scen_keep <- scen_raw[grepl("^S\\d+$", scen_raw)]
    scenarios <- unique(scen_keep)
    scenarios <- scenarios[order(as.integer(sub("^S(\\d+)$", "\\1", scenarios)))]
    if (length(scenarios) == 0) stop("Could not auto-detect scenarios from names(model_list).")
  }

  ## ---- Precompute wide summaries for every scenario (once) ----
  wlist_SDM <- list(); wlist_GLM <- list()
  get_mod <- function(scen, model, type) {
    nm <- paste0(scen, model, ".", type)
    if (!nm %in% names(model_list)) stop("Model not found: ", nm)
    model_list[[nm]]
  }

  for (scen in scenarios) {
    if (verbose) message("Collecting summaries for ", scen, " ...")
    models_SDM <- list(Fox = get_mod(scen, "F", "SDM"),
                       Schaefer = get_mod(scen, "S", "SDM"),
                       Pella = get_mod(scen, "P", "SDM"))
    models_GLM <- list(Fox = get_mod(scen, "F", "GLM"),
                       Schaefer = get_mod(scen, "S", "GLM"),
                       Pella = get_mod(scen, "P", "GLM"))

    df_SDM_long <- make_spict_summary_df1_nopipe(models_SDM, names(models_SDM))
    df_GLM_long <- make_spict_summary_df1_nopipe(models_GLM, names(models_GLM))

    wlist_SDM[[scen]] <- pivot_to_wide(df_SDM_long)
    wlist_GLM[[scen]] <- pivot_to_wide(df_GLM_long)
  }

  ## ---- Build one big data.frame: rows = scenarios; sections = parameters ----
  rows <- list()
  section_index <- integer(0)
  current_row <- 0L

  block_names <- names(block_def)
  for (b in seq_along(block_def)) {
    params <- block_def[[b]]
    for (param in params) {
      # section header row count start
      start_row <- current_row + 1L
      for (scen in scenarios) {
        wS <- wlist_SDM[[scen]]; wG <- wlist_GLM[[scen]]

        iS <- match(param, wS$Parameters); iG <- match(param, wG$Parameters)
        if (is.na(iS) || is.na(iG)) {
          stop("Parameter ", param, " missing in scenario ", scen, ".")
        }

        row <- data.frame(
          Section  = param,              # will be used to group_rows
          Scenario = scen,
          Fox_SDM       = format_cell(wS$Fox_estimate[iS],      wS$Fox_low[iS],      wS$Fox_up[iS]),
          Fox_GLM       = format_cell(wG$Fox_estimate[iG],      wG$Fox_low[iG],      wG$Fox_up[iG]),
          Schaefer_SDM  = format_cell(wS$Schaefer_estimate[iS], wS$Schaefer_low[iS], wS$Schaefer_up[iS]),
          Schaefer_GLM  = format_cell(wG$Schaefer_estimate[iG], wG$Schaefer_low[iG], wG$Schaefer_up[iG]),
          Pella_SDM     = format_cell(wS$Pella_estimate[iS],    wS$Pella_low[iS],    wS$Pella_up[iS]),
          Pella_GLM     = format_cell(wG$Pella_estimate[iG],    wG$Pella_low[iG],    wG$Pella_up[iG]),
          stringsAsFactors = FALSE
        )
        rows[[length(rows) + 1L]] <- row
        current_row <- current_row + 1L
      }
      # section header row count end
      section_index <- c(section_index, current_row - length(scenarios) + 1L, current_row)
    }
  }

  big <- do.call(rbind, rows)

  ## ---- Build LaTeX table ----
  if (is.null(caption)) {
    caption <- "\\captionsetup{width=\\textwidth}Comparison of parameter estimates (SDM & GLM; Fox, Schaefer, Pella) across scenarios. Cells show bold estimate with 95\\% CI on a new line."
  }

  out_file <- file.path(output_dir, "table_compare_all_scenarios_SDM_GLM.tex")

  kbl_out <- kableExtra::kbl(
    big[, c("Scenario","Fox_SDM","Fox_GLM","Schaefer_SDM","Schaefer_GLM","Pella_SDM","Pella_GLM")],
    format   = "latex",
    escape   = FALSE,
    booktabs = TRUE,
    align    = "lllllll",
    col.names = c("\\textbf{Scenario}",
                  "SDM","GLM","SDM","GLM","SDM","GLM"),
    caption  = caption
  )
  kbl_out <- kableExtra::add_header_above(
    kbl_out,
    header = c(" " = 1, "Fox" = 2, "Schaefer" = 2, "Pella-Tomlinson" = 2),
    escape = FALSE
  )

  ## group_rows for each parameter section, using the index we tracked
  to_readable <- function(x) {
    # customize if you want pretty names like "alpha1" etc.
    x
  }

  idx_ptr <- 1L
  for (b in seq_along(block_def)) {
    params <- block_def[[b]]
    for (param in params) {
      start <- section_index[idx_ptr]
      end   <- section_index[idx_ptr + 1L]
      idx_ptr <- idx_ptr + 2L
      kbl_out <- kableExtra::group_rows(kbl_out, label = to_readable(param),
                                        start_row = start, end_row = end,
                                        escape = FALSE, bold = TRUE,
                                        latex_gap_space = "0.25em")
    }
    # optional: thin rule between blocks
    kbl_out <- kableExtra::row_spec(kbl_out, row = 0, extra_latex_after = "\\\addlinespace")
  }

  kbl_out <- kableExtra::kable_styling(
    kbl_out, latex_options = "hold_position", full_width = FALSE, font_size = 8, position = "center"
  )
  kbl_out <- kableExtra::column_spec(kbl_out, 2:7, width = "3cm", latex_column_spec = "p{3cm}")
  kbl_out <- kableExtra::footnote(
    kbl_out,
    general = "‘Fixed’ indicates the parameter was not estimated and was assigned a fixed value.",
    general_title = "Note: ",
    footnote_as_chunk = TRUE,
    escape = FALSE
  )

  writeLines(as.character(kbl_out), out_file)
  if (verbose) message("Written: ", out_file)
  invisible(kbl_out)
}
