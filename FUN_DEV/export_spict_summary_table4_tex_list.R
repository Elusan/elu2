#' Export All SPiCT Summary Tables for All Scenarios and Types in a Flat Named List
#'
#' @description
#' Given a flat named list (e.g., S1P.SDM, S1S.GLM, ...), auto-group by scenario and type,
#' then export a LaTeX summary table for each using `make_spict_summary_df1_nopipe()`
#' and `export_spict_summary_table4_tex()`.
#'
#' @param model_list Flat named list of SPiCT models (see above).
#' @param output_dir Directory to save .tex files.
#' @param verbose Print progress.
#' @return Invisibly, a list of LaTeX table objects (named by scenario/type).
#' @export
export_spict_summary_table4_tex_list <- function(
    model_list,
    output_dir = ".",
    verbose = TRUE
) {
  stopifnot(dir.exists(output_dir))
  stopifnot(is.list(model_list))

  # Parse names to extract scenario, model, type
  info <- do.call(rbind, lapply(names(model_list), function(nm) {
    # Extract scenario (e.g., S1), model (P/S/F), type (SDM/GLM)
    m <- regexec("^([A-Za-z0-9]+)([PSF])\\.([A-Za-z0-9]+)$", nm)
    parts <- regmatches(nm, m)[[1]]
    if (length(parts) == 4) {
      data.frame(
        full_name = nm,
        scenario = parts[2],  # S1
        model = parts[3],     # P/S/F
        type = parts[4],      # SDM/GLM
        stringsAsFactors = FALSE
      )
    } else {
      stop(sprintf("Name '%s' does not match expected pattern 'S1P.SDM'", nm))
    }
  }))

  # Find all unique scenario/type combos
  combos <- unique(info[, c("scenario", "type")])

  latex_outputs <- list()
  for (i in seq_len(nrow(combos))) {
    scen <- combos$scenario[i]
    typ  <- combos$type[i]
    if (verbose) message(sprintf("Exporting summary for scenario: %s, type: %s", scen, typ))

    # Build named model list in correct order for summary table
    make_model <- function(letter, label) {
      row <- info$scenario == scen & info$type == typ & info$model == letter
      if (sum(row) != 1) stop(sprintf("Missing model %s for scenario %s, type %s", letter, scen, typ))
      m <- model_list[[info$full_name[row]]]
      nm <- paste0(scen, "_", label)
      setNames(list(m), nm)
    }
    models <- c(
      make_model("F", "Fox"),
      make_model("S", "Schaefer"),
      make_model("P", "Pella")
    )
    model_names <- names(models)

    # Create summary data frame and export table
    df_sum <- make_spict_summary_df1_nopipe(models, model_names)
    scen_name <- paste0(scen, "_", typ)
    latex_outputs[[scen_name]] <- export_spict_summary_table4_tex(
      df_sum,
      scenario_name = scen_name,
      output_dir = output_dir
    )
  }

  invisible(latex_outputs)
}
