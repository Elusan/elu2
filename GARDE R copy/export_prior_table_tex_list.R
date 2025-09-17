#' Export Prior Tables for All Scenarios and Types in a Flat Named List
#'
#' @description
#' Given a flat named list of models (e.g., S1P.SDM, S1S.GLM, ...), auto-group
#' by scenario (e.g., S1) and type (e.g., SDM), and export a LaTeX table for each.
#'
#' @param model_list Flat named list of models (see above).
#' @param output_dir Directory to save .tex files.
#' @param verbose Print progress.
#' @return Invisibly, a list of LaTeX table strings (named by scenario/type).
#' @export
export_prior_table_tex_list <- function(
    model_list,
    output_dir = ".",
    verbose = TRUE
) {
  stopifnot(dir.exists(output_dir))
  stopifnot(is.list(model_list))

  # Parse all names: extract scenario, model, type
  # Name format assumed: "S1P.SDM" or "S2F.GLM", etc
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

  # Find all unique scenario/type combinations
  combos <- unique(info[, c("scenario", "type")])

  latex_outputs <- list()
  for (i in seq_len(nrow(combos))) {
    scen <- combos$scenario[i]
    typ  <- combos$type[i]
    if (verbose) message(sprintf("Exporting scenario: %s, type: %s", scen, typ))

    # Find the model names for this scenario/type
    get_model <- function(letter) {
      row <- info$scenario == scen & info$type == typ & info$model == letter
      if (sum(row) != 1) stop(sprintf("Missing model %s for scenario %s, type %s", letter, scen, typ))
      model_list[[info$full_name[row]]]
    }
    S_Pella    <- get_model("P")
    S_Schaefer <- get_model("S")
    S_Fox      <- get_model("F")

    scen_name <- paste0(scen, "_", typ)
    latex_outputs[[scen_name]] <- export_prior_table_tex(
      S_Fox = S_Fox,
      S_Schaefer = S_Schaefer,
      S_Pella = S_Pella,
      scenario_name = scen_name,
      output_dir = output_dir
    )
  }

  invisible(latex_outputs)
}
