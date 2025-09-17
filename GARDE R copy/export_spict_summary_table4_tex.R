#' Export SPiCT Summary Table to LaTeX (Formatted for Beamer)
#'
#' @description This function takes a summary data frame of SPiCT model estimates
#' (as produced by `make_spict_summary_df1_nopipe`) and creates a LaTeX-formatted
#' table with appropriate model headers and footnotes, saving it to disk.
#'
#' @param df_summary A data frame of summary estimates for Fox, Schaefer, and Pella models.
#'                   Must contain columns like: Scenario, Parameters, Fox_estimate, Fox_low, etc.
#' @param scenario_name A character string used to name the output file (e.g., `"1"` or `"S1"`).
#' @param output_dir A directory path to save the LaTeX `.tex` file (default: current directory).
#'
#' @return Invisibly returns the LaTeX table object as a string.
#'
#' @importFrom kableExtra kbl add_header_above kable_styling footnote cell_spec
#' @export
export_spict_summary_table4_tex <- function(df_summary, scenario_name = "1", output_dir = ".") {
  require(kableExtra)

  # Define all groups
  group1 <- c("$\\alpha_{1}$", "$\\alpha_{2}$", "$\\beta$", "$r$", "$m$", "$K$",
              "$q_{1}$", "$q_{2}$", "$n$", "$\\sigma_B$", "$\\sigma_F$",
              "$\\sigma_{I_{1}}$", "$\\sigma_{I_{2}}$", "$\\sigma_C$")
  group2 <- c("$B_{\\mathrm{MSY}}$", "$F_{\\mathrm{MSY}}$", "$\\mathrm{MSY}$")
  group3 <- c("$B_{2022}$", "$F_{2022}$", "$B_{2022}/B_{\\mathrm{MSY}}$", "$F_{2022}/F_{\\mathrm{MSY}}$")
  group4 <- c("$B_{2024}$", "$F_{2024}$", "$B_{2024}/B_{\\mathrm{MSY}}$", "$F_{2024}/F_{\\mathrm{MSY}}$")

  # Strip \\textbf{} from Parameters to compare
  clean_param <- function(x) gsub("^\\\\textbf\\{(.+)\\}$", "\\1", x)

  # Combine into a full sequence
  group_sequence <- list(group1, group2, group3, group4)
  plain_parameters <- vapply(df_summary$Parameters, clean_param, character(1))

  # Tag each parameter row with a group number
  df_summary$group <- NA_integer_
  for (i in seq_along(group_sequence)) {
    df_summary$group[plain_parameters %in% group_sequence[[i]]] <- i
  }

  # Add \\hline to every row and \\addlinespace only *after group transition*
  df_summary$add_linespace <- c(
    diff(df_summary$group) == 1,  # TRUE when group changes
    FALSE  # last row no space
  )

  # Define grouped parameter labels (must match LaTeX-formatted names from param_labels)
  # (kept as comments from your original)
  # group1 <- ...
  # group2 <- ...
  # group3 <- ...
  # group4 <- ...

  df_summary$Parameters <- kableExtra::cell_spec(df_summary$Parameters, format = "latex", bold = TRUE, escape = FALSE)

  n_label <- kableExtra::cell_spec("$n$", format = "latex", bold = TRUE, escape = FALSE)

  # Replace values for "n" rows
  for (col in c("Fox", "Schaefer")) {
    idx <- which(df_summary$Parameters == n_label)
    if (length(idx)) {
      df_summary[[paste0(col, "_estimate")]][idx] <- "Fixed"
      df_summary[[paste0(col, "_low")]][idx] <- ""
      df_summary[[paste0(col, "_up")]][idx] <- ""
    }
  }

  # Identify rows to add spacing (only at the end of groups 1–3)
  df_summary$add_linespace <- plain_parameters %in% c(
    tail(group1, 1),  # "$\\sigma_C$"
    tail(group2, 1),  # "$\\mathrm{MSY}$"
    tail(group3, 1)   # "$F_{2022}/F_{\\mathrm{MSY}}$"
  )

  out_file <- file.path(output_dir, paste0("table_scenario", scenario_name, "_estimates.tex"))

  df_print <- df_summary[, c("Scenario", "Parameters",
                             "Fox_estimate", "Fox_low", "Fox_up",
                             "Schaefer_estimate", "Schaefer_low", "Schaefer_up",
                             "Pella_estimate", "Pella_low", "Pella_up")]

  kbl_out <- kableExtra::kbl(
    df_print,
    format = "latex",
    escape = FALSE,
    linesep = "\\addlinespace",
    booktabs = TRUE,
    col.names = c(
      "\\textbf{Scenario}", "\\textbf{Parameters}",
      "\\textbf{estimate}", "\\textbf{low}", "\\textbf{up}",
      "\\textbf{estimate}", "\\textbf{low}", "\\textbf{up}",
      "\\textbf{estimate}", "\\textbf{low}", "\\textbf{up}"
    )
  )

  kbl_out <- kableExtra::add_header_above(
    kbl_out,
    header = c(
      " " = 2,
      "Fox" = 3,
      "Schaefer" = 3,
      "Pella-Tomlinson" = 3
    ),
    escape = FALSE
  )

  kbl_out <- kableExtra::kable_styling(
    kbl_out,
    latex_options = c("hold_position"),
    full_width = FALSE,
    position = "center",
    font_size = 8
  )

  kbl_out <- kableExtra::footnote(
    kbl_out,
    general = "‘Fixed’ indicates that the parameter was not estimated but assigned a fixed value.",
    general_title = "Note: ",
    footnote_as_chunk = TRUE,
    escape = FALSE
  )

  for (i in seq_len(nrow(df_summary))) {
    extra_code <- ""

    # Add \\addlinespace only *between groups*
    if (df_summary$add_linespace[i]) {
      extra_code <- paste0(extra_code, "\\addlinespace\n")
    }

    # Always add \\hline
    extra_code <- paste0(extra_code, "\\hline\n")

    # Apply formatting
    kbl_out <- kableExtra::row_spec(kbl_out, row = i, extra_latex_after = extra_code)
  }

  writeLines(as.character(kbl_out), out_file)

  invisible(kbl_out)
}
