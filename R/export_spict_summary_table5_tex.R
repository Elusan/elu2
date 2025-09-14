#' Export SPiCT Summary Table to LaTeX (Formatted for Beamer)
#'
#' @description This function takes a summary data frame of SPiCT model estimates
#' (as produced by `make_spict_summary_df1_nopipe`) and creates a LaTeX-formatted
#' table with appropriate model headers and footnotes, saving it to disk.
#'
#' @param df_summary A data frame of summary estimates for Fox, Schaefer, and Pella models.
#'                   Must contain columns like: Scenario, Parameters, Fox_estimate, Fox_low, etc.
#' @param scenario_name A character string used to name the output file (e.g., "1" or "S1").
#' @param output_dir A directory path to save the LaTeX `.tex` file (default: current directory).
#'
#' @return Invisibly returns the LaTeX table object as a string.
#'
#' @importFrom kableExtra kbl add_header_above kable_styling footnote cell_spec
#' @export
export_spict_summary_table5_tex <- function(df_summary, scenario_name = "1", output_dir = ".") {
  require(kableExtra)

  # Group definitions
  group1 <- c("$\\alpha_{1}$", "$\\alpha_{2}$", "$\\beta$", "$r$", "$m$", "$K$",
              "$q_{1}$", "$q_{2}$", "$n$", "$\\sigma_B$", "$\\sigma_F$",
              "$\\sigma_{I_{1}}$", "$\\sigma_{I_{2}}$", "$\\sigma_C$")
  group2 <- c("$B_{\\mathrm{MSY}}$", "$F_{\\mathrm{MSY}}$", "$\\mathrm{MSY}$")
  group3 <- c("$B_{2022}$", "$F_{2022}$", "$B_{2022}/B_{\\mathrm{MSY}}$", "$F_{2022}/F_{\\mathrm{MSY}}$")
  group4 <- c("$B_{2024}$", "$F_{2024}$", "$B_{2024}/B_{\\mathrm{MSY}}$", "$F_{2024}/F_{\\mathrm{MSY}}$")

  clean_param <- function(x) gsub("^\\\\textbf\\{(.+)\\}$", "\\1", x)
  plain_parameters <- vapply(df_summary$Parameters, clean_param, character(1))

  group_sequence <- list(group1, group2, group3, group4)
  df_summary$group <- NA_integer_
  for (i in seq_along(group_sequence)) {
    df_summary$group[plain_parameters %in% group_sequence[[i]]] <- i
  }

  df_summary$add_linespace <- plain_parameters %in% c(
    tail(group1, 1), tail(group2, 1), tail(group3, 1)
  )

  df_summary$Parameters <- kableExtra::cell_spec(df_summary$Parameters, format = "latex", bold = TRUE, escape = FALSE)
  n_label <- kableExtra::cell_spec("$n$", format = "latex", bold = TRUE, escape = FALSE)

  for (col in c("Fox", "Schaefer")) {
    idx <- which(df_summary$Parameters == n_label)
    if (length(idx)) {
      df_summary[[paste0(col, "_estimate")]][idx] <- "Fixed"
      df_summary[[paste0(col, "_low")]][idx] <- ""
      df_summary[[paste0(col, "_up")]][idx] <- ""
    }
  }

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
  ) %>%
    add_header_above(
      header = c(" " = 2, "Fox" = 3, "Schaefer" = 3, "Pella-Tomlinson" = 3),
      escape = FALSE
    ) %>%
    kable_styling(
      latex_options = c("hold_position"),
      full_width = FALSE,
      position = "center",
      font_size = 8
    ) %>%
    footnote(
      general = "‘Fixed’ indicates that the parameter was not estimated but assigned a fixed value.",
      general_title = "Note: ",
      footnote_as_chunk = TRUE,
      escape = FALSE
    )

  # --- Manual LaTeX modification to insert \\hline after every row ---
  latex_lines <- as.character(kbl_out)
  lines_split <- strsplit(latex_lines, "\n")[[1]]

  # Correctly escape \\begin and \\end for matching in R strings
  start_line <- which(grepl("^\\\\begin\\{tabular", lines_split))
  end_line   <- which(grepl("^\\\\end\\{tabular", lines_split))

  new_body <- c()
  nrows_df <- nrow(df_summary)
  row_idx <- 1  # Index for df_summary rows

  for (i in seq((start_line + 1), (end_line - 1))) {
    line <- lines_split[i]

    if (nzchar(trimws(line))) {
      # Append row content with line break and hline
      new_body <- c(new_body, paste0(line, " \\\\ \\\\hline"))

      # Optionally add linespace after specific rows
      if (!is.na(df_summary$add_linespace[row_idx]) && df_summary$add_linespace[row_idx]) {
        new_body <- c(new_body, "\\\\addlinespace")
      }

      row_idx <- row_idx + 1
    } else {
      new_body <- c(new_body, line)
    }
  }

  final_lines <- c(
    lines_split[1:start_line],
    new_body,
    lines_split[end_line:length(lines_split)]
  )

  writeLines(final_lines, out_file)

  invisible(kbl_out)
}
