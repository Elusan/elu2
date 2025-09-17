#' Export Prior Table for a Scenario to LaTeX
#'
#' @description This function extracts the priors from SPiCT fits (Fox, Schaefer, Pella),
#' formats them in LaTeX, and saves a `table_scenarioX.tex` file for Beamer inclusion.
#'
#' @param S_Fox A SPiCT model fit for the Fox model.
#' @param S_Schaefer A SPiCT model fit for the Schaefer model.
#' @param S_Pella A SPiCT model fit for the Pella-Tomlinson model.
#' @param scenario_name A string used to name the output file (e.g., `"1"` or `"S1"`).
#' @param output_dir A path to the folder where the .tex file should be saved. Default is current directory.
#'
#' @return Invisibly returns the LaTeX string for further use if needed.
#' @export
export_prior_table_tex <- function(S_Fox, S_Schaefer, S_Pella,
                                   scenario_name = "1",
                                   output_dir = ".") {
  stopifnot(dir.exists(output_dir))

  # Format Prior column
  format_prior <- function(df) {
    df$Prior <- gsub("(\\d+(\\.\\d+)?)(\\^2)", "(\\1)^{2}", df$Prior)
    df$Prior <- paste0("\\( ", df$Prior, " \\)")
    df
  }

  df_fox      <- format_prior(sumspict.priors.df(S_Fox))
  df_schaefer <- format_prior(sumspict.priors.df(S_Schaefer))
  df_pella    <- format_prior(sumspict.priors.df(S_Pella))

  df_fox$Model      <- "Fox"
  df_schaefer$Model <- "Schaefer"
  df_pella$Model    <- "Pella-Tomlinson"

  blank_row <- data.frame(Model = "", Parameter = "", Prior = "")

  df_all <- dplyr::bind_rows(
    df_pella,
    blank_row,
    df_schaefer,
    blank_row,
    df_fox
  )

  # Filename
  fname <- file.path(output_dir, paste0("table_scenario", scenario_name, ".tex"))

  # Generate LaTeX table
  tmp_kbl <- kableExtra::kbl(
    df_all,
    format = "latex",
    booktabs = FALSE,
    escape = FALSE,
    linesep = "\addlinespace",
    col.names = c("\\textbf{Parameter}", "\\textbf{Prior}", "\\textbf{Model}")
  )

  tab_latex <- kableExtra::kable_styling(
    tmp_kbl,
    latex_options = c("hold_position"),
    stripe_color = c("blue!20", "blue!10"),
    font_size = 8,
    full_width = FALSE
  )


  # Wrap in LaTeX table environment
  custom_table <- paste0(
    "\\begin{footnotesize}\n",
    "\\rowcolors{1}{blue!20}{blue!10}\n",
    "\\begin{table}[htbp]\n",
    tab_latex,
    "\\end{table}\n",
    "\\end{footnotesize}\n"
  )

  # Save to .tex
  writeLines(custom_table, con = fname)

  invisible(custom_table)
}
