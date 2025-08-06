#' Export SPiCT Summary Table (Base R, LaTeX)
#'
#' Creates and exports a summary parameter table from multiple SPiCT models
#' (Pella, Fox, Schaefer), arranged by scenario and parameter, in LaTeX-ready format.
#'
#' @param models A named list of SPiCT model objects. Names should contain scenario and model type (e.g., S1_Pella).
#' @return A formatted LaTeX table printed to console (use `cat()` to capture) or to file using capture.output().
#' @importFrom tibble tibble
#' @importFrom knitr kable is_latex_output
#' @importFrom kableExtra kable_styling add_header_above group_rows
#' @export
export_spict_summary_table5_baed_R_tex <- function(models) {
  library(tibble)
  library(knitr)
  library(kableExtra)

  cleaned_results <- list()

  for (i in seq_along(models)) {
    model_name <- names(models)[i]
    model <- models[[i]]

    result <- tryCatch({
      df1 <- sumspict.parest(model)
      df2 <- sumspict.srefpoints(model)[, 1:4]
      df3 <- sumspict.states(model)[, 1:4]
      df4 <- sumspict.predictions(model)[, 1:4]

      df <- rbind(df1, df2, df3, df4)
      df <- as.data.frame(df)
      df$Parameters <- rownames(df)
      rownames(df) <- NULL

      drop_rows <- c(5, 6, 28, 29)
      drop_rows <- drop_rows[drop_rows <= nrow(df)]
      if (length(drop_rows) > 0) {
        df <- df[-drop_rows, ]
      }

      if (ncol(df) > 4) {
        df <- df[, 1:4]
      }

      colnames(df)[3:4] <- c("low", "up")

      for (j in seq_along(df)) {
        if (is.numeric(df[[j]])) {
          df[[j]] <- round(df[[j]], 3)
        }
      }

      scenario <- sub("_(Pella|Schaefer|Fox)", "", model_name)
      df$Scenario <- scenario
      df <- df[, c("Scenario", "Parameters", "estimate", "low", "up")]
      df$Model_Name <- model_name
      df
    }, error = function(e) {
      message(paste("Error in model", model_name, ":", e$message))
      return(NULL)
    })

    if (!is.null(result)) {
      cleaned_results[[model_name]] <- result
    }
  }

  all_models_df <- do.call(rbind, cleaned_results)
  model_vec <- character(nrow(all_models_df))
  for (i in seq_len(nrow(all_models_df))) {
    model_vec[i] <- if (grepl("Fox", all_models_df$Model_Name[i])) {
      "Fox"
    } else if (grepl("Schaefer", all_models_df$Model_Name[i])) {
      "Schaefer"
    } else {
      "Pella"
    }
  }
  all_models_df$Model <- model_vec

  unique_keys <- unique(all_models_df[, c("Scenario", "Parameters")])
  wide_df <- data.frame()

  for (i in seq_len(nrow(unique_keys))) {
    row_key <- unique_keys[i, ]
    scenario_i <- row_key$Scenario
    param_i <- row_key$Parameters

    sub_df <- all_models_df[all_models_df$Scenario == scenario_i &
                              all_models_df$Parameters == param_i, ]

    new_row <- list(Scenario = scenario_i, Parameters = param_i)

    for (j in seq_len(nrow(sub_df))) {
      model_name <- sub_df$Model[j]
      new_row[[paste0(model_name, "_estimate")]] <- sub_df$estimate[j]
      new_row[[paste0(model_name, "_low")]]     <- sub_df$low[j]
      new_row[[paste0(model_name, "_up")]]      <- sub_df$up[j]
    }

    wide_df <- rbind(wide_df, as.data.frame(new_row, stringsAsFactors = FALSE))
  }

  param_levels <- c(
    "alpha1", "alpha2", "beta", "r  ", "m", "K", "q1", "q2", "n",
    "sdb", "sdf", "sdi1", "sdi2", "sdc", "Bmsys", "Fmsys", "MSYs",
    "B_2022.94", "F_2022.94", "B_2022.94/Bmsy", "F_2022.94/Fmsy",
    "B_2024.00", "F_2024.00", "B_2024.00/Bmsy", "F_2024.00/Fmsy"
  )

  wide_df$Parameters <- factor(wide_df$Parameters, levels = param_levels)
  wide_df <- wide_df[!is.na(wide_df$Parameters), ]
  wide_df <- wide_df[order(wide_df$Parameters, wide_df$Scenario), ]

  latex_labels <- c(
    alpha1 = "$\\alpha_{1}$", alpha2 = "$\\alpha_{2}$", beta = "$\\beta$",
    "r  " = "$r$", m = "$m$", K = "$K$",
    q1 = "$q_{1}$", q2 = "$q_{2}$", n = "$n$",
    sdb = "$\\sigma_B$", sdf = "$\\sigma_F$",
    sdi1 = "$\\sigma_{I_{1}}$", sdi2 = "$\\sigma_{I_{2}}$", sdc = "$\\sigma_C$",
    Bmsys = "$B_{\\mathrm{MSY}}$", Fmsys = "$F_{\\mathrm{MSY}}$", MSYs = "$\\mathrm{MSY}$",
    "B_2022.94" = "$B_{2022}$", "F_2022.94" = "$F_{2022}$",
    "B_2022.94/Bmsy" = "$B_{2022}/B_{\\mathrm{MSY}}$",
    "F_2022.94/Fmsy" = "$F_{2022}/F_{\\mathrm{MSY}}$",
    "B_2024.00" = "$B_{2024}$", "F_2024.00" = "$F_{2024}$",
    "B_2024.00/Bmsy" = "$B_{2024}/B_{\\mathrm{MSY}}$",
    "F_2024.00/Fmsy" = "$F_{2024}/F_{\\mathrm{MSY}}$"
  )

  latex_vec <- character(nrow(wide_df))
  for (i in seq_len(nrow(wide_df))) {
    param_raw <- as.character(wide_df$Parameters[i])
    latex_vec[i] <- if (param_raw %in% names(latex_labels)) {
      latex_labels[[param_raw]]
    } else {
      param_raw
    }
  }
  wide_df$Parameters <- latex_vec

  group1 <- c("$\\alpha_{1}$", "$\\alpha_{2}$", "$\\beta$", "$r$", "$m$", "$K$",
              "$q_{1}$", "$q_{2}$", "$n$", "$\\sigma_B$", "$\\sigma_F$",
              "$\\sigma_{I_{1}}$", "$\\sigma_{I_{2}}$", "$\\sigma_C$")
  group2 <- c("$B_{\\mathrm{MSY}}$", "$F_{\\mathrm{MSY}}$", "$\\mathrm{MSY}$")
  group3 <- c("$B_{2022}$", "$F_{2022}$", "$B_{2022}/B_{\\mathrm{MSY}}$", "$F_{2022}/F_{\\mathrm{MSY}}$")
  group4 <- c("$B_{2024}$", "$F_{2024}$", "$B_{2024}/B_{\\mathrm{MSY}}$", "$F_{2024}/F_{\\mathrm{MSY}}$")
  all_groups <- list(group1, group2, group3, group4)

  group_id <- rep(NA_integer_, nrow(wide_df))
  for (g in seq_along(all_groups)) {
    group_id[wide_df$Parameters %in% all_groups[[g]]] <- g
  }
  wide_df$group <- group_id

  param_counts <- table(wide_df$Parameters)
  param_counts <- as.list(param_counts)

  col_names <- c("Scenario", "Parameter",
                 "Fox Estimate", "Fox Low", "Fox Up",
                 "Schaefer Estimate", "Schaefer Low", "Schaefer Up",
                 "Pella Estimate", "Pella Low", "Pella Up")

  out <- kbl(
    wide_df,
    format = ifelse(knitr::is_latex_output(), "latex", "html"),
    booktabs = TRUE,
    escape = FALSE,
    row.names = FALSE,
    col.names = col_names
  )

  out <- add_header_above(out, c(" " = 2,
                                 "Fox" = 3,
                                 "Schaefer" = 3,
                                 "Pella-Tomlinson" = 3))

  out <- group_rows(out, index = param_counts)

  out <- kable_styling(
    out,
    latex_options = c("hold_position", "scale_down"),
    bootstrap_options = c("striped", "hover"),
    full_width = FALSE,
    position = "center"
  )

  return(out)
}
