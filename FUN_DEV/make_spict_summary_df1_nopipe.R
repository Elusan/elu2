#' Create a Summary Data Frame of SPiCT Model Estimates (No Pipe Version)
#'
#' @description This function extracts and formats parameter estimates, reference points,
#' state variables, and predictions from multiple SPiCT model fits. It returns a wide-format
#' data frame suitable for table export, with LaTeX-formatted parameter names and rounded estimates.
#'
#' @param model_list A list of SPiCT model objects, such as output from `fit.spict()`.
#' @param model_names A character vector of model names corresponding to `model_list`.
#'                   Each name should follow the format "Scenario_Model" (e.g., "S1_Fox").
#' @param exclude_parameters A character vector of parameter names to exclude from the output
#'                           (default includes "rc", "rold", "Catch_2023.00", "E(B_inf)").
#'
#' @return A wide-format data frame with columns for each model's estimate, lower CI, and upper CI,
#'         and LaTeX-formatted parameter names.
#'
#' @importFrom dplyr filter mutate select relocate bind_rows
#' @importFrom tidyr pivot_wider
#' @export
make_spict_summary_df1_nopipe <- function(model_list, model_names,
                                          exclude_parameters = c("rc", "rold", "Catch_2023.00", "E(B_inf)")) {
  require(dplyr)
  require(tidyr)

  param_labels <- c(
    "alpha1" = "$\\alpha_{1}$", "alpha2" = "$\\alpha_{2}$", "beta" = "$\\beta$",
    "r" = "$r$", "m" = "$m$", "K" = "$K$", "q1" = "$q_{1}$", "q2" = "$q_{2}$",
    "n" = "$n$", "sdb" = "$\\sigma_B$", "sdf" = "$\\sigma_F$", "sdi1" = "$\\sigma_{I_{1}}$",
    "sdi2" = "$\\sigma_{I_{2}}$", "sdc" = "$\\sigma_C$", "Bmsys" = "$B_{\\mathrm{MSY}}$",
    "Fmsys" = "$F_{\\mathrm{MSY}}$", "MSYs" = "$\\mathrm{MSY}$",
    "B_2022.94" = "$B_{2022}$", "F_2022.94" = "$F_{2022}$",
    "B_2022.94/Bmsy" = "$B_{2022}/B_{\\mathrm{MSY}}$", "F_2022.94/Fmsy" = "$F_{2022}/F_{\\mathrm{MSY}}$",
    "B_2024.00" = "$B_{2024}$", "F_2024.00" = "$F_{2024}$",
    "B_2024.00/Bmsy" = "$B_{2024}/B_{\\mathrm{MSY}}$", "F_2024.00/Fmsy" = "$F_{2024}/F_{\\mathrm{MSY}}$"
  )

  cleaned_list <- list()
  for (i in seq_along(model_list)) {
    model <- model_list[[i]]
    model_name <- model_names[[i]]

    summary_data <- tryCatch({
      rbind(
        sumspict.parest(model),
        sumspict.srefpoints(model)[, 1:4],
        sumspict.states(model)[, 1:4],
        sumspict.predictions(model)[, 1:4]
      )
    }, error = function(e) {
      message(paste("Error in model", model_name, ":", e$message))
      return(NULL)
    })

    if (!is.null(summary_data)) {
      df <- as.data.frame(summary_data)
      df$Parameters <- trimws(rownames(df))
      df <- df[, c("Parameters", "estimate", "cilow", "ciupp")]
      df <- df[!df$Parameters %in% exclude_parameters, , drop = FALSE]
      df$estimate <- round(df$estimate, 3)
      df$cilow <- round(df$cilow, 3)
      df$ciupp <- round(df$ciupp, 3)
      df$Scenario <- sub("_(Pella|Schaefer|Fox)", "", model_name)
      df$Model <- sub("^.*_(Pella|Schaefer|Fox)$", "\\1", model_name)
      df$Parameters <- ifelse(df$Parameters %in% names(param_labels),
                              param_labels[df$Parameters],
                              df$Parameters)
      cleaned_list[[i]] <- df[, c("Scenario", "Parameters", "Model", "estimate", "cilow", "ciupp")]
    }
  }

  all_df <- do.call(rbind, cleaned_list)

  df_wide <- tidyr::pivot_wider(
    data = all_df,
    names_from = "Model",
    values_from = c("estimate", "cilow", "ciupp"),
    names_glue = "{Model}_{.value}"
  )

  df_final <- df_wide[, c("Scenario", "Parameters",
                          "Fox_estimate", "Fox_cilow", "Fox_ciupp",
                          "Schaefer_estimate", "Schaefer_cilow", "Schaefer_ciupp",
                          "Pella_estimate", "Pella_cilow", "Pella_ciupp")]

  names(df_final) <- c(
    "Scenario", "Parameters",
    "Fox_estimate", "Fox_low", "Fox_up",
    "Schaefer_estimate", "Schaefer_low", "Schaefer_up",
    "Pella_estimate", "Pella_low", "Pella_up"
  )

  # Reorder Parameters by original param_labels vector
  ordered_levels <- unname(param_labels)

  df_final$Parameters <- factor(df_final$Parameters, levels = ordered_levels)
  df_final <- df_final[order(df_final$Parameters), ]


  return(df_final)
}
