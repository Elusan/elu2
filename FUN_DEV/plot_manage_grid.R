#' Plot SPiCT Management Scenario Outputs for Multiple Models
#'
#' @param managed_models A named list of SPiCT objects (e.g., list(Fox = ..., Schaefer = ..., Pella = ...))
#' @param include_abs Logical, include absolute values for B and F? Default FALSE.
#' @param include_unc Logical, include 95% CI ribbons? Default TRUE.
#' @param caption Optional string for plot title.
#'
#' @return A ggplot2 object combining B/Bmsy and F/Fmsy across models and scenarios.
#' @export
plot_manage_grid <- function(managed_models,
                             include_abs = FALSE,
                             include_unc = TRUE,
                             caption = NULL) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)

  # Helper to extract and reshape data from one model
  extract_manage_df <- function(model_obj, model_name) {
    res <- sumspict.manage(model_obj, include.abs = include_abs, include.unc = include_unc, verbose = FALSE)
    est_df <- res$est %>%
      as.data.frame() %>%
      tibble::rownames_to_column("Scenario") %>%
      mutate(Model = model_name)

    if (include_unc) {
      unc_df <- res$unc %>%
        as.data.frame() %>%
        tibble::rownames_to_column("Scenario") %>%
        mutate(Model = model_name)

      # Combine for ribbon plotting
      est_long <- est_df %>%
        select(Model, Scenario, `B/Bmsy`, `F/Fmsy`) %>%
        pivot_longer(cols = c(`B/Bmsy`, `F/Fmsy`), names_to = "Metric", values_to = "Estimate")

      unc_long <- unc_df %>%
        select(Model, Scenario, starts_with("B/Bmsy"), starts_with("F/Fmsy")) %>%
        pivot_longer(cols = -c(Model, Scenario), names_to = "Metric", values_to = "Value") %>%
        mutate(Quantile = ifelse(grepl("\\.lo$", Metric), "Lower",
                                 ifelse(grepl("\\.hi$", Metric), "Upper", "Median")),
               Metric = gsub("\\.(lo|hi)$", "", Metric)) %>%
        pivot_wider(names_from = Quantile, values_from = Value)

      return(left_join(est_long, unc_long, by = c("Model", "Scenario", "Metric")))
    } else {
      return(est_df %>%
               select(Model, Scenario, `B/Bmsy`, `F/Fmsy`) %>%
               pivot_longer(cols = c(`B/Bmsy`, `F/Fmsy`), names_to = "Metric", values_to = "Estimate"))
    }
  }

  # Combine all model results
  df_all <- imap_dfr(managed_models, extract_manage_df)

  # Plot
  p <- ggplot(df_all, aes(x = Scenario, y = Estimate, fill = Metric, group = interaction(Metric, Model))) +
    geom_col(position = position_dodge(width = 0.8), width = 0.6, color = "black") +
    facet_wrap(~Model, nrow = 1) +
    scale_fill_manual(values = c("B/Bmsy" = "steelblue", "F/Fmsy" = "firebrick")) +
    coord_cartesian(ylim = c(0, NA)) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    ) +
    labs(x = NULL, y = "Relative Status", fill = NULL,
         title = caption %||% "Management Scenarios Across Models")

  if (include_unc) {
    p <- p +
      geom_errorbar(
        aes(ymin = Lower, ymax = Upper),
        position = position_dodge(width = 0.8),
        width = 0.2
      )
  }

  return(p)
}
