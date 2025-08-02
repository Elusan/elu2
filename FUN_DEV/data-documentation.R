#' sdmTMB Input List
#'
#' List of inputs for sdmTMB modeling, pre-processed for use in spatio-temporal stock assessment of shrimp.
#'
#' @format A list.
#' @source Created from original SPICT_DATA/inp1.sdmTMB.rds file.
"sdmTMB_inp"

#' GLM Input List
#'
#' List of inputs for delta-GLM modeling, prepared for stock assessment.
#'
#' @format A list.
#' @source Created from original SPICT_DATA/inp1.glm.rds file.
"GLM_inp"

#' List of sdmTMB Model Outputs
#'
#' List of sdmTMB model output objects.
#'
#' @format A list.
#' @source From SPICT_DATA/Long_sdmTMB_list.rds
"Long_sdmTMB_list"

#' List of delta-GLM Model Outputs
#'
#' List of delta-GLM model output objects.
#'
#' @format A list.
#' @source From SPICT_DATA/Long_delta_GLM_list.rds
"Long_delta_GLM_list"

#' Delta-GLM Short Index Data (1992–2022)
#'
#' This dataset contains standardized delta-GLM index values and catch records
#' from 1992 to 2022. It is used for visualizing long and short index trends
#' (pre-2017 and post-2017) alongside catch values.
#'
#' @format A data frame with 31 rows and 3 columns:
#' \describe{
#'   \item{Year}{Year from 1992 to 2022}
#'   \item{Catch}{Catch in metric tons}
#'   \item{Index_stand}{Standardized delta-GLM index}
#' }
#' @source Derived from full_data1992_2022_catch_and_Index.rds
"Data_glm_short"

#' Delta-GLM Full Index Data (1992–2022)
#'
#' This dataset contains an alternative version of the delta-GLM standardized index
#' and catch values over the full period 1992–2022.
#'
#' @format A data frame with 31 rows and 3 columns:
#' \describe{
#'   \item{Year}{Year from 1992 to 2022}
#'   \item{Catch}{Catch in metric tons}
#'   \item{Index_stand}{Alternative delta-GLM index values}
#' }
#' @source Derived from Dataglm_full_catch_and_index.rds
"Data_glm_full"
#' sdmTMB Index Data (1992–2022)
#'
#' Standardized index from the sdmTMB delta model for shrimp distribution,
#' including catch values. Used in post-2017 time-series visualization.
#'
#' @format A data frame with 31 rows and 3 columns:
#' \describe{
#'   \item{Year}{Year from 1992 to 2022}
#'   \item{Catch}{Catch in metric tons}
#'   \item{Index}{Standardized sdmTMB index}
#' }
#' @source Derived from gg_sdmTMB_long_short_data.rds
"Data_sdm"




