#' Delta GLM Model List
#'
#' A list of delta-GLM fitted models used for shrimp CPUE standardization.
#'
#' @format A named list of fitted model objects. Each element corresponds to a year-season combination.
#'   The models are generally of class `glm` or `lm`, depending on the family used.
#'
#' @details This object was created during index standardization prior to inclusion in the SPiCT model.
#'   It is used to extract time series of standardized CPUE for short/long series comparison.
#'
#' @seealso [inp1.glm], [Long_delta_GLM_list]
#'
#' @keywords datasets
"Delta_GLM_list"
