#' Generate Scenario Definitions with Flexible Prior Overrides
#'
#' Constructs a structured list of model scenario definitions for SPiCT fitting,
#' with the ability to override any prior parameter per scenario.
#'
#' @param scenario_names Character vector of scenario names (e.g., \code{c("S1", "S2", ...)}).
#' @param prior_list_per_scenario A named list of prior overrides per scenario. Each element can contain
#'   custom values for any of \code{logr}, \code{logn}, \code{logalpha}, \code{logsdi}, etc.
#' @param base_priors A named list of default prior values to be used unless overridden (default: \code{base_priors}).
#' @param ini_vals A named list with initial values for Fox (\code{"F"}) and Schaefer (\code{"S"}).
#' @param logn_def Default prior for \code{logn} if not overridden (default: \code{c(1, 1, 0)}).
#' @param phase_def Default phase value for \code{logn} (default: -1).
#'
#' @return A named list of scenario definitions, each with submodels \code{SxP}, \code{SxF}, and \code{SxS}.
#'
#' @examples
#' prior_overrides <- list(
#'   S1 = list(),
#'   S2 = list(logr = c(log(0.6), 0.4, 1)),
#'   S3 = list(logr = c(log(0.9), 0.2, 1), logalpha = c(0, 1, 1))
#' )
#' scenario_definitions <- generate_scenario_definitions(
#'   scenario_names = names(prior_overrides),
#'   prior_list_per_scenario = prior_overrides
#' )
#'
#' @export
generate_scenario_definitions <- function(
    scenario_names,
    prior_list_per_scenario = list(),
    base_priors = base_priors,
    ini_vals = list(F = log(1.001), S = log(2)),
    logn_def = c(1, 1, 0),
    phase_def = -1
) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  out <- list()

  for (scen in scenario_names) {
    priors <- prior_list_per_scenario[[scen]]
    if (is.null(priors)) priors <- list()

    # Merge with base_priors
    p <- modifyList(base_priors, priors)

    scen_list <- list()

    # Pella
    scen_list[[paste0(scen, "P")]] <- list()
    if (!is.null(p$logr)) scen_list[[paste0(scen, "P")]]$logr <- p$logr

    # Fox
    scen_list[[paste0(scen, "F")]] <- list(
      logr  = p$logr,
      logn  = p$logn %||% logn_def,
      phase = phase_def,
      ini   = ini_vals$F
    )

    # Schaefer
    scen_list[[paste0(scen, "S")]] <- list(
      logr  = p$logr,
      logn  = p$logn %||% logn_def,
      phase = phase_def,
      ini   = ini_vals$S
    )

    out[[scen]] <- scen_list
  }

  return(out)
}
