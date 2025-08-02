#' Extract Equilibrium Production Curve from SPiCT Model Fit
#'
#' Computes and returns equilibrium surplus production values across a biomass gradient for a fitted SPiCT model (e.g., Fox, Schaefer, or Pella-Tomlinson) for plotting production curves.
#'
#' @param rep A fitted SPiCT model object (usually the output of `fit.spict()` or a compatible function). Must contain model parameter estimates and their names.
#' @param model_name Character. Name for the model to be included in the returned data frame (default: "Model").
#'
#' @details
#' For each fitted model, this function:
#' \itemize{
#'   \item Retrieves the most recent estimated parameters: carrying capacity (K), intrinsic rate of increase (m), shape parameter (n), and scaling coefficient (`gamma`).
#'   \item Computes the equilibrium production across a grid of relative biomass values (`B/K`), using the generic surplus production formula.
#'   \item Returns a data frame suitable for use in ggplot2 or further aggregation.
#' }
#' If the input model object contains a `"sderr"` entry (indicating a stochastic or error model, rather than a deterministic fit), the function returns `NULL`.
#'
#' @return
#' A data frame with columns:
#'   \describe{
#'     \item{B_K}{Relative biomass (B/K).}
#'     \item{Production}{Estimated equilibrium production at each biomass.}
#'     \item{Model}{Model name, for use in faceting or coloring in plots.}
#'   }
#' Returns `NULL` if the input is not compatible.
#'
#' @examples
#' \dontrun{
#'   df <- get_production_data(fit.spict_result, model_name = "Pella-Tomlinson")
#'   ggplot(df, aes(x = B_K, y = Production, color = Model)) + geom_line()
#' }
#'
#' @export
get_production_data <- function(rep, model_name = "Model") {
  if (!"sderr" %in% names(rep)) {
    inp <- rep$inp
    tvgflag <- rep$inp$timevaryinggrowth | rep$inp$logmcovflag
    Kest <- get.par("logK", rep, exp = TRUE)
    mest <- get.par("logm", rep, exp = TRUE)
    gamma <- get.par("gamma", rep)
    n <- get.par("logn", rep, exp = TRUE)
    nr <- dim(mest)[1]

    Bplot <- seq(0.5 * 1e-08, Kest[2], length.out = 200)
    pfun <- function(gamma, m, K, n, B) gamma * m/K * B * (1 - (B/K)^(n - 1))

    prod <- pfun(gamma[2], mest[nr, 2], Kest[2], n[2], Bplot)

    data.frame(
      B_K = Bplot / Kest[2],
      Production = prod,
      Model = model_name
    )
  } else {
    NULL
  }
}
