#' Estimate MSY from catch and index overlap
#'
#' This function provides a heuristic estimate of Maximum Sustainable Yield (MSY) using
#' the overlap between catch and index time series. If sufficient overlap exists,
#' a linear model is fitted to approximate surplus production dynamics. Otherwise,
#' the mean catch is returned.
#'
#' @param inp A SPiCT input list containing at least \code{obsC} (observed catch) and,
#' optionally, one or more index series \code{obsI}.
#' @param all.return Logical; if \code{TRUE}, a list of detailed regression outputs is returned
#'                   instead of just the estimated MSY. Default is \code{FALSE}.
#'
#' @details
#' When both catch and index data overlap sufficiently in time, the function:
#' \itemize{
#'   \item Extracts overlapping values of catch (\code{y}) and index (\code{z}),
#'   \item Computes the ratio \code{x = y/z},
#'   \item Fits a linear model \code{lm(z ~ x)},
#'   \item Uses the fitted coefficients \code{a, b} to compute:
#'     \deqn{MSY = -0.25 \cdot \frac{a^2}{b}}
#'     \deqn{E_{MSY} = -0.5 \cdot \frac{a}{b}}
#'   \item If \code{MSY <= 0} or \code{R² < 0.3}, returns the mean catch instead.
#' }
#' If no index data is available or the overlap is insufficient, the mean catch is returned.
#'
#' @return If \code{all.return = FALSE}, a single numeric value (MSY estimate).
#' If \code{all.return = TRUE}, a list with the following components:
#' \describe{
#'   \item{MSY}{Estimated Maximum Sustainable Yield}
#'   \item{Emsy}{Effort at MSY, from linear approximation}
#'   \item{a, b}{Intercept and slope of the linear model}
#'   \item{x}{Computed x values (catch/index)}
#'   \item{y}{Overlapping catch values}
#'   \item{z}{Overlapping index values}
#'   \item{mod0}{Fitted \code{lm} object}
#' }
#'
#' @seealso \code{\link{get.catchindexoverlap}}, \code{\link{lm}}
#'
#' @examples
#' \dontrun{
#' inp <- make.spict.data(...)  # or prepare your SPiCT input manually
#' guess.m(inp)                 # simple MSY guess
#' guess.m(inp, all.return = TRUE)  # full regression output
#' }
#'
#' @export
guess.m <- function (inp, all.return = FALSE)
{
  meancatch <- mean(unlist(inp$obsC))
  flag <- FALSE
  if ("obsI" %in% names(inp)) {
    if (length(inp$obsI) > 0) {
      out <- get.catchindexoverlap(inp)
      if (!is.null(out)) {
        ty <- out$ty
        y <- out$y
        tz <- out$tz
        z <- out$z
        if (length(y) == length(z)) {
          flag <- TRUE
        }
      }
    }
  }
  if (flag) {
    x <- y/z
    mod0 <- lm(z ~ x)
    a <- mod0$coefficients[1]
    b <- mod0$coefficients[2]
    MSY <- -0.25 * a^2/b
    if (MSY <= 0 | summary(mod0)$r.squared < 0.3) {
      MSY <- meancatch
    }
    if (all.return) {
      Emsy <- -0.5 * a/b
      return(list(MSY = MSY, Emsy = Emsy, a = a, b = b,
                  x = x, y = y, z = z, mod0 = mod0))
    }
    else {
      return(MSY)
    }
  }
  else {
    return(meancatch)
  }
}
