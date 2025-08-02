#' @title Extract and Format Priors from SPiCT Model
#'
#' @description This function extracts the active prior distributions used in a fitted SPiCT model object
#' (including custom priors and random effects priors), and formats them into a data frame suitable
#' for LaTeX export or reporting. The output includes the parameter name and its associated prior distribution
#' in a readable mathematical format (e.g., `dnorm`, `dgamma`).
#'
#' @param rep A fitted SPiCT model object, typically the result of `fit.spict()` or a custom equivalent
#' like `fit.elu2()`. This object must contain a valid `rep$inp$priors` and optionally `rep$inp$matrixpriors`.
#' @param ndigits Number of digits to use for rounding values in the output (currently not used directly in this function).
#'
#' @return A `data.frame` with columns:
#' \describe{
#'   \item{Parameter}{The name of the parameter with an active prior.}
#'   \item{Prior}{A string describing the prior distribution (e.g., `~ dnorm[log(0.7), 0.2^2]`).}
#' }
#' If no priors are used, the returned data frame will contain a single row stating `"No priors are used"`.
#'
#' @details This function:
#' \itemize{
#'   \item Identifies which priors are active using `priorsuseflag`.
#'   \item Handles standard priors, matrix priors, and random effect priors (e.g., `logB`, `logF`).
#'   \item Converts gamma priors to mean/sd format via `shaperate2meanvar`.
#'   \item Outputs clean formatted strings that can be easily exported using `kable()` or `xtable()`.
#' }
#'
#' @seealso \code{\link{sumspict.priors}} for a console-printing version of this function.
#'
#' @examples
#' \dontrun{
#'   df_priors <- sumspict.priors.df(fit.spict(my_data))
#'   kableExtra::kbl(df_priors, format = "latex")
#' }
#'
#' @export
sumspict.priors.df <- function(rep, ndigits = 8) {
  indso <- which(rep$inp$priorsuseflag == 1)
  if (length(indso) > 0) {
    priors <- rep$inp$priors[indso]
    usepriors <- names(priors)
    usepriors <- gsub('gamma', '', usepriors)
    repriors <- c('logB', 'logF', 'logBBmsy', 'logFFmsy')

    if (any(repriors %in% usepriors)) {
      inds <- na.omit(match(repriors, usepriors))
      for (i in seq_along(inds)) {
        names(priors)[inds[i]] <- paste0(usepriors[inds[i]], fd(priors[[inds[i]]][4]))
        priors[[inds[i]]] <- priors[[inds[i]]][1:3]
      }
    }

    matpriors <- rep$inp$matrixpriors
    nmmatpriors <- names(matpriors)
    priorsmat <- do.call(rbind, priors[!names(priors) %in% nmmatpriors])
    priorsmat <- rbind(priorsmat, do.call(rbind, matpriors))
    storage.mode(priorsmat) <- "double"
    priorsmat <- priorsmat[priorsmat[, 3] == 1, , drop = FALSE]

    gammainds <- grep("gamma", rownames(priorsmat))
    npriors <- nrow(priorsmat)
    usepriors <- rownames(priorsmat)

    Prior <- character(npriors)

    for (i in 1:npriors) {
      if (i %in% gammainds) {
        shape <- priorsmat[i, 1]
        rate <- priorsmat[i, 2]
        vec <- shaperate2meanvar(shape, rate)
        Prior[i] <- paste0(
          "~ dgamma[", shape, ", ", rate, "] ",
          "(mean=", vec[1], ", sd=", vec[3], ")"
        )
      } else {
        Prior[i] <- paste0(
          "~ dnorm[log(", exp(priorsmat[i, 1]), "), ",
          priorsmat[i, 2], "^2]",
          ifelse(priorsmat[i, 2] <= 1e-3, " (fixed)", "")
        )
      }
    }

    df <- data.frame(
      Parameter = usepriors,
      Prior = Prior,
      stringsAsFactors = FALSE
    )
    return(df)
  } else {
    return(data.frame(Parameter = "None", Prior = "No priors are used"))
  }
}
