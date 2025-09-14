#' Save Prior–Posterior Plots for All Scenarios and Models
#'
#' Iterates over a nested \code{models_all} object and saves one PNG per fitted
#' model using a user-supplied plotting function (defaults to \code{priors.elu6}).
#' This is designed for objects produced by your scenario runs, where each model
#' entry may be either a direct fitted object or a small list like
#' \code{list(result = <fit>, error = <condition_or_NULL>)}.
#'
#' @param models_all A nested named list of the form
#'   \code{models_all[[scenario]][[model_tag]]}, where each leaf is either
#'   (a) a fitted ELU/SPiCT-like object, or
#'   (b) a list with elements \code{result} (the fitted object) and optional
#'   \code{error}. If \code{error} is non-\code{NULL}, that model is skipped.
#' @param priors_fun Function that returns a ggplot/patchwork object for a given
#'   fitted model. Default: \code{priors.elu6}. If \code{priors_fun} has a formal
#'   argument named \code{model_id}, it will be supplied automatically.
#' @param outdir Base output directory. One subfolder per scenario is created
#'   inside this directory. Default: \code{"FIG/Priors"}.
#' @param width,height Numeric; size in inches for the saved PNGs. Defaults: 8, 6.
#' @param dpi Numeric; dots per inch for the saved PNGs. Default: 300.
#'
#' @return (Invisibly) a character vector of file paths that were successfully saved.
#' @examples
#' \dontrun{
#'   paths <- elu_nice_priors_all(models_all)
#' }
#' @export
elu_nice_priors_all <- function(models_all,
                                priors_fun = priors.elu6,
                                outdir = "FIG/Priors",
                                width = 8, height = 6, dpi = 300) {
  if (!is.list(models_all) || length(models_all) == 0L) {
    stop("`models_all` must be a non-empty named list of scenarios.")
  }

  # Create base outdir (recursively, quietly)
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  }

  # Helper: extract a "fit" object from a leaf that might be a list(result=..., error=...)
  get_fit_obj <- function(x) {
    if (is.list(x) && !is.null(x$error)) return(NULL)   # flagged error → skip
    if (is.list(x) && !is.null(x$result)) return(x$result)
    x  # otherwise assume x itself is the fitted object
  }

  # Helper: sanitize file names (no slashes, spaces, etc.)
  sanitize <- function(x) {
    x <- gsub("[/\\:*?\"<>|]", "_", x)
    x <- gsub("\\s+", "_", x)
    x
  }

  # Helper: call priors_fun robustly with or without model_id
  call_priors_fun <- function(fun, fit, model_id) {
    fmls <- try(formals(fun), silent = TRUE)
    use_id <- FALSE
    if (!inherits(fmls, "try-error") && !is.null(fmls)) {
      use_id <- "model_id" %in% names(fmls)
    }
    if (use_id) {
      fun(fit, model_id = model_id)
    } else {
      fun(fit)
    }
  }

  saved <- character(0L)

  scenarios <- names(models_all)
  if (is.null(scenarios)) scenarios <- as.character(seq_along(models_all))

  for (scenario in scenarios) {
    scen_list <- models_all[[scenario]]
    if (is.null(scen_list)) next
    if (!is.list(scen_list) || length(scen_list) == 0L) next

    # per-scenario output folder
    scen_dir <- file.path(outdir, sanitize(as.character(scenario)))
    if (!dir.exists(scen_dir)) dir.create(scen_dir, recursive = TRUE, showWarnings = FALSE)

    model_tags <- names(scen_list)
    if (is.null(model_tags)) model_tags <- as.character(seq_along(scen_list))

    message("Processing scenario: ", scenario)

    for (model_tag in model_tags) {
      leaf <- scen_list[[model_tag]]

      # skip if marked as error
      if (is.list(leaf) && !is.null(leaf$error)) {
        message(sprintf("  - Skipping %s (error in fit)", model_tag))
        next
      }

      fit <- get_fit_obj(leaf)
      if (is.null(fit)) {
        message(sprintf("  - Skipping %s (no result)", model_tag))
        next
      }

      # Build plot (robustly)
      message(sprintf("  - Saving priors for model: %s", model_tag))
      plt <- try(call_priors_fun(priors_fun, fit, model_tag), silent = TRUE)
      if (inherits(plt, "try-error") || is.null(plt)) {
        message(sprintf("    ! priors_fun failed for %s: %s", model_tag, as.character(plt)))
        next
      }

      # If the function returned a list instead of a gg object, try to wrap
      if (is.list(plt) && !inherits(plt, "gg") && !"ggplot" %in% class(plt)) {
        # try to assemble with patchwork if available
        if ("package:patchwork" %in% search() || requireNamespace("patchwork", quietly = TRUE)) {
          plt <- try(patchwork::wrap_plots(plt), silent = TRUE)
          if (inherits(plt, "try-error") || is.null(plt)) {
            message(sprintf("    ! Could not assemble list plot for %s", model_tag))
            next
          }
        } else {
          message(sprintf("    ! priors_fun returned a list and patchwork is unavailable for %s", model_tag))
          next
        }
      }

      # final path
      fn <- file.path(scen_dir,
                      paste0("priors_",
                             sanitize(as.character(scenario)), "_",
                             sanitize(as.character(model_tag)), ".png"))

      ok <- try({
        ggplot2::ggsave(filename = fn, plot = plt,
                        width = width, height = height, dpi = dpi)
      }, silent = TRUE)

      if (inherits(ok, "try-error")) {
        message(sprintf("    ! ggsave failed for %s: %s", model_tag, as.character(ok)))
        next
      } else {
        saved <- c(saved, fn)
        message(sprintf("    ✓ Saved: %s", fn))
      }
    }
  }

  invisible(saved)
}
