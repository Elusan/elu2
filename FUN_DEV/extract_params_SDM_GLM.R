#' Extract K, r, MSY, B/Bmsy, F/Fmsy (with CIs) from SDM+GLM SPiCT Models
#'
#' Builds a tidy data frame of selected parameters and their 95% confidence
#' intervals for **both** SDM and GLM model types across scenarios and Fox /
#' Schaefer / Pella families. The function relies on your helper
#' \code{make_spict_summary_df1_nopipe()} to create per-scenario summary tables
#' and then reshapes and filters to the requested parameters.
#'
#' @section Required input format:
#' \itemize{
#'   \item \strong{Named} \code{model_list} with names that match
#'         \code{"<scenario><family>.<type>"}; e.g.:
#'         \code{"S1F.SDM"}, \code{"S1S.GLM"}, \code{"S1P.SDM"}, … where
#'         \code{scenario = S\\d+}, \code{family in \{F,S,P\}} (Fox, Schaefer, Pella),
#'         and \code{type in \{SDM, GLM\}}.
#'   \item \code{make_spict_summary_df1_nopipe()} must be available in the search path.
#'         It should return a wide table with columns:
#'         \code{Scenario}, \code{Parameters},
#'         \code{Fox_estimate, Fox_low, Fox_up},
#'         \code{Schaefer_estimate, Schaefer_low, Schaefer_up},
#'         \code{Pella_estimate, Pella_low, Pella_up}.
#' }
#'
#' @details
#' The function parses LaTeX-style parameter labels produced by
#' \code{make_spict_summary_df1_nopipe()} and maps them to the keys
#' \code{"K"}, \code{"r"}, \code{"MSY"}, \code{"B/Bmsy"}, \code{"F/Fmsy"}.
#' If a ratio label embeds a year (e.g., \code{$B_{2022}/B_{MSY}$}), the year
#' is extracted into the \code{year} column; otherwise \code{year} is left
#' as \code{NA_integer_}. Parameters not in the requested set are discarded.
#'
#' @param model_list A \strong{named} flat list of SPiCT fits whose names follow
#'   the pattern \code{<scenario><family>.<type>} (e.g., \code{"S3F.SDM"}).
#' @param scenarios Optional character vector of scenarios to include
#'   (e.g., \code{c("S1","S2")}). If \code{NULL}, scenarios are auto-detected
#'   from \code{names(model_list)}.
#' @param verbose Logical; print progress messages per scenario. Default \code{TRUE}.
#'
#' @return A tidy \code{data.frame} with columns:
#' \describe{
#'   \item{scenario}{Scenario ID (e.g., \code{"S1"}).}
#'   \item{type}{Model type (\code{"SDM"} or \code{"GLM"}).}
#'   \item{family}{\code{"Fox"}, \code{"Schaefer"}, or \code{"Pella"}.}
#'   \item{param}{One of \code{"K"}, \code{"r"}, \code{"MSY"}, \code{"B/Bmsy"}, \code{"F/Fmsy"}.}
#'   \item{year}{Year embedded in a ratio label (if any), otherwise \code{NA_integer_}.}
#'   \item{est}{Point estimate (numeric).}
#'   \item{lo}{Lower 95\% CI (numeric).}
#'   \item{hi}{Upper 95\% CI (numeric).}
#' }
#' Factor levels are set to \code{family = c("Fox","Schaefer","Pella")},
#' \code{type = c("SDM","GLM")}, and
#' \code{param = c("K","r","MSY","B/Bmsy","F/Fmsy")}.
#'
#' @examples
#' \dontrun{
#' # Sensitivity_SDM_GLM_models <- list(
#' #   S1P.SDM = S1_Pella.SDM, S1S.SDM = S1_Schaefer.SDM, S1F.SDM = S1_Fox.SDM,
#' #   S1P.GLM = S1_Pella.GLM,  S1S.GLM = S1_Schaefer.GLM,  S1F.GLM = S1_Fox.GLM,
#' #   ... # up to S7*
#' # )
#' #
#' # Ensure make_spict_summary_df1_nopipe() is defined/loaded:
#' # source("make_spict_summary_df1_nopipe.R")
#' #
#' # Extract tidy table:
#' tab <- extract_params_SDM_GLM(Sensitivity_SDM_GLM_models, scenarios = paste0("S", 1:7))
#' head(tab)
#' }
#'
#' @seealso \code{\link{make_spict_summary_df1_nopipe}}
#' @export
extract_params_SDM_GLM <- function(model_list, scenarios = NULL, verbose = TRUE) {
  if (is.null(names(model_list))) stop("`model_list` must be a named list.")

  # ---- helpers (aligned to your naming/layout) ----
  family_from_code <- function(code) {
    z <- toupper(trimws(as.character(code)))
    if (nchar(z) > 1L) z <- substr(z, 1L, 1L)
    if (z == "F") "Fox" else if (z == "S") "Schaefer" else if (z == "P") "Pella" else NA_character_
  }
  parse_name <- function(nm) {
    base <- sub("\\..*$", "", nm)               # "S1F"
    scen <- sub("^((S\\d+)).*$", "\\1", base)   # "S1"
    type <- toupper(sub("^.*\\.", "", nm))      # "SDM"/"GLM"
    code <- substr(base, nchar(base), nchar(base))
    list(scenario = scen, type = type, family = family_from_code(code))
  }

  # Parameter label normalizer -> (key, year)
  normalize_param <- function(label) {
    if (is.null(label) || is.na(label)) return(list(key = NA_character_, year = NA_integer_))
    lab <- gsub("^\\s+|\\s+$", "", as.character(label))

    if (lab == "K" || lab == "$K$")  return(list(key = "K", year = NA_integer_))
    if (lab == "r" || lab == "$r$")  return(list(key = "r", year = NA_integer_))

    # MSY: ensure it's not Bmsy or Fmsy ratio (no slash, has MSY)
    if (grepl("MSY", lab) && !grepl("/", lab) &&
        !grepl("B[_\\{].*MSY", lab) && !grepl("F[_\\{].*MSY", lab)) {
      return(list(key = "MSY", year = NA_integer_))
    }

    # B/Bmsy (optionally with year: B_{2022}/B_{MSY})
    if (grepl("B.*/.*MSY", lab)) {
      yr <- NA_integer_
      m <- regexpr("\\{(\\d{4})\\}", lab, perl = TRUE)
      if (m > 0) {
        yrtxt <- regmatches(lab, m)
        yrnum <- gsub("[^0-9]", "", yrtxt)
        if (nchar(yrnum)) yr <- as.integer(yrnum)
      }
      return(list(key = "B/Bmsy", year = yr))
    }

    # F/Fmsy (optionally with year)
    if (grepl("F.*/.*MSY", lab)) {
      yr <- NA_integer_
      m <- regexpr("\\{(\\d{4})\\}", lab, perl = TRUE)
      if (m > 0) {
        yrtxt <- regmatches(lab, m)
        yrnum <- gsub("[^0-9]", "", yrtxt)
        if (nchar(yrnum)) yr <- as.integer(yrnum)
      }
      return(list(key = "F/Fmsy", year = yr))
    }

    list(key = NA_character_, year = NA_integer_)
  }

  # ---------- Auto-detect scenarios from names(model_list) ----------
  if (is.null(scenarios) || !length(scenarios)) {
    scen_raw  <- gsub("^((S\\d+)).*$", "\\1", sub("\\..*$", "", names(model_list)))
    scen_keep <- scen_raw[grepl("^S\\d+$", scen_raw)]
    scenarios <- unique(scen_keep)
    scenarios <- scenarios[order(as.integer(sub("^S(\\d+)$", "\\1", scenarios)))]
    if (!length(scenarios)) stop("Could not auto-detect scenarios.")
  }

  wanted_keys <- c("K","r","MSY","B/Bmsy","F/Fmsy")
  out <- list()

  for (scen in scenarios) {
    if (isTRUE(verbose)) message("Extracting: ", scen)

    # Fetch a model by family code + type
    get_mod <- function(fam_code, type) {
      nm <- paste0(scen, fam_code, ".", type) # e.g., "S1F.SDM"
      if (!nm %in% names(model_list)) return(NULL)
      model_list[[nm]]
    }

    # Triplets (Fox, Schaefer, Pella) for SDM & GLM
    models_SDM <- list(Fox = get_mod("F", "SDM"),
                       Schaefer = get_mod("S", "SDM"),
                       Pella = get_mod("P", "SDM"))
    models_GLM <- list(Fox = get_mod("F", "GLM"),
                       Schaefer = get_mod("S", "GLM"),
                       Pella = get_mod("P", "GLM"))

    have_any_sdm <- any(vapply(models_SDM, function(x) !is.null(x), FALSE))
    have_any_glm <- any(vapply(models_GLM, function(x) !is.null(x), FALSE))
    if (!have_any_sdm && !have_any_glm) next

    # ---- Call YOUR helper with scenario-aware model_names ----
    pivot_triplet_to_long <- function(models_triplet, scenario_id) {
      fams <- c("Fox","Schaefer","Pella")
      # model_names must be like "S1_Fox", "S1_Schaefer", "S1_Pella"
      model_names <- paste0(scenario_id, "_", fams)

      dfw <- make_spict_summary_df1_nopipe(models_triplet, model_names)
      # dfw columns: Scenario, Parameters, Fox_estimate/low/up, Schaefer_*, Pella_*

      # Long format (without pipes)
      rows <- list(); k <- 0L
      for (i in seq_len(nrow(dfw))) {
        par_lab <- as.character(dfw$Parameters[i])
        for (fam in fams) {
          e_col <- paste0(fam, "_estimate")
          l_col <- paste0(fam, "_low")
          u_col <- paste0(fam, "_up")
          if (!all(c(e_col, l_col, u_col) %in% names(dfw))) next
          est <- dfw[i, e_col]; lo <- dfw[i, l_col]; hi <- dfw[i, u_col]
          if (is.na(est) && is.na(lo) && is.na(hi)) next
          k <- k + 1L
          rows[[k]] <- data.frame(
            Scenario  = as.character(dfw$Scenario[i]),
            Parameters= par_lab,
            family    = fam,
            est       = as.numeric(est),
            lo        = as.numeric(lo),
            hi        = as.numeric(hi),
            stringsAsFactors = FALSE
          )
        }
      }
      if (k == 0L) {
        return(data.frame(Scenario=character(0), Parameters=character(0),
                          family=character(0), est=numeric(0), lo=numeric(0), hi=numeric(0)))
      }
      do.call(rbind, rows)
    }

    long_SDM <- if (have_any_sdm) pivot_triplet_to_long(models_SDM, scen) else NULL
    long_GLM <- if (have_any_glm) pivot_triplet_to_long(models_GLM, scen) else NULL

    # Attach type, normalize param labels, keep wanted
    add_meta_and_filter <- function(dd, type_lab) {
      if (is.null(dd) || !nrow(dd)) {
        return(data.frame(scenario=character(0), type=character(0), family=character(0),
                          param=character(0), year=integer(0), est=numeric(0), lo=numeric(0), hi=numeric(0)))
      }
      n <- nrow(dd)
      param_key <- character(n); param_year <- rep(NA_integer_, n)
      for (i in seq_len(n)) {
        inf <- normalize_param(dd$Parameters[i])
        param_key[i]  <- inf$key
        param_year[i] <- inf$year
      }
      keep <- param_key %in% wanted_keys
      if (!any(keep)) {
        return(data.frame(scenario=character(0), type=character(0), family=character(0),
                          param=character(0), year=integer(0), est=numeric(0), lo=numeric(0), hi=numeric(0)))
      }
      data.frame(
        scenario = as.character(dd$Scenario[keep]),
        type     = rep(type_lab, sum(keep)),
        family   = dd$family[keep],
        param    = param_key[keep],
        year     = param_year[keep],
        est      = dd$est[keep],
        lo       = dd$lo[keep],
        hi       = dd$hi[keep],
        stringsAsFactors = FALSE
      )
    }

    part_sdm <- add_meta_and_filter(long_SDM, "SDM")
    part_glm <- add_meta_and_filter(long_GLM, "GLM")

    chunk <- rbind(part_sdm, part_glm)
    if (nrow(chunk)) out[[length(out) + 1L]] <- chunk
  }

  if (!length(out)) {
    if (isTRUE(verbose)) message("No parameters extracted.")
    return(data.frame(scenario=character(0), type=character(0), family=character(0),
                      param=character(0), year=integer(0), est=numeric(0), lo=numeric(0), hi=numeric(0)))
  }

  res <- do.call(rbind, out)
  res$family <- factor(res$family, levels = c("Fox","Schaefer","Pella"))
  res$type   <- factor(res$type, levels = c("SDM","GLM"))
  res$param  <- factor(res$param, levels = c("K","r","MSY","B/Bmsy","F/Fmsy"))
  rownames(res) <- NULL
  res
}
