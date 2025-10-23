#' Extract K, r, MSY, B/Bmsy, F/Fmsy (+ CIs) from SDM+GLM SPiCT models
#'
#' Builds a tidy table of selected parameters and relative status metrics for
#' **both** SDM and GLM model types across scenarios and Fox / Schaefer / Pella
#' families. This wrapper relies on your helper `make_spict_summary_df1_nopipe()`
#' to produce wide per-scenario summaries, then reshapes and filters them.
#'
#' @section Required input layout:
#' The function expects a **named** list `model_list` whose names follow
#' `"<scenario><family>.<type>"`, for example:
#' `"S1F.SDM"`, `"S1S.GLM"`, `"S1P.SDM"`, … where
#' - `scenario = S\\d+` (e.g., `"S1"`, `"S2"`, …),
#' - `family ∈ {F,S,P}` → Fox, Schaefer, Pella,
#' - `type ∈ {SDM, GLM}`.
#'
#' Your helper `make_spict_summary_df1_nopipe(models_triplet, model_names)` should
#' return a **wide** data frame with columns like:
#' `Scenario`, `Parameters`,
#' `Fox_estimate`, `Fox_low`, `Fox_up`,
#' `Schaefer_estimate`, `Schaefer_low`, `Schaefer_up`,
#' `Pella_estimate`, `Pella_low`, `Pella_up`.
#'
#' @details
#' Parameter labels in the wide table are parsed and mapped to the keys
#' `"K"`, `"r"`, `"MSY"`, `"B/Bmsy"`, `"F/Fmsy"`. If a ratio label embeds a year
#' (e.g., `$B_{2022}/B_{MSY}$`), that year is extracted into the `year` column;
#' otherwise `year` is `NA_integer_`. Non-target parameters are discarded.
#'
#' Scenarios are auto-detected from `names(model_list)` when `scenarios` is `NULL`.
#' Missing triplet members (e.g., no `"S1F.GLM"`) are skipped silently.
#'
#' @param model_list A **named** list of SPiCT fits named like
#'   `"<scenario><family>.<type>"` (e.g., `"S3F.SDM"`).
#' @param scenarios Optional character vector of scenarios to include
#'   (e.g., `c("S1","S2")`). If `NULL`, scenarios are inferred from names.
#' @param verbose Logical; print progress by scenario? Default `TRUE`.
#'
#' @return
#' A tidy `data.frame` with columns:
#' - `scenario` — e.g., `"S1"`,
#' - `type` — `"SDM"` or `"GLM"`,
#' - `family` — `"Fox"`, `"Schaefer"`, or `"Pella"`,
#' - `param` — one of `"K"`, `"r"`, `"MSY"`, `"B/Bmsy"`, `"F/Fmsy"`,
#' - `year` — embedded year for ratio labels if present, else `NA_integer_`,
#' - `est`, `lo`, `hi` — estimate and bounds.
#'
#' Factors are set to:
#' `family = c("Fox","Schaefer","Pella")`,
#' `type   = c("SDM","GLM")`,
#' `param  = c("K","r","MSY","B/Bmsy","F/Fmsy")`.
#'
#' @examples
#' \dontrun{
#' # Example structure:
#' fits <- list(
#'   S1F.SDM = fit_fox_sdm,     S1S.SDM = fit_schaefer_sdm, S1P.SDM = fit_pella_sdm,
#'   S1F.GLM = fit_fox_glm,     S1S.GLM = fit_schaefer_glm, S1P.GLM = fit_pella_glm,
#'   S2F.SDM = fit2_fox_sdm     # etc…
#' )
#'
#' # Ensure your helper is available:
#' # source("make_spict_summary_df1_nopipe.R")
#'
#' df <- extract_params_SDM_GLM(fits, scenarios = c("S1","S2"))
#' subset(df, param %in% c("B/Bmsy","F/Fmsy"))
#' }
#'
#' @seealso
#'  \code{\link[=compute_param_pct_change]{compute_param_pct_change}} for
#'  computing % changes vs. a base scenario from this tidy output.
#'
#' @export
extract_params_SDM_GLM555 <- function(model_list, scenarios = NULL, verbose = TRUE) {
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
