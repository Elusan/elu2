#' Select models by WAIC within each (scenario × type)
#'
#' @description
#' Computes WAIC-based model selection for your nested structure
#' `all_models$Sx$Sx{P|S|F}.{SDM|GLM}`. It:
#'
#' * extracts `WAIC_plug_total` robustly (from `model$report$WAIC_plug_total`
#'   or `model$report()$WAIC_plug_total`);
#' * keeps only converged fits (`opt$convergence == 0`) with finite WAIC;
#' * computes ΔWAIC and Akaike-like weights **within each (scenario × type)**;
#' * identifies the best model per (scenario × type);
#' * returns a wide table with WAIC/ΔWAIC/weights per family (Pella, Schaefer, Fox);
#' * emits notes when data counts differ across models (compare WAIC with caution).
#'
#' @param all_models Nested list like `all_models$S1$S1P.SDM`, `all_models$S1$S1F.GLM`, ...
#'
#' @return A list with:
#' \describe{
#'   \item{`full`}{Long data.frame of all retained fits with columns:
#'     `scenario`, `set_name`, `type`, `family`, `waic`, `convergence`,
#'     `nC`, `nI1`, `nI2`, `delta_st`, `weight_st`.}
#'   \item{`best`}{One row per (scenario × type) — the smallest WAIC.}
#'   \item{`wide`}{One row per (scenario × type); columns `WAIC_*`, `dWAIC_*`, `w_*`
#'     for families Pella, Schaefer, Fox.}
#'   \item{`notes`}{Character vector of diagnostic messages (possibly length 0).}
#' }
#'
#' @examples
#' \dontrun{
#'   sel <- select_models_by_WAIC(all_models)
#'   sel$best            # winner per (scenario × type)
#'   sel$wide            # compact comparison table
#'   sel$full            # full long table
#'   sel$notes           # any comparability warnings
#' }
#'
#' @export
select_models_by_WAIC <- function(all_models) {
  # Extract WAIC robustly (works if $report is a list or a function)
  extract_waic <- function(o) {
    v <- NA_real_
    if (!is.null(o$report) && !is.null(o$report$WAIC_plug_total)) {
      v <- o$report$WAIC_plug_total
    } else if (is.function(o$report)) {
      rr <- o$report()
      if (!is.null(rr$WAIC_plug_total)) v <- rr$WAIC_plug_total
    }
    v
  }

  # Family from code letter
  fam_from_code <- function(code) {
    if (code == "P") "Pella" else if (code == "S") "Schaefer" else if (code == "F") "Fox" else NA_character_
  }

  rows <- list(); k <- 0L
  scen_names <- names(all_models)

  for (i in seq_along(all_models)) {
    sc <- scen_names[i]                     # "S1", "S2", ...
    models_i <- all_models[[i]]
    mdl_names <- names(models_i)

    for (j in seq_along(models_i)) {
      mdl <- models_i[[j]]
      nm  <- mdl_names[j]                   # e.g., "S1P.SDM"

      parts <- strsplit(nm, "\\.")[[1]]
      short <- if (length(parts) >= 1) parts[1] else NA_character_   # "S1P"
      type  <- if (length(parts) >= 2) parts[2] else NA_character_   # "SDM" or "GLM"

      fam_cd <- if (!is.na(short)) substring(short, nchar(short), nchar(short)) else NA_character_
      family <- fam_from_code(fam_cd)

      waic <- extract_waic(mdl)
      conv <- if (!is.null(mdl$opt$convergence)) mdl$opt$convergence else NA_integer_

      # Optional comparability checks (counts)
      nC  <- if (!is.null(mdl$inp$timeC)) length(mdl$inp$timeC) else NA_integer_
      nI1 <- NA_integer_; nI2 <- NA_integer_
      if (!is.null(mdl$inp$obsI)) {
        if (is.list(mdl$inp$obsI)) {
          if (length(mdl$inp$obsI) >= 1) nI1 <- length(mdl$inp$obsI[[1]])
          if (length(mdl$inp$obsI) >= 2) nI2 <- length(mdl$inp$obsI[[2]])
        } else {
          nI1 <- length(mdl$inp$obsI)
        }
      }

      k <- k + 1L
      rows[[k]] <- data.frame(
        scenario = sc,
        set_name = nm,        # "S1P.SDM"
        type     = type,      # "SDM"/"GLM"
        family   = family,    # "Pella"/"Fox"/"Schaefer"
        waic     = waic,
        convergence = conv,
        nC = nC, nI1 = nI1, nI2 = nI2,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) == 0L) {
    return(list(full = data.frame(), best = data.frame(), wide = data.frame(), notes = character(0)))
  }

  df <- do.call(rbind, rows)

  # Keep only converged + finite WAIC
  ok <- is.finite(df$waic) & (is.na(df$convergence) | df$convergence == 0)
  df_ok <- df[ok, , drop = FALSE]
  if (!nrow(df_ok)) {
    return(list(full = df_ok, best = data.frame(), wide = data.frame(), notes = "No converged models with finite WAIC."))
  }

  # Compute ΔWAIC and weights within each (scenario × type)
  df_ok$delta_st <- NA_real_
  df_ok$weight_st <- NA_real_

  key <- paste(df_ok$scenario, df_ok$type, sep = "||")
  ukey <- unique(key)
  for (kk in ukey) {
    idx <- which(key == kk)
    wv  <- df_ok$waic[idx]
    d   <- wv - min(wv)
    ww  <- exp(-0.5 * d); ww <- ww / sum(ww)
    df_ok$delta_st[idx]  <- d
    df_ok$weight_st[idx] <- ww
  }

  # Winner per (scenario × type)
  winner_list <- list(); h <- 0L
  for (kk in ukey) {
    idx <- which(key == kk)
    sub <- df_ok[idx, , drop = FALSE]
    sub <- sub[order(sub$waic), ]
    h <- h + 1L
    winner_list[[h]] <- data.frame(
      scenario = sub$scenario[1],
      type     = sub$type[1],
      best_set = sub$set_name[1],
      best_family = sub$family[1],
      waic     = sub$waic[1],
      delta    = sub$delta_st[1],
      weight   = sub$weight_st[1],
      stringsAsFactors = FALSE
    )
  }
  best_by_set <- do.call(rbind, winner_list)
  best_by_set <- best_by_set[order(best_by_set$scenario, best_by_set$type), ]

  # Wide table per (scenario × type): WAIC_* and weight_* per family
  fams <- c("Pella", "Schaefer", "Fox")  # expected 3 models
  make_row <- function(sc, ty) {
    idx <- which(df_ok$scenario == sc & df_ok$type == ty)
    sub <- df_ok[idx, , drop = FALSE]
    out <- data.frame(scenario = sc, type = ty, stringsAsFactors = FALSE)
    for (f in fams) {
      ii <- which(sub$family == f)
      out[[paste0("WAIC_", f)]]   <- if (length(ii) == 1) sub$waic[ii] else NA_real_
      out[[paste0("dWAIC_", f)]]  <- if (length(ii) == 1) sub$delta_st[ii] else NA_real_
      out[[paste0("w_", f)]]      <- if (length(ii) == 1) sub$weight_st[ii] else NA_real_
    }
    out
  }

  wide_rows <- list(); wj <- 0L
  for (sc in unique(df_ok$scenario)) {
    for (ty in unique(df_ok$type[df_ok$scenario == sc])) {
      wj <- wj + 1L
      wide_rows[[wj]] <- make_row(sc, ty)
    }
  }
  wide <- do.call(rbind, wide_rows)
  wide <- wide[order(wide$scenario, wide$type), ]

  # Diagnostics: flag sets with unequal data counts
  messages <- character(0)
  for (kk in ukey) {
    idx <- which(key == kk)
    sub <- df_ok[idx, , drop = FALSE]
    if (length(unique(sub$nC))  > 1 ||
        length(unique(sub$nI1)) > 1 ||
        length(unique(sub$nI2)) > 1) {
      parts <- strsplit(kk, "\\|\\|")[[1]]
      messages <- c(messages, paste("Note:", parts[1], parts[2], "have unequal data counts across models; compare WAIC with caution."))
    }
  }

  # Nice sort for full table
  df_ok <- df_ok[order(df_ok$scenario, df_ok$type, df_ok$waic), ]

  list(
    full   = df_ok,        # long table with WAIC, Δ and weights per (scenario×type)
    best   = best_by_set,  # winner per (scenario×type)
    wide   = wide,         # one row per (scenario×type), columns for each family
    notes  = messages
  )
}
