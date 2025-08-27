#' Exact Mean Absolute Percentage Error (MAPE)
#'
#' Implements the textbook/paper definition:
#' \deqn{\mathrm{MAPE}=\frac{100}{n}\sum_{t=1}^{n}\left|\frac{\hat y_t-y_t}{y_t}\right|}{}
#' Requires finite inputs and forbids zero observations.
#'
#' @param obs Numeric vector of observed values (finite, non-zero).
#' @param pred Numeric vector of predicted values (finite).
#'
#' @return Numeric scalar: MAPE in percent.
#' @export
mape_exact <- function(obs, pred){
  stopifnot(length(obs) == length(pred))
  if (any(!is.finite(obs)) || any(!is.finite(pred)))
    stop("Non-finite values present.")
  if (any(obs == 0))
    stop("MAPE undefined when any observed value is zero.")
  mean(abs(pred - obs) / abs(obs)) * 100
}

# ---- Build per-observation index labels -------------------------------------
#' @keywords internal
#' @noRd
.get_obs_index_labels <- function(inp){
  if (!("obsI" %in% names(inp))) stop(".get_obs_index_labels: inp$obsI missing.")
  nI   <- length(inp$obsI)
  if (!("nobsI" %in% names(inp))) stop(".get_obs_index_labels: inp$nobsI missing.")
  nobs <- as.integer(inp$nobsI)
  if (length(nobs) != nI) stop(".get_obs_index_labels: length(nobsI) != length(obsI).")
  totI <- sum(nobs)

  # Helper: validate a candidate label vector against nobsI
  check_labels <- function(lab){
    if (is.null(lab)) return(NULL)
    if (is.list(lab)) lab <- unlist(lab, use.names = FALSE)
    lab <- as.integer(lab)
    if (length(lab) != totI) return(NULL)
    if (min(lab, na.rm = TRUE) == 0L && max(lab, na.rm = TRUE) == (nI - 1L)) lab <- lab + 1L
    if (any(!is.finite(lab))) return(NULL)
    if (any(lab < 1L | lab > nI)) return(NULL)
    cnt <- tabulate(lab, nbins = nI)
    if (!identical(cnt, nobs)) return(NULL)
    lab
  }

  if ("iiin" %in% names(inp)) {
    lab <- check_labels(inp$iiin)
    if (!is.null(lab)) return(lab)
  }
  if ("ii" %in% names(inp)) {
    lab <- check_labels(inp$ii)
    if (!is.null(lab)) return(lab)
  }

  # Fallback: construct from nobsI (block labels)
  lab <- integer(totI)
  pos <- 1L
  i   <- 1L
  while (i <= nI){
    cnt <- nobs[i]
    if (cnt > 0L) {
      lab[pos:(pos + cnt - 1L)] <- i
      pos <- pos + cnt
    }
    i <- i + 1L
  }
  lab
}

# ---- Map helper: observation -> prediction indices ---------------------------
#' @keywords internal
#' @noRd
.get_obs_pred_map <- function(inp){
  if (!("iqin" %in% names(inp))) return(NULL)
  iq <- inp$iqin
  if (is.null(iq)) return(NULL)
  if (is.list(iq)) iq <- unlist(iq, use.names = FALSE)
  iq <- as.integer(iq)
  if (!length(iq)) return(NULL)
  if (any(!is.finite(iq))) return(NULL)
  if (min(iq, na.rm = TRUE) == 0L) iq <- iq + 1L
  iq
}

# ---- Safe subsetting by iq (grid -> obs) -------------------------------------
#' @keywords internal
#' @noRd
.subset_by_iq <- function(v, iq, totI, is_log = FALSE){
  if (is.matrix(v) || is.data.frame(v)) v <- as.numeric(v)
  if (!is.numeric(v)) stop(".subset_by_iq: 'v' must be numeric after coercion.")
  if (length(iq) != totI)
    warning(".subset_by_iq: length(iq) != sum(nobsI); proceeding by available mapping.")
  if (any(!is.finite(iq))) stop(".subset_by_iq: non-finite indices in iq.")
  if (max(iq, na.rm = TRUE) > length(v) || min(iq, na.rm = TRUE) < 1L)
    stop(".subset_by_iq: iq contains out-of-range indices for the candidate vector.")
  res <- v[iq]
  if (is_log) res <- exp(as.numeric(res)) else res <- as.numeric(res)
  res
}

# ---- Expand per-index vector (length = nI) to per-observation ----------------
#' @keywords internal
#' @noRd
.expand_per_index_to_obs <- function(v, inp, is_log = FALSE){
  nI   <- length(inp$obsI)
  nobs <- as.integer(inp$nobsI)
  if (length(v) != nI) stop(".expand_per_index_to_obs: length(v) != number of indices.")
  if (is_log) v <- exp(as.numeric(v)) else v <- as.numeric(v)
  out <- numeric(sum(nobs))
  pos <- 1L
  i   <- 1L
  while (i <= nI){
    cnt <- nobs[i]
    if (cnt > 0L){
      out[pos:(pos + cnt - 1L)] <- v[i]
      pos <- pos + cnt
    }
    i <- i + 1L
  }
  out
}

# ---- Locate predicted INDEX vector in fit$report -----------------------------
#' @keywords internal
#' @noRd
.locate_pred_I_vector <- function(fit, cand_pred_names = NULL){
  if (is.null(fit$report)) stop(".locate_pred_I_vector: fit$report is NULL.")
  if (is.null(fit$inp))    stop(".locate_pred_I_vector: fit$inp is NULL.")

  totI <- sum(fit$inp$nobsI)
  rep  <- fit$report
  nm   <- names(rep); if (is.null(nm)) nm <- character()
  iq   <- .get_obs_pred_map(fit$inp)  # may be NULL

  # 0) OSAR often contains Ihat exactly at observation times
  if ("osar" %in% names(fit) && "Ihat" %in% names(fit$osar)) {
    v <- fit$osar$Ihat
    if (is.numeric(v) && length(v) == totI) return(as.numeric(v))
    if (is.numeric(v) && !is.null(iq) && max(iq, na.rm = TRUE) <= length(v))
      return(.subset_by_iq(v, iq, totI, is_log = FALSE))
  }

  # 1) Preferred explicit candidates (try 'logIp' first given your build)
  if (is.null(cand_pred_names)) {
    cand_pred_names <- c(
      "logIp", "logIhat", "logIpred", "logIp_obs", "logIhat_obs", "logIpred_obs",
      "Ip", "Ihat", "Ipred", "Ip_obs", "Ihat_obs", "Ipred_obs"
    )
  }

  i <- 1L
  while (i <= length(cand_pred_names)){
    cn <- cand_pred_names[i]
    if (cn %in% nm) {
      v <- rep[[cn]]
      if (is.matrix(v) || is.data.frame(v)) v <- as.numeric(v)
      if (is.numeric(v)) {
        if (length(v) == totI) {
          return(if (startsWith(cn, "log")) exp(as.numeric(v)) else as.numeric(v))
        }
        if (!is.null(iq) && max(iq, na.rm = TRUE) <= length(v)) {
          return(.subset_by_iq(v, iq, totI, is_log = startsWith(cn, "log")))
        }
        nI <- length(fit$inp$obsI)
        if (length(v) == nI) {
          return(.expand_per_index_to_obs(v, fit$inp, is_log = startsWith(cn, "log")))
        }
      }
    }
    i <- i + 1L
  }

  # 2) Pattern fallback (avoid Cp/Ep etc.)
  likeI <- nm[grepl("^(log)?I(p|hat|pred)(_|$)|^(log)?I$|^Ip$", nm, perl = TRUE)]
  likeI <- setdiff(likeI, nm[grepl("^(C|E)p", nm)])  # drop Cp/Ep if present

  i <- 1L
  while (i <= length(likeI)){
    cn <- likeI[i]
    v <- rep[[cn]]
    if (is.matrix(v) || is.data.frame(v)) v <- as.numeric(v)
    if (is.numeric(v)) {
      if (length(v) == totI) {
        return(if (startsWith(cn, "log")) exp(as.numeric(v)) else as.numeric(v))
      }
      if (!is.null(iq) && max(iq, na.rm = TRUE) <= length(v)) {
        return(.subset_by_iq(v, iq, totI, is_log = startsWith(cn, "log")))
      }
      nI <- length(fit$inp$obsI)
      if (length(v) == nI) {
        return(.expand_per_index_to_obs(v, fit$inp, is_log = startsWith(cn, "log")))
      }
    }
    i <- i + 1L
  }

  stop(
    "Could not locate predicted index in fit$report.\n",
    "Tried: per-observation, iq-mapped, and per-index expansion for Ip/Ihat/Ipred/logIp ",
    "(obs len=", totI, "). Available (first 30): ",
    paste(utils::head(names(rep), 30), collapse = ", "), " ..."
  )
}

# ---- Public helper: list candidates of correct length ------------------------
#' Peek report items with the right length for index predictions
#'
#' Lists report entries whose lengths match \code{sum(nobsI)}.
#'
#' @param fit A fitted elu2/SPiCT object.
#' @return A data.frame with columns \code{name}, \code{length}.
#' @export
peek_index_candidates <- function(fit){
  totI <- sum(fit$inp$nobsI)
  nm <- names(fit$report)
  lens <- vapply(fit$report, function(x) if (is.numeric(x)) length(x) else NA_integer_, 1L)
  hits <- which(lens == totI)
  data.frame(name = nm[hits], length = totI, stringsAsFactors = FALSE)
}

# ---- Extract observed vs predicted for index k -------------------------------
#' @keywords internal
#' @noRd
extract_obs_pred_index_elu2 <- function(fit, idx, cand_pred_names = NULL){
  if (idx < 1L || idx > length(fit$inp$obsI)) stop("extract_obs_pred_index_elu2: invalid 'idx'.")
  pred_all <- .locate_pred_I_vector(fit, cand_pred_names = cand_pred_names)
  lab <- .get_obs_index_labels(fit$inp)
  if (length(lab) != length(pred_all))
    stop("extract_obs_pred_index_elu2: length(iiin/ii/nobsI labels) != length(pred_all).")

  sel  <- lab == idx
  obs  <- fit$inp$obsI[[idx]]
  yrs  <- fit$inp$timeI[[idx]]
  pred <- pred_all[sel]

  if (length(obs) != length(pred))
    stop("Index ", idx, ": obs and pred lengths differ (",
         length(obs), " vs ", length(pred), ").")

  data.frame(
    year  = as.numeric(yrs),
    obs   = as.numeric(obs),
    pred  = as.numeric(pred),
    index = as.integer(idx),
    stringsAsFactors = FALSE
  )
}

# ---- Per-fit MAPE (per-index + pooled) --------------------------------------
#' @keywords internal
#' @noRd
mape_for_fit_elu2 <- function(fit, model_name = NA_character_,
                              scenario   = NA_character_,
                              type       = NA_character_,
                              cand_pred_names = NULL){
  nI <- length(fit$inp$obsI)
  rows <- vector("list", nI + 1L)
  pooled_num <- 0
  pooled_den <- 0

  k <- 1L
  while (k <= nI){
    d <- extract_obs_pred_index_elu2(fit, k, cand_pred_names = cand_pred_names)

    # Exact MAPE per index (errors if zeros or non-finite present)
    Ms <- mape_exact(d$obs, d$pred)

    rows[[k]] <- data.frame(
      scenario = scenario,
      model    = model_name,
      type     = type,
      index    = k,
      n        = length(d$obs),
      MAPE     = Ms,
      stringsAsFactors = FALSE
    )

    # Exact pooled accumulation
    pooled_num <- pooled_num + sum(abs(d$pred - d$obs) / abs(d$obs))
    pooled_den <- pooled_den + length(d$obs)
    k <- k + 1L
  }

  pooled_mape <- if (pooled_den > 0) 100 * (pooled_num / pooled_den) else NA_real_
  rows[[nI + 1L]] <- data.frame(
    scenario = scenario,
    model    = model_name,
    type     = type,
    index    = "AllIndices",
    n        = pooled_den,
    MAPE     = pooled_mape,
    stringsAsFactors = FALSE
  )

  do.call(rbind, rows)
}

# ---- Main entry: nested all_models walker -----------------------------------
#' Compute MAPE table for nested elu2 model lists by scenario (exact MAPE)
#'
#' Walks the usual \code{all_models$Sx$SxP.SDM} structure, computing per-index
#' and pooled MAPE (exact definition) for each fitted model.
#'
#' @param all_models Nested list like \code{all_models$S1$S1P.SDM}, ...
#' @param cand_pred_names Optional character vector of report names to try first.
#'
#' @return data.frame with columns: \code{scenario}, \code{model}, \code{type},
#'   \code{index}, \code{n}, \code{MAPE}.
#' @export
mape_by_scenario <- function(all_models, cand_pred_names = NULL){
  if (!is.list(all_models)) stop("mape_by_scenario: 'all_models' must be a list.")
  out <- list(); oi <- 1L
  scen_names <- names(all_models)
  if (is.null(scen_names)) scen_names <- rep(NA_character_, length(all_models))

  si <- 1L
  while (si <= length(all_models)){
    scen <- all_models[[si]]
    scen_label_outer <- scen_names[si]

    if (is.list(scen)){
      inner_names <- names(scen); if (is.null(inner_names)) inner_names <- rep(NA_character_, length(scen))
      j <- 1L
      while (j <= length(scen)){
        fit <- scen[[j]]
        if (is.null(fit) || !is.list(fit) || is.null(fit$report)) { j <- j + 1L; next }

        nm <- inner_names[j]
        scenario <- if (!is.na(scen_label_outer) && grepl("^S\\d+$", scen_label_outer)) scen_label_outer else NA_character_
        base     <- if (!is.na(nm)) sub("\\..*$", "", nm) else NA_character_
        if (is.na(scenario) && !is.na(base) && grepl("^S\\d+", base)) scenario <- sub("^((S\\d+)).*$", "\\1", base)
        model <- if (!is.na(base)) sub("^S\\d+", "", base) else NA_character_
        type  <- if (!is.na(nm) && grepl("\\.", nm)) sub("^.*\\.", "", nm) else NA_character_
        if (is.na(scenario)) scenario <- "Unknown"
        if (is.na(model))    model    <- "Unknown"
        if (is.na(type))     type     <- "Unknown"

        out[[oi]] <- mape_for_fit_elu2(
          fit,
          model_name = model,
          scenario   = scenario,
          type       = type,
          cand_pred_names = cand_pred_names
        )
        oi <- oi + 1L
        j  <- j + 1L
      }
    } else {
      # Flat list fallback
      fit <- scen
      nm  <- scen_label_outer
      scenario <- if (!is.na(nm) && grepl("^S\\d+", nm)) sub("^((S\\d+)).*$", "\\1", nm) else "Unknown"
      base     <- if (!is.na(nm)) sub("\\..*$", "", nm) else NA_character_
      model    <- if (!is.na(base)) sub("^S\\d+", "", base) else "Unknown"
      type     <- if (!is.na(nm) && grepl("\\.", nm)) sub("^.*\\.", "", nm) else "Unknown"

      out[[oi]] <- mape_for_fit_elu2(
        fit,
        model_name = model,
        scenario   = scenario,
        type       = type,
        cand_pred_names = cand_pred_names
      )
      oi <- oi + 1L
    }

    si <- si + 1L
  }

  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}
