#' Compute % changes vs base scenario (default "S4")
#' @param param_df ...
#' @param base_scenario ...
#' @param keep_ratio_year ...
#' @param include_ci ...
#' @param families_order ...
#' @param params_order ...
#' @return ...
#' @export
# add export above (normal comment)
compute_param_pct_change <- function(
    param_df,
    base_scenario   = "S4",
    keep_ratio_year = 2022,
    include_ci      = FALSE,
    families_order  = c("Pella","Fox","Schaefer"),
    params_order    = c("K","r","MSY","B/Bmsy","F/Fmsy")
) {
  # ---- checks
  # Ensure we have the minimum columns to work with.
  need <- c("scenario","type","family","param","est")
  if (is.null(param_df) || !is.data.frame(param_df) || !all(need %in% names(param_df))) {
    stop("`param_df` must contain columns: scenario, type, family, param, est (lo, hi and year are recommended).", call. = FALSE)
  }

  # ---- normalize types / add missing columns
  # Convert to character/factor-friendly columns and ensure `lo`/`hi` exist.
  df <- param_df
  df$scenario <- as.character(df$scenario)
  df$type     <- as.character(df$type)
  df$family   <- as.character(df$family)
  df$param    <- as.character(df$param)
  if ("year" %in% names(df)) df$year <- suppressWarnings(as.integer(df$year))
  if (!("lo" %in% names(df))) df$lo <- NA_real_
  if (!("hi" %in% names(df))) df$hi <- NA_real_

  # ---- keep one year in ratio panels (optional)
  # If your helper exists, reduce B/Bmsy and F/Fmsy to a single target year.
  if (exists(".keep_ratio_year", mode = "function")) {
    df <- .keep_ratio_year(df, keep_year = keep_ratio_year)
  }

  # ---- drop rows without finite estimates
  # % change requires finite numerators; remove invalid rows early.
  ok <- is.finite(df$est)
  df <- df[ok, , drop = FALSE]
  if (!nrow(df)) stop("No finite rows remain in `param_df`.", call. = FALSE)

  # ---- focus on requested parameters
  # If user-provided order has items not present, intersect to those that exist.
  params_present <- intersect(params_order, unique(df$param))
  if (!length(params_present)) {
    warning("None of the requested parameters found. Using all present params.")
    params_present <- sort(unique(df$param))
  }
  df <- df[df$param %in% params_present, , drop = FALSE]

  # ---- family ordering
  # Factorize family to a stable plotting/table order, but keep unseen families too.
  fams_present <- intersect(families_order, unique(df$family))
  fams_present <- unique(c(fams_present, setdiff(unique(df$family), fams_present)))
  df$family <- factor(df$family, levels = fams_present)

  # ---- scenario ordering (base first)
  # Put the base scenario first for auditing; others follow in original order.
  sc_all <- unique(df$scenario)
  if (!base_scenario %in% sc_all) {
    stop(paste0("Base scenario '", base_scenario, "' not found in `param_df$scenario`."))
  }
  sc_ord <- c(base_scenario, setdiff(sc_all, base_scenario))
  df$scenario <- factor(df$scenario, levels = sc_ord)

  # ---- split by type (SDM / GLM)
  types <- sort(unique(as.character(df$type)))

  # Small helper to create a unique key for a (type, family, param) cell.
  keyfun <- function(typ, fam, par) paste(typ, as.character(fam), par, sep = "||")

  # ---- identify base rows
  # Take only rows belonging to the base scenario (after filtering).
  base_rows <- df[as.character(df$scenario) == base_scenario, , drop = FALSE]
  if (!nrow(base_rows)) stop("No rows for base scenario after filtering.", call. = FALSE)

  # Index base rows by (type, family, param); if multiple, keep the first per key.
  base_key <- keyfun(as.character(base_rows$type), base_rows$family, base_rows$param)
  keep_first <- !duplicated(base_key)
  base_rows  <- base_rows[keep_first, , drop = FALSE]
  base_key   <- base_key[keep_first]

  # Named maps for quick lookup of base denominators (and optional CIs).
  base_est_map <- base_rows$est; names(base_est_map) <- base_key
  base_lo_map  <- base_rows$lo ; names(base_lo_map)  <- base_key
  base_hi_map  <- base_rows$hi ; names(base_hi_map)  <- base_key

  # ---- compute % changes vs. base
  n <- nrow(df)
  pct_est <- rep(NA_real_, n)
  pct_lo  <- rep(NA_real_, n)
  pct_hi  <- rep(NA_real_, n)

  for (i in seq_len(n)) {
    # Build the lookup key for this row and fetch the base denominator.
    kk <- keyfun(as.character(df$type[i]), df$family[i], df$param[i])
    if (!kk %in% names(base_est_map)) next            # skip if no matching base
    b0 <- base_est_map[[kk]]
    if (!is.finite(b0) || b0 == 0) next               # avoid divide-by-zero/NaN

    # % change on point estimate
    pct_est[i] <- 100 * (df$est[i] - b0) / b0

    # Optional % change on CI bounds (only if requested and finite)
    if (isTRUE(include_ci)) {
      if (is.finite(df$lo[i])) pct_lo[i] <- 100 * (df$lo[i] - b0) / b0
      if (is.finite(df$hi[i])) pct_hi[i] <- 100 * (df$hi[i] - b0) / b0
    }
  }

  # ---- build tidy output
  tidy <- data.frame(
    scenario = as.character(df$scenario),
    type     = as.character(df$type),
    family   = as.character(df$family),
    param    = as.character(df$param),
    pct_est  = pct_est,
    pct_lo   = if (isTRUE(include_ci)) pct_lo else NA_real_,  # keep columns stable
    pct_hi   = if (isTRUE(include_ci)) pct_hi else NA_real_,
    stringsAsFactors = FALSE
  )

  # Drop base scenario from deltas (always 0; table focuses on changes only).
  tidy <- tidy[tidy$scenario != base_scenario, , drop = FALSE]
  rownames(tidy) <- NULL

  # ---- build wide tables per type
  # Rows = Family × Param; Columns = scenarios excluding base; Cells = formatted "%".
  make_wide <- function(tidy_sub, type_label) {
    # Column order: all scenarios except base, keeping the factor order.
    sc_others <- setdiff(levels(df$scenario), base_scenario)
    sc_others <- sc_others[sc_others %in% unique(tidy_sub$scenario)]

    # Row order: requested families × requested params (limited to those present).
    families_use <- fams_present
    params_use   <- params_present

    # Initialize the skeleton wide table.
    nr <- length(families_use) * length(params_use)
    wide <- data.frame(
      Family = rep(as.character(families_use), each = length(params_use)),
      Param  = rep(as.character(params_use), times = length(families_use)),
      stringsAsFactors = FALSE
    )
    # Add one character column per scenario (fill later).
    if (length(sc_others)) {
      for (sc in sc_others) {
        wide[[sc]] <- NA_character_
      }
    }

    # Pretty formatter: "+3.1%" / "-0.8%" / "" (if NA).
    fmt <- function(x) {
      if (!is.finite(x)) return("")
      sprintf("%+.1f%%", x)
    }

    # Fill each cell with the first matching row for (family, param, scenario).
    for (i in seq_len(nrow(wide))) {
      fa <- wide$Family[i]
      pa <- wide$Param[i]
      # All rows for this (type, family, param)
      subset_rows <- tidy_sub[ tidy_sub$family == fa & tidy_sub$param == pa, , drop = FALSE]
      if (!nrow(subset_rows)) next
      # Write per-scenario, if present
      for (sc in sc_others) {
        j <- which(subset_rows$scenario == sc)
        if (length(j) >= 1L && is.finite(subset_rows$pct_est[j[1L]])) {
          wide[i, sc] <- fmt(subset_rows$pct_est[j[1L]])
        } else {
          wide[i, sc] <- ""
        }
      }
    }

    wide
  }

  # Assemble `wide` as a list split by `type` (e.g., SDM and GLM).
  wide_list <- list()
  for (tt in types) {
    sub <- tidy[ tidy$type == tt, , drop = FALSE]
    if (!nrow(sub)) next
    wide_list[[tt]] <- make_wide(sub, tt)
  }

  # Return tidy deltas, the per-type wide tables, and the exact base rows used.
  list(tidy = tidy, wide = wide_list, base = base_rows)
}
