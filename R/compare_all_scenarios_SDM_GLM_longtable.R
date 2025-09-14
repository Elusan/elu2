#' Export one LaTeX LONGTABLE comparing parameter values across ALL scenarios
#' (Beamer-friendly; shows a visible Parameter column)
#'
#' @param model_list Named list like "S1F.SDM", "S1F.GLM", "S1S.SDM", ...
#' @param output_dir Directory for the .tex file (must exist)
#' @param scenarios Character vector (auto-detected from names(model_list) if NULL)
#' @param block_def Named list: names = block headers; values = vector of parameter labels
#'                   exactly as in the summary `Parameters` column (e.g., "$r$", "$\\alpha_1$").
#' @param caption LaTeX caption (plain; ampersands auto-escaped)
#' @param verbose logical
#' @return (invisible) kableExtra object
#' @export
compare_all_scenarios_SDM_GLM_longtable <- function(
    model_list,
    output_dir = ".",
    scenarios = NULL,  # auto-detect by default
    block_def = NULL,
    caption = NULL,
    verbose = TRUE
) {
  if (!dir.exists(output_dir)) stop("Output directory does not exist!")

  # ---- Auto-detect scenarios from names(model_list) if not provided ----
  if (is.null(scenarios) || length(scenarios) == 0) {
    nm <- names(model_list)
    if (is.null(nm)) stop("`model_list` must be a named list with elements like 'S1F.SDM'.")
    scen_raw  <- gsub("^((S\\d+)).*$", "\\1", nm)
    scen_keep <- scen_raw[grepl("^S\\d+$", scen_raw)]
    scenarios <- unique(scen_keep)
    scenarios <- scenarios[order(as.integer(sub("^S(\\d+)$", "\\1", scenarios)))]
    if (length(scenarios) == 0) stop("Could not auto-detect scenarios from names(model_list).")
  }

  # ---------- helpers ----------
  pivot_to_wide <- function(df) {
    cols_needed <- c("Parameters",
                     "Fox_estimate", "Fox_low", "Fox_up",
                     "Schaefer_estimate", "Schaefer_low", "Schaefer_up",
                     "Pella_estimate", "Pella_low", "Pella_up")
    df_clean <- df[, cols_needed]
    keep <- !apply(df_clean, 1, function(row) all(is.na(row[-1])))
    df_clean <- df_clean[keep, , drop = FALSE]
    unique_params <- unique(as.character(df_clean$Parameters))

    wide <- data.frame(Parameters = unique_params, stringsAsFactors = FALSE)
    for (nm in c("Fox","Schaefer","Pella")) {
      wide[[paste0(nm,"_estimate")]] <- NA_real_
      wide[[paste0(nm,"_low")]]      <- NA_real_
      wide[[paste0(nm,"_up")]]       <- NA_real_
    }
    for (param in unique_params) {
      rows <- df_clean[df_clean$Parameters == param, , drop = FALSE]
      for (mod in c("Fox","Schaefer","Pella")) {
        est <- rows[[paste0(mod,"_estimate")]]
        low <- rows[[paste0(mod,"_low")]]
        up  <- rows[[paste0(mod,"_up")]]
        idx <- which(!is.na(est))
        if (length(idx) == 0) {
          est_val <- low_val <- up_val <- NA_real_
        } else {
          est_val <- est[idx[1]]; low_val <- low[idx[1]]; up_val <- up[idx[1]]
        }
        wide[wide$Parameters == param, paste0(mod,"_estimate")] <- est_val
        wide[wide$Parameters == param, paste0(mod,"_low")]      <- low_val
        wide[wide$Parameters == param, paste0(mod,"_up")]       <- up_val
      }
    }
    wide
  }

  # Two-line cell: bold estimate + small CI below (no nested tabular)
  format_cell <- function(est, low, up) {
    if (is.na(est)) return("\\shortstack[l]{~ \\\\ ~}")
    est_str <- sprintf("\\textbf{%.3f}", as.numeric(est))
    ci_str  <- sprintf("\\scriptsize(%.1f--%.1f)", as.numeric(low), as.numeric(up))
    paste0("\\shortstack[l]{", est_str, " \\\\ ", ci_str, "}")
  }

  # Escape bare ampersands in caption
  esc_amp <- function(x) gsub("(?<!\\\\)&", "\\\\&", x, perl = TRUE)

  # ---- SPiCT-driven block helpers (order + labels) ----

  # Grab the first spictcls object from model_list
  .first_spict_rep <- function(model_list) {
    for (nm in names(model_list)) {
      if (inherits(model_list[[nm]], "spictcls")) return(model_list[[nm]])
    }
    stop("No 'spictcls' objects found in model_list.")
  }

  # Build a default block_def using SPiCT summaries, then map to your LaTeX labels
  .spict_block_def <- function(rep_obj, pool_labels) {
    # last observed calendar year (rounded down)
    yr <- tryCatch({
      tt <- rep_obj$inp$timerangeObs[2]
      if (length(tt) == 1 && is.finite(tt)) floor(tt) else NA_integer_
    }, error = function(e) NA_integer_)

    # Fetch SPiCT summaries (best-effort)
    parest <- tryCatch(sumspict.parest(rep_obj),       error = function(e) NULL)
    sref   <- tryCatch(sumspict.srefpoints(rep_obj),   error = function(e) NULL)
    states <- tryCatch(sumspict.states(rep_obj),       error = function(e) NULL)

    # Extract raw names in the order SPiCT prints them
    get_names <- function(df) {
      if (is.null(df)) return(character(0))
      cand <- intersect(names(df), c("par","Par","parameter","Parameter","name","Name"))
      if (length(cand)) as.character(df[[cand[1]]]) else character(0)
    }
    par_names   <- get_names(parest)
    sref_names  <- get_names(sref)
    state_names <- get_names(states)

    # Map SPiCT short names -> LaTeX labels used by your table
    label_map <- new.env(parent = emptyenv())
    label_map[["r"]]    <- "$r$"
    label_map[["m"]]    <- "$m$"
    label_map[["n"]]    <- "$n$"
    label_map[["K"]]    <- "$K$"

    label_map[["Bmsy"]] <- "$B_{\\mathrm{MSY}}$"
    label_map[["Fmsy"]] <- "$F_{\\mathrm{MSY}}$"
    label_map[["MSY"]]  <- "$\\mathrm{MSY}$"

    label_map[["sdb"]]    <- "$\\sigma_B$"
    label_map[["logsdb"]] <- "$\\sigma_B$"

    label_map[["sdf"]]    <- "$\\sigma_F$"
    label_map[["logsdf"]] <- "$\\sigma_F$"

    label_map[["sdi"]]     <- "$\\sigma_{I_{1}}$"
    label_map[["logsdi"]]  <- "$\\sigma_{I_{1}}$"
    label_map[["sdi1"]]    <- "$\\sigma_{I_{1}}$"
    label_map[["logsdi1"]] <- "$\\sigma_{I_{1}}$"
    label_map[["sdi2"]]    <- "$\\sigma_{I_{2}}$"
    label_map[["logsdi2"]] <- "$\\sigma_{I_{2}}$"

    label_map[["sdc"]]    <- "$\\sigma_C$"
    label_map[["logsdc"]] <- "$\\sigma_C$"

    if (is.finite(yr)) {
      label_map[[paste0("B_", yr)]]     <- sprintf("$B_{%d}$", yr)
      label_map[[paste0("F_", yr)]]     <- sprintf("$F_{%d}$", yr)
      label_map[[paste0("BBmsy_", yr)]] <- sprintf("$B_{%d}/B_{\\mathrm{MSY}}$", yr)
      label_map[[paste0("FFmsy_", yr)]] <- sprintf("$F_{%d}/F_{\\mathrm{MSY}}$", yr)
      label_map[["Bnow"]]   <- sprintf("$B_{%d}$", yr)
      label_map[["Fnow"]]   <- sprintf("$F_{%d}$", yr)
      label_map[["BBmsy"]]  <- sprintf("$B_{%d}/B_{\\mathrm{MSY}}$", yr)
      label_map[["FFmsy"]]  <- sprintf("$F_{%d}/F_{\\mathrm{MSY}}$", yr)
    }

    map_keep <- function(keys) {
      lbls <- vapply(keys, function(k) {
        v <- label_map[[k]]
        if (is.null(v)) NA_character_ else v
      }, character(1))
      lbls <- lbls[!is.na(lbls)]
      lbls[lbls %in% pool_labels]
    }

    dynamics_keys  <- par_names[par_names %in% c("r","m","n","K")]
    variance_keys  <- par_names[grepl("^logsd", par_names) |
                                  grepl("^sd",    par_names) |
                                  par_names %in% c("sdb","sdf","sdi","sdi1","sdi2","sdc",
                                                   "logsdb","logsdf","logsdi","logsdi1","logsdi2","logsdc")]

    msy_keys    <- sref_names[sref_names %in% c("Bmsy","Fmsy","MSY")]
    status_keys <- state_names[state_names %in% c(
      paste0("B_", yr), paste0("F_", yr), paste0("BBmsy_", yr), paste0("FFmsy_", yr),
      "Bnow","Fnow","BBmsy","FFmsy"
    )]

    dynamics_lbls <- map_keep(dynamics_keys)
    variances_lbls<- map_keep(variance_keys)
    msy_lbls      <- map_keep(msy_keys)
    status_lbls   <- map_keep(status_keys)

    out <- list()
    if (length(variances_lbls)) out[["Process & obs. variances"]] <- variances_lbls
    if (length(msy_lbls))       out[["MSY quantities"]]           <- msy_lbls
    if (length(status_lbls)) {
      nm <- if (is.finite(yr)) sprintf("Status (%d)", yr) else "Status"
      out[[nm]] <- status_lbls
    }
    if (length(dynamics_lbls))  out[["Dynamics"]]                 <- dynamics_lbls
    if (!length(out)) out[["All parameters"]] <- pool_labels
    out
  }

  # ---------- summaries per scenario ----------
  wlist_SDM <- list(); wlist_GLM <- list()
  get_mod <- function(scen, model, type) {
    nm <- paste0(scen, model, ".", type)
    if (!nm %in% names(model_list)) stop("Model not found: ", nm)
    model_list[[nm]]
  }
  for (scen in scenarios) {
    if (verbose) message("Collecting summaries for ", scen, " ...")
    models_SDM <- list(Fox = get_mod(scen, "F", "SDM"),
                       Schaefer = get_mod(scen, "S", "SDM"),
                       Pella = get_mod(scen, "P", "SDM"))
    models_GLM <- list(Fox = get_mod(scen, "F", "GLM"),
                       Schaefer = get_mod(scen, "S", "GLM"),
                       Pella = get_mod(scen, "P", "GLM"))
    df_SDM_long <- make_spict_summary_df1_nopipe(models_SDM, names(models_SDM))
    df_GLM_long <- make_spict_summary_df1_nopipe(models_GLM, names(models_GLM))
    wlist_SDM[[scen]] <- pivot_to_wide(df_SDM_long)
    wlist_GLM[[scen]] <- pivot_to_wide(df_GLM_long)
  }

  # ---------- SPiCT-driven default blocks (order taken from SPiCT) ----------
  if (is.null(block_def)) {
    # pool of LaTeX labels actually present in your wide tables, across all scenarios & both datasets
    pool_labels <- unique(c(
      unlist(lapply(wlist_SDM, function(w) if (is.null(w)) character(0) else as.character(w$Parameters))),
      unlist(lapply(wlist_GLM, function(w) if (is.null(w)) character(0) else as.character(w$Parameters)))
    ))
    pool_labels <- pool_labels[is.finite(match(pool_labels, pool_labels))]

    rep0 <- .first_spict_rep(model_list)
    block_def <- .spict_block_def(rep0, pool_labels)
  }

  # ---------- assemble big table ----------
  rows <- list()
  block_ranges <- list()  # for group_rows by block
  current_row <- 0L
  block_names <- names(block_def)

  for (b in seq_along(block_def)) {
    params <- block_def[[b]]
    block_start <- current_row + 1L

    for (param in params) {
      for (scen in scenarios) {
        wS <- wlist_SDM[[scen]]; wG <- wlist_GLM[[scen]]
        iS <- match(param, wS$Parameters); iG <- match(param, wG$Parameters)
        rows[[length(rows) + 1L]] <- data.frame(
          Block         = block_names[b],
          Parameter     = param,            # <-- visible in first column
          Scenario      = scen,
          Fox_SDM       = if (!is.na(iS)) format_cell(wS$Fox_estimate[iS],      wS$Fox_low[iS],      wS$Fox_up[iS]) else format_cell(NA,NA,NA),
          Fox_GLM       = if (!is.na(iG)) format_cell(wG$Fox_estimate[iG],      wG$Fox_low[iG],      wG$Fox_up[iG]) else format_cell(NA,NA,NA),
          Schaefer_SDM  = if (!is.na(iS)) format_cell(wS$Schaefer_estimate[iS], wS$Schaefer_low[iS], wS$Schaefer_up[iS]) else format_cell(NA,NA,NA),
          Schaefer_GLM  = if (!is.na(iG)) format_cell(wG$Schaefer_estimate[iG], wG$Schaefer_low[iG], wG$Schaefer_up[iG]) else format_cell(NA,NA,NA),
          Pella_SDM     = if (!is.na(iS)) format_cell(wS$Pella_estimate[iS],    wS$Pella_low[iS],    wS$Pella_up[iS]) else format_cell(NA,NA,NA),
          Pella_GLM     = if (!is.na(iG)) format_cell(wG$Pella_estimate[iG],    wG$Pella_low[iG],    wG$Pella_up[iG]) else format_cell(NA,NA,NA),
          stringsAsFactors = FALSE
        )
        current_row <- current_row + 1L
      }
    }

    block_ranges[[block_names[b]]] <- c(block_start, current_row)
  }

  big <- do.call(rbind, rows)

  # ---------- caption ----------
  if (is.null(caption)) {
    caption <- "Comparison of parameter estimates (SDM & GLM; Fox, Schaefer, Pella) across scenarios. Cells show bold estimate with 95\\% CI on a new line."
  }
  caption <- esc_amp(caption)  # SDM \& GLM

  # ---------- build longtable (now includes Parameter column) ----------
  out_file <- file.path(output_dir, "table_compare_all_scenarios_SDM_GLM_longtable.tex")

  kbl_out <- kableExtra::kbl(
    big[, c("Parameter","Scenario","Fox_SDM","Fox_GLM","Schaefer_SDM","Schaefer_GLM","Pella_SDM","Pella_GLM")],
    format    = "latex",
    booktabs  = TRUE,
    longtable = TRUE,
    escape    = FALSE,
    linesep = "\addlinespace",
    align     = "llllllll",   # 8 columns
    col.names = c("\\textbf{Para.}","\\textbf{Scena.}","SDM","GLM","SDM","GLM","SDM","GLM"),
    caption   = caption
  )

  kbl_out <- kableExtra::add_header_above(
    kbl_out,
    header = c(" " = 2, "Fox" = 2, "Schaefer" = 2, "Pella" = 2),  # 2 leading cols
    escape = FALSE
  )

  # Group by block (whole sections)
  for (bn in names(block_ranges)) {
    rng <- block_ranges[[bn]]
    kbl_out <- kableExtra::group_rows(
      kbl_out, label = bn, start_row = rng[1], end_row = rng[2],
      escape = FALSE, bold = TRUE, latex_gap_space = "0.35em"
    )
  }

  # Styling + repeated header on continuation pages
  kbl_out <- kableExtra::kable_styling(
    kbl_out,
    latex_options = c("repeat_header"),
    repeat_header_text = "(continued)",
    full_width = FALSE,
    position = "center",
    font_size = 8
  )

  # Widths for the six value columns only (cols 3..8)
  kbl_out <- kableExtra::column_spec(kbl_out, 3:8, width = "2cm", latex_column_spec = "p{2cm}")

  kbl_out <- kableExtra::footnote(
    kbl_out,
    general = "‘Fixed’ indicates the parameter was not estimated and was assigned a fixed value.",
    general_title = "Note: ",
    footnote_as_chunk = TRUE,
    escape = FALSE
  )

  # Add longtable layout tweaks
  tex_body <- as.character(kbl_out)
  prefix <- paste(
    "\\setlength{\\LTleft}{0pt}\\setlength{\\LTright}{0pt}",
    "\\captionsetup[table]{singlelinecheck=false,justification=justified,width=\\linewidth}",
    sep = "\n"
  )
  tex_body <- paste0(prefix, "\n", tex_body)

  writeLines(tex_body, out_file)
  if (verbose) message("Written: ", out_file)
  invisible(kbl_out)
}
