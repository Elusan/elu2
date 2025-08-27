#' Export one LaTeX LONGTABLE per Parameter (spliced from the combined table)
#' (Beamer-friendly; identical styling to your combined longtable)
#'
#' @param model_list Named list like "S1F.SDM", "S1F.GLM", "S1S.SDM", ...
#' @param output_dir Directory for the .tex files (must exist)
#' @param scenarios Character vector (auto-detected from names(model_list) if NULL)
#' @param block_def Named list for parameter order. If NULL, uses SPiCT order
#'                  restricted to your selected parameters/groups below.
#' @param caption LaTeX caption (plain; ampersands auto-escaped). If NULL, uses default.
#' @param verbose logical
#' @return (invisible) named character vector of written file paths (names = parameters)
#' @export
compare_all_scenarios_SDM_GLM_longtable_split2 <- function(
    model_list,
    output_dir = ".",
    scenarios = NULL,  # CHANGED: auto-detect by default
    block_def = NULL,
    caption = NULL,
    verbose = TRUE
) {
  if (!dir.exists(output_dir)) stop("Output directory does not exist!")
  if (!requireNamespace("kableExtra", quietly = TRUE)) stop("kableExtra required.")

  # CHANGED: Auto-detect scenarios from names(model_list) if not provided
  if (is.null(scenarios) || length(scenarios) == 0) {
    nm <- names(model_list)
    if (is.null(nm)) stop("`model_list` must be a named list with elements like 'S1F.SDM'.")
    scen_raw  <- gsub("^((S\\d+)).*$", "\\1", nm)
    scen_keep <- scen_raw[grepl("^S\\d+$", scen_raw)]
    scenarios <- unique(scen_keep)
    scenarios <- scenarios[order(as.integer(sub("^S(\\d+)$", "\\1", scenarios)))]
    if (length(scenarios) == 0) stop("Could not auto-detect scenarios from names(model_list).")
  }

  ## ---------- helpers ----------
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

  # normalization utils for robust matching
  .unlatex <- function(x) {
    x <- gsub("\\\\mathrm\\{([^}]*)\\}", "\\1", x)
    x <- gsub("\\\\textit\\{([^}]*)\\}", "\\1", x)
    x <- gsub("\\$","", x); x <- gsub("\\\\","", x)
    x <- gsub("[{}]","", x); x <- gsub("\\s+","", x)
    x <- gsub("_\\{","_", x, perl = TRUE); x <- gsub("\\}","", x)
    tolower(x)
  }
  .norm_key <- function(x) gsub("_","", .unlatex(x))
  .match_idx <- function(lbl, vec) {
    nk <- .norm_key(lbl); vv <- .norm_key(vec)
    idx <- which(vv == nk)
    if (length(idx)) idx[1] else NA_integer_
  }

  esc_amp <- function(x) gsub("(?<!\\\\)&", "\\\\&", x, perl = TRUE)
  .safe_stem <- function(lbl) {
    s <- gsub("\\\\mathrm\\{([^}]*)\\}", "\\1", lbl)
    s <- gsub("\\\\", "", s); s <- gsub("\\$", "", s)
    s <- gsub("[{}]", "", s); s <- gsub("\\^", "", s)
    s <- gsub("\\s+", "_", s); s <- gsub("[^A-Za-z0-9_\\-]+", "", s)
    s <- gsub("_+", "_", s); tolower(s)
  }
  .first_spict_rep <- function(model_list) {
    for (nm in names(model_list)) if (inherits(model_list[[nm]], "spictcls")) return(model_list[[nm]])
    stop("No 'spictcls' objects found in model_list.")
  }

  # ------- exact q1/q2 override from sumspict.parest() (adds rows if missing) ------
  .get_parest_row <- function(df, key) {
    if (is.null(df)) return(c(NA_real_, NA_real_, NA_real_))
    rn <- rownames(df)
    cols <- intersect(colnames(df), c("estimate","cilow","ciupp"))
    if (length(cols) < 3) return(c(NA_real_, NA_real_, NA_real_))
    if (length(rn) && key %in% rn) {
      as.numeric(df[key, cols])
    } else if ("par" %in% names(df) && key %in% df$par) {
      i <- which(df$par == key)[1]
      as.numeric(df[i, cols])
    } else {
      c(NA_real_, NA_real_, NA_real_)
    }
  }
  .ensure_row <- function(wide, label) {
    if (!any(.norm_key(wide$Parameters) == .norm_key(label))) {
      add <- wide[0, , drop = FALSE]
      add[1, ] <- NA
      add$Parameters <- label
      wide <- rbind(wide, add)
      rownames(wide) <- NULL
    }
    wide
  }
  .override_q_from_spict <- function(wide, models) {
    for (mod in c("Fox","Schaefer","Pella")) {
      rep <- models[[mod]]
      if (is.null(rep)) next
      pr <- tryCatch(sumspict.parest(rep), error = function(e) NULL)
      if (is.null(pr)) next

      # q1
      wide <- .ensure_row(wide, "$q_{1}$")
      vals <- .get_parest_row(pr, "q1")
      if (!all(is.na(vals))) {
        ridx <- .match_idx("$q_{1}$", wide$Parameters)
        wide[ridx, paste0(mod,"_estimate")] <- vals[1]
        wide[ridx, paste0(mod,"_low")]      <- vals[2]
        wide[ridx, paste0(mod,"_up")]       <- vals[3]
      }
      # q2
      wide <- .ensure_row(wide, "$q_{2}$")
      vals <- .get_parest_row(pr, "q2")
      if (!all(is.na(vals))) {
        ridx <- .match_idx("$q_{2}$", wide$Parameters)
        wide[ridx, paste0(mod,"_estimate")] <- vals[1]
        wide[ridx, paste0(mod,"_low")]      <- vals[2]
        wide[ridx, paste0(mod,"_up")]       <- vals[3]
      }
    }
    wide
  }

  # two-line cell; q’s use SPiCT-like exact scientific printing (%.6e)
  #.fmt_q <- function(x) if (is.na(x)) NA_character_ else sprintf("%.6e", as.numeric(x))
  .fmt_q <- function(x) if (is.na(x)) NA_character_ else sprintf("%.2e", as.numeric(x))

  format_cell <- function(est, low, up, param_label) {
    if (is.na(est)) return("\\shortstack[l]{~ \\\\ ~}")
    is_q <- .norm_key(param_label) %in% c("q1","q2")
    if (is_q) {
      est_str <- paste0("\\textbf{", .fmt_q(est), "}")
      ci_str  <- paste0("\\CIstyle{(", .fmt_q(low), "--", .fmt_q(up), ")}")
    } else {
      est_str <- sprintf("\\textbf{%.3f}", as.numeric(est))
      ci_str  <- sprintf("\\CIstyle{(%.1f--%.1f)}", as.numeric(low), as.numeric(up))
    }
    paste0("\\shortstack[l]{", est_str, " \\\\ ", ci_str, "}")
  }

  # ---------- SPiCT-driven block order, restricted to your selections ----------
  .spict_block_def <- function(rep_obj, pool_labels) {
    group1_sel <- c("$\\alpha_{1}$", "$\\alpha_{2}$", "$\\beta$", "$r$", "$m$", "$K$",
                    "$q_{1}$", "$q_{2}$", "$n$", "$\\sigma_B$", "$\\sigma_F$",
                    "$\\sigma_{I_{1}}$", "$\\sigma_{I_{2}}$", "$\\sigma_C$")
    group2_sel <- c("$B_{\\mathrm{MSY}}$", "$F_{\\mathrm{MSY}}$", "$\\mathrm{MSY}$")
    group3_sel <- c("$B_{2022}$", "$F_{2022}$",
                    "$B_{2022}/B_{\\mathrm{MSY}}$", "$F_{2022}/F_{\\mathrm{MSY}}$")
    group4_sel <- c("$B_{2024}$", "$F_{2024}$",
                    "$B_{2024}/B_{\\mathrm{MSY}}$", "$F_{2024}/F_{\\mathrm{MSY}}$")

    parest <- tryCatch(sumspict.parest(rep_obj),     error = function(e) NULL)
    sref   <- tryCatch(sumspict.srefpoints(rep_obj), error = function(e) NULL)
    states <- tryCatch(sumspict.states(rep_obj),     error = function(e) NULL)

    get_names <- function(df) {
      if (is.null(df)) return(character(0))
      cand <- intersect(names(df), c("par","Par","parameter","Parameter","name","Name","row.names"))
      if (length(cand)) as.character(df[[cand[1]]]) else { rn <- rownames(df); if (length(rn)) rn else character(0) }
    }
    par_names   <- get_names(parest)
    sref_names  <- get_names(sref)
    state_names <- get_names(states)

    label_map <- new.env(parent = emptyenv())
    label_map[["alpha1"]] <- "$\\alpha_{1}$"; label_map[["alpha2"]] <- "$\\alpha_{2}$"
    label_map[["beta"]]   <- "$\\beta$";      label_map[["r"]]      <- "$r$"
    label_map[["m"]]      <- "$m$";           label_map[["K"]]      <- "$K$"
    label_map[["q1"]]     <- "$q_{1}$";       label_map[["q2"]]     <- "$q_{2}$"
    label_map[["n"]]      <- "$n$"
    label_map[["sdb"]]    <- "$\\sigma_B$";   label_map[["logsdb"]] <- "$\\sigma_B$"
    label_map[["sdf"]]    <- "$\\sigma_F$";   label_map[["logsdf"]] <- "$\\sigma_F$"
    label_map[["sdi1"]]   <- "$\\sigma_{I_{1}}$"; label_map[["logsdi1"]] <- "$\\sigma_{I_{1}}$"
    label_map[["sdi2"]]   <- "$\\sigma_{I_{2}}$"; label_map[["logsdi2"]] <- "$\\sigma_{I_{2}}$"
    label_map[["sdi"]]    <- "$\\sigma_{I_{1}}$"; label_map[["logsdi"]]  <- "$\\sigma_{I_{1}}$"
    label_map[["sdc"]]    <- "$\\sigma_C$";   label_map[["logsdc"]] <- "$\\sigma_C$"

    label_map[["Bmsy"]]   <- "$B_{\\mathrm{MSY}}$"; label_map[["Bmsys"]] <- "$B_{\\mathrm{MSY}}$"
    label_map[["Fmsy"]]   <- "$F_{\\mathrm{MSY}}$"; label_map[["Fmsys"]] <- "$F_{\\mathrm{MSY}}$"
    label_map[["MSY"]]    <- "$\\mathrm{MSY}$";     label_map[["MSYs"]]  <- "$\\mathrm{MSY}$"

    to_state_label <- function(key) {
      m <- regexec("^([BF])_(\\d{4})(?:\\.[0-9]+)?(?:/(B|F)msy)?$", key)
      mm <- regmatches(key, m)[[1]]
      if (length(mm) == 0) return(NA_character_)
      BF <- mm[2]; Y <- mm[3]; ratio <- mm[4]
      if (is.na(ratio) || ratio == "") {
        if (BF == "B") return(sprintf("$B_{%s}$", Y))
        if (BF == "F") return(sprintf("$F_{%s}$", Y))
      } else {
        if (BF == "B" && ratio == "B") return(sprintf("$B_{%s}/B_{\\mathrm{MSY}}$", Y))
        if (BF == "F" && ratio == "F") return(sprintf("$F_{%s}/F_{\\mathrm{MSY}}$", Y))
      }
      NA_character_
    }
    map_to_labels <- function(keys) {
      vapply(keys, function(k) {
        v <- label_map[[k]]
        if (!is.null(v)) return(v)
        s <- to_state_label(k); if (!is.na(s)) return(s)
        NA_character_
      }, character(1))
    }

    pool_map <- function(labels_wanted, pool_labels) {
      wanted <- labels_wanted[!is.na(labels_wanted)]
      # match by normalized key to the ACTUAL labels present in pool_labels
      vapply(wanted, function(lbl) {
        idx <- .match_idx(lbl, pool_labels)
        if (is.na(idx)) NA_character_ else pool_labels[idx]
      }, character(1)) -> out
      out[!is.na(out)]
    }

    sel1 <- c("$\\alpha_{1}$", "$\\alpha_{2}$", "$\\beta$", "$r$", "$m$", "$K$",
              "$q_{1}$", "$q_{2}$", "$n$", "$\\sigma_B$", "$\\sigma_F$",
              "$\\sigma_{I_{1}}$", "$\\sigma_{I_{2}}$", "$\\sigma_C$")
    g1 <- pool_map(map_to_labels(par_names), pool_labels)
    g1 <- g1[.norm_key(g1) %in% .norm_key(sel1)]

    sel2 <- c("$B_{\\mathrm{MSY}}$", "$F_{\\mathrm{MSY}}$", "$\\mathrm{MSY}$")
    g2 <- pool_map(map_to_labels(sref_names), pool_labels)
    g2 <- g2[.norm_key(g2) %in% .norm_key(sel2)]

    g3 <- pool_map(c("$B_{2022}$", "$F_{2022}$",
                     "$B_{2022}/B_{\\mathrm{MSY}}$", "$F_{2022}/F_{\\mathrm{MSY}}$"), pool_labels)
    g4 <- pool_map(c("$B_{2024}$", "$F_{2024}$",
                     "$B_{2024}/B_{\\mathrm{MSY}}$", "$F_{2024}/F_{\\mathrm{MSY}}$"), pool_labels)

    out <- list()
    if (length(g1)) out[["Process, dynamics & obs. (selected)"]] <- g1
    if (length(g2)) out[["MSY quantities"]] <- g2
    if (length(g3)) out[["Status (2022)"]]  <- g3
    if (length(g4)) out[["Status (2024)"]]  <- g4
    if (!length(out)) out[["Selected parameters"]] <- intersect(sel1, pool_labels)
    out
  }

  ## ---------- summaries per scenario ----------
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

    wS <- pivot_to_wide(df_SDM_long)
    wG <- pivot_to_wide(df_GLM_long)

    # **Critical**: force q1/q2 from raw parest; add rows if missing
    wS <- .override_q_from_spict(wS, models_SDM)
    wG <- .override_q_from_spict(wG, models_GLM)

    wlist_SDM[[scen]] <- wS
    wlist_GLM[[scen]] <- wG
  }

  ## ---------- default blocks ----------
  if (is.null(block_def)) {
    pool_labels <- unique(c(
      unlist(lapply(wlist_SDM, function(w) if (is.null(w)) character(0) else as.character(w$Parameters))),
      unlist(lapply(wlist_GLM, function(w) if (is.null(w)) character(0) else as.character(w$Parameters)))
    ))
    pool_labels <- pool_labels[!is.na(pool_labels)]
    rep0 <- .first_spict_rep(model_list)
    block_def <- .spict_block_def(rep0, pool_labels)
  }

  ## ---------- build table rows ----------
  rows <- list()
  block_names <- names(block_def)
  for (b in seq_along(block_def)) {
    params <- block_def[[b]]
    for (param in params) {
      for (scen in scenarios) {
        wS <- wlist_SDM[[scen]]; wG <- wlist_GLM[[scen]]
        iS <- .match_idx(param, wS$Parameters)
        iG <- .match_idx(param, wG$Parameters)
        rows[[length(rows) + 1L]] <- data.frame(
          Block         = block_names[b],
          Parameter     = param,
          Scenario      = scen,
          Fox_SDM       = if (!is.na(iS)) format_cell(wS$Fox_estimate[iS],      wS$Fox_low[iS],      wS$Fox_up[iS],      param) else format_cell(NA,NA,NA,param),
          Fox_GLM       = if (!is.na(iG)) format_cell(wG$Fox_estimate[iG],      wG$Fox_low[iG],      wG$Fox_up[iG],      param) else format_cell(NA,NA,NA,param),
          Schaefer_SDM  = if (!is.na(iS)) format_cell(wS$Schaefer_estimate[iS], wS$Schaefer_low[iS], wS$Schaefer_up[iS], param) else format_cell(NA,NA,NA,param),
          Schaefer_GLM  = if (!is.na(iG)) format_cell(wG$Schaefer_estimate[iG], wG$Schaefer_low[iG], wG$Schaefer_up[iG], param) else format_cell(NA,NA,NA,param),
          Pella_SDM     = if (!is.na(iS)) format_cell(wS$Pella_estimate[iS],    wS$Pella_low[iS],    wS$Pella_up[iS],    param) else format_cell(NA,NA,NA,param),
          Pella_GLM     = if (!is.na(iG)) format_cell(wG$Pella_estimate[iG],    wG$Pella_low[iG],    wG$Pella_up[iG],    param) else format_cell(NA,NA,NA,param),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  big <- do.call(rbind, rows)

  ## ---------- caption ----------
  if (is.null(caption)) {
    caption <- "Comparison of parameter estimates (SDM & GLM; Fox, Schaefer, Pella) across scenarios. Cells show bold estimate with 95\\% CI on a new line."
  }
  caption <- esc_amp(caption)

  ## ---------- write one .tex per parameter ----------
  param_order <- unlist(block_def, use.names = FALSE)
  param_order <- param_order[!is.na(match(.norm_key(param_order), .norm_key(unique(big$Parameter))))]
  written <- character(0)

  for (param in param_order) {
    sub <- big[.norm_key(big$Parameter) == .norm_key(param), , drop = FALSE]
    rownames(sub) <- NULL

    kbl_out <- kableExtra::kbl(
      sub[, c("Parameter","Scenario","Fox_SDM","Fox_GLM","Schaefer_SDM","Schaefer_GLM","Pella_SDM","Pella_GLM")],
      format    = "latex",
      booktabs  = TRUE,
      longtable = TRUE,
      escape    = FALSE,
      linesep   = "",
      align     = "llllllll",
      col.names = c("\\textbf{Para.}","\\textbf{Scena.}","SDM","GLM","SDM","GLM","SDM","GLM"),
      caption   = caption,
      row.names = FALSE
    )

    kbl_out <- kableExtra::add_header_above(
      kbl_out,
      header = c(" " = 2, "Fox" = 2, "Schaefer" = 2, "Pella" = 2),
      escape = FALSE
    )

    kbl_out <- kableExtra::kable_styling(
      kbl_out,
      latex_options = c("repeat_header"),
      repeat_header_text = "(continued)",
      full_width = FALSE,
      position = "center",
      font_size = 6
    )

    kbl_out <- kableExtra::column_spec(kbl_out, 3:8, width = "1.4cm", latex_column_spec = "p{1.4cm}")

    kbl_out <- kableExtra::footnote(
      kbl_out,
      general = "‘Fixed’ indicates the parameter was not estimated and was assigned a fixed value.",
      general_title = "Note: ",
      footnote_as_chunk = TRUE,
      escape = FALSE
    )

    tex_body <- as.character(kbl_out)
    prefix <- paste(
      "\\setlength{\\LTleft}{\\fill}\\setlength{\\LTright}{\\fill}",
      "\\captionsetup[table]{singlelinecheck=false,justification=justified,width=\\linewidth}",
      sep = "\n"
    )
    tex_body <- paste0(prefix, "\n", tex_body)

    out_file <- file.path(output_dir, paste0("table_compare_param_", .safe_stem(param), "_SDM_GLM_longtable.tex"))
    writeLines(tex_body, out_file)
    if (verbose) message("Written: ", out_file)
    written[param] <- out_file
  }

  invisible(written)
}
