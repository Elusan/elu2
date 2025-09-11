# ---- Theme (compact, crisp border, tidy legend) -----------------------------
theme_minimal_compact2_for_PANELS <- function(base_size = 10, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = element_text(face = "bold", size = 12),
      axis.text  = element_text(size = 12, face = "bold"),
      legend.position = c(0.88, 0.98),
      legend.justification = c("left", "top"),
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      legend.title = element_blank(),
      legend.text  = element_text(size = 12),
      legend.key.size = unit(1, "lines"),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, colour = "grey35", linewidth = 2),
      axis.ticks = element_line(linewidth = 0.7, color = "grey35"),
      axis.ticks.length = unit(3, "pt"),
      # keep default strip bg here; we will override to blank + inside in the plot call
      strip.background = element_rect(fill = "grey35", color = "grey35", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = rel(1)),  # base; overridden in plot
      text = element_text(face = "bold", size = 11),
      plot.margin = margin(3, 3, 3, 3)
    )
}

# ---- Plot function ----------------------------------------------------------
plot_param_panels_by_model <- function(model_list,
                                       model_name = c("Pella","Schaefer","Fox"),
                                       scenario_levels = c("S1","S2","S3","S4"),
                                       year = NULL,
                                       font_size_base = 11,
                                       save_path = NULL,
                                       dodge_width = 0.36,            # tighter SDM–GLM gap
                                       x_expand_mult = c(0.02, 0.02)  # trim outer padding
) {
  model_name <- match.arg(model_name)
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
  if (!exists("make_spict_summary_df1_nopipe")) stop("make_spict_summary_df1_nopipe() must be available.")

  .first_spict_rep <- function(x) { for (i in seq_along(x)) if (inherits(x[[i]], "spictcls")) return(x[[i]]); stop("No 'spictcls' in model_list.") }
  .infer_year <- function(x) { rep0 <- .first_spict_rep(x); yy <- tryCatch(rep0$inp$timerangeObs[2], error=function(e) NA_real_); if (!is.finite(yy)) stop("Cannot infer terminal year."); as.integer(floor(yy)) }
  .tex_label_map <- function(y) c(
    BBmsy_year = sprintf("$B_{%d}/B_{\\mathrm{MSY}}$", y),
    Bmsy="$B_{\\mathrm{MSY}}$", B_year=sprintf("$B_{%d}$", y), K="$K$", n="$n$",
    FFmsy_year=sprintf("$F_{%d}/F_{\\mathrm{MSY}}$", y), F_year=sprintf("$F_{%d}$", y),
    Fmsy="$F_{\\mathrm{MSY}}$", MSY="$\\mathrm{MSY}$", r="$r$"
  )
  .plotmath_string_map <- function(y) {
    tex <- .tex_label_map(y)
    setNames(c(sprintf("B[%d]/B[MSY]", y), "B[MSY]", sprintf("B[%d]", y), "K","n",
               sprintf("F[%d]/F[MSY]", y), sprintf("F[%d]", y), "F[MSY]", "MSY", "r"), tex)
  }
  .get_models <- function(scen, type) list(
    Fox=model_list[[paste0(scen,"F.",type)]],
    Schaefer=model_list[[paste0(scen,"S.",type)]],
    Pella=model_list[[paste0(scen,"P.",type)]]
  )
  .model_column_names <- function(model_name) paste0(model_name, c("_estimate","_low","_up"))

  if (is.null(year)) year <- .infer_year(model_list)
  tex_map <- .tex_label_map(year)
  target_order_keys <- c("BBmsy_year","Bmsy","B_year","K","n","FFmsy_year","F_year","Fmsy","MSY","r")
  target_labels_tex <- unname(tex_map[target_order_keys])

  est_col <- .model_column_names(model_name)[1]
  lo_col  <- .model_column_names(model_name)[2]
  up_col  <- .model_column_names(model_name)[3]

  all_rows <- list(); k <- 0L
  for (scen in scenario_levels) {
    mods_SDM <- .get_models(scen,"SDM")
    df_SDM <- make_spict_summary_df1_nopipe(mods_SDM, names(mods_SDM))
    sub_SDM <- df_SDM[!is.na(match(df_SDM$Parameters, target_labels_tex)),
                      c("Parameters", est_col, lo_col, up_col), drop=FALSE]
    if (nrow(sub_SDM)) for (i in seq_len(nrow(sub_SDM))) {
      k <- k+1L; all_rows[[k]] <- data.frame(
        scenario=scen, dtype="SDM", parameter_tex=as.character(sub_SDM$Parameters[i]),
        est=as.numeric(sub_SDM[i,est_col]), lwr=as.numeric(sub_SDM[i,lo_col]),
        upr=as.numeric(sub_SDM[i,up_col]), stringsAsFactors=FALSE)
    }
    mods_GLM <- .get_models(scen,"GLM")
    df_GLM <- make_spict_summary_df1_nopipe(mods_GLM, names(mods_GLM))
    sub_GLM <- df_GLM[!is.na(match(df_GLM$Parameters, target_labels_tex)),
                      c("Parameters", est_col, lo_col, up_col), drop=FALSE]
    if (nrow(sub_GLM)) for (i in seq_len(nrow(sub_GLM))) {
      k <- k+1L; all_rows[[k]] <- data.frame(
        scenario=scen, dtype="GLM", parameter_tex=as.character(sub_GLM$Parameters[i]),
        est=as.numeric(sub_GLM[i,est_col]), lwr=as.numeric(sub_GLM[i,lo_col]),
        upr=as.numeric(sub_GLM[i,up_col]), stringsAsFactors=FALSE)
    }
  }
  if (!length(all_rows)) stop("No data found for the requested model / scenarios.")
  dat <- do.call(rbind, all_rows)

  dat$scenario      <- factor(dat$scenario, levels=scenario_levels, ordered=TRUE)
  dat$dtype         <- factor(dat$dtype, levels=c("SDM","GLM"), ordered=TRUE)
  dat$parameter_tex <- factor(dat$parameter_tex, levels=target_labels_tex, ordered=TRUE)

  keep <- is.finite(dat$est) & is.finite(dat$lwr) & is.finite(dat$upr)
  dat  <- dat[keep, , drop=FALSE]
  if (!nrow(dat)) stop("All rows had missing estimate/CI.")

  facet_map_strings <- .plotmath_string_map(year)
  pd <- ggplot2::position_dodge(width = dodge_width)

  gg <- ggplot2::ggplot(dat, ggplot2::aes(x = scenario, group = dtype)) +
    ggplot2::geom_linerange(
      ggplot2::aes(ymin = lwr, ymax = upr, linetype = dtype),
      position = pd, linewidth = 0.6
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = est, shape = dtype),
      position = pd, size = 2.5, colour = "red"
    ) +
    ggplot2::facet_wrap(
      ~ parameter_tex,
      scales = "free_y",
      ncol = 2,
      labeller = ggplot2::as_labeller(facet_map_strings, ggplot2::label_parsed)
    ) +
    ggplot2::labs(
      x = NULL, y = NULL,
      title = paste0(model_name, " model — Parameters with 95% CI (SDM vs GLM, year ", year, ")"),
      subtitle = paste("Scenarios:", paste(scenario_levels, collapse = ", "))
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.14))) +  # a bit more top padding
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = x_expand_mult)) +
    ggplot2::scale_linetype_manual(values = c(SDM = "solid", GLM = "dashed")) +
    ggplot2::scale_shape_manual(values = c(SDM = 16, GLM = 17)) +
    ggplot2::guides(linetype = ggplot2::guide_legend(title = NULL),
                    shape    = ggplot2::guide_legend(title = NULL)) +
    # Base theme
    theme_minimal_compact2_for_PANELS(base_size = font_size_base) +
    # Put facet titles INSIDE, bold, larger, no gray band
    ggplot2::theme(
      strip.placement = "inside",
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = rel(1.25),
                                         margin = margin(t = 2, r = 3, b = 2, l = 3)),
      plot.subtitle = ggplot2::element_text(hjust = 0.5)
    )

  if (!is.null(save_path)) {
    ext <- tolower(tools::file_ext(save_path))
    if (ext %in% c("pdf","png","jpeg","jpg","tiff")) {
      ggplot2::ggsave(filename = save_path, plot = gg, width = 20, height = 12, units = "in", dpi = 300)
    } else {
      warning("Unknown file extension for save_path; skipping save.")
    }
  }
  invisible(gg)
}
