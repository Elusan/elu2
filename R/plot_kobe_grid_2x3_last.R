
#' 2×3 Kobe grid (rows = data type; cols = model) + in-panel tags + optional save
#'
#' @description
#' Builds a 2×3 grid of Kobe panels using your existing
#' \code{my_plot_kobe_all_management_scenario()}:
#'   Row 1 = SDM (Pella, Schaefer, Fox)
#'   Row 2 = GLM (Pella, Schaefer, Fox)
#'
#' Each cell shows a small bold tag at the top-left with the original list name
#' (e.g., "S1P.SDM"). No plot titles are added.
#'
#' @param models_list Named list for one scenario containing up to six model fits:
#'   SxP.SDM, SxS.SDM, SxF.SDM, SxP.GLM, SxS.GLM, SxF.GLM (names are used).
#' @param man.legend Logical; passed to \code{my_plot_kobe_all_management_scenario()}.
#' @param text_size Numeric; size of the in-panel tag text.
#' @param save Logical; if TRUE, saves a PNG to \code{out_dir}.
#' @param out_dir Character; output directory for the PNG. Created if missing.
#' @param width,height Numeric; device size in inches for saving.
#' @param dpi Numeric; resolution for saving.
#' @param file_basename Optional base name (without extension). If NULL, derived
#'   from the first list name (e.g., "S1_KobeGrid_2x3").
#' @param ... Extra args forwarded to \code{my_plot_kobe_all_management_scenario()}.
#'
#' @return Invisibly returns a list with elements:
#'   \item{plot}{the patchwork object}
#'   \item{file}{path to the saved file (or NULL if not saved)}
#' @export
plot_kobe_grid_2x3_last <- function(models_list,
                                    man.legend   = FALSE,
                                    text_size    = 4,
                                    save         = TRUE,
                                    out_dir      = "FIG/KobePhasesNEW",
                                    width        = 20,
                                    height       = 9,
                                    dpi          = 300,
                                    file_basename = NULL,
                                    ...) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  data_types <- c("SDM", "GLM")      # rows
  mod_keys   <- c("P", "S", "F")     # cols
  mod_names  <- c(P = "Pella", S = "Schaefer", F = "Fox")

  .first_or_null <- function(x) { if (length(x) >= 1) x[[1]] else NULL }

  .find_entry_name <- function(models_list, key, dtype) {
    pat <- paste0(key, "\\.", dtype, "$")   # e.g., "P.SDM$"
    nm  <- names(models_list)
    hit <- nm[grep(pat, nm)]
    .first_or_null(hit)
  }

  .blank_panel <- function(tag) {
    p0 <- ggplot2::ggplot() +
      ggplot2::geom_blank() +
      ggplot2::expand_limits(x = c(0, 1), y = c(0, 1)) +
      ggplot2::theme_void() +
      ggplot2::theme(
        panel.border = ggplot2::element_rect(color = "grey80", fill = NA, linewidth = 0.6),
        plot.margin  = ggplot2::margin(4, 4, 4, 4)
      )
    p0 <- p0 +
      ggplot2::annotate("text",
                        x = -Inf, y =  Inf, label = tag,
                        hjust = -0.1, vjust = 2,
                        fontface = "bold", size = text_size) +
      ggplot2::coord_cartesian(clip = "off")
    p0
  }

  .add_inpanel_tag <- function(p, tag) {
    p +
      ggplot2::annotate("text",
                        x = -Inf, y =  Inf, label = tag,
                        hjust = -0.1, vjust = 2,
                        fontface = "bold", size = text_size) +
      ggplot2::coord_cartesian(clip = "off")
  }

  plots <- vector("list", length(data_types) * length(mod_keys))
  idx <- 0
  r <- 1
  while (r <= length(data_types)) {
    dtype <- data_types[r]
    c <- 1
    while (c <= length(mod_keys)) {
      key <- mod_keys[c]
      idx <- idx + 1

      nm <- .find_entry_name(models_list, key, dtype)

      if (!is.null(nm)) {
        rep_obj <- models_list[[nm]]
        p <- my_plot_kobe_all_management_scenario(rep = rep_obj,
                                                  man.legend = man.legend, ...)
        p <- .add_inpanel_tag(p, nm)
      } else {
        fallback <- paste0("(", dtype, " – ", mod_names[[key]], ")")
        p <- .blank_panel(fallback)
      }

      plots[[idx]] <- p
      c <- c + 1
    }
    r <- r + 1
  }

  g_row1 <- plots[[1]] | plots[[2]] | plots[[3]]
  g_row2 <- plots[[4]] | plots[[5]] | plots[[6]]
  g <- g_row1 / g_row2

  out_path <- NULL
  if (isTRUE(save)) {
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    # Derive a sensible basename when not provided.
    if (is.null(file_basename)) {
      nm_all <- names(models_list)
      base <- if (length(nm_all) >= 1) {
        # Take the common scenario prefix up to the first dot/group, e.g., "S1" from "S1P.SDM"
        sub("([PSF]).*$", "", nm_all[1])
      } else {
        "Scenario"
      }
      if (!nzchar(base)) base <- "Scenario"
      file_basename <- paste0(base, "KobeGrid_2x3")
    }
    out_path <- file.path(out_dir, paste0(file_basename, ".png"))
    ggplot2::ggsave(filename = out_path, plot = g,
                    width = width, height = height, dpi = dpi, units = "in")
  }

  print(g)
  invisible(list(plot = g, file = out_path))
}
