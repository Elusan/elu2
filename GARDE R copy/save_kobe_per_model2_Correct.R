#' Save individual Kobe plots per model (robust; no blank PNGs)
#'
#' @description
#' For each requested scenario in `all_models`, builds a Kobe ggplot for every
#' fitted model and saves to `FIG/per_model/<scenario>/<model>.png`.
#'
#' Robustness:
#' * Builder cascade: `kobe_safe()` → `Original_kobe_all_in_gg()` → `kobe_all_in_one_gg()`
#'   → minimal placeholder ggplot (never NULL).
#' * Saver cascade: **ragg::agg_png** → **ggplot2::ggsave(device = "png")**
#'   → **grDevices::png** with `print()`. No silent white images.
#' * Filenames sanitized, dirs auto-created, errors downgraded to warnings.
#'
#' @param all_models Named list of scenarios; each scenario is a named list of
#'   fitted SPiCT objects (e.g. `list(S1 = list(S1P.SDM = rep1, S1F.SDM = rep2))`).
#' @param scenario_names Which scenarios to process. Default: `names(all_models)`.
#'   Missing names are skipped with a warning (never stops the run).
#' @param out_dir Root directory. Per-scenario subdirs created automatically.
#'   Default `file.path("FIG","per_model")`.
#' @param width,height PNG size in inches. Default `6, 5`.
#' @param dpi PNG resolution. Default `300`.
#' @param bg Background color. Default `"white"`.
#' @param add_id Add a top-left ID (model name) inside the panel when possible.
#'   Default `TRUE`.
#' @param id_size Size for the ID label (ggplot text size). Default `4`.
#' @param overwrite Overwrite existing files. Default `TRUE`.
#' @param builder One of `"auto"`, `"kobe_safe"`, `"Original_kobe_all_in_gg"`,
#'   `"kobe_all_in_one_gg"`. `"auto"` tries in that priority order. Default `"auto"`.
#' @param filename_fun Optional function `(scenario, model_name) -> filename`
#'   (without extension) to customize file names.
#' @param verbose Chatty progress. Default `TRUE`.
#' @param ... Extra arguments forwarded to the builder (e.g., `rel.axes`, `CI`, `logax`).
#'
#' @return Invisibly a named list of character vectors (file paths) per scenario.
#'   A run log is attached as `attr(result, "log")` with columns:
#'   `scenario, model, path, ok, note`.
#'
#' @examples
#' \dontrun{
#' res <- save_kobe_per_model(
#'   all_models,
#'   scenario_names = names(all_models),
#'   rel.axes = FALSE, CI = 0.95
#' )
#' attr(res, "log")
#' }
#' @import ggplot2
#' @importFrom grDevices png dev.off
#' @importFrom utils getAnywhere
#' @export
save_kobe_per_model2_Correct <- function(all_models,
                                scenario_names = names(all_models),
                                out_dir = file.path("FIG", "per_model"),
                                width = 6, height = 5, dpi = 300, bg = "white",
                                add_id = TRUE, id_size = 4,
                                overwrite = TRUE,
                                builder = c("auto","kobe_safe","Original_kobe_all_in_gg","kobe_all_in_one_gg"),
                                filename_fun = NULL,
                                verbose = TRUE,
                                ...) {
  builder <- match.arg(builder)

  # ----------------- tiny utilities (no throws) -----------------
  have_gg   <- function() isTRUE(requireNamespace("ggplot2", quietly = TRUE))
  have_ragg <- function() isTRUE(requireNamespace("ragg",     quietly = TRUE))
  say <- function(...) if (isTRUE(verbose)) message(...)
  nz1 <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  ensure_dir <- function(path) { if (!dir.exists(path)) dir.create(path, TRUE, FALSE); dir.exists(path) }
  sanitize <- function(x) { x <- ifelse(is.na(x) | !nzchar(x), "unnamed", x); gsub("[^[:alnum:]_.-]+", "_", x, perl = TRUE) }
  exists_fun <- function(fname) exists(fname, mode = "function", inherits = TRUE) || !is.null(utils::getAnywhere(fname)$objs)
  get_fun <- function(fname) if (exists_fun(fname)) get(fname, inherits = TRUE) else NULL
  mk_filename <- function(scen, model) {
    base <- if (is.function(filename_fun)) {
      out <- try(filename_fun(scen, model), silent = TRUE)
      if (inherits(out, "try-error") || !nz1(out)) sanitize(model) else sanitize(as.character(out))
    } else sanitize(model)
    paste0(base, ".png")
  }

  # ----------------- plot builders (cascade) --------------------
  build_with_kobe_safe <- function(rep_obj, ...) {
    f <- get_fun("kobe_safe"); if (is.null(f)) return(NULL)
    tryCatch(f(rep_obj, ...), error = function(e) NULL)
  }
  build_with_original <- function(rep_obj, ...) {
    f <- get_fun("Original_kobe_all_in_gg"); if (is.null(f)) return(NULL)
    tryCatch({ p <- f(rep = rep_obj, ...); p }, error = function(e) NULL)
  }
  build_with_modern <- function(rep_obj, ...) {
    f <- get_fun("kobe_all_in_one_gg"); if (is.null(f)) return(NULL)
    tryCatch({ p <- f(rep = rep_obj, ...); p }, error = function(e) NULL)
  }
  placeholder_plot <- function(msg = "Kobe plot unavailable") {
    if (!have_gg()) return(NULL)
    ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = msg,
                        colour = "red", fontface = "bold", size = 4) +
      ggplot2::xlim(-1, 1) + ggplot2::ylim(-1, 1) +
      ggplot2::theme_void()
  }
  add_id_if_possible <- function(p, id_text, id_size) {
    if (!inherits(p, "ggplot")) return(p)
    af <- get_fun("add_id_label")
    if (!is.null(af)) {
      return(tryCatch(af(p, id_text, id_size = id_size), error = function(e) p))
    }
    # tiny inline annotate if user helper not available
    tryCatch({
      p + ggplot2::annotate("text", x = -Inf, y = Inf, label = id_text,
                            hjust = -0.2, vjust = 2, fontface = "bold", size = id_size) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::theme(plot.margin = ggplot2::margin(6, 6, 6, 6))
    }, error = function(e) p)
  }

  build_plot <- function(rep_obj, ...) {
    order <- switch(builder,
                    "kobe_safe" = "kobe_safe",
                    "Original_kobe_all_in_gg" = "Original_kobe_all_in_gg",
                    "kobe_all_in_one_gg" = "kobe_all_in_one_gg",
                    "auto" = c("kobe_safe","Original_kobe_all_in_gg","kobe_all_in_one_gg")
    )
    p <- NULL
    for (nm in order) {
      p <- switch(nm,
                  "kobe_safe"              = build_with_kobe_safe(rep_obj, ...),
                  "Original_kobe_all_in_gg"= build_with_original(rep_obj, ...),
                  "kobe_all_in_one_gg"     = build_with_modern(rep_obj, ...),
                  NULL
      )
      if (inherits(p, "ggplot")) break
    }
    # guarantee a ggplot (prevents white device with nothing drawn)
    if (!inherits(p, "ggplot")) p <- placeholder_plot()
    p
  }

  # ----------------- file saving (no blank outputs) -------------
  save_plot_png <- function(plot_obj, path, width, height, dpi, bg) {
    # 1) ragg is best and avoids Quartz/Cairo quirks on macOS
    if (have_ragg() && inherits(plot_obj, "ggplot")) {
      ok <- try({
        dev <- ragg::agg_png(filename = path, width = width, height = height, units = "in", res = dpi, background = bg)
        on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
        print(plot_obj) # IMPORTANT: draw content to device
        TRUE
      }, silent = TRUE)
      if (!inherits(ok, "try-error")) return(TRUE)
    }

    # 2) ggsave (device = "png")
    if (have_gg() && inherits(plot_obj, "ggplot")) {
      ok <- try({
        ggplot2::ggsave(filename = path, plot = plot_obj, device = "png",
                        width = width, height = height, units = "in",
                        dpi = dpi, bg = bg)
        TRUE
      }, silent = TRUE)
      if (!inherits(ok, "try-error")) return(TRUE)
    }

    # 3) last resort: base device + print()
    ok <- try({
      grDevices::png(filename = path, width = width, height = height,
                     units = "in", res = dpi, type = if (.Platform$OS.type == "windows") "windows" else "quartz", bg = bg)
      on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
      if (inherits(plot_obj, "ggplot")) {
        print(plot_obj)  # CRITICAL: ensure plot is drawn
      } else {
        par(mar = c(0, 0, 0, 0)); plot.new(); rect(0,0,1,1, col = bg, border = NA)
        text(0.5, 0.55, "Kobe plot unavailable", cex = 1.1, font = 2)
      }
      TRUE
    }, silent = TRUE)
    !inherits(ok, "try-error")
  }

  # ----------------- validate inputs & prepare ------------------
  if (!is.list(all_models) || !length(all_models)) {
    warning("`all_models` is empty or not a list; nothing to do.")
    return(invisible(structure(list(), log = data.frame())))
  }
  scenario_names <- unique(as.character(stats::na.omit(scenario_names)))
  have <- intersect(scenario_names, names(all_models))
  miss <- setdiff(scenario_names, names(all_models))
  if (length(miss)) warning("Skipping missing scenarios: ", paste(miss, collapse = ", "))

  if (!ensure_dir(out_dir)) warning("Could not create/confirm out_dir: ", out_dir)

  # ----------------- main loop ---------------------------------
  res_list <- vector("list", length(have)); names(res_list) <- have
  log_rows <- list()

  for (scen in have) {
    sc_dir <- file.path(out_dir, sanitize(scen))
    if (!ensure_dir(sc_dir)) {
      warning("Cannot create scenario dir: ", sc_dir, ". Using out_dir instead.")
      sc_dir <- out_dir
    }

    scen_list <- all_models[[scen]]
    if (!is.list(scen_list) || !length(scen_list)) {
      warning("Empty scenario '", scen, "'. Skipping.")
      res_list[[scen]] <- character(0)
      next
    }

    fpaths <- character(length(scen_list))
    names(fpaths) <- names(scen_list)

    for (i in seq_along(scen_list)) {
      mname <- names(scen_list)[i]; if (!nz1(mname)) mname <- paste0("model_", i)
      rep_obj <- scen_list[[i]]

      p <- build_plot(rep_obj, ...)
      if (isTRUE(add_id) && inherits(p, "ggplot")) {
        p <- add_id_if_possible(p, mname, id_size)
      }

      fname <- mk_filename(scen, mname)
      fpath <- file.path(sc_dir, fname)

      if (file.exists(fpath) && !overwrite) {
        say("Skip existing: ", fpath)
        fpaths[i] <- fpath
        log_rows[[length(log_rows) + 1L]] <- data.frame(
          scenario = scen, model = mname, path = fpath, ok = TRUE,
          note = "exists (kept)", stringsAsFactors = FALSE
        )
        next
      }

      ok <- save_plot_png(p, fpath, width, height, dpi, bg)
      if (!ok) {
        # absolute last-ditch: try to at least write a placeholder
        ok2 <- try({
          grDevices::png(filename = fpath, width = width, height = height,
                         units = "in", res = dpi, bg = bg)
          on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
          par(mar = c(0, 0, 0, 0)); plot.new(); rect(0,0,1,1, col = bg, border = NA)
          text(0.5, 0.62, "Kobe plot could not be built", cex = 1.05, font = 2)
          text(0.5, 0.48, paste0("Scenario: ", scen), cex = 0.9)
          text(0.5, 0.40, paste0("Model: ", mname), cex = 0.9)
          TRUE
        }, silent = TRUE)
        ok <- !inherits(ok2, "try-error")
      }

      fpaths[i] <- fpath
      log_rows[[length(log_rows) + 1L]] <- data.frame(
        scenario = scen, model = mname, path = fpath, ok = isTRUE(ok),
        note = if (isTRUE(ok)) "ok" else "failed (placeholder)", stringsAsFactors = FALSE
      )
      if (isTRUE(ok)) say("Wrote: ", fpath) else warning("Problem writing: ", fpath)
    }

    res_list[[scen]] <- fpaths
  }

  # build run log
  log_df <- if (length(log_rows)) do.call(rbind, log_rows) else
    data.frame(scenario = character(), model = character(), path = character(),
               ok = logical(), note = character(), stringsAsFactors = FALSE)

  attr(res_list, "log") <- log_df
  invisible(res_list)
}
