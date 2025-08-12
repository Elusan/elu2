#' Save individuals and grids together (one-liner helper)
#' @export
save_kobe_everything <- function(all_models,
                                 out_root_individuals = file.path("FIG","KobeIndividuals"),
                                 out_root_grids = file.path("FIG","KobePhases"),
                                 ...) {
  log_df <- save_kobe_individuals(all_models, out_root = out_root_individuals, ...)
  invisible(build_kobe_grids_scenario(all_models,
                                      layout = "3x2",
                                      save = TRUE,
                                      out_dir = out_root_grids,
                                      ...))
  return(log_df)
}
