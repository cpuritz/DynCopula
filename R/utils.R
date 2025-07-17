#' Load Python module
#'
#' @description Utility function to lazily load a Python module.
#'
#' @param module Name of module to load.
py_load <- function(module = c("dynamic_gaussian", "dynamic_t")) {
    module <- match.arg(module)
    reticulate::import_from_path(
        module = module,
        path = system.file("python", package = "DynCopula"),
        delay_load = TRUE
    )

}

#' Build Python virtual environment
#'
#' @description Utility function to build a Python virtual environment.
#'
#' @export
setup <- function() {
    envname <- file.path(rappdirs::user_cache_dir("DynCopula"), "venv",
                         "r-dyncopula")
    reticulate::virtualenv_create(
        envname = envname,
        force = TRUE,
        packages = c("numpy", "torch", "botorch")
    )
    reticulate::use_virtualenv(envname, required = TRUE)

    config_dir <- rappdirs::user_config_dir("DynCopula")
    if (!dir.exists(config_dir)) {
        dir.create(config_dir)
    }
    jsonlite::write_json(
        x = list(envname = jsonlite::unbox(envname)),
        path = file.path(config_dir, "config.json")
    )
    message("Setup complete!")
}
