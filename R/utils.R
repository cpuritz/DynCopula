#' Load python functions
#'
#' Utility function to lazily load necessary python functions.
py_load <- function() {
    reticulate::import_from_path(
        "dynamic_fit",
        path = system.file("python", package = "DynCopula"),
        delay_load = TRUE
    )
}

#' Builds Python virtual environment
#'
#' Builds Python virtual environment
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
