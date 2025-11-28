###############################################################################

#' Build Python virtual environment
#'
#' @description Build a Python virtual environment for the package to use. Only
#' needs to be called once after package installation.
#'
#' @param python The path to a Python interpreter.
#'
#' @details The virtual environment is located at
#' \code{rappdirs::user_cache_dir("DynCopula")}. A JSON file is also written to
#' \code{rappdirs::user_config_dir("DynCopula")}.
#'
#' @returns Nothing, called for side effects.
#'
#' @export
pkg_setup <- function(python = reticulate::virtualenv_starter(NULL)) {
    #python <- "/opt/homebrew/opt/python@3.9/bin/python3.9"
    reticulate::use_python(python, required = TRUE)

    # Path to virtual environment
    cache_dir <- rappdirs::user_cache_dir("DynCopula")
    envname <- file.path(cache_dir, "venv", "r-dyncopula")

    # Create virtual environment
    reticulate::virtualenv_create(
        envname = envname,
        python = python,
        force = TRUE,
        packages = c("numpy", "torch")
    )

    # Link headers
    # inc <- reticulate::import("sysconfig")$get_config_var("INCLUDEPY")
    # unlink(file.path(envname, "Headers"), recursive = TRUE, force = TRUE)
    # file.symlink(inc, file.path(envname, "Headers"))

    # Record path to virtual environment in config file
    config_dir <- rappdirs::user_config_dir("DynCopula")
    if (!dir.exists(config_dir)) {
        dir.create(config_dir)
    }
    jsonlite::write_json(
        x = list(envname = jsonlite::unbox(envname)),
        path = file.path(config_dir, "config.json")
    )
    message("Setup complete! Please restart R now.")
}

###############################################################################
