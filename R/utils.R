###############################################################################

#' Setup the package
#'
#' @description Setup the package by checking that necessary Python modules are
#' installed, and then storing the path to the Python interpreter.
#'
#' @param python_path The path to a Python interpreter.
#'
#' @details The NumPy and PyTorch libraries are required for computations. This
#' function writes a JSON file to \code{rappdirs::user_config_dir("DynCopula")}.
#'
#' @returns Nothing, called for side effects.
#'
#' @export
pkg_setup <- function(python_path) {
    message("DynCopula package setup")
    message("    Loading Python interpreter...")
    reticulate::use_python(python = python_path, required = TRUE)

    # Check that necessary modules are installed
    message("    Verifying that necessary modules are installed...")
    mods <- c("numpy", "torch")
    for (m in mods) {
        if (!reticulate::py_module_available(m)) {
            msg <- paste0("'", m, "' import failed. Ensure that '", m,
                          "' is installed.")
            stop(msg)
        }
    }

    # Record the path to the Python interpreter
    message("    Saving config...")
    config_dir <- rappdirs::user_config_dir("DynCopula")
    if (!dir.exists(config_dir)) {
        dir.create(config_dir)
    }
    jsonlite::write_json(
        x = list(python_path = jsonlite::unbox(python_path)),
        path = file.path(config_dir, "config.json")
    )
    message("Setup complete! Please restart R now.")
}

###############################################################################
