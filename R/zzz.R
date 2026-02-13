###############################################################################

.onLoad <- function(libname, pkgname) {
    # Select the correct Python interpreter
    config_path <- file.path(rappdirs::user_config_dir(pkgname), "config.json")
    if (file.exists(config_path)) {
        config <- jsonlite::read_json(config_path)
        py_path <- config$python_path
        if (!is.null(py_path) && file.exists(py_path)) {
            reticulate::use_python(py_path, required = TRUE)
        }
    }
}

###############################################################################

.onAttach <- function(libname, pkgname) {
    # Don't print startup messages during development
    if ("devtools" %in% loadedNamespaces()) {
        return()
    }

    # Print message if a Python interpreter has not been specified yet
    config_path <- file.path(rappdirs::user_config_dir(pkgname), "config.json")
    msg <- paste0("Python interpreter not yet specified. Run '", pkgname,
                  "::pkg_setup()'.")
    if (!file.exists(config_path)) {
        packageStartupMessage(msg)
    } else {
        config <- jsonlite::read_json(config_path)
        # Should never get here, but check to be safe
        if ((is.null(config$path) || !file.exists(config$path))) {
            packageStartupMessage(msg)
        }
    }
}

###############################################################################
