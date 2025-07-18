.onLoad <- function(libname, pkgname) {
    # Load the package Python environment
    config_path <- file.path(rappdirs::user_config_dir(pkgname), "config.json")
    if (file.exists(config_path)) {
        config <- jsonlite::read_json(config_path)
        if (!is.null(config$envname) && file.exists(config$envname)) {
            reticulate::use_virtualenv(config$envname, required = TRUE)
        }
    }
}

.onAttach <- function(libname, pkgname) {
    # Don't print startup messages during development
    if ("devtools" %in% loadedNamespaces()) return()
    # Print message if the package Python environment cannot be found
    config_path <- file.path(rappdirs::user_config_dir(pkgname), "config.json")
    msg <- paste0("Python environment not yet configured. ", "Run '", pkgname,
                  "::setup()' to initialize.")
    if (!file.exists(config_path)) {
        packageStartupMessage(msg)
    } else {
        config <- jsonlite::read_json(config_path)
        if ((is.null(config$envname) || !file.exists(config$envname))) {
            packageStartupMessage(msg)
        }
    }
}
