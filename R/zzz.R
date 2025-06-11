.onLoad <- function(libname, pkgname) {
    config_path <- file.path(rappdirs::user_config_dir(pkgname), "config.json")
    if (file.exists(config_path)) {
        config <- jsonlite::read_json(config_path)
        if (!is.null(config$envname) && file.exists(config$envname)) {
            reticulate::use_virtualenv(config$envname, required = TRUE)
        }
    }
}

.onAttach <- function(libname, pkgname) {
    if ("devtools" %in% loadedNamespaces()) return()
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
