###############################################################################

.onLoad <- function(libname, pkgname) {
    config_path <- file.path(
        rappdirs::user_config_dir(pkgname),
        "config.json"
    )
    if (file.exists(config_path)) {
        envname <- jsonlite::read_json(config_path)$envname
        if (!is.null(envname)) {
            reticulate::use_virtualenv(envname, required = FALSE)
        }
    }
}

###############################################################################

.onAttach <- function(libname, pkgname) {
    # Don't print startup messages during development
    if ("devtools" %in% loadedNamespaces()) {
        return()
    }

    config_path <- file.path(
        rappdirs::user_config_dir(pkgname),
        "config.json"
    )
    msg <- paste0(pkgname, " is not configured. Run pkg_setup().")
    if (!file.exists(config_path)) {
        packageStartupMessage(msg)
        return()
    }

    envname <- jsonlite::read_json(config_path)$envname
    if (is.null(envname)) {
        # Will only trigger if config file got corrupted
        packageStartupMessage(msg)
    }
}

###############################################################################
