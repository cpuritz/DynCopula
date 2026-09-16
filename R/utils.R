###############################################################################
#' Setup the package
#'
#' @description Creates a dedicated Python virtual environment for gamgc,
#' installs the required Python dependencies (NumPy and PyTorch), and saves the
#' configuration for future sessions.
#'
#' @param envname Name of the virtual environment to create. Default is
#' \code{"r-gamgc"}. An existing environment with this name will be removed
#' and recreated.
#' @param torch_version PyTorch version constraint passed to pip, e.g.
#' \code{"torch==2.3.0"} or \code{"torch>=2.0"}. Defaults to \code{"torch"}
#' (latest stable).
#'
#' @details This function:
#' \enumerate{
#'   \item Creates a new Python virtual environment managed by reticulate.
#'   \item Installs NumPy and PyTorch (CPU build) into that environment via pip.
#'   \item Verifies the installations by importing both modules.
#'   \item Writes the environment name to a JSON config file in
#'     \code{rappdirs::user_config_dir("gamgc")} so that future sessions
#'     activate it automatically.
#' }
#'
#' Call this function once after installing gamgc. Restart R afterwards.
#' The environment is self-contained and will not affect other Python
#' installations on your system.
#'
#' @returns Nothing, called for side effects.
#'
#' @export
pkg_setup <- function(envname = "r-gamgc", torch_version = "torch") {
    message("gamgc package setup")

    # Create the virtual environment
    message("  Creating Python virtual environment '", envname, "'...")
    if (reticulate::virtualenv_exists(envname)) {
        message("  Existing environment found; removing it before reinstalling.")
        reticulate::virtualenv_remove(envname, confirm = FALSE)
    }
    reticulate::virtualenv_create(envname)

    # Install packages
    message("  Installing Python dependencies (this may take a few minutes)...")
    reticulate::virtualenv_install(
        envname = envname,
        packages = c("numpy", torch_version),
        pip_options = c("--quiet", "--index-url",
                        "https://download.pytorch.org/whl/cpu")
    )

    # -- 3. Verify installations -------------------------------------------------
    message("  Verifying installations...")
    reticulate::use_virtualenv(envname, required = TRUE)

    for (mod in c("numpy", "torch")) {
        if (!reticulate::py_module_available(mod)) {
            stop(
                "Installation of '", mod, "' appears to have failed. Try ",
                "running pkg_setup() again or install manually with:\n",
                "  reticulate::virtualenv_install('", envname, "', '", mod, "')"
            )
        }
    }

    torch <- reticulate::import("torch")
    message("  PyTorch ", torch$`__version__`, " installed successfully.")

    # Save config file
    message("  Saving config...")
    config_dir <- rappdirs::user_config_dir("gamgc")
    if (!dir.exists(config_dir)) {
        dir.create(config_dir, recursive = TRUE)
    }
    jsonlite::write_json(
        x = list(envname = jsonlite::unbox(envname)),
        path = file.path(config_dir, "config.json")
    )

    message("Setup complete! Please restart R now.")
    invisible(NULL)
}
