###############################################################################

#' Fit a dynamic Gaussian copula model
#'
#' @description Fit a dynamic Gaussian copula model.
#'
#' @param FX Matrix of pseudo-observations at covariate values.
#' @param x Vector of covariate values corresponding to \code{FX}. Must be
#' sorted and have no duplicates.
#' @param x0 Covariate values to estimate copula parameters at.
#' @param h Kernel bandwidth. Must satisfy \code{0 < h < 1}.
#' @param degree Degree of local polynomial approximation. Default is \code{0}.
#' @param control A \code{list} of control parameters for optimization.
#' @param cores Number of cores to use for parallel optimization. Parallelized
#' over \code{x0}. Default is \code{1}.
#'
#' @details Optimization is performed using L-BFGS. The \code{control} argument
#' is a list that supplies control parameters for optimization. The following
#' parameters can be supplied:
#' \itemize{
#'   \item \code{max_outer} Maximum number of outer iterations. Default is
#'   \code{1}.
#'   \item \code{max_itr} Maximum number of inner iterations. Default is
#'   \code{100}.
#'   \item \code{history_size} History size. Default is \code{30}.
#'   \item \code{tolerance_grad} Termination tolerance for gradient. Default is
#'   \code{1e-7}.
#'   \item \code{tolerance_change} Termination tolerance for log-likelihood.
#'   Default is \code{1e-9}.
#' }
#'
#' @return A list with the following components:
#' \itemize{
#'   \item \code{x}: The input argument \code{x}.
#'   \item \code{x0}: The input argument \code{x0}.
#'   \item \code{h}: The input argument \code{h}.
#'   \item \code{NX}: Normal-transformed pseudo-observations.
#'   \item \code{eta}: Matrix of estimated calibration coefficients.
#'   \item \code{rho}: Matrix of estimated pairwise correlation coefficients.
#' }
#'
#' @export
fit_dynamic_gaussian <- function(FX,
                                 x,
                                 x0,
                                 h,
                                 degree = 0L,
                                 control = list(),
                                 cores = 1L) {
    # Basic checks
    assert_that(
        is.vector(x, mode = "numeric"),
        is.matrix(FX),
        is.numeric(FX),
        dim(FX)[1] == length(x),
        dim(FX)[2] > 1L,
        is.numeric(h) && h > 0 && h < 1,
        is.numeric(cores) && cores >= 1,
        is.list(control),
        !anyDuplicated(x) && !is.unsorted(x),
        is.numeric(x0),
        is.numeric(degree) && degree >= 0
    )
    cores <- as.integer(cores)
    degree <- as.integer(degree)

    # Default control parameters
    defaults <- list(
        max_outer = 1L,
        max_itr = 100L,
        history_size = 30L,
        tolerance_grad = 1e-7,
        tolerance_change = 1e-9
    )
    control <- utils::modifyList(defaults, control)
    assert_that(all(names(control) %in% names(defaults)))

    # Verify control parameters
    assert_that(
        all(sapply(control, is.numeric)),
        control$max_outer >= 1,
        control$max_itr >= 1,
        control$history_size >= 1,
        control$tolerance_grad > 0,
        control$tolerance_change > 0
    )
    control$max_outer <- as.integer(control$max_outer)
    control$max_itr <- as.integer(control$max_itr)
    control$history_size <- as.integer(control$history_size)

    # Scale covariates to [0, 1]
    min_x <- x[1]
    dx <- x[length(x)] - min_x
    x <- (x - min_x) / dx
    x0 <- (x0 - min_x) / dx

    # Normal-transform pseudo-observations
    NX <- stats::qnorm(FX)

    if (cores > 1L) {
        # Set up futures plan
        cl <- parallel::makeCluster(cores)
        future::plan(future::cluster, workers = cl)
        on.exit({ future::plan(future::sequential); parallel::stopCluster(cl) },
                add = TRUE)

        # Since Python functions are not serializable, we can't load the
        # optimization function in the global environment. Instead, the function
        # needs to be loaded on each worker. This environment stores the
        # function once it has been loaded to avoid having to do it multiple
        # times on the same worker.
        .fit_env <- new.env(parent = emptyenv())
        fit_gaussian <- function(par0, x, NX, h, control, x0) {
            # Load the module if it hasn't been loaded yet
            if (!exists("fit_fun", envir = .fit_env, inherits = FALSE)) {
                .fit_env$fit_fun <- reticulate::import_from_path(
                    module = "dynamic_gaussian",
                    path = system.file("python", package = "DynCopula"),
                    delay_load = FALSE
                )$fit_gaussian
            }
            .fit_env$fit_fun(
                par0 = par0,
                x = x,
                NX = NX,
                h = h,
                control = control,
                x0 = x0,
                degree = degree
            )
        }

        # Parallelized with progress bar
        progressr::with_progress({
            pbar <- progressr::progressor(along = x0)
            eta_est <- future.apply::future_lapply(
                X = x0,
                FUN = function(t0) {
                    y <- fit_gaussian(
                        par0 = .init_par(t0, x, NX, h),
                        x = x,
                        NX = NX,
                        h = h,
                        control = control,
                        x0 = t0
                    )
                    pbar()
                    return(y)
                },
                future.seed = TRUE,
                future.globals = TRUE
            )
        })
    } else {
        # Load Python module
        fit_fun <- reticulate::import_from_path(
            module = "dynamic_gaussian",
            path = system.file("python", package = "DynCopula"),
            delay_load = FALSE
        )$fit_gaussian

        progressr::with_progress({
            pbar <- progressr::progressor(along = x0)
            eta_est <- lapply(
                X = x0,
                FUN = function(t0) {
                    y <- fit_fun(
                        par0 = .init_par(t0, x, NX, h),
                        x = x,
                        NX = NX,
                        h = h,
                        control = control,
                        x0 = t0,
                        degree = degree
                    )
                    pbar()
                    return(y)
                }
            )
        })
    }

    # Matrix of estimated calibration coefficients
    Hhat <- do.call(rbind, eta_est)
    # Matrix of estimated correlation coefficients
    Rhat <- t(apply(Hhat, 1, function(v) {
        copula::P2p(vec2cor(v))
    }))

    # Ensure consistent shape of Rhat across all dimensions
    d <- dim(FX)[2]
    if (d == 2) {
        Rhat <- t(Rhat)
    }

    # Add numbered eta/rho labels
    colnames(Hhat) <- paste0("eta", seq_len(choose(d, 2)))
    ix_lab <- apply(utils::combn(seq_len(d), 2), 2, function(x) {
        paste(x, collapse = '_')
    })
    colnames(Rhat) <- paste0("rho", ix_lab)

    # Rescale covariates back to their original scale
    x <- x * dx + min_x
    x0 <- x0 * dx + min_x

    return(list(
        x = x,
        x0 = x0,
        h = h,
        NX = NX,
        eta = Hhat,
        rho = Rhat
    ))
}

###############################################################################
