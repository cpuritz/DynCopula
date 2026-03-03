###############################################################################

#' Spline estimation of dynamic Gaussian copula model
#'
#' @description Fit a dynamic Gaussian copula model using penalized splines.
#'
#' @param FX Matrix of pseudo-observations.
#' @param design Design matrix. Rows correspond to rows in \code{FX}. A single
#' continuous covariate can be included in a column named \code{"time"}. All
#' other columns are treated as discrete covariates. If no covariates should be
#' included in the model, pass a \code{data.frame} with a column of all ones.
#' @param formula Formula for discrete covariates. Default is \code{~1}.
#' @param lambda Vector of smoothing parameters to test. Default is
#' \code{10^(seq(-5, 5, length.out = 7))}.
#' @param K Dimension of the spline basis matrix. Default is \code{30}.
#' @param nfold Number of folds for cross-validation. Default is \code{5}.
#' @param cores Number of cores to use. Default is \code{1}.
#' @param control A \code{list} of control parameters for optimization.
#'
#' @details If a continuous covariate is included in the design matrix, it is
#' modeled using penalized splines with the smoothing parameter selected via
#' cross-validation. If no continuous covariate is included, then the arguments
#' \code{lambda}, \code{K}, \code{nfold}, and \code{cores} have no effect.
#'
#' Optimization is performed using L-BFGS. The \code{control} argument
#' is a list that supplies control parameters for optimization. The following
#' parameters can be supplied:
#' \itemize{
#'   \item \code{max_itr} Maximum number of iterations. Default is \code{100}.
#'   \item \code{history_size} History size. Default is \code{30}.
#'   \item \code{tolerance_grad} Termination tolerance for gradient. Default is
#'   \code{1e-7}.
#'   \item \code{tolerance_change} Termination tolerance for log-likelihood.
#'   Default is \code{1e-9}.
#'   \item \code{precision}: Floating point precision for calculations. Either
#'   \code{"float32"} or \code{"float64"}. Default is \code{"float64"}.
#'   \item \code{boundary}: Boundary extension for spline knot boundaries.
#'   Default is \code{0.05}.
#' }
#' Any control parameters not specified are replaced by their default values.
#'
#' @return A \code{gamGaussianCopula} object with the following components:
#' \itemize{
#'   \item \code{nobs}: Number of observations.
#'   \item \code{dim}: Dimension of data.
#'   \item \code{time}: Continuous covariate values at which the model was fit.
#'   \item \code{beta}: Matrix of estimated model coefficients.
#'   \item \code{lambda}: The optimal smoothing parameter.
#'   \item \code{cv}: Cross-validation results.
#'   \item \code{B}: Spline basis matrix.
#'   \item \code{continuous}: Whether a continuous covariate was modeled.
#'   \item \code{formula}: The formula for discrete covariates.
#' }
#'
#' @export
fit_dyn_gc <- function(FX,
                       design,
                       formula = ~1,
                       lambda = 10^(seq(-5, 5, length.out = 7)),
                       K = 30,
                       nfold = 5,
                       cores = 1,
                       control = list()) {
    # Basic argument checks
    assert_that(
        is.numeric(FX) && is.matrix(FX),
        is.data.frame(design),
        dim(FX)[1] == dim(design)[1],
        methods::is(formula, "formula"),
        is.vector(lambda, mode = "numeric") && all(lambda > 0),
        is.numeric(K) && K >= 3,
        is.numeric(nfold) && nfold >= 2,
        is.numeric(cores) && cores >= 1,
        is.list(control)
    )

    # Default control parameters
    defaults <- list(
        max_itr = 100L,
        history_size = 30L,
        tolerance_grad = 1e-7,
        tolerance_change = 1e-9,
        precision = "float64",
        boundary = 0.05
    )
    control <- utils::modifyList(defaults, control)
    assert_that(all(names(control) %in% names(defaults)))

    # Verify control parameters
    assert_that(
        all(sapply(control[names(control) != "precision"], is.numeric)),
        control$max_itr >= 1,
        control$history_size >= 1,
        control$tolerance_grad > 0,
        control$tolerance_change > 0,
        control$precision %in% c("float32", "float64"),
        control$boundary >= 0
    )
    control$max_itr <- as.integer(control$max_itr)
    control$history_size <- as.integer(control$history_size)

    # Normal-transform pseudo-observations
    NX <- stats::qnorm(FX)

    # Ensure correct Python interpreter is used
    conf_dir <- rappdirs::user_config_dir("DynCopula")
    conf_path <- file.path(conf_dir, "config.json")
    py_intr <- jsonlite::read_json(conf_path)$python_path
    Sys.setenv(
        RETICULATE_PYTHON = py_intr,
        RETICULATE_AUTOCONFIGURE = "FALSE"
    )

    # Location of Python files
    py_path <- system.file("python", package = "DynCopula")

    # Construct discrete design matrix
    design_disc <- design[, colnames(design) != "time", drop = FALSE]
    missing_vars <- setdiff(all.vars(formula), colnames(design_disc))
    if (length(missing_vars) > 0L) {
        stop("The following variables are missing from 'design': ",
             paste(missing_vars, collapse = ", "))
    }
    Z <- stats::model.matrix(formula, design_disc)

    # If no continuous covariate is specified, fit a GLM and return
    if (!"time" %in% colnames(design)) {
        # Force Python initialization
        reticulate::py_config()

        # Load fit function
        fit_fun <- reticulate::import_from_path(
            module = "gaussian_fit",
            path = py_path,
            delay_load = FALSE
        )$fit_gaussian_linear

        # Fit model
        beta_hat <- fit_fun(NX, Z, control)

        # Return model
        res <- list(
            nobs = dim(FX)[1],
            dim = dim(FX)[2],
            beta = beta_hat,
            continuous = FALSE,
            formula = formula
        )
        class(res) <- "gamGaussianCopula"
        return(res)
    }

    # Whether parallelization is required. CV can be run with only one core,
    # but parallelization is over lambda so no parallelization is done if CV
    # is not going to be run.
    cores <- as.integer(cores)
    run_cv <- (length(lambda) > 1)
    run_parallel <- (cores > 1L) && run_cv

    # Sort continuous covariate and scale to [0, 1]
    x <- design[["time"]]
    if (is.unsorted(x)) {
        ord <- order(x)
        x <- x[ord]
        FX <- FX[ord, , drop = FALSE]
        design <- design[ord, , drop = FALSE]
    }
    min_x <- x[1]
    dx <- x[length(x)] - min_x
    x <- (x - min_x) / dx

    # Spline basis matrix
    B <- splines::ns(
        x = x,
        df = K,
        intercept = FALSE,
        Boundary.knots = c(-control$boundary, 1 + control$boundary)
    )

    # Since Python functions are not serializable, we can't load the
    # optimization functions in the global environment. Instead, the functions
    # needs to be loaded on each worker. This environment stores the
    # Python module once it has been loaded to avoid having to do it multiple
    # times on the same worker.
    .fit_env <- new.env(parent = emptyenv())
    fit_fun <- function(lam) {
        # Load the module if it hasn't been loaded yet
        if (!exists("module", envir = .fit_env, inherits = FALSE)) {
            .fit_env$module <- reticulate::import_from_path(
                module = "gaussian_fit",
                path = py_path,
                delay_load = FALSE
            )
        }
        # Everything but lam already exported to workers
        .fit_env$module$fit_gaussian_spline(
            NX = NX,
            B = B,
            Z = Z,
            lam = lam,
            control = control
        )
    }

    cv_fun <- function(lambda_ix, min_test_ix, max_test_ix) {
        # Load the module if it hasn't been loaded yet
        if (!exists("module", envir = .fit_env, inherits = FALSE)) {
            .fit_env$module <- reticulate::import_from_path(
                module = "gaussian_fit",
                path = py_path,
                delay_load = FALSE
            )
        }
        # NX, B, Z, control already exported to workers
        .fit_env$module$gaussian_spline_cv(
            NX = NX,
            B = B,
            Z = Z,
            lam = lambda[lambda_ix],
            control = control,
            min_test_ix = min_test_ix - 1L,  # convert to 0-indexing
            max_test_ix = max_test_ix - 1L
        )
    }

    if (run_cv) {
        # N-fold cross validation with equal-sized contiguous blocks
        nfold <- as.integer(nfold)
        fold_ids <- cut(seq_along(x), breaks = nfold, labels = FALSE)
        # All combinations of smoothing parameter and test fold ID
        combs <- expand.grid(
            lambda_ix = seq_along(lambda),
            fold = seq_len(nfold)
        )
        ncomb <- dim(combs)[1]
        if (ncomb < cores) {
            message("NOTE: ", cores, " cores have been requested, but there ",
                    "are only ", ncomb, " tasks to run. Only using ", ncomb,
                    " cores.")
            cores <- ncomb
        }
    }

    if (!run_parallel) {
        # Force Python initialization
        reticulate::py_config()
    } else {
        # Create cluster
        cl <- parallel::makeCluster(cores)
        # These variables never change and will be needed for all Python
        # function calls, so we'll export them now.
        parallel::clusterExport(
            cl = cl,
            varlist = c("NX", "B", "lambda", "Z", "fold_ids", "combs",
                        "control", "py_path", ".fit_env"),
            envir = environment()
        )
        # Initialize cluster
        parallel::clusterEvalQ(cl, {
            library(reticulate)

            # Force Python initialization on workers
            reticulate::py_config()

            # Disable torch multithreading to prevent oversubscription
            torch <- reticulate::import("torch", delay_load = FALSE)
            torch$set_num_interop_threads(1L)
            torch$set_num_threads(1L)

            NULL
        })

        # Set up futures plan
        future::plan(future::cluster, workers = cl)
        on.exit({
            future::plan(future::sequential)
            parallel::stopCluster(cl)
        }, add = TRUE)
    }

    if (run_cv) {
        # Avoid setting up futures if no parallelization is requested
        lfun <- ifelse(run_parallel, future.apply::future_lapply, lapply)
        largs <- list(X = seq_len(ncomb))

        progressr::with_progress({
            # Add an extra step to reflect the final model fitting after CV
            pbar <- progressr::progressor(along = seq_len(ncomb + 1L))

            if (run_parallel) {
                largs <- c(largs, list(future.seed = TRUE))
            }

            # Model selection via cross-validation
            ll_fun <- function(i) {
                # Since the folds are contiguous blocks, we can save resources
                # by only passing the start/end points for the testing block
                test_ix <- which(fold_ids == combs$fold[i])
                ll <- cv_fun(
                    lambda_ix = combs$lambda_ix[i],
                    min_test_ix = min(test_ix),
                    max_test_ix = max(test_ix)
                )
                pbar()
                return(ll)
            }
            cv <- do.call(what = lfun, args = c(largs, list(FUN = ll_fun)))
            cv <- unlist(cv)

            # Sum likelihoods across folds
            ll_sum <- tapply(
                X = cv,
                INDEX = list(combs$lambda_ix),
                FUN = sum
            )
            cv_df <- as.data.frame(as.table(ll_sum))
            names(cv_df) <- c("lambda", "ll")

            # Convert factor IDs to indices and then to values
            cv_df$lambda <- lambda[as.integer(as.character(cv_df$lambda))]

            # Optimal value of lambda
            lambda_opt <- cv_df$lambda[which.max(cv_df$ll)]
        })
    } else {
        # Only one lambda value passed, so no cross-validation needed
        lambda_opt <- lambda
        cv_df <- NULL
    }

    # Fit model using the optimal smoothing parameter
    beta_hat <- fit_fun(lambda_opt)

    # Rescale covariates back to their original scale
    x <- x * dx + min_x

    res <- list(
        nobs = dim(FX)[1],
        dim = dim(FX)[2],
        time = x,
        beta = beta_hat,
        lambda = lambda_opt,
        cv = cv_df,
        B = B,
        continuous = TRUE,
        formula = formula
    )
    class(res) <- "gamGaussianCopula"
    return(res)
}

###############################################################################

#' Print method for GAM Gaussian copula models
#'
#' @param x Object of class \code{gamGaussianCopula}.
#' @param ... Additional arguments.
#'
#' @export
#' @method print gamGaussianCopula
print.gamGaussianCopula <- function(x, ...) {
    l1 <- paste0(x$dim, "-dimensional GAM Gaussian copula fit")
    l2 <- paste("nobs =", x$nobs)
    cat(paste(l1, l2, sep = "\n"))
}

###############################################################################
