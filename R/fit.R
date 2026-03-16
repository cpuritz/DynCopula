###############################################################################

#' Generalized Additive Model for a Gaussian Copula
#'
#' @description Fit a generalized additive model for a Gaussian copula.
#'
#' @param FX Matrix of pseudo-observations.
#' @param design A \code{data.frame} specifying the design matrix. Rows
#' correspond to rows in \code{FX}. If no covariates should be included in the
#' model, pass a \code{data.frame} with a column of all ones.
#' @param formula Formula for covariates.
#' @param lambda Vector of penalty parameters to test. Default is
#' \code{10^(seq(-5, 5, length.out = 7))}.
#' @param K Dimension of the spline basis matrix. Default is \code{30}.
#' @param nfold Number of folds for cross-validation. Default is \code{5}.
#' @param cores Number of cores to use. Default is \code{1}.
#' @param control A \code{list} of control parameters for optimization.
#'
#' @details The formula can include a single smooth covariate and any number of
#' linear covariates. Specify the smooth covariate by \code{s(t)}, replacing
#' \code{t} with the name of the smooth covariate in the design matrix. If
#' interaction terms between linear covariates and the smooth covariate are
#' included, a baseline (intercept) smooth function will be included in the
#' model.
#'
#' If a smooth covariate is included in the formula, it is modeled
#' using penalized splines with a roughness penalty. The linear covariates are
#' penalized by an L2 penalty. The penalty parameters are selected via
#' cross-validation. A separate penalty parameter is used to for each term
#' which includes the smooth covariate. This includes a baseline smooth
#' function as well as any interactions between the smooth covariate and
#' linear covariates.
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
#' @returns A \code{gamGaussianCopula} object with the following components:
#' \itemize{
#'   \item \code{nobs}: Number of observations.
#'   \item \code{dim}: Dimension of data.
#'   \item \code{time}: Smooth covariate values at which the model was fit.
#'   \item \code{beta}: Matrix of estimated model coefficients.
#'   \item \code{lambda}: The optimal smoothing parameter.
#'   \item \code{cv}: Cross-validation results.
#'   \item \code{B}: Spline basis matrix.
#'   \item \code{smooth_name}: The name of the smooth covariate.
#'   \item \code{Z_names}: Column names of the design matrix for linear
#'   covariates.
#'   \item \code{M_names}: Column names of the design matrix for interactions
#'   between the smooth
#'   covariate and linear covariates.
#'   \item \code{lin_formula}: Formula for linear covariates.
#'   \item \code{int_formula}: Formula for interactions between the smooth
#'   covariate and linear covariates.
#'   \item \code{colnames}: Column names of the pseudo-observation matrix.
#' }
#'
#' @examples
#' \dontrun{
#' library(DynCopula)
#' library(copula)
#' set.seed(0)
#'
#' N <- 1000
#' # Smooth covariate
#' t <- runif(N)
#' # Linear covariates
#' x1 <- sample(c("a", "b"), N, replace = TRUE)
#' x2 <- sample(seq(5), N, replace = TRUE)
#' # Design matrix
#' design <- data.frame(t = t, x1 = x1, x2 = x2)
#'
#' # Covariate-dependent correlation function
#' rho <- function(t, x1, x2) {
#'     0.7 * cos(4 * pi * t) + (x1 == "b") * 0.1 - 0.01 * x2
#' }
#' # Sample from copula
#' U <- t(sapply(seq(N), function(i) {
#'     rho_i <- rho(t[i], x1[i], x2[i])
#'     cop <- normalCopula(param = rho_i, dim = 2, dispstr = "un")
#'     return(rCopula(1L, cop))
#' }))
#'
#' # Fit with no covariates
#' gc1 <- fit_gamgc(U, design = data.frame(rep(1, N)), formula = ~1)
#' # Fit with a smooth covariate
#' gc2 <- fit_gamgc(U, design, formula = ~s(t))
#' # Fit with a smooth covariate and two linear covariates
#' gc3 <- fit_gamgc(U, design, formula = ~x1 + x2 + s(t))
#' # Fit with a smooth covariate, two linear covariates, and an interaction term
#' gc4 <- fit_gamgc(U, design, formula = ~x1 + x2*s(t))
#' }
#'
#' @export
fit_gamgc <- function(FX,
                      design,
                      formula,
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

    cores <- as.integer(cores)

    # Parse formula
    formulas <- .parse_formula(formula)
    lin_formula <- formulas$linear
    int_formula <- formulas$smooth_linear_int
    smooth_name <- formulas$smooth_name

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

    # Ensure correct Python interpreter is used
    conf_dir <- rappdirs::user_config_dir("DynCopula")
    conf_path <- file.path(conf_dir, "config.json")
    py_intr <- jsonlite::read_json(conf_path)$python_path
    Sys.setenv(
        RETICULATE_PYTHON = py_intr,
        RETICULATE_AUTOCONFIGURE = "FALSE"
    )

    if (is.null(int_formula)) {
        .glm_fit(
            FX = FX,
            design = design,
            lin_formula = lin_formula,
            lambda = lambda,
            nfold = nfold,
            cores = cores,
            control = control
        )
    } else {
        .gam_fit(
            FX = FX,
            design = design,
            lin_formula = lin_formula,
            int_formula = int_formula,
            smooth_name = smooth_name,
            lambda = lambda,
            K = K,
            nfold = nfold,
            cores = cores,
            control = control
        )
    }
}

###############################################################################

#' Internal function to fit a GLM Gaussian copula model
#'
#' @inheritParams fit_gamgc
#' @param lin_formula Formula for linear covariates.
.glm_fit <- function(FX,
                     design,
                     lin_formula,
                     lambda,
                     nfold,
                     cores,
                     control) {
    # Normal-transform pseudo-observations
    NX <- stats::qnorm(FX)

    # Location of Python files
    py_path <- system.file("python", package = "DynCopula")

    # Construct linear design matrix
    missing_vars <- setdiff(all.vars(lin_formula), colnames(design))
    if (length(missing_vars) > 0L) {
        stop("The following variables are missing from 'design': ",
             paste(missing_vars, collapse = ", "))
    }
    Z <- stats::model.matrix(lin_formula, design)

    N <- dim(NX)[1]
    d <- dim(NX)[2]
    L <- dim(Z)[2]

    # Whether parallelization is required. CV can be run with only one core,
    # but parallelization is over lambda so no parallelization is done if CV
    # is not going to be run.
    run_cv <- (length(lambda) > 1)
    run_parallel <- (cores > 1L) && run_cv

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
                module = "glm_fit",
                path = py_path,
                delay_load = FALSE
            )
        }
        # Everything but lam already exported to workers
        .fit_env$module$fit_gaussian_glm(
            NX = NX,
            Z = Z,
            lam = lam,
            control = control
        )
    }

    cv_fun <- function(lambda_ix, min_test_ix, max_test_ix) {
        # Load the module if it hasn't been loaded yet
        if (!exists("module", envir = .fit_env, inherits = FALSE)) {
            .fit_env$module <- reticulate::import_from_path(
                module = "glm_fit",
                path = py_path,
                delay_load = FALSE
            )
        }
        # NX, Z, control already exported to workers
        .fit_env$module$gaussian_glm_cv(
            NX = NX,
            Z = Z,
            lam = lambda[lambda_ix],
            control = control,
            min_test_ix = min_test_ix - 1L,  # convert to 0-indexing
            max_test_ix = max_test_ix - 1L
        )
    }

    if (run_cv) {
        # N-fold cross validation
        nfold <- as.integer(nfold)
        fold_ids <- cut(seq_len(N), breaks = nfold, labels = FALSE)
        folds_ids <- sample(fold_ids)

        # All combinations of penalty parameter and test fold ID
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
            varlist = c("NX", "lambda", "Z", "fold_ids", "combs",
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
        })

        # Sum likelihoods across folds
        ll_sum <- tapply(
            X = unlist(cv),
            INDEX = list(combs$lambda_ix),
            FUN = sum
        )
        cv_df <- as.data.frame(as.table(ll_sum))
        names(cv_df) <- c("lambda_ix", "ll")
        cv_df$lambda <- lambda[cv_df$lambda_ix]

        # Optimal value of lambda
        lambda_opt <- cv_df$lambda[which.max(cv_df$ll)]
    } else {
        # Only one lambda value passed, so no cross-validation needed
        lambda_opt <- lambda
        cv_df <- NULL
    }

    # Fit model using the optimal smoothing parameter
    beta_hat <- fit_fun(lambda_opt)

    # Save column names
    if (is.null(colnames(FX))) {
        cnames <- as.character(seq_len(dim(FX)[2]))
    } else {
        cnames <- colnames(FX)
    }

    res <- list(
        nobs = N,
        dim = d,
        beta = beta_hat,
        lambda = lambda_opt,
        cv = cv_df,
        Z_names = colnames(Z),
        lin_formula = lin_formula,
        colnames = cnames
    )
    class(res) <- "gamGaussianCopula"
    return(res)
}

###############################################################################

#' Internal function to fit a GAM Gaussian copula model
#'
#' @inheritParams fit_gamgc
#' @param lin_formula Formula for linear covariates.
#' @param int_formula Formula for interactions between linear covariates and
#' the smooth covariate.
#' @param smooth_name Name of the smooth covariate.
.gam_fit <- function(FX,
                     design,
                     lin_formula,
                     int_formula,
                     smooth_name,
                     lambda,
                     K,
                     nfold,
                     cores,
                     control) {
    # Normal-transform pseudo-observations
    NX <- stats::qnorm(FX)

    # Location of Python files
    py_path <- system.file("python", package = "DynCopula")

    # Sort by smooth covariate
    x <- design[[smooth_name]]
    if (is.unsorted(x)) {
        ord <- order(x)
        x <- x[ord]
        NX <- NX[ord, , drop = FALSE]
        design <- design[ord, , drop = FALSE]
    }

    # Construct linear design matrix
    all_vars <- c(all.vars(lin_formula), all.vars(int_formula))
    missing_vars <- setdiff(all_vars, colnames(design))
    if (length(missing_vars) > 0L) {
        stop("The following variables are missing from 'design': ",
             paste(missing_vars, collapse = ", "))
    }
    M <- stats::model.matrix(int_formula, design)
    Z <- stats::model.matrix(lin_formula, design)

    N <- dim(NX)[1]
    d <- dim(NX)[2]
    L1 <- dim(M)[2]
    L2 <- dim(Z)[2]

    # Whether parallelization is required. CV can be run with only one core,
    # but parallelization is over lambda so no parallelization is done if CV
    # is not going to be run.
    run_cv <- (length(lambda) > 1)
    run_parallel <- (cores > 1L) && run_cv

    # Scale smooth covariate to [0, 1]
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
                module = "gam_fit",
                path = py_path,
                delay_load = FALSE
            )
        }
        # Everything but lam already exported to workers
        .fit_env$module$fit_gaussian_gam(
            NX = NX,
            B = B,
            Z = Z,
            M = M,
            lam = lam,
            control = control
        )
    }

    cv_fun <- function(lambda_ix, min_test_ix, max_test_ix) {
        # Load the module if it hasn't been loaded yet
        if (!exists("module", envir = .fit_env, inherits = FALSE)) {
            .fit_env$module <- reticulate::import_from_path(
                module = "gam_fit",
                path = py_path,
                delay_load = FALSE
            )
        }
        # NX, B, Z, control already exported to workers
        .fit_env$module$gaussian_gam_cv(
            NX = NX,
            B = B,
            Z = Z,
            M = M,
            lam = unlist(lambda_grid[lambda_ix, ]),
            control = control,
            min_test_ix = min_test_ix - 1L,  # convert to 0-indexing
            max_test_ix = max_test_ix - 1L
        )
    }

    if (run_cv) {
        # N-fold cross validation with equal-sized contiguous blocks
        nfold <- as.integer(nfold)
        fold_ids <- cut(seq_len(N), breaks = nfold, labels = FALSE)

        # Grid search for smoothing parameters
        if (L2 > 1) {
            # If there are linear covariates besides the intercept, then
            # we include an extra lambda for the L2 penalty. But if the only
            # linear covariate is the intercept, then the penalty would amount
            # to a rescaling of the intercept, so we set the penalty to 0.
            lambda_grid <- expand.grid(rep(list(lambda), L1 + 1L))
        } else {
            lambda_grid <- expand.grid(rep(list(lambda), L1))
            lambda_grid <- cbind(rep(0, dim(lambda_grid)[1]), lambda_grid)
            colnames(lambda_grid)[1] <- "Var0"
        }

        # All combinations of smoothing parameter and test fold ID
        combs <- expand.grid(
            lambda_ix = seq_len(dim(lambda_grid)[1]),
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
            varlist = c("NX", "B", "lambda_grid", "Z", "M", "fold_ids",
                        "combs", "control", "py_path", ".fit_env"),
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
            names(cv_df) <- c("lambda_ix", "ll")
            cv_df$lambda <- lambda_grid[cv_df$lambda_ix, , drop = FALSE]

            # Optimal value of lambda
            lambda_opt <- unlist(cv_df$lambda[which.max(cv_df$ll), ])
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

    # Save column names
    if (is.null(colnames(FX))) {
        cnames <- as.character(seq_len(dim(FX)[2]))
    } else {
        cnames <- colnames(FX)
    }

    res <- list(
        nobs = N,
        dim = d,
        time = x,
        beta = beta_hat,
        lambda = lambda_opt,
        cv = cv_df,
        B = B,
        Z_names = colnames(Z),
        M_names = colnames(M),
        smooth_name = smooth_name,
        lin_formula = lin_formula,
        int_formula = int_formula,
        colnames = colnames(FX)
    )
    class(res) <- "gamGaussianCopula"
    return(res)
}

###############################################################################

#' Print method for class gamGaussianCopula
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
