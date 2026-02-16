###############################################################################

#' Spline estimation of dynamic Gaussian copula model
#'
#' @description Fit a dynamic Gaussian copula model using penalized splines.
#'
#' @param FX Matrix of pseudo-observations at covariate values.
#' @param x Vector of covariate values corresponding to \code{FX}. Must have no
#' duplicates.
#' @param lambda Vector of smoothing parameters to test. Default is
#' \code{10^(seq(-5, 5, length.out = 7))}.
#' @param df Vector of degrees of freedom for the spline basis matrix to test.
#' Default is \code{c(10, 50, 100)}.
#' @param model_select Method for model selection. Either \code{"aic"} (Akaike
#' information criterion) or \code{"cv"} (cross-validation).
#' @param nfold Number of folds for cross-validation. Default is \code{5}.
#' Only used if \code{model_select = "cv"}.
#' @param cores Number of cores to use. Default is \code{1}.
#' @param control A \code{list} of control parameters for optimization.
#'
#' @details Optimization is performed using L-BFGS. The \code{control} argument
#' is a list that supplies control parameters for optimization. The following
#' parameters can be supplied:
#' \itemize{
#'   \item \code{max_outer} Number of outer iterations. Default is \code{1}.
#'   \item \code{max_itr} Maximum number of inner iterations. Default is
#'   \code{100}.
#'   \item \code{history_size} History size. Default is \code{30}.
#'   \item \code{tolerance_grad} Termination tolerance for gradient. Default is
#'   \code{1e-7}.
#'   \item \code{tolerance_change} Termination tolerance for log-likelihood.
#'   Default is \code{1e-9}.
#'   \item \code{precision}: Floating precision. Either \code{"float32"} or
#'   \code{"float64"}. Default is \code{"float32"}.
#'   \item \code{boundary}: Boundary extension for spline knot boundaries.
#'   Default is \code{0.05}.
#' }
#'
#' The model hyperparameters are selected from the values specified by
#' \code{lambda} and \code{df}. The model selection method is specified by
#' \code{model_select}. Cross-validation will take approximately \code{nfold}
#' times longer than AIC, but may be more accurate.
#'
#' @return A list with the following components:
#' \itemize{
#'   \item \code{x}: Covariate values at which correlation coefficients were
#'   estimated.
#'   \item \code{NX}: Normal-transformed pseudo-observations.
#'   \item \code{eta}: Matrix of estimated calibration coefficients.
#'   \item \code{rho}: Matrix of estimated pairwise correlation coefficients.
#'   \item \code{lambda}: The optimal smoothing parameter.
#'   \item \code{df}: The optimal degrees of freedom.
#'   \item \code{model_select}: Model selection results.
#' }
#'
#' @export
fit_spline_gaussian <- function(FX,
                                x,
                                lambda = 10^(seq(-5, 5, length.out = 7)),
                                df = c(10, 50, 100),
                                model_select = c("aic", "cv"),
                                nfold = 5,
                                cores = 1,
                                control = list()) {
    assert_that(
        is.numeric(FX) && is.matrix(FX),
        is.vector(x, mode = "numeric"),
        !anyDuplicated(x),
        dim(FX)[1] == length(x),
        is.vector(lambda, mode = "numeric") && all(lambda > 0),
        is.vector(df, mode = "numeric") && all(df >= 1),
        is.character(model_select),
        is.numeric(cores) && cores >= 1,
        is.list(control)
    )

    # Model selection method
    model_select <- match.arg(model_select)
    if (model_select == "cv") {
        assert_that(is.numeric(nfold) && nfold >= 2)
        nfold <- as.integer(nfold)
    }

    # Whether parallelization is required
    cores <- as.integer(cores)
    run_parallel <- (cores > 1L)

    # Default control parameters
    defaults <- list(
        max_outer = 1L,
        max_itr = 100L,
        history_size = 30L,
        tolerance_grad = 1e-7,
        tolerance_change = 1e-9,
        precision = "float32",
        boundary = 0.05
    )
    control <- utils::modifyList(defaults, control)
    assert_that(all(names(control) %in% names(defaults)))

    # Verify control parameters
    assert_that(
        all(sapply(control[names(control) != "precision"], is.numeric)),
        control$max_outer >= 1,
        control$max_itr >= 1,
        control$history_size >= 1,
        control$tolerance_grad > 0,
        control$tolerance_change > 0,
        control$precision %in% c("float32", "float64"),
        control$boundary >= 0
    )
    control$max_outer <- as.integer(control$max_outer)
    control$max_itr <- as.integer(control$max_itr)
    control$history_size <- as.integer(control$history_size)

    # Dimension and number of parameters
    d <- dim(FX)[2]
    npar <- choose(d, 2)

    # Sort covariates if needed
    if (is.unsorted(x)) {
        ord <- order(x)
        x <- x[ord]
        FX <- FX[ord, ]
    }

    # Scale covariates to [0, 1]
    min_x <- x[1]
    dx <- x[length(x)] - min_x
    x <- (x - min_x) / dx

    # Normal-transform pseudo-observations
    NX <- stats::qnorm(FX)

    # Spline basis matrix with knot boundary extension
    knot_eps <- control$boundary
    get_basis <- function(x, df) {
        splines::ns(
            x = x,
            df = df,
            intercept = TRUE,
            Boundary.knots = c(-knot_eps, 1 + knot_eps)
        )
    }

    ## Python configuration ##
    # Ensure correct Python interpreter is used
    conf_dir <- rappdirs::user_config_dir("DynCopula")
    conf_path <- file.path(conf_dir, "config.json")
    py_intr <- jsonlite::read_json(conf_path)$python_path
    Sys.setenv(
        RETICULATE_PYTHON = py_intr,
        RETICULATE_AUTOCONFIGURE = "FALSE"
    )

    if (!run_parallel) {
        reticulate::py_config()
    } else {
        # Set up cluster
        cl <- parallel::makeCluster(cores)

        parallel::clusterEvalQ(cl, {
            library(reticulate)

            # Force Python initialization on workers
            reticulate::py_config()

            # Disable multithreading to prevent oversubscription
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

    # Since Python functions are not serializable, we can't load the
    # optimization function in the global environment. Instead, the function
    # needs to be loaded on each worker. This environment stores the
    # function once it has been loaded to avoid having to do it multiple
    # times on the same worker.
    .fit_env <- new.env(parent = emptyenv())
    py_path <- system.file("python", package = "DynCopula")
    fit_fun <- function(par0, x, NX, B, lam, control, compute_ll, compute_edf,
                        py_path) {
        # Load the module if it hasn't been loaded yet
        if (!exists("fit_fun", envir = .fit_env, inherits = FALSE)) {
            .fit_env$fit_fun <- reticulate::import_from_path(
                module = "spline_gaussian",
                path = py_path,
                delay_load = FALSE
            )$fit_gaussian_spline
        }
        .fit_env$fit_fun(
            par0 = par0,
            x = x,
            NX = NX,
            B = B,
            lam = lam,
            control = control,
            compute_ll = compute_ll,
            compute_edf = compute_edf
        )
    }

    if (model_select == "aic") {
        # All combinations of lambda and df
        combs <- expand.grid(
            lambda_ix = seq_along(lambda),
            df_ix = seq_along(df)
        )
    } else {
        # K-fold cross validation with equal-sized contiguous blocks
        fold_ids <- cut(seq_along(x), breaks = nfold, labels = FALSE)
        # All combinations of lambda df, and folds
        combs <- expand.grid(
            lambda_ix = seq_along(lambda),
            df_ix = seq_along(df),
            fold = seq_len(nfold)
        )
    }
    ncomb <- dim(combs)[1]

    # Avoid setting up futures if no parallelization is requested
    lfun <- ifelse(run_parallel, future.apply::future_lapply, lapply)
    largs <- list(X = seq_len(ncomb))

    message("Starting")

    progressr::with_progress({
        pbar <- progressr::progressor(along = seq_len(ncomb + 1L))
        if (model_select == "aic") {
            ## Model selection via AIC ##

            if (run_parallel) {
                largs <- c(largs, list(
                    future.seed = TRUE,
                    future.globals = list(
                        x = x, NX = NX, lambda = lambda, df = df, combs = combs,
                        pbar = pbar, get_basis = get_basis, npar = npar,
                        py_path = py_path, control = control,
                        .fit_env = .fit_env
                    ),
                    future.packages = c("splines")
                ))
            }

            aic_fun <- function(i) {
                # Spline basis matrix
                B <- get_basis(x, df[combs$df_ix[i]])

                # Initial coefficient estimates
                par0 <- matrix(0, nrow = dim(B)[2], ncol = npar)

                # Fit model
                model_fit <- fit_fun(
                    par0 = par0,
                    x = x,
                    NX = NX,
                    B = B,
                    lam = lambda[combs$lambda_ix[i]],
                    control = control,
                    compute_ll = TRUE,
                    compute_edf = TRUE,
                    py_path = py_path
                )
                pbar()
                return(model_fit[c("ll", "edf")])
            }
            largs <- c(largs, list(FUN = aic_fun))
            res <- do.call(what = lfun, args = largs)
            edf <- sapply(res, '[[', "edf")
            ll <- sapply(res, '[[', "ll")

            # Convert indices for lambda and df to values
            model_df <- data.frame(
                lambda = lambda[combs$lambda_ix],
                df = df[combs$df_ix],
                edf = edf,
                ll = ll,
                aic = -2 * ll + 2 * edf
            )

            # Index for optimal hyperparameters
            ix_opt <- which.min(model_df$aic)
        } else {
            ## Model selection via cross-validation ##
            if (run_parallel) {
                largs <- c(largs, list(
                    future.seed = TRUE,
                    future.globals = list(
                        x = x, NX = NX, lambda = lambda, df = df, combs = combs,
                        pbar = pbar, get_basis = get_basis, py_path = py_path,
                        fold_ids = fold_ids, npar = npar, control = control,
                        .fit_env = .fit_env
                    ),
                    future.packages = c("splines", "mvtnorm", "copula")
                ))
            }

            cv_fun <- function(i) {
                # Train and test folds
                test_ix <- which(fold_ids == combs$fold[i])
                train_ix <- which(fold_ids != combs$fold[i])

                # Spline basis matrices
                B <- get_basis(x, df[combs$df_ix[i]])
                B_train <- B[train_ix, , drop = FALSE]
                B_test <- B[test_ix, , drop = FALSE]

                # Initial coefficient estimates
                par0 <- matrix(0, nrow = dim(B)[2], ncol = npar)

                # Fit using training data
                beta_est <- fit_fun(
                    par0 = par0,
                    x = x[train_ix],
                    NX = NX[train_ix, , drop = FALSE],
                    B = B_train,
                    lam = lambda[combs$lambda_ix[i]],
                    control = control,
                    compute_ll = FALSE,
                    compute_edf = FALSE,
                    py_path = py_path
                )$beta

                # Predictions for test data
                H_test <- B_test %*% beta_est

                # Cross-validated likelihood criterion
                margin_ll <- apply(
                    X = stats::dnorm(NX[test_ix, , drop = FALSE], log = TRUE),
                    MARGIN = 1,
                    FUN = sum
                )
                copula_ll <- sapply(seq_along(test_ix), function(j) {
                    eta <- H_test[j, ]
                    if (any(is.na(eta) | is.nan(eta))) {
                        return(-1e12)
                    }
                    mvtnorm::dmvnorm(
                        x = NX[test_ix[j], ],
                        sigma = vec2cor(eta),
                        log = TRUE,
                        checkSymmetry = FALSE
                    )
                })
                ll <- sum(copula_ll - margin_ll)

                pbar()
                return(ll)
            }
            largs <- c(largs, list(FUN = cv_fun))
            cv <- do.call(what = lfun, args = largs)

            # Sum likelihoods across folds
            ll_sum <- tapply(
                X = unlist(cv),
                INDEX = list(combs$lambda_ix, combs$df_ix),
                FUN = sum
            )
            cv_df <- as.data.frame(as.table(ll_sum))
            names(cv_df) <- c("lambda_ix", "df_ix", "ll")

            # Convert factor IDs to indices
            cv_df$lambda_ix <- as.integer(as.character(cv_df$lambda_ix))
            cv_df$df_ix <- as.integer(as.character(cv_df$df_ix))

            # Convert indices for lambda and df to values
            model_df <- data.frame(
                lambda = lambda[cv_df$lambda_ix],
                df = df[cv_df$df_ix],
                ll = cv_df$ll
            )

            # Index for optimal hyperparameters
            ix_opt <- which.max(model_df$ll)
        }

        # Optimal hyperparameters
        lambda_opt <- model_df$lambda[ix_opt]
        df_opt <- model_df$df[ix_opt]

        # Estimation using the selected hyperparameters
        B <- get_basis(x, df_opt)
        par0 <- matrix(data = 0, nrow = df_opt, ncol = npar)
        beta_opt <- fit_fun(
            par0 = par0,
            x = x,
            NX = NX,
            B = B,
            lam = lambda_opt,
            control = control,
            compute_ll = FALSE,
            compute_edf = FALSE,
            py_path = py_path
        )$beta
        pbar()
    })

    # Matrix of estimated calibration coefficients
    Hhat <- B %*% beta_opt
    # Matrix of estimated correlation coefficients
    Rhat <- t(apply(Hhat, 1, function(v) {
        copula::P2p(vec2cor(v))
    }))

    # Ensure consistent shape of Rhat across all dimensions
    if (d == 2L) {
        Rhat <- t(Rhat)
    }

    # Add numbered eta/rho labels
    colnames(Hhat) <- paste0("eta", seq_len(npar))
    ix_lab <- apply(utils::combn(seq_len(d), 2), 2, function(x) {
        paste(x, collapse = '_')
    })
    colnames(Rhat) <- paste0("rho", ix_lab)

    # Rescale covariates back to their original scale
    x <- x * dx + min_x

    return(list(
        x = x,
        NX = NX,
        eta = Hhat,
        rho = Rhat,
        lambda = lambda_opt,
        df = df_opt,
        model_select = model_df
    ))
}

###############################################################################
