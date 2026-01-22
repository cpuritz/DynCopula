###############################################################################

#' Spline estimation of dynamic Gaussian copula model
#'
#' @description Fit a dynamic Gaussian copula model using B-splines.
#'
#' @param FX Matrix of pseudo-observations at covariate values.
#' @param x Vector of covariate values corresponding to \code{FX}. Must be
#' sorted and have no duplicates.
#' @param lambda Vector of smoothing parameters.
#' @param lambda_blocks Either an integer specifying the number of smoothing
#' blocks, or a vector of integers specifying block assignments. Each block
#' has its own smoothing parameter. Default is \code{1}.
#' @param df Vector of degrees of freedom.
#' @param nfold Number of folds for cross-validation.
#' @param cores Number of cores to use. Default is \code{1}.
#' @param control A \code{list} of control parameters for optimization.
#'
#' @details Cross-validation is used to select the values of \code{lambda} and
#' \code{df}. The number of sets of parameters evaluated is of size
#' \code{len(lambda_blocks)^(lambda_blocks) * len(df)}.
#'
#' Optimization is performed using L-BFGS. The \code{control} argument
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
#'   \item \code{h}: The bandwidth used at each value in \code{x0}.
#'   \item \code{eta}: Matrix of estimated calibration coefficients.
#'   \item \code{rho}: Matrix of estimated pairwise correlation coefficients.
#'   \item \code{cv}: Data frame of cross-validation results.
#' }
#'
#' @export
fit_spline_gaussian <- function(FX,
                                x,
                                lambda,
                                lambda_blocks = 1L,
                                df,
                                nfold = 5L,
                                cores = 1L,
                                control = list()) {
    assert_that(
        is.vector(x, mode = "numeric"),
        !anyDuplicated(x),
        is.matrix(FX),
        is.numeric(FX),
        dim(FX)[1] == length(x),
        is.numeric(lambda),
        is.numeric(lambda_blocks),
        is.numeric(df) && all(df >= 1),
        is.numeric(nfold) && nfold >= 2,
        is.numeric(cores) && cores >= 1,
        is.list(control)
    )
    cores <- as.integer(cores)

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

    # Dimension
    d <- dim(FX)[2]
    npar <- choose(d, 2)

    # Blocks of smoothing parameters
    if (length(lambda_blocks) == 1) {
        lambda_blocks <- as.integer(lambda_blocks)
        assert_that(lambda_blocks >= 1L)
        # Break into equal sized blocks
        nrep <- ceiling(npar / lambda_blocks)
        lambda_blocks <- sort(rep(seq_len(lambda_blocks), nrep)[1:npar])
    }
    assert_that(length(lambda_blocks) == npar)
    # Convert to sequential integers starting at 1
    lambda_blocks <- as.integer(factor(lambda_blocks))
    nblock <- length(unique(lambda_blocks))

    # Scale covariates to [0, 1]
    min_x <- x[1]
    dx <- x[length(x)] - min_x
    x <- (x - min_x) / dx

    # Normal-transform pseudo-observations
    NX <- stats::qnorm(FX)

    # Gaussian copula log likelihood
    loglik <- function(NX, eta) {
        copula_ll <- mvtnorm::dmvnorm(
            x = NX,
            sigma = vec2cor(eta),
            log = TRUE,
            checkSymmetry = FALSE
        )
        margin_ll <- sum(stats::dnorm(NX, log = TRUE))
        return(copula_ll - margin_ll)
    }

    get_basis <- function(x, K) {
        splines::ns(x, df = K, intercept = TRUE)
    }

    # K-fold cross validation
    folds <- seq_len(nfold)
    fold_ids <- rep(folds, ceiling(length(x) / nfold))[1:length(x)]

    # All combinations of lambda and df
    combs <- expand.grid(c(rep(list(lambda), nblock), list(df = df)))
    colnames(combs)[seq_len(nblock)] <- paste0("lambda", seq_len(nblock))
    ncomb <- dim(combs)[1]

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
        fit_fun <- function(par0, x, NX, B, lam, control) {
            # Load the module if it hasn't been loaded yet
            if (!exists("fit_fun", envir = .fit_env, inherits = FALSE)) {
                .fit_env$fit_fun <- reticulate::import_from_path(
                    module = "spline_gaussian",
                    path = system.file("python", package = "DynCopula"),
                    delay_load = FALSE
                )$fit_gaussian_spline
            }
            .fit_env$fit_fun(
                par0 = par0,
                x = x,
                NX = NX,
                B = B,
                lam = lam,
                control = control
            )
        }

        progressr::with_progress({
            pbar <- progressr::progressor(along = seq_len(ncomb))
            cv <- future.apply::future_lapply(
                X = seq_len(ncomb),
                FUN = function(i) {
                    pars_i <- unlist(combs[i, ])
                    # Lambda value for each block
                    group_lambdas <- pars_i[seq_len(nblock)]
                    # Lambda value for each component of eta
                    lambda_i <- group_lambdas[lambda_blocks]
                    # Basis size
                    df_i <- pars_i[length(pars_i)]

                    cv_i <- 0
                    for (fold in folds) {
                        test_ix <- which(fold_ids == fold)
                        train_ix <- setdiff(seq_along(x), test_ix)

                        # Spline basis matrices
                        B_train <- get_basis(x[train_ix], df_i)
                        B_test <- get_basis(x[test_ix], df_i)

                        # Initial estimate
                        par0 <- matrix(0, nrow = df_i, ncol = npar)

                        # Fit using training data
                        beta_est <- fit_fun(
                            par0 = par0,
                            x = x[train_ix],
                            NX = NX[train_ix, , drop = FALSE],
                            B = B_train,
                            lam = lambda_i,
                            control = control
                        )

                        # Predictions for test data
                        H_test <- B_test %*% beta_est

                        # Cross-validated likelihood criterion
                        ll <- sum(sapply(seq_along(test_ix), function(j) {
                            loglik(NX[test_ix[j], ], H_test[j, ])
                        }))
                        cv_i <- cv_i + ll
                    }
                    pbar()
                    return(cv_i)
                },
                future.seed = TRUE,
                future.globals = TRUE
            )
        })
        cv <- unlist(cv)
    } else {
        # Load Python module
        fit_fun <- reticulate::import_from_path(
            module = "spline_gaussian",
            path = system.file("python", package = "DynCopula"),
            delay_load = FALSE
        )$fit_gaussian_spline

        progressr::with_progress({
            pbar <- progressr::progressor(along = seq_len(ncomb))

            cv <- sapply(seq_len(ncomb), function(i) {
                pars_i <- unlist(combs[i, ])
                # Lambda value for each block
                group_lambdas <- pars_i[seq_len(nblock)]
                # Lambda value for each component of eta
                lambda_i <- group_lambdas[lambda_blocks]
                # Basis size
                df_i <- pars_i[length(pars_i)]

                cv_i <- 0
                for (fold in folds) {
                    test_ix <- which(fold_ids == fold)
                    train_ix <- setdiff(seq_along(x), test_ix)

                    # Spline basis matrices
                    B_train <- get_basis(x[train_ix], df_i)
                    B_test <- get_basis(x[test_ix], df_i)

                    # Initial estimate
                    par0 <- matrix(0, nrow = df_i, ncol = npar)

                    # Fit using training data
                    beta_est <- fit_fun(
                        par0 = par0,
                        x = x[train_ix],
                        NX = NX[train_ix, , drop = FALSE],
                        B = B_train,
                        lam = lambda_i,
                        control = control
                    )

                    # Predictions for test data
                    H_test <- B_test %*% beta_est

                    # Cross-validated likelihood criterion
                    ll <- sum(sapply(seq_along(test_ix), function(j) {
                        loglik(NX[test_ix[j], ], H_test[j, ])
                    }))
                    cv_i <- cv_i + ll
                }
                pbar()
                return(cv_i)
            })
        })
    }

    ix_opt <- which.max(cv)
    # Optimal lambda values for each component of eta
    lambda_opt <- unlist(combs[ix_opt, seq_len(nblock)])[lambda_blocks]
    # Optimal basis size
    df_opt <- combs[ix_opt, dim(combs)[2]]
    # Estimation using the CV optimal lambda and df
    B <- get_basis(x, df_opt)
    par0 <- matrix(0, nrow = df_opt, ncol = npar)
    beta_opt <- fit_fun(
        par0 = par0,
        x = x,
        NX = NX,
        B = B,
        lam = lambda_opt,
        control = control
    )

    # Matrix of estimated calibration coefficients
    Hhat <- B %*% beta_opt
    # Matrix of estimated correlation coefficients
    Rhat <- t(apply(Hhat, 1, function(v) {
        copula::P2p(vec2cor(v))
    }))

    # Ensure consistent shape of Rhat across all dimensions
    if (d == 2) {
        Rhat <- t(Rhat)
    }

    # Estimate roughness
    x_unif <- seq(0, 1, length.out = 1e3)
    roughness <- apply(Hhat, 2, function(y) {
        h_unif <- stats::predict(stats::smooth.spline(x = x, y = y), x_unif)$y
        sum(diff(diff(h_unif))^2)
    })
    roughness <- roughness / min(roughness)

    # Add numbered eta/rho labels
    colnames(Hhat) <- paste0("eta", seq_len(npar))
    ix_lab <- apply(utils::combn(seq_len(d), 2), 2, function(x) {
        paste(x, collapse = '_')
    })
    colnames(Rhat) <- paste0("rho", ix_lab)

    # Rescale covariates back to their original scale
    x <- x * dx + min_x

    # Save cross validation results
    combs$ll <- cv

    return(list(
        x = x,
        NX = NX,
        eta = Hhat,
        rho = Rhat,
        cv = combs,
        roughness = roughness
    ))
}

###############################################################################
