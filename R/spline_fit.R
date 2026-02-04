###############################################################################

#' Spline estimation of dynamic Gaussian copula model
#'
#' @description Fit a dynamic Gaussian copula model using penalized splines.
#'
#' @param FX Matrix of pseudo-observations at covariate values.
#' @param x Vector of covariate values corresponding to \code{FX}. Must be
#' sorted and have no duplicates.
#' @param lambda Vector of smoothing parameters. Default is
#' \code{10^(seq(-5, 5, length.out = 7))}.
#' @param lambda_blocks Number of smoothing blocks. Default is \code{1}.
#' @param df Vector of degrees of freedom. Default is \code{c(10, 50, 100)}.
#' @param nfold Number of folds for cross-validation. Default is \code{10}.
#' @param cores Number of cores to use. Default is \code{1}.
#' @param control A \code{list} of control parameters for optimization.
#'
#' @details Cross-validation is used to select the values of \code{lambda} and
#' \code{df}. An initial estimate is first computed by cross-validation over all
#' values of \code{lambda} and \code{df}. If \code{lambda_blocks > 1}, a second
#' sweep is performed to choose block-specific smoothing parameters, with the
#' parameter set of size \code{|lambda|^(lambda_blocks)}.
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
#'   \item \code{NX}: Normal-transformed pseudo-observations.
#'   \item \code{eta}: Matrix of estimated calibration coefficients.
#'   \item \code{rho}: Matrix of estimated pairwise correlation coefficients.
#'   \item \code{lambda}: The optimal smoothing parameter.
#'   \item \code{df}: The optimal degrees of freedom.
#'   \item \code{cv}: Cross-validation results
#' }
#'
#' @export
fit_spline_gaussian <- function(FX,
                                x,
                                lambda = 10^(seq(-5, 5, length.out = 7)),
                                lambda_blocks = 1,
                                df = c(10, 50, 100),
                                nfold = 10,
                                cores = 1,
                                control = list()) {
    assert_that(
        is.numeric(FX) && is.matrix(FX),
        is.vector(x, mode = "numeric"),
        !anyDuplicated(x),
        !is.unsorted(x),
        dim(FX)[1] == length(x),
        is.vector(lambda, mode = "numeric"),
        is.numeric(lambda_blocks) && lambda_blocks >= 1,
        is.vector(df, mode = "numeric") && all(df >= 1),
        is.numeric(nfold) && nfold >= 2,
        is.numeric(cores) && cores >= 1,
        is.list(control)
    )
    lambda_blocks <- as.integer(lambda_blocks)
    nfold <- as.integer(nfold)

    # Set up futures plan
    cl <- parallel::makeCluster(as.integer(cores))
    future::plan(future::cluster, workers = cl)
    on.exit({ future::plan(future::sequential); parallel::stopCluster(cl) },
            add = TRUE)

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

    # Normal-transform pseudo-observations
    NX <- stats::qnorm(FX)

    if (lambda_blocks != 1L) {
        message("Fitting pilot estimate")
    }

    # First fit pilot estimate
    res_pilot <- .spline_fit(
        NX = NX,
        x = x,
        lambda = lambda,
        lambda_blocks = 1L,
        df = df,
        nfold = nfold,
        control = control,
        cl = cl
    )

    if (lambda_blocks == 1L) {
        res_pilot$beta <- NULL
        return(res_pilot)
    }

    message("Fitting multiple smoothing blocks")

    # Quantify roughness of pilot curve estimates
    D2 <- diff(diff(diag(dim(res_pilot$beta)[1])))
    S <- t(D2) %*% D2
    roughness <- diag(t(res_pilot$beta) %*% S %*% res_pilot$beta)

    # Identify clusters of curves by roughness
    sink <- utils::capture.output(
        groups <- mclust::Mclust(roughness, G = lambda_blocks)$classification
    )

    # Fit using multiple smoothing blocks using pilot df
    res_spline <- .spline_fit(
        NX = NX,
        x = x,
        lambda = lambda,
        lambda_blocks = groups,
        df = res_pilot$df,
        nfold = nfold,
        control = control,
        cl = cl
    )
    res_spline$beta <- NULL
    return(res_spline)
}

###############################################################################

#' Spline estimation of dynamic Gaussian copula model
#'
#' @description Fit a dynamic Gaussian copula model using penalized splines.
#' Internal function.
#'
#' @param NX Matrix of normal-transformed pseudo-observations at covariate
#' values.
#' @param x Vector of covariate values corresponding to \code{FX}. Must be
#' sorted and have no duplicates.
#' @param lambda Vector of smoothing parameters.
#' @param lambda_blocks A vector of integers specifying block assignments. Each
#' block has its own smoothing parameter.
#' @param df Vector of degrees of freedom.
#' @param nfold Number of folds for cross-validation.
#' @param control A \code{list} of control parameters for optimization.
#' @param cl A cluster for parallel computations.
#'
#' @return A list with the following components:
#' \itemize{
#'   \item \code{x}: The input argument \code{x}.
#'   \item \code{NX}: Normal-transformed pseudo-observations.
#'   \item \code{eta}: Matrix of estimated calibration coefficients.
#'   \item \code{rho}: Matrix of estimated pairwise correlation coefficients.
#'   \item \code{lambda}: The optimal smoothing parameters.
#'   \item \code{df}: The optimal degrees of freedom.
#'   \item \code{beta}: Estimated basis coefficients.
#'   \item \code{cv}: Cross-validation results.
#' }
.spline_fit <- function(NX,
                        x,
                        lambda,
                        lambda_blocks,
                        df,
                        nfold,
                        control,
                        cl) {
    # Dimension
    d <- dim(NX)[2]
    npar <- choose(d, 2)

    # Scale covariates to [0, 1]
    min_x <- x[1]
    dx <- x[length(x)] - min_x
    x <- (x - min_x) / dx

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

    # K-fold cross validation
    fold_ids <- cut(seq_along(x), breaks = nfold, labels = FALSE)

    # Knot boundary extension
    knot_eps <- 0.05
    boundary_knot <- c(-knot_eps, 1 + knot_eps)

    # Number of smoothing blocks
    nblock <- length(unique(lambda_blocks))

    # All combinations of lambda and df
    combs <- expand.grid(c(rep(list(lambda), nblock), list(df = df)))
    colnames(combs)[seq_len(nblock)] <- paste0("lambda", seq_len(nblock))
    ncomb <- dim(combs)[1]

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
                for (fold in seq_len(nfold)) {
                    test_ix <- which(fold_ids == fold)
                    train_ix <- which(fold_ids != fold)

                    # Spline basis matrices
                    B <- splines::ns(
                        x = x,
                        df = df_i,
                        intercept = TRUE,
                        Boundary.knots = boundary_knot
                    )
                    B_train <- B[train_ix, , drop = FALSE]
                    B_test <- B[test_ix, , drop = FALSE]

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
    cv_df <- cbind(combs, data.frame(ll = cv))
    comb_opt <- unlist(combs[which.max(cv), ])

    # Optimal lambda values for each component of eta
    lambda_opt <- comb_opt[seq_len(nblock)][lambda_blocks]

    # Optimal basis size
    df_opt <- comb_opt[nblock + 1L]

    # Estimation using the CV optimal lambda and df
    B <- splines::ns(
        x = x,
        df = df_opt,
        intercept = TRUE,
        Boundary.knots = boundary_knot
    )
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
        beta = beta_opt,
        cv = cv_df
    ))
}

###############################################################################
