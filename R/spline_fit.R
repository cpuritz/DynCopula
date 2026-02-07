###############################################################################

#' Spline estimation of dynamic Gaussian copula model
#'
#' @description Fit a dynamic Gaussian copula model using penalized splines.
#'
#' @param FX Matrix of pseudo-observations at covariate values.
#' @param x Vector of covariate values corresponding to \code{FX}. Must have no
#' duplicates.
#' @param lambda Vector of smoothing parameters. Default is
#' \code{10^(seq(-5, 5, length.out = 7))}.
#' @param df Vector of degrees of freedom. Default is \code{c(10, 50, 100)}.
#' @param nfold Number of folds for cross-validation. Default is \code{5}.
#' @param cores Number of cores to use. Default is \code{1}.
#' @param control A \code{list} of control parameters for optimization.
#'
#' @details Cross-validation is used to select the values of \code{lambda} and
#' \code{df}.
#'
#' Optimization is performed using L-BFGS. The \code{control} argument
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
#' }
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
#'   \item \code{cv}: Cross-validation results.
#' }
#'
#' @export
fit_spline_gaussian <- function(FX,
                                x,
                                lambda = 10^(seq(-5, 5, length.out = 7)),
                                df = c(10, 50, 100),
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
        is.numeric(nfold) && nfold >= 2,
        is.numeric(cores) && cores >= 1,
        is.list(control)
    )
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

    # Gaussian copula log likelihood
    loglik <- function(NU, eta) {
        if (any(is.na(eta) | is.nan(eta))) {
            return(-1e12)
        }
        copula_ll <- mvtnorm::dmvnorm(
            x = NU,
            sigma = vec2cor(eta),
            log = TRUE,
            checkSymmetry = FALSE
        )
        margin_ll <- sum(stats::dnorm(NU, log = TRUE))
        return(copula_ll - margin_ll)
    }

    # K-fold cross validation with contiguous blocks
    fold_ids <- cut(seq_along(x), breaks = nfold, labels = FALSE)

    # Spline basis matrix
    get_basis <- function(x, df) {
        # Knot boundary extension
        knot_eps <- 0.05
        splines::ns(
            x = x,
            df = df,
            intercept = TRUE,
            Boundary.knots = c(-knot_eps, 1 + knot_eps)
        )
    }

    # All combinations of lambda and df
    combs <- expand.grid(
        lambda_ix = seq_along(lambda),
        df_ix = seq_along(df),
        fold = seq_len(nfold)
    )
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
                    control = control
                )

                # Predictions for test data
                H_test <- B_test %*% beta_est

                # Cross-validated likelihood criterion
                ll <- sum(sapply(seq_along(test_ix), function(j) {
                    loglik(NX[test_ix[j], ], H_test[j, ])
                }))
                pbar()
                return(ll)
            },
            future.seed = TRUE,
            future.globals = TRUE
        )
    })

    # Sum likelihoods across folds
    ll_sum <- tapply(
        X = unlist(cv),
        INDEX = list(combs$lambda_ix, combs$df_ix),
        FUN = sum
    )
    cv_df <- as.data.frame(as.table(ll_sum))
    names(cv_df) <- c("lambda_ix", "df_ix", "ll")
    cv_df$lambda_ix <- as.integer(as.character(cv_df$lambda_ix))
    cv_df$df_ix <- as.integer(as.character(cv_df$df_ix))
    cv_df <- data.frame(
        lambda = lambda[cv_df$lambda_ix],
        df = df[cv_df$df_ix],
        ll = cv_df$ll
    )

    # Optimal hyperparameters
    lambda_opt <- cv_df[which.max(cv_df$ll), "lambda"]
    df_opt <- cv_df[which.max(cv_df$ll), "df"]

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
        cv = cv_df
    ))
}

###############################################################################
