###############################################################################

#' Local likelihood estimation of conditional Gaussian copula model
#'
#' @description Fit a conditional Gaussian copula model using local likelihood.
#'
#' @param FX Matrix of pseudo-observations at covariate values.
#' @param x Vector of covariate values corresponding to \code{FX}. Must be
#' sorted and have no duplicates.
#' @param bandwidths Vector of global kernel bandwidths to test.
#' @param variable Whether to use a variable bandwidth. Default is \code{TRUE}.
#' @param alpha Powers to test for the variable bandwidth function. Ignored if
#' \code{variable = FALSE}. Must be between \code{0} and \code{1}. Default is
#' \code{seq(0, 1, 0.1)}.
#' @param beta Scale factors to test for the variable bandwidth function.
#' Ignored if \code{variable = FALSE}. Must be greater than \code{0}. Default is
#' \code{seq(0.75, 1.25, 0.05)}.
#' @param ncv Number of covariate values to use for LOOCV. Default is
#' \code{length(x)}.
#' @param degree Degree of local polynomial approximation. Default is \code{0}.
#' @param control A \code{list} of control parameters for optimization.
#' @param cores Number of cores to use. Default is \code{1}.
#'
#' @details Optimization is performed using L-BFGS. The \code{control} argument
#' is a list that supplies control parameters for optimization. The following
#' parameters can be supplied:
#' \itemize{
#'   \item \code{max_outer}: Maximum number of outer iterations. Default is
#'   \code{1}.
#'   \item \code{max_itr}: Maximum number of inner iterations. Default is
#'   \code{100}.
#'   \item \code{history_size}: History size. Default is \code{30}.
#'   \item \code{tolerance_grad}: Termination tolerance for gradient. Default is
#'   \code{1e-7}.
#'   \item \code{tolerance_change}: Termination tolerance for log-likelihood.
#'   Default is \code{1e-9}.
#' }
#' Any control parameters not specified are assigned their default values.
#'
#' If only one global kernel bandwidth is specified, the model is fit using this
#' value and returned. Otherwise, leave-one-out cross-validation (LOOCV) is used
#' to select the optimal bandwidth out of the bandwidths specified. The optimal
#' bandwidth is the one that maximizes the cross-validated likelihood criterion.
#' The model returned uses this optimal bandwidth.
#'
#' If a variable bandwidth should be used (\code{variable = TRUE}), a global
#' pilot bandwidth is first selected as discussed above. LOOCV is then used
#' to select parameters for the variable bandwidth function from the values in
#' \code{alpha} and \code{beta}. The variable bandwidth function is defined as
#' \deqn{
#' h(x;\alpha,\beta)=\beta h_{0}\big(\hat{f}_{x}(x)/G\big)^{-\alpha}
#' }
#' where \eqn{h_{0}} is the global pilot bandwidth, \eqn{\hat{f}_{x}} is a
#' kernel density estimator for the covariate values, and \eqn{G} is the
#' geometric mean of \eqn{\hat{f}_{x}(x)}.
#'
#' The argument \code{ncv} specifies the number of covariate values to use for
#' LOOCV. Full LOOCV corresponds to \code{ncv = length(x)}. If
#' \code{ncv < length(x)}, then only a subset of the covariates are used to
#' reduce the run time. The indices of the chosen covariates are equally spaced
#' along \code{seq_along(x)}. This is only an estimate to full LOOCV and may be
#' quite inaccurate for \code{ncv << length(x)}.
#'
#' @return A list with the following components:
#' \itemize{
#'   \item \code{x}: The input argument \code{x}.
#'   \item \code{x0}: The input argument \code{x0}.
#'   \item \code{bandwidths}: The bandwidths used to fit the final model.
#'   \item \code{NX}: Normal-transformed pseudo-observations.
#'   \item \code{eta}: Matrix of estimated calibration coefficients.
#'   \item \code{rho}: Matrix of estimated pairwise correlation coefficients.
#'   \item \code{global_cv}: Cross-validation results for the global bandwidth
#'   (if performed).
#'   \item \code{var_cv}: Cross-validation results for the variable bandwidth
#'   parameters (if performed).
#' }
#'
#' @export
fit_local_gaussian <- function(FX,
                               x,
                               bandwidths,
                               variable = TRUE,
                               alpha = seq(0, 1, 0.1),
                               beta = seq(0.75, 1.25, 0.05),
                               ncv = length(x),
                               degree = 0,
                               control = list(),
                               cores = 1) {

    # Set up futures plan
    assert_that(is.numeric(cores) && cores >= 1)
    cl <- parallel::makeCluster(as.integer(cores))
    future::plan(future::cluster, workers = cl)
    on.exit({ future::plan(future::sequential); parallel::stopCluster(cl) },
            add = TRUE)

    if (length(bandwidths) == 1L && !variable) {
        # No bandwidth selection, just need to fit model. Error checking
        # is handled by .local_fit.
        return(.local_fit(
            FX = FX,
            x = x,
            x0 = x,
            h = bandwidths,
            degree = degree,
            control = control,
            cl = cl
        ))
    }

    assert_that(
        is.vector(x, mode = "numeric"),
        is.matrix(FX),
        is.numeric(FX),
        dim(FX)[1] == length(x),
        is.numeric(bandwidths) && all(bandwidths > 0) && all(bandwidths < 1),
        !anyDuplicated(x) && !is.unsorted(x),
        is.numeric(alpha) && all(alpha >= 0) && all(alpha <= 1),
        is.numeric(beta) && all(beta > 0),
        is.numeric(ncv) && ncv > 1,
        is.numeric(degree) && degree >= 0
    )

    # Require that global model (alpha = 0, beta = 1) is included in the set of
    # variable bandwidth parameters
    alpha <- unique(alpha)
    beta <- unique(beta)
    if (!any(alpha < 1e-16)) {
        alpha <- c(0, alpha)
    }
    if (!any(abs(beta - 1.0) < 1e-16)) {
        beta <- c(1.0, beta)
    }

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

    # Use ncv equally spaced covariate values
    ncv <- as.integer(ncv)
    if (ncv == length(x)) {
        cv_ix <- seq_along(x)
    } else {
        cv_ix <- unique(floor(seq(1, length(x), length.out = ncv)))
    }

    # Normal-transform pseudo-observations
    NX <- stats::qnorm(FX)

    # Scale covariates to [0, 1]
    min_x <- min(x)
    x <- (x - min_x) / (max(x) - min_x)

    # Log-likelihood of eta given X under a Gaussian copula model
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

    # Since Python functions are not serializable, we can't load the
    # optimization function in the global environment. Instead, the function
    # needs to be loaded on each worker. This environment stores the
    # function once it has been loaded to avoid having to do it multiple
    # times on the same worker.
    .fit_env <- new.env(parent = emptyenv())
    compute_one_eta <- function(NX, x, x0, h) {
        # Load the module if it hasn't been loaded yet
        if (!exists("fit_fun", envir = .fit_env, inherits = FALSE)) {
            .fit_env$fit_fun <- reticulate::import_from_path(
                module = "local_gaussian",
                path = system.file("python", package = "DynCopula"),
                delay_load = FALSE
            )$fit_gaussian
        }
        eta <- .fit_env$fit_fun(
            par0 = .init_par(x0, x, NX, h),
            x = x,
            NX = NX,
            h = h,
            control = control,
            x0 = x0,
            degree = degree
        )
        return(as.vector(eta))
    }

    if (length(bandwidths) > 1L) {
        message("Selecting global bandwidth")
        progressr::with_progress({
            steps <- expand.grid(bandwidth = bandwidths, cv_ix = cv_ix)
            pbar <- progressr::progressor(steps = dim(steps)[1])
            ll_all <- future.apply::future_apply(
                X = steps,
                MARGIN = 1,
                FUN = function(r) {
                    h <- r[1]
                    j <- as.integer(r[2])
                    # Estimate calibration coefficients when leaving out
                    # observation at j
                    eta_j <- compute_one_eta(
                        NX = NX[-j, ],
                        x = x[-j],
                        x0 = x[j],
                        h = h
                    )
                    # Log-likelihood of estimated calibration coefficients
                    # at observation j
                    print(paste("j =", j))
                    print(paste("h =", h))
                    print(paste("x0 =", x[j]))
                    print(paste("NAs =", sum(is.na(eta_j))))
                    print(" ")
                    print(" ")
                    ll <- loglik(NX[j, ], eta_j)
                    pbar()
                    return(c(ll, h))
                },
                future.seed = TRUE,
                future.globals = TRUE,
                future.scheduling = 0
            )
            ll_all <- sapply(bandwidths, function(h) {
                sum(ll_all[1, ][ll_all[2, ] == h])
            })
        })
        h_opt <- bandwidths[which.max(ll_all)]
    } else {
        # No global bandwidth selection
        h_opt <- bandwidths[1]
    }

    if (variable) {
        message("Selecting variable bandwidth parameters")

        # Estimate density function of covariates
        x_dens <- kde1d::dkde1d(x, kde1d::kde1d(x))
        # Geometric mean across all covariate values
        G <- exp(mean(log(x_dens)))

        # Compute adaptive bandwidths
        get_h_adapt <- function(alpha, beta, eps = 1e-6) {
            h_adapt <- beta * h_opt * (x_dens / G)^(-alpha)
            # Restrict bandwidths to (0, 1)
            h_adapt <- pmin(pmax(h_adapt, eps), 1 - eps)
            return(h_adapt)
        }

        # Generally, we find that the for any given alpha, the behavior of the
        # cross-validation likelihood is similar across beta values. So instead
        # of performing expensive nested CV, we first select alpha using
        # beta = 1. We then use this value of alpha to refine beta. Note that
        # this set of models includes the global model (alpha = 0, beta = 1).
        progressr::with_progress({
            # We only need to do CV once for (alpha, beta) = (alpha_opt, 1)
            total_steps <- (length(alpha) + (length(beta) - 1)) * length(cv_ix)
            pbar <- progressr::progressor(steps = total_steps)

            # First select alpha with beta = 1
            steps <- expand.grid(alpha = alpha, cv_ix = cv_ix)
            ll_alpha <- future.apply::future_apply(
                X = steps,
                MARGIN = 1,
                FUN = function(r) {
                    a_val <- r[1]
                    j <- r[2]
                    h_adapt <- get_h_adapt(a_val, 1)
                    eta_j <- compute_one_eta(
                        NX = NX[-j, ],
                        x = x[-j],
                        x0 = x[j],
                        h = h_adapt[j]
                    )
                    ll <- loglik(NX[j, ], eta_j)
                    pbar()
                    return(c(ll, a_val))
                },
                future.seed = TRUE,
                future.globals = TRUE,
                future.scheduling = 0
            )
            ll_alpha <- sapply(alpha, function(a) {
                sum(ll_alpha[1, ][ll_alpha[2, ] == a])
            })
            ix_opt <- which.max(ll_alpha)
            alpha_opt <- alpha[ix_opt]
            # Likelihood for (alpha, beta) = (alpha_opt, 1)
            ll_alpha_opt <- ll_alpha[ix_opt]

            # Next select beta with alpha = alpha_opt, skipping beta = 1
            beta_no_one <- beta[abs(beta - 1) > 1e-16]
            steps <- expand.grid(beta = beta_no_one, cv_ix = cv_ix)
            ll_beta <- future.apply::future_apply(
                X = steps,
                MARGIN = 1,
                FUN = function(r) {
                    b_val <- r[1]
                    j <- r[2]
                    h_adapt <- get_h_adapt(alpha_opt, b_val)
                    eta_j <- compute_one_eta(
                        NX = NX[-j, ],
                        x = x[-j],
                        x0 = x[j],
                        h = h_adapt[j]
                    )
                    ll <- loglik(NX[j, ], eta_j)
                    pbar()
                    return(c(ll, b_val))
                },
                future.seed = TRUE,
                future.globals = TRUE,
                future.scheduling = 0
            )
            ll_beta <- sapply(beta_no_one, function(b) {
                sum(ll_beta[1, ][ll_beta[2, ] == b])
            })

            # Add back in likelihood for beta = 1
            beta_opt <- c(1.0, beta_no_one)[which.max(c(ll_alpha_opt, ll_beta))]
        })
        h_final <- get_h_adapt(alpha_opt, beta_opt)
    } else {
        h_final <- h_opt
    }

    # Fit final model
    message("Fitting model")
    res_opt <- .local_fit(
        FX = FX,
        x = x,
        x0 = x,
        h = h_final,
        degree = degree,
        control = control,
        cl = cl
    )
    if (length(bandwidths) > 1L) {
        res_opt$global_cv <- data.frame(bandwidth = bandwidths, loglik = ll_all)
    }
    if (variable) {
        a_df <- data.frame(alpha = alpha, beta = 1.0, ll = ll_alpha)
        b_df <- data.frame(alpha = alpha_opt, beta = beta_no_one, ll = ll_beta)
        res_opt$var_cv <- rbind(a_df, b_df)
    }
    return(res_opt)
}

###############################################################################

#' Local likelihood estimation of conditional Gaussian copula model
#'
#' @description Fit a conditional Gaussian copula model using local likelihood.
#' Internal function.
#'
#' @param FX Matrix of pseudo-observations at covariate values.
#' @param x Vector of covariate values corresponding to \code{FX}. Must be
#' sorted and have no duplicates.
#' @param x0 Covariate values to estimate copula parameters at.
#' @param h Kernel bandwidth. Either a single value specifying a global
#' bandwidth, or a vector of length \code{length(x0)} specifying the bandwidth
#' to use at each value in \code{x0}.
#' @param degree Degree of local polynomial approximation.
#' @param control A \code{list} of control parameters for optimization.
#' @param cl A cluster for parallel computations.
#'
#' @return A list with the following components:
#' \itemize{
#'   \item \code{x}: The input argument \code{x}.
#'   \item \code{x0}: The input argument \code{x0}.
#'   \item \code{bandwidths}: The bandwidth used at each value in \code{x0}.
#'   \item \code{NX}: Normal-transformed pseudo-observations.
#'   \item \code{eta}: Matrix of estimated calibration coefficients.
#'   \item \code{rho}: Matrix of estimated pairwise correlation coefficients.
#' }
.local_fit <- function(FX,
                       x,
                       x0,
                       h,
                       degree,
                       control,
                       cl) {
    # Basic checks
    assert_that(
        is.vector(x, mode = "numeric"),
        is.matrix(FX),
        is.numeric(FX),
        dim(FX)[1] == length(x),
        dim(FX)[2] > 1L,
        length(h) == 1L || length(h) == length(x0),
        is.numeric(h) && all(h > 0) && all(h < 1),
        is.list(control),
        !anyDuplicated(x) && !is.unsorted(x),
        is.numeric(x0),
        is.numeric(degree) && degree >= 0
    )
    degree <- as.integer(degree)

    if (length(h) == 1L) {
        h <- rep(h, length(x0))
    }

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

    # Since Python functions are not serializable, we can't load the
    # optimization function in the global environment. Instead, the function
    # needs to be loaded on each worker. This environment stores the
    # function once it has been loaded to avoid having to do it multiple
    # times on the same worker.
    .fit_env <- new.env(parent = emptyenv())
    fit_gaussian <- function(par0, x, NX, h, x0) {
        # Load the module if it hasn't been loaded yet
        if (!exists("fit_fun", envir = .fit_env, inherits = FALSE)) {
            .fit_env$fit_fun <- reticulate::import_from_path(
                module = "local_gaussian",
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
            X = seq_along(x0),
            FUN = function(i) {
                y <- fit_gaussian(
                    par0 = .init_par(x0[i], x, NX, h[i]),
                    x = x,
                    NX = NX,
                    h = h[i],
                    x0 = x0[i]
                )
                pbar()
                return(y)
            },
            future.seed = TRUE,
            future.globals = TRUE,
            future.scheduling = 0
        )
    })

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

    # Scale covariates back to their original scale
    x <- x * dx + min_x
    x0 <- x0 * dx + min_x

    return(list(
        x = x,
        x0 = x0,
        bandwidths = h,
        NX = NX,
        eta = Hhat,
        rho = Rhat
    ))
}

###############################################################################
