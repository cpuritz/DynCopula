###############################################################################

#' Fitting dynamic Gaussian copula model
#'
#' @description Fit a dynamic Gaussian copula to time series data.
#'
#' @param FX Matrix of pseudo-observations F_i(X_i) at time points.
#' @param FXm Matrix of left limits of distribution functions F_i(X_i - 1) at
#' time points. If \code{NULL}, it is assumed that all margins are continuous.
#' @param x Vector of time points corresponding to \code{FX} and \code{FXm}.
#' @param x0 Time points to estimate copula parameters at.
#' @param band Kernel bandwidth. Default is \code{0.10}. Must satisfy
#' \code{0 < band < 1}.
#' @param optMethod Optimization method. Either \code{"SGD"} or \code{"L-BFGS"}.
#' \code{"SGD"} implements stochastic gradient descent with gradients computed
#' using automatic differentiation. \code{"L-BFGS"} implements the quasi-Newton
#' limited memory BFGS with gradients estimated numerically. Default is
#' \code{"SGD"}.
#' @param control A \code{list} of control parameters for optimization.
#' @param cores Number of cores to use for parallel optimization. Parallelized
#' over \code{x0}. Default is \code{1}.
#' @param cores2 If using \code{"L-BFGS"}, each optimization step can be
#' parallelized, in addition to parallelization over \code{x0}. The total
#' number of cores used is \code{cores * cores2}. Default is \code{1}.
#'
#' @details The \code{control} argument is a list that supplies control
#' parameters for optimization. For SGD, the following parameters can be
#' supplied:
#' \itemize{
#'   \item \code{maxit} Maximum number of iterations. Default is \code{100}.
#'   \item \code{lr} Learning rate for RMSProp. Default is \code{1e-2}.
#'   \item \code{reltol} Relative convergence tolerance.
#'   \item \code{patience} Optimization stops if the relative log-likelihood
#'   has not decreased by a factor of \code{reltol} within the last
#'   \code{patience} iterations. Default is \code{5}.
#'   \item \code{weight_decay} Weight decay for RMSProp. Default is \code{0}.
#'   \item \code{R0} Fraction of neighbors to use to estimate the initial
#'   correlation matrix. Default is \code{0.10}. Must satisfy \code{0 < R0 <= 1}.
#' }
#' For L-BFGS, the available control parameters and default values are the same
#' as those of the \code{\link[stats]{optim}} function from the \strong{stats}
#' package.
#'
#' @return \itemize{
#'   \item \code{x}: The input argument \code{x}.
#'   \item \code{x0}: The input argument \code{x0}.
#'   \item \code{eta}: Matrix of estimated coefficients in the unconstrained
#'   space.
#'   \item \code{rho}: Matrix of estimated pairwise correlation coefficients.
#'   \item \code{convergence}: Convergence codes for each coefficient.
#'   \code{0} indicates successful completion. \code{1} indicates that
#'   iteration limit had been reached.
#'   \item \code{loss}: Loss history for each coefficient. Only available if
#'   \code{optMethod} is \code{"SGD"}.
#' }
#'
#' @export
fit_dynamic_gaussian <- function(FX,
                                 FXm = NULL,
                                 x,
                                 x0,
                                 band = 0.10,
                                 optMethod = c("SGD", "L-BFGS"),
                                 control = list(),
                                 cores = 1L,
                                 cores2 = 1L) {
    optMethod <- match.arg(optMethod)
    if (is.null(FXm)) {
        optMethod <- "SGD"
    }

    # Basic checks
    assertthat::assert_that(
        dim(FX)[1] == length(x),
        is.null(FXm) || all(dim(FX) == dim(FXm)),
        min(x0) >= min(x),
        max(x0) <= max(x),
        is.character(optMethod),
        is.numeric(band) && band > 0 && band < 1,
        is.numeric(cores) && cores >= 1,
        is.numeric(cores2) && cores2 >= 1
    )

    cores <- as.integer(cores)
    cores2 <- as.integer(cores2)
    if (cores2 > 1L && optMethod == "SGD") {
        message("'cores2' is not used for SGD optimization.")
    }

    ## Control parameters
    # Control parameters common to both optimization methods
    assertthat::assert_that(all(sapply(control, is.numeric)))
    if ("maxit" %in% names(control)) {
        control$maxit <- as.integer(control$maxit)
    } else {
        control <- c(list(maxit = 100L), control)
    }
    if (!"R0" %in% names(control)) {
        control$R0 <- 0.10
    }
    assertthat::assert_that(
        control$maxit >= 1L,
        control$R0 > 0 && control$R0 <= 1
    )

    # SGD specific control parameters
    if (optMethod == "SGD") {
        if (!"lr" %in% names(control)) {
            control$lr <- 1e-2
        }
        if (!"reltol" %in% names(control)) {
            control$reltol <- 1e-5
        }
        if (!"patience" %in% names(control)) {
            control$patience <- 5L
        } else {
            control$patience <- as.integer(control$patience)
        }
        if (!"weight_decay" %in% names(control)) {
            control$weight_decay <- 0
        }
        assertthat::assert_that(
            control$lr > 0,
            control$reltol > 0,
            control$patience >= 1L,
            control$weight_decay >= 0
        )
    }

    # Sort covariate values and scale to [0, 1]
    ord <- order(x)
    x <- x[ord]
    FX <- FX[ord, ]
    if (!is.null(FXm)) {
        FXm <- FXm[ord, ]
    }
    min_x <- x[1]
    dx <- x[length(x)] - min_x
    x <- (x - min_x) / dx
    x0 <- (x0 - min_x) / dx

    args <- as.list(environment())
    if (is.null(FXm)) {
        fun <- .fit_dynamic_gaussian_cts
    } else {
        fun <- .fit_dynamic_gaussian_discrete
    }
    res <- do.call(fun, args[names(formals(fun))])
    eta_vals <- res[["eta_vals"]]

    # Estimated eta matrix
    Hhat <- do.call(rbind, eta_vals)
    # Estimated correlation matrix
    Rhat <- t(apply(Hhat, 1, function(v) {
        copula::P2p(vec2cor(v))
    }))

    ix_lab <- apply(utils::combn(seq(dim(FX)[2]), 2), 2, function(x) {
        paste(x, collapse = '')
    })
    colnames(Hhat) <- paste0("eta", ix_lab)
    colnames(Rhat) <- paste0("rho", ix_lab)

    res <- c(res[setdiff(names(res), "eta_vals")],
             list(x = x * dx + min_x,
                  x0 = x0 * dx + min_x,
                  eta = data.frame(Hhat),
                  rho = data.frame(Rhat)))
    return(res)
}

###############################################################################

#' Fitting dynamic t copula model
#'
#' @description Fit a dynamic t copula to time series data.
#'
#' @param FX Matrix of pseudo-observations F_i(X_i) at time points.
#' @param nu Degrees of freedom.
#' @param x Vector of time points corresponding to \code{FX} and \code{FXm}.
#' @param x0 Time points to estimate copula parameters at.
#' @param band Kernel bandwidth. Default is \code{0.10}. Must satisfy
#' \code{0 < band < 1}.
#' @param control A \code{list} of control parameters for optimization.
#' @param cores Number of cores to use for parallel optimization. Parallelized
#' over \code{x0}. Default is \code{1}.
#'
#' @details The \code{control} argument is a list that supplies control
#' parameters for optimization. The following parameters can be
#' supplied:
#' \itemize{
#'   \item \code{maxit} Maximum number of iterations. Default is \code{100}.
#'   \item \code{lr} Learning rate for RMSProp. Default is \code{1e-2}.
#'   \item \code{reltol} Relative convergence tolerance.
#'   \item \code{patience} Optimization stops if the relative log-likelihood
#'   has not decreased by a factor of \code{reltol} within the last
#'   \code{patience} iterations. Default is \code{5}.
#'   \item \code{weight_decay} Weight decay for RMSProp. Default is \code{0}.
#'   \item \code{R0} Fraction of neighbors to use to estimate the initial
#'   correlation matrix. Default is \code{0.10}. Must satisfy \code{0 < R0 <= 1}.
#' }
#'
#' @return \itemize{
#'   \item \code{x}: The input argument \code{x}.
#'   \item \code{x0}: The input argument \code{x0}.
#'   \item \code{eta}: Matrix of estimated coefficients in the unconstrained
#'   space.
#'   \item \code{rho}: Matrix of estimated pairwise correlation coefficients.
#'   \item \code{convergence}: Convergence codes for each coefficient.
#'   \code{0} indicates successful completion. \code{1} indicates that
#'   iteration limit had been reached.
#'   \item \code{loss}: Loss history for each coefficient.
#' }
#'
#' @export
fit_dynamic_t <- function(FX,
                          nu,
                          x,
                          x0,
                          band = 0.10,
                          control = list(),
                          cores = 1L) {
    # Basic checks
    assertthat::assert_that(
        dim(FX)[1] == length(x),
        nu >= 1,
        min(x0) >= min(x),
        max(x0) <= max(x),
        is.numeric(band) && band > 0 && band < 1,
        is.numeric(cores) && cores >= 1
    )
    cores <- as.integer(cores)

    # Control parameters
    assertthat::assert_that(all(sapply(control, is.numeric)))
    if ("maxit" %in% names(control)) {
        control$maxit <- as.integer(control$maxit)
    } else {
        control <- c(list(maxit = 100L), control)
    }
    if (!"R0" %in% names(control)) {
        control$R0 <- 0.10
    }
    if (!"lr" %in% names(control)) {
        control <- c(list(lr = 1e-2), control)
    }
    if (!"reltol" %in% names(control)) {
        control <- c(list(reltol = 1e-5), control)
    }
    if (!"patience" %in% names(control)) {
        control <- c(list(patience = 5L), control)
    } else {
        control$patience <- as.integer(control$patience)
    }
    if (!"weight_decay" %in% names(control)) {
        control <- c(list(weight_decay = 0), control)
    }
    assertthat::assert_that(
        control$maxit >= 1L,
        control$lr > 0,
        control$reltol > 0,
        control$patience >= 1L,
        control$weight_decay >= 0,
        control$R0 > 0 && control$R0 <= 1
    )

    # Sort covariate values and scale to [0, 1]
    ord <- order(x)
    x <- x[ord]
    FX <- FX[ord, ]
    min_x <- x[1]
    dx <- x[length(x)] - min_x
    x <- (x - min_x) / dx
    x0 <- (x0 - min_x) / dx

    args <- as.list(environment())
    fun <- .fit_dynamic_t_cts
    res <- do.call(fun, args[names(formals(fun))])
    eta_vals <- res[["eta_vals"]]

    # Estimated eta matrix
    Hhat <- do.call(rbind, eta_vals)
    # Estimated correlation matrix
    Rhat <- t(apply(Hhat, 1, function(v) {
        copula::P2p(vec2cor(v))
    }))

    ix_lab <- apply(utils::combn(seq(dim(FX)[2]), 2), 2, function(x) {
        paste(x, collapse = '')
    })
    colnames(Hhat) <- paste0("eta", ix_lab)
    colnames(Rhat) <- paste0("rho", ix_lab)

    res <- c(res[setdiff(names(res), "eta_vals")],
             list(x = x * dx + min_x,
                  x0 = x0 * dx + min_x,
                  eta = data.frame(Hhat),
                  rho = data.frame(Rhat)))
    return(res)
}

###############################################################################

#' Fitting dynamic Gaussian copula model for discrete data
#'
#' @inheritParams fit_dynamic_gaussian
#'
#' @return Estimated coefficients and information about optimization
.fit_dynamic_gaussian_discrete <- function(FX,
                                           FXm,
                                           x,
                                           x0,
                                           band,
                                           optMethod,
                                           control,
                                           cores,
                                           cores2) {
    # Set up futures plan
    cl <- parallel::makeCluster(cores)
    future::plan(future::cluster, workers = cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    # Convert to standard normal margins
    NX <- stats::qnorm(FX)
    NXm <- stats::qnorm(FXm)

    R0 <- control$R0
    control$R0 <- NULL

    lbfgs_optim <- function(x0i) {
        # Use points nearby to estimate initial correlation matrix
        dx <- x - x0i
        max_thr <- sort(abs(dx))[min(length(dx), dim(NX)[2])]
        thr <- max(stats::quantile(abs(dx), R0), max_thr)
        NX_loc <- NX[which(abs(dx) <= thr), ]

        eta0 <- cor2vec(stats::cor(NX_loc, method = "pearson"))
        par0 <- c(eta0, rep(0, length(eta0)))

        # Negative local log likelihood function
        loglik <- function(par) {
            par <- as.numeric(par)

            # Epanechnikov kernel weights
            wgt <- 3/(4 * band) * pmax(1 - (dx / band)^2, 0)

            # Only compute densities when weight is nonzero
            pos_ix <- which(wgt > 0)
            P <- sapply(pos_ix, function(j) {
                # Reconstruct correlation matrix
                npar <- as.integer(length(par) / 2)
                eta0 <- par[1:npar]
                eta1 <- par[(npar + 1):(2 * npar)]
                R <- vec2cor(eta0 + eta1 * dx[j])

                if (min(eigen(R, only.values = TRUE)$values) <= 0) {
                    # Matrix is theoretically PD but may not be
                    # numerically PD
                    return(0)
                } else {
                    # Numerical evaluation of Gaussian CDF
                    return(tryCatch(
                        TruncatedNormal::mvNcdf(
                            l = NXm[j, ],
                            u = NX[j, ],
                            Sig = R,
                            n = 1e3
                        )$prob,
                        error = function(e) {
                            # Error often thrown for poorly conditioned
                            # matrices, return probability of 0
                            0
                        }
                    ))
                }
            })
            if (any(P <= 0)) {
                # Probabilities are theoretically nonnegative but may
                # be numerically negative. Return large number instead
                # of infinity.
                return(1e20)
            } else {
                return(-sum(wgt[pos_ix] * log(P)))
            }
        }

        # Create inner cluster
        cl2 <- parallel::makeCluster(cores2)
        on.exit(parallel::stopCluster(cl2), add = TRUE)
        parallel::clusterExport(
            cl = cl2,
            varlist = c("vec2cor", "band", "NX", "NXm")
        )

        # Parallelized L-BFGS optimization
        res <- optimParallel::optimParallel(
            par = par0,
            fn = loglik,
            parallel = list(cl = cl2),
            control = control
        )
        pbar()
        return(list(
            par = res$par[1:length(eta0)],
            convergence = res$convergence
        ))
    }

    sgd_safe <- function(x0i) {
        # Use points nearby to estimate initial correlation matrix
        dx <- x - x0i
        max_thr <- sort(abs(dx))[min(length(dx), dim(NX)[2])]
        thr <- max(stats::quantile(abs(dx), R0), max_thr)
        NX_loc <- NX[which(abs(dx) <= thr), ]

        eta0 <- cor2vec(stats::cor(NX_loc, method = "pearson"))
        par0 <- c(eta0, rep(0, length(eta0)))

        dyn_fit <- py_load("dynamic_gaussian")$fit_discrete_gaussian
        res <- dyn_fit(par0, dx, NX, NXm, band, control)
        loss <- res$hist

        if (res$convergence == 2) {
            message("Failed prematurely")
            par1 <- res$par[1:length(eta0)]
            n_itr <- max(control$maxit - length(res$hist), 1L)
            message("Running ", n_itr, " additional iterations")

            loglik_lbfgs <- function(par) {
                par <- as.numeric(par)

                # Epanechnikov kernel weights
                wgt <- 3/(4 * band) * pmax(1 - (dx / band)^2, 0)

                # Only compute densities when weight is nonzero
                pos_ix <- which(wgt > 0)
                P <- sapply(pos_ix, function(j) {
                    # Reconstruct correlation matrix
                    npar <- as.integer(length(par) / 2)
                    eta0 <- par[1:npar]
                    eta1 <- par[(npar + 1):(2 * npar)]
                    R <- vec2cor(eta0 + eta1 * dx[j])

                    if (min(eigen(R, only.values = TRUE)$values) <= 0) {
                        # Matrix is theoretically PD but may not be
                        # numerically PD
                        return(0)
                    } else {
                        # Numerical evaluation of Gaussian CDF
                        return(tryCatch(
                            TruncatedNormal::mvNcdf(
                                l = NXm[j, ],
                                u = NX[j, ],
                                Sig = R,
                                n = 1e3
                            )$prob,
                            error = function(e) {
                                # Error often thrown for poorly conditioned
                                # matrices, return probability of 0
                                0
                            }
                        ))
                    }
                })
                if (any(P <= 0)) {
                    # Probabilities are theoretically nonnegative but may
                    # be numerically negative. Return large number instead
                    # of infinity.
                    return(1e20)
                } else {
                    return(-sum(wgt[pos_ix] * log(P)))
                }
            }
            res <- stats::optim(
                par = par1,
                fn = loglik_lbfgs,
                method = "L-BFGS-B",
                control = list(maxit = n_itr)
            )
        }
        pbar()

        return(list(
            par = res$par[1:length(eta0)],
            convergence = res$convergence,
            loss = loss
        ))
    }

    progressr::with_progress({
        pbar <- progressr::progressor(along = x0)
        if (optMethod == "L-BFGS") {
            opt_res <- future.apply::future_lapply(
                X = x0,
                FUN = lbfgs_optim,
                future.seed = TRUE,
                future.packages = c("parallel", "optimParallel")
            )
            return(list(
                eta_vals = lapply(opt_res, '[[', "par"),
                convergence = sapply(opt_res, '[[', "convergence")
            ))
        } else if (optMethod == "SGD") {
            opt_res <- future.apply::future_lapply(
                X = x0,
                FUN = sgd_safe,
                future.seed = TRUE
            )
            return(list(
                eta_vals = lapply(opt_res, '[[', "par"),
                convergence = sapply(opt_res, '[[', "convergence"),
                loss = lapply(opt_res, '[[', "loss")
            ))
        }
    })
}

###############################################################################

#' Fitting dynamic Gaussian copula model for continuous data
#'
#' @inheritParams fit_dynamic_gaussian
#'
#' @return Estimated coefficients and information about optimization
.fit_dynamic_gaussian_cts <- function(FX,
                                      x,
                                      x0,
                                      band,
                                      control,
                                      cores) {
    # Convert to standard normal margins
    NX <- stats::qnorm(FX)

    # Set up futures plan
    cl <- parallel::makeCluster(cores)
    future::plan(future::cluster, workers = cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    progressr::with_progress({
        pbar <- progressr::progressor(along = x0)

        sgd_optim <- function(x0i) {
            # Use points nearby to estimate initial correlation matrix
            dx <- x - x0i
            max_thr <- sort(abs(dx))[min(length(dx), dim(NX)[2])]
            thr <- max(stats::quantile(abs(dx), control$R0), max_thr)
            NX_loc <- NX[which(abs(dx) <= thr), ]

            eta0 <- cor2vec(stats::cor(NX_loc, method = "pearson"))
            par0 <- c(eta0, rep(0, length(eta0)))

            dyn_fit <- py_load("dynamic_gaussian")$fit_continuous_gaussian
            res <- dyn_fit(par0, dx, NX, band, control)
            pbar()

            return(list(
                par = res$par[1:length(eta0)],
                loss = res$hist,
                convergence = res$convergence
            ))
        }

        res <- future.apply::future_lapply(
            X = x0,
            FUN = sgd_optim,
            future.seed = TRUE
        )
    })

    return(list(
        eta_vals = lapply(res, '[[', "par"),
        convergence = sapply(res, '[[', "convergence"),
        loss = lapply(res, '[[', "loss")
    ))
}

###############################################################################

#' Fitting dynamic t copula model for continuous data
#'
#' @inheritParams fit_dynamic_t
#'
#' @return Estimated coefficients and information about optimization
.fit_dynamic_t_cts <- function(FX,
                               nu,
                               x,
                               x0,
                               band,
                               control,
                               cores) {
    # Convert to t margins
    TX <- stats::qt(FX, df = nu)

    # Set up futures plan
    cl <- parallel::makeCluster(cores)
    future::plan(future::cluster, workers = cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    progressr::with_progress({
        pbar <- progressr::progressor(along = x0)

        optimize <- function(x0i) {
            # Use points nearby to estimate initial correlation matrix
            dx <- x - x0i
            max_thr <- sort(abs(dx))[min(length(dx), dim(TX)[2])]
            thr <- max(stats::quantile(abs(dx), control$R0), max_thr)
            TX_loc <- TX[which(abs(dx) <= thr), ]
            eta0 <- cor2vec(stats::cor(TX_loc, method = "pearson"))
            par0 <- c(eta0, rep(0, length(eta0)))

            dyn_fit <- py_load("dynamic_t")$fit_continuous_t
            res <- dyn_fit(par0, nu, dx, TX, band, control)
            pbar()

            return(list(
                par = res$par[1:length(eta0)],
                loss = res$hist,
                convergence = res$convergence
            ))
        }

        res <- future.apply::future_lapply(
            X = x0,
            FUN = optimize,
            future.seed = TRUE
        )
    })

    return(list(
        eta_vals = lapply(res, '[[', "par"),
        convergence = sapply(res, '[[', "convergence"),
        loss = lapply(res, '[[', "loss")
    ))
}

###############################################################################

#' Log-likelihood
#'
#' Compute log-likelihood
#'
#' @param FX Matrix of pseudo-observations F_i(X_i).
#' @param FXm Matrix of left limits of distribution functions F_i(X_i - 1).
#' If \code{NULL}, it is assumed that all margins are continuous.
#' @param R Correlation matrix.
#'
#' @return Log-likelihood
#'
#' @export
loglik <- function(FX, FXm = NULL, R)  {
    NX <- stats::qnorm(FX)
    if (is.null(FXm)) {
        ll <- apply(NX, 1, function(x) {
            mvtnorm::dmvnorm(x = x, sigma = R)
        })
    } else {
        NXm <- stats::qnorm(FXm)
        ll <- sapply(seq(dim(FX)[1]), function(i) {
            TruncatedNormal::mvNcdf(
                l = NXm[i, ],
                u = NX[i, ],
                Sig = R,
                n = 1e3
            )$prob
        })
    }
    return(sum(ll))
}

###############################################################################
