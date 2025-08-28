###############################################################################

#' Fitting dynamic Gaussian copula model
#'
#' @description Fit a dynamic Gaussian copula to a time series.
#'
#' @param FX Matrix of pseudo-observations F_i(X_i) at time points.
#' @param FXm Matrix of left limits of pseudo-observations at time points. If
#' \code{NULL}, data is assumed to be continuous. Otherwise, data is assumed
#' to be count-valued.
#' @param x Vector of time points corresponding to \code{FX} and \code{FXm}.
#' @param x0 Time points to estimate copula parameters at.
#' @param x0_ix Indices of time points to estimate copula parameters at.
#' @param h Kernel bandwidth. Must satisfy \code{0 < h < 1}.
#' @param control A \code{list} of control parameters for optimization.
#' @param cores Number of cores to use for parallel optimization. Parallelized
#' over \code{x0}. Default is \code{1}.
#'
#' @details Optimization is performed using gradient descent. The \code{control}
#' argument is a list that supplies control parameters for optimization. The
#' following parameters can be supplied:
#' \itemize{
#'   \item \code{max_itr} Maximum number of iterations. Default is \code{100}.
#'   \item \code{reltol} Relative convergence tolerance. Default is \code{1e-4}.
#'   \item \code{lr} Learning rate. Default is \code{1e-5}.
#'   \item \code{patience} Optimization stops if the relative log-likelihood
#'   has not decreased by a factor of \code{reltol} within the last
#'   \code{patience} iterations. Default is \code{3}.
#'   \item \code{momentum} Momentum factor. Default is \code{0.9}.
#'   \item \code{max_grad} Gradients with an L-infinity norm above this value
#'   are clipped. Default is \code{1e3}.
#' }
#' If no improvement is made in the first \code{patience} iterations, the
#' learning rate is increased by a factor of \code{10}.
#'
#' @return A list with the following components:
#' \itemize{
#'   \item \code{x}: The input argument \code{x}.
#'   \item \code{x0}: The input argument \code{x0}.
#'   \item \code{h}: The input argument \code{h}.
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
fit_dynamic_gaussian <- function(FX,
                                 FXm = NULL,
                                 x,
                                 x0 = NULL,
                                 x0_ix = NULL,
                                 h,
                                 control = list(),
                                 cores = 1L) {
    # Basic checks
    assertthat::assert_that(
        is.vector(x, mode = "numeric"),
        is.numeric(FX) && is.matrix(FX),
        is.null(FXm) || (is.numeric(FXm) && is.matrix(FXm)),
        dim(FX)[1] == length(x),
        is.numeric(h) && h > 0 && h < 1,
        is.numeric(cores) && cores >= 1,
        is.list(control)
    )
    cores <- as.integer(cores)

    if (!is.null(FXm)) {
        # Data is count-valued, jitter to produce continuous data
        V <- matrix(stats::runif(prod(dim(FX))), ncol = ncol(FX))
        FX <- FXm + V * (FX - FXm)
    }

    # Control parameters
    defaults <- list(
        max_itr = 100L,
        lr = 1e-5,
        reltol = 1e-4,
        patience = 3L,
        momentum = 0.9,
        max_grad = 1e3
    )
    control <- utils::modifyList(defaults, control)
    assertthat::assert_that(all(names(control) %in% names(defaults)))

    assertthat::assert_that(
        all(sapply(control, is.numeric)),
        control$max_itr >= 1L,
        control$lr > 0,
        control$reltol > 0,
        control$patience >= 1L,
        control$momentum >= 0 && control$momentum < 1,
        control$max_grad > 0
    )

    # Sort covariate values and scale to [0, 1]
    ord <- order(x)
    x <- x[ord]
    FX <- FX[ord, ]
    min_x <- x[1]
    dx <- x[length(x)] - min_x
    x <- (x - min_x) / dx
    #x0 <- (x0 - min_x) / dx

    # Convert to standard normal margins
    NX <- stats::qnorm(FX)

    # Set up futures plan
    cl <- parallel::makeCluster(cores)
    future::plan(future::cluster, workers = cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    if (is.null(x0)) {
        # Time points to perform inference at are time points at which the time
        # series was sampled. The only difference is that AIC will be computed.
        vec <- x0_ix
        optim <- function(t0) {
            py_load("dynamic_gaussian")$fit_gaussian(
                par0 = .init_par(x[t0], x, NX, h),
                x = x,
                NX = NX,
                h = h,
                control = control,
                i = t0 - 1
            )
        }
    } else {
        # Time points to perform inference at are not necessarily time points
        # at which the time series was sampled. AIC will not be computed.
        vec <- x0
        optim <- function(t0) {
            py_load("dynamic_gaussian")$fit_gaussian(
                par0 = .init_par(t0, x, NX, h),
                x = x,
                NX = NX,
                h = h,
                control = control,
                x0 = t0
            )
        }
    }

    # Parallelized with progress bar
    progressr::with_progress({
        pbar <- progressr::progressor(along = vec)
        res <- future.apply::future_lapply(
            X = vec,
            FUN = function(i) {
                y <- optim(i)
                pbar()
                return(y)
            },
            future.seed = TRUE
        )
    })

    output <- list(
        eta_vals = lapply(res, '[[', "par"),
        convergence = sapply(res, '[[', "convergence"),
        loss = lapply(res, '[[', "loss_hist"),
        hist = lapply(res, '[[', "eta_hist")
    )
    if ("aic" %in% names(res[[1]])) {
        output <- c(output, list(aic = sapply(res, '[[', "aic")))
    }

    # Estimated eta matrix
    Hhat <- do.call(rbind, output[["eta_vals"]])
    # Estimated correlation matrix
    Rhat <- t(apply(Hhat, 1, function(v) {
        copula::P2p(vec2cor(v))
    }))

    d <- dim(FX)[2]
    if (d == 2) {
        Rhat <- t(Rhat)
    }

    # Add numbered eta/rho labels
    ix_lab <- apply(utils::combn(seq_len(d), 2), 2, function(x) {
        paste(x, collapse = '')
    })
    colnames(Hhat) <- paste0("eta", ix_lab)
    colnames(Rhat) <- paste0("rho", ix_lab)

    if ("hist" %in% names(output) && is.matrix(output$hist)) {
        labs <- paste0("eta", ix_lab)
        output$hist <- t(lapply(output$hist, function(x) {
            rownames(x) <- labs
            colnames(x) <- seq_len(dim(x)[2])
            return(x)
        }))
    }

    output <- c(
        output[setdiff(names(output), "eta_vals")],
        list(x = x * dx + min_x,
             x0 = x0, # * dx + min_x,
             h = h,
             eta = data.frame(Hhat),
             rho = data.frame(Rhat))
    )
    return(output)
}

###############################################################################

#' Initial parameter guess
#'
#' @param x0 Time point to center around.
#' @param x Vector of time points.
#' @param NX Matrix of normal-transformed pseudo-observations.
#' @param h Bandwidth.
#'
#' @returns Parameter vector
.init_par <- function(x0, x, NX, h) {
    # Use points nearby to estimate initial correlation matrix
    dx <- abs(x0 - x)
    # Ensure that at least d + 1 points are used to avoid a singular
    # correlation matrix
    d <- dim(NX)[2]
    thr <- max((d + 1) / length(x), h)
    NX_loc <- NX[which(abs(dx) <= thr), ]
    cor_loc <- stats::cor(NX_loc, method = "pearson")

    # If computing Cholesky decomposition would fail (cor_loc is probably
    # not numerically positive definite), use identity matrix instead
    if (inherits(try(chol(cor_loc), silent = TRUE), "try-error")) {
        cor_loc <- diag(d)
    }

    # Convert to vector
    par0 <- cor2vec(cor_loc)
    return(par0)
}

###############################################################################
