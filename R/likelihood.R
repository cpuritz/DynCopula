###############################################################################

#' Time-varying Gaussian log-likelihood
#'
#' @description Compute log-likelihood for a time-varying Gaussian copula. All
#' arguments should correspond to the same time points.
#'
#' @param FX Matrix of pseudo-observations F_i(X_i).
#' @param FXm Matrix of left limits of distribution functions F_i(X_i - 1).
#' If \code{NULL}, it is assumed that all margins are continuous.
#' @param rho Correlation coefficients.
#'
#' @return Log-likelihood
#'
#' @export
loglik_gaussian <- function(FX, FXm = NULL, rho)  {
    assertthat::assert_that(
        dim(FX)[1] == dim(rho)[1],
        is.null(FXm) || all(dim(FX) == dim(FXm))
    )

    NX <- stats::qnorm(FX)
    if (is.null(FXm)) {
        ll <- sapply(seq(dim(FX)[1]), function(i) {
            mvtnorm::dmvnorm(x = NX[i, ], sigma = copula::p2P(rho[i, ]))
        })
    } else {
        NXm <- stats::qnorm(FXm)
        # ll <- sapply(seq(dim(FX)[1]), function(i) {
        #     TruncatedNormal::mvNcdf(
        #         l = NXm[i, ],
        #         u = NX[i, ],
        #         Sig = R,
        #         n = 1e3
        #     )$prob
        # })
        stop("Not implemented")
    }
    return(sum(log(ll)))
}

###############################################################################

#' t log-likelihood
#'
#' @description Compute log-likelihood for a t copula.
#'
#' @param FX Matrix of pseudo-observations F_i(X_i).
#' @param FXm Matrix of left limits of distribution functions F_i(X_i - 1).
#' If \code{NULL}, it is assumed that all margins are continuous.
#' @param nu Degrees of freedom.
#' @param R Correlation matrix.
#'
#' @return Log-likelihood
#'
#' @export
loglik_t <- function(FX, FXm = NULL, nu, R)  {
    assertthat::assert_that(
        all(dim(R) == dim(FX)[2]),
        is.null(FXm) || all(dim(FX) == dim(FXm)),
        min(eigen(R, only.values = TRUE)$values) > 0,
        nu > 0
    )

    TX <- stats::qt(FX, df = nu)
    if (is.null(FXm)) {
        ll <- apply(TX, 1, function(x) {
            mvtnorm::dmvt(x, sigma = R, df = nu)
        })
    } else {
        TXm <- stats::qt(FXm, df = nu)
        stop("Not implemented")
    }
    return(sum(log(ll)))
}

###############################################################################
