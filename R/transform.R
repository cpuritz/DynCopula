###############################################################################

#' Correlation matrix to vector
#'
#' @description Convert a correlation matrix to an unconstrained vector.
#'
#' @param R A \code{d}x\code{d} correlation matrix.
#'
#' @returns A vector of length \code{choose(d, 2)}.
#'
#' @details For numerical stability, the maximum magnitude of correlation
#' coefficients is capped at \code{1 - 1e-4}.
#'
#' @export
cor2vec <- function(R) {
    scale <- 0.5
    rho_max <- 1 - 1e-4
    eps <- 1e-12

    d <- dim(R)[1]

    # Compute Cholesky factor
    L <- t(chol(R))

    # Parametrize the Cholesky factor to have a real and positive diagonal
    H <- matrix(0, nrow = d, ncol = d)
    diag(H) <- 1
    H[, 1] <- L[, 1]
    if (d > 2) {
        for (i in 3:d) {
            H[i, 2:(i - 1)] <- L[i, 2:(i - 1)] / sqrt(1 - cumsum(L[i, 1:(i - 2)]^2))
        }
    }

    # Map to an unconstrained vector in R^(choose(d, 2))
    vH <- H[lower.tri(H)]
    vH_overflow <- (abs(vH) > rho_max)
    vH[vH_overflow] <- sign(vH[vH_overflow]) * (rho_max - eps)
    vH <- (1 / scale) * atanh(vH / rho_max)
    return(vH)
}

###############################################################################

#' Vector to correlation matrix
#'
#' @description Convert an unconstrained vector to a correlation matrix.
#'
#' @param v A vector of length \code{choose(d, 2)}.
#'
#' @returns A \code{d}x\code{d} correlation matrix.
#'
#' @details For numerical stability, the maximum magnitude of correlation
#' coefficients is capped at \code{1 - 1e-4}.
#'
#' @export
vec2cor <- function(v) {
    scale <- 0.5
    rho_max <- 1 - 1e-4

    d <- as.integer((1 + sqrt(1 + 8 * length(v))) / 2)

    # Map unconstrained vector to d x d matrix with entries in (-1, 1)
    H <- matrix(0, nrow = d, ncol = d)
    diag(H) <- 1
    H[lower.tri(H)] <- rho_max * tanh(scale * v)

    # Map back to Cholesky factor space
    for (i in 2:d) {
        H[i, 2:i] <- H[i, 2:i] * sqrt(cumprod(1 - H[i, 1:(i - 1)]^2))
    }

    # Return correlation matrix
    return(tcrossprod(H, H))
}

###############################################################################
