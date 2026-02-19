###############################################################################

#' Correlation matrix to vector
#'
#' @description Convert a correlation matrix to an unconstrained vector. This is
#' the inverse of \link[DynCopula]{vec2cor}.
#'
#' @param R A \code{d}x\code{d} correlation matrix.
#'
#' @returns A vector of length \code{choose(d, 2)}.
#'
#' @export
cor2vec <- function(R) {
    scale <- 0.5

    d <- dim(R)[1]

    # Compute lower-triangular Cholesky factor
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
    vH <- (1 / scale) * atanh(vH)
    return(vH)
}

###############################################################################

#' Vector to correlation matrix
#'
#' @description Convert an unconstrained vector to a correlation matrix. This is
#' the inverse of \link[DynCopula]{cor2vec}.
#'
#' @param v A vector of length \code{choose(d, 2)}.
#'
#' @returns A \code{d}x\code{d} correlation matrix.
#'
#' @export
vec2cor <- function(v) {
    scale <- 0.5

    d <- as.integer((1 + sqrt(1 + 8 * length(v))) / 2)

    # Map unconstrained vector to d x d matrix with entries in (-1, 1)
    H <- matrix(0, nrow = d, ncol = d)
    diag(H) <- 1
    H[lower.tri(H)] <- tanh(scale * v)

    # Map back to Cholesky factor space
    for (i in 2:d) {
        H[i, 2:i] <- H[i, 2:i] * sqrt(cumprod(1 - H[i, 1:(i - 1)]^2))
    }

    # Compute correlation matrix
    R <- tcrossprod(H, H)

    return(R)
}

###############################################################################
