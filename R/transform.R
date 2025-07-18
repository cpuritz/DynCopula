###############################################################################

#' Correlation matrix to vector
#'
#' @description Convert a correlation matrix to vector.
#'
#' @param R A correlation matrix.
#'
#' @returns A vector.
#'
#' @export
cor2vec <- function(R) {
    d <- dim(R)[1]

    # Compute Cholesky factor
    L <- t(chol(R))

    # Parametrize the Cholesky factor to have a real and positive diagonal
    H <- matrix(0, nrow = d, ncol = d)
    diag(H) <- 1
    H[, 1] <- L[, 1]
    for (i in 3:d) {
        H[i, 2:(i - 1)] <- L[i, 2:(i - 1)] / sqrt(1 - cumsum(L[i, 1:(i - 2)]^2))
    }

    # Map to an unconstrained vector in R^(d choose 2)
    vH <- H[lower.tri(H)]
    vH <- atanh(vH)
    return(vH)
}

###############################################################################

#' Vector to correlation matrix
#'
#' @description Convert a vector to a correlation matrix.
#'
#' @param v A vector.
#'
#' @returns A correlation matrix.
#'
#' @export
vec2cor <- function(v) {
    d <- as.integer(round((1 + sqrt(1 + 8 * length(v))) / 2))

    # Map unconstrained vector to d x d matrix with entries in (-1, 1)
    H <- matrix(0, nrow = d, ncol = d)
    diag(H) <- 1
    H[lower.tri(H)] <- tanh(v)

    # Map back to Cholesky factor space
    for (i in 2:d) {
        H[i, 2:i] <- H[i, 2:i] * sqrt(cumprod(1 - H[i, 1:(i - 1)]^2))
    }

    # Return correlation matrix
    return(tcrossprod(H, H))
}

###############################################################################
