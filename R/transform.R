###############################################################################

#' cor2vec
#'
#' cor2vec
#'
#' @param R R
#' @param scale scale
#'
#' @returns value
#'
#' @export
cor2vec <- function(R, scale) {
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
    vH <- scale * atanh(vH)
    return(vH)
}

###############################################################################

#' vec2cor
#'
#' vec2cor
#'
#' @param v v
#' @param scale scale
#'
#' @returns value
#'
#' @export
vec2cor <- function(v, scale) {
    d <- as.integer(round((1 + sqrt(1 + 8 * length(v))) / 2))

    # Map unconstrained vector to d x d matrix with entries in (-1, 1)
    H <- matrix(0, nrow = d, ncol = d)
    diag(H) <- 1
    H[lower.tri(H)] <- tanh(v / scale)

    # Map back to Cholesky factor space
    for (i in 2:d) {
        H[i, 2:i] <- H[i, 2:i] * sqrt(cumprod(1 - H[i, 1:(i - 1)]^2))
    }

    # Return correlation matrix
    return(tcrossprod(H, H))
}

###############################################################################
