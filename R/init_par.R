###############################################################################

#' Initial parameter guess
#'
#' @param x0 Time point to center around.
#' @param x Vector of time points.
#' @param NX Matrix of normal-transformed pseudo-observations.
#' @param h Kernel bandwidth.
#'
#' @returns Parameter vector
.init_par <- function(x0, x, NX, h) {
    # Use nearby points to estimate initial correlation matrix
    dx <- abs(x0 - x)
    # Ensure that at least d + 1 points are used to avoid a singular
    # correlation matrix
    d <- dim(NX)[2]
    thr <- max((d + 1) / length(x), h)
    NX_loc <- NX[which(abs(dx) <= thr), ]
    cor_loc <- stats::cor(NX_loc, method = "pearson")

    # Sometimes the Cholesky decomposition fails. This generally occurs if
    # cor_loc is numerically not positive definite, even though it theoretically
    # is. If this happens, fall back to the identity matrix.
    if (inherits(try(chol(cor_loc), silent = TRUE), "try-error")) {
        cor_loc <- diag(d)
    }

    # Convert to vector
    par0 <- cor2vec(cor_loc)
    return(par0)
}

###############################################################################
