###############################################################################

#' Initial estimate of calibration coefficients
#'
#' @param x0 Time point to center around.
#' @param x Vector of time points.
#' @param NX Matrix of normal-transformed pseudo-observations.
#' @param h Kernel bandwidth.
#'
#' @returns Vector of calibration coefficients.
.init_par <- function(x0, x, NX, h) {
    d <- dim(NX)[2]
    dx <- abs(x0 - x)

    # Use points within smoothing window to estimate initial correlation matrix.
    # Ensure that at least max(5, d + 1) points are used to avoid a singular
    # correlation matrix (the d + 1) and to ensure that the estimate is not just
    # noise (the 5).
    min_pt <- min(max(5, d + 1), length(dx))
    thr <- max(quantile(dx, min_pt / length(dx)), h)
    NX_loc <- NX[dx <= thr, ]
    cor_loc <- stats::cor(NX_loc, method = "pearson")

    # Avoid near perfect initial correlations
    rho_max <- 0.95
    cor_loc_vec <- copula::P2p(cor_loc)
    is_large <- abs(cor_loc_vec) > rho_max
    cor_loc_vec[is_large] <- sign(cor_loc_vec[is_large]) * rho_max
    cor_loc <- copula::p2P(cor_loc_vec)

    # Ensure final matrix is a valid correlation matrix. If nearPD does not
    # converge (which throws a warning), fall back to identity matrix.
    cor_loc <- tryCatch({
        Matrix::nearPD(cor_loc, corr = TRUE, base.matrix = TRUE, maxit = 200)$mat
    }, warning = function(w) {
        return(diag(d))
    })

    # Convert to vector
    par0 <- cor2vec(cor_loc)
    return(par0)
}

###############################################################################
