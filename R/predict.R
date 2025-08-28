###############################################################################

#' Estimate correlation coefficients
#'
#' @description Estimate correlation coefficients using smoothing splines fit
#' in the unconstrained space.
#'
#' @param res Output of \code{fit_dynamic_gaussian}.
#' @param x Times at which to estimate coefficients. Default is \code{res$x}.
#' @param df Number of degrees of freedom for spline fitting. Default is
#' \code{NULL}. If \code{NULL}, generalized cross validation is used instead.
#'
#' @details If \code{df} is not set to \code{length(res$x0)}, the splines will
#' not necessarily pass through the coefficient values originally predicted.
#'
#' @return A list containing
#' \itemize{
#'    \item \code{eta} A matrix of predicted unconstrained coefficients.
#'    \item \code{rho} A matrix of predicted correlation coefficients.
#' }
#'
#' @export
predict_rho <- function(res,
                        x = res$x,
                        df = NULL) {
    assertthat::assert_that(length(res$x0) >= 4)

    # Predict eta values using smooth splines
    eta_pred <- apply(res$eta, 2, function(y) {
        if (is.null(df)) {
            fit <- stats::smooth.spline(res$x0, y)
        } else {
            fit <- stats::smooth.spline(res$x0, y, df = df)
        }
        return(stats::predict(fit, x)$y)
    })

    # Convert to correlation matrices
    R_pred <- apply(eta_pred, 1, vec2cor, simplify = FALSE)

    # Extract time series for each coefficient
    ix <- which(lower.tri(R_pred[[1]]), arr.ind = TRUE)
    rho_pred <- apply(ix, 1, function(v) {
        sapply(R_pred, function(R) { R[v[1], v[2]] })
    })
    return(list(eta = eta_pred, rho = rho_pred))
}

###############################################################################
