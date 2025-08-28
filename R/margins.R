###############################################################################

#' Construct Empirical CDF
#'
#' @description Construct a function that evaluates the empirical CDF of the
#' input vector.
#'
#' @param x A numeric vector.
#'
#' @returns A function.
#'
#' @export
empcdf <- function(x) {
    assertthat::assert_that(
        is.numeric(x) && length(x) > 0
    )

    x <- sort(x)
    n <- length(x)
    vals <- unique(x)
    y <- cumsum(tabulate(match(x, vals))) / (n + 1)

    fun <- stats::approxfun(
        x = vals,
        y = y,
        method = "constant",
        yleft = 0,
        yright = y[length(y)],
        f = 0,
        ties = "ordered"
    )
    return(fun)
}

###############################################################################

#' Compute pseudo-observations
#'
#' @description Compute pseudo-observations and left-limits of
#' pseudo-observations for count-valued data.
#'
#' @param X A \code{matrix} or \code{data.frame}.
#'
#' @returns A list containing
#' \itemize{
#'    \item \code{FX} Pseudo-observations
#'    \item \code{FXm} Left limits of pseudo-observations
#' }
#'
#' @export
pseudo_obs <- function(X) {
    assertthat::assert_that(
        is.matrix(X) || is.data.frame(X)
    )
    pX <- apply(X, 2, empcdf)
    FX <- do.call(cbind, lapply(seq_along(pX), function(i) {
        pX[[i]](X[, i])
    }))
    FXm <- do.call(cbind, lapply(seq_along(pX), function(i) {
        pX[[i]](X[, i] - 1)
    }))
    return(list(FX = FX, FXm = FXm))
}

###############################################################################
