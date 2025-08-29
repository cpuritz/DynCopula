###############################################################################

#' Construct Empirical CDF
#'
#' @description Construct a function that evaluates the empirical CDF of the
#' input vector.
#'
#' @param x A numeric vector.
#'
#' @returns A function.
.empcdf <- function(x) {
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
#' @param cores Number of cores to use. Default is \code{1}.
#'
#' @returns A list containing
#' \itemize{
#'    \item \code{FX} Pseudo-observations
#'    \item \code{FXm} Left limits of pseudo-observations
#' }
#'
#' @export
pseudo_obs <- function(X, cores = 1L) {
    assertthat::assert_that(
        is.matrix(X) || is.data.frame(X) || methods::is(X, "Matrix"),
        is.numeric(cores) && cores >= 1L
    )
    cores <- as.integer(cores)

    # Construct empirical CDF functions
    pX <- apply(X, 2, .empcdf)

    # Parallel computation of pseudo-observations with progress bar
    cl <- parallel::makeCluster(cores)
    future::plan(future::cluster, workers = cl)
    on.exit({ future::plan(future::sequential); parallel::stopCluster(cl) },
            add = TRUE)
    progressr::with_progress({
        pbar <- progressr::progressor(along = pX)
        res <- future.apply::future_lapply(
            X = seq_along(pX),
            FUN = function(i) {
                y <- list(FX = pX[[i]](X[, i]),
                          FXm = pX[[i]](X[, i] - 1))
                pbar()
                return(y)
            },
            future.packages = "Matrix"
        )
    })
    FX <- do.call(cbind, lapply(res, '[[', "FX"))
    FXm <- do.call(cbind, lapply(res, '[[', "FXm"))

    return(list(FX = FX, FXm = FXm))
}

###############################################################################
