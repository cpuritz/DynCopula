###############################################################################

#' Predict method for GAM Gaussian copula models
#'
#' @description Predict coefficients for new data from a previously fit GAM
#' Gaussian copula model.
#'
#' @param object A fitted object of class \code{gamGaussianCopula}.
#' @param FX_new Matrix of new pseudo-observations.
#' @param design_new New design matrix. Rows correspond to rows in
#' \code{FX_new}. Must include the same columns as in the design matrix in the
#' original \link[DynCopula]{fit_dyn_gc} call.
#' @param type The type of prediction required. The default is correlation
#' coefficients (\code{"response"}). The \code{"link"} option returns
#' predictions of the calibration coefficients.
#' @param ... Additional arguments.
#'
#' @return A matrix of coefficients.
#'
#' @export
#'
#' @method predict gamGaussianCopula
predict.gamGaussianCopula <- function(object,
                                      FX_new,
                                      design_new,
                                      type = c("response", "link"),
                                      ...) {
    assert_that(
        methods::is(object, "gamGaussianCopula"),
        is.numeric(FX_new) && is.matrix(FX_new),
        is.data.frame(design_new),
        dim(FX_new)[1] == dim(design_new)[1],
        is.character(type)
    )

    type <- match.arg(type)

    design_disc <- design_new[, colnames(design_new) != "time", drop = FALSE]
    missing_vars <- setdiff(all.vars(object$formula), colnames(design_disc))
    if (length(missing_vars) > 0L) {
        stop("The following columns are missing in the design matrix: ",
             paste(missing_vars, collapse = ", "))
    }
    Z_new <- stats::model.matrix(object$formula, design_disc)

    # Construct matrix of estimated calibration coefficients
    if (!object$continuous) {
        Hhat <- Z_new %*% object$beta
    } else {
        # Extract new continuous covariate values
        x_new <- design_new[["time"]]

        # Standard scale continuous covariate on original scale
        min_x <- min(object$time)
        max_x <- max(object$time)
        x_new <- (x_new - min_x) / (max_x - min_x)

        # Construct basis matrix for new covariates
        B_attr <- attributes(object$B)
        B_new <- splines::ns(
            x = x_new,
            knots = B_attr$knots,
            Boundary.knots = B_attr$Boundary.knots,
            intercept = B_attr$intercept
        )

        # Matrix of estimated calibration coefficients
        beta <- object$beta
        K <- dim(object$B)[2]
        cts_coef <- beta[1:K, , drop = FALSE]
        disc_coef <- beta[(K + 1):dim(beta)[1], , drop = FALSE]
        Hhat <- B_new %*% cts_coef + Z_new %*% disc_coef
    }

    # Add column labels
    colnames(Hhat) <- paste0("eta", seq_len(choose(d, 2)))
    if (type == "link") {
        return(Hhat)
    }

    # Matrix of estimated correlation coefficients
    Rhat <- t(apply(Hhat, 1, function(v) {
        copula::P2p(vec2cor(v))
    }))

    # Ensure consistent shape of Rhat across all dimensions
    d <- dim(FX_new)[2]
    if (d == 2L) {
        Rhat <- t(Rhat)
    }

    # Add column labels
    ix_lab <- apply(utils::combn(seq_len(d), 2), 2, function(x) {
        paste(x, collapse = '_')
    })
    colnames(Rhat) <- paste0("rho", ix_lab)

    return(Rhat)
}

###############################################################################
