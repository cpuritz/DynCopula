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
#' original \link[DynCopula]{fit_gamgc} call.
#' @param type The type of prediction required. The default is correlation
#' coefficients (\code{"response"}). The \code{"link"} option returns
#' predictions of the calibration coefficients.
#' @param ... Additional arguments.
#'
#' @returns A matrix of coefficients.
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

    all_vars <- c(all.vars(object$disc_formula), all.vars(object$int_formula))
    design_disc <- design_new[, colnames(design_new) != "time", drop = FALSE]
    missing_vars <- setdiff(all_vars, colnames(design_disc))
    if (length(missing_vars) > 0L) {
        stop("The following columns are missing in the design matrix: ",
             paste(missing_vars, collapse = ", "))
    }
    Z_new <- stats::model.matrix(object$disc_formula, design_disc)

    # Construct matrix of estimated calibration coefficients
    if (!object$smooth) {
        Hhat <- Z_new %*% object$beta
    } else {
        M_new <- stats::model.matrix(object$int_formula, design_disc)

        # Extract new smooth covariate values
        x_new <- design_new[["time"]]

        # Standard scale smooth covariate on original scale
        min_x <- min(object$time)
        max_x <- max(object$time)
        x_new <- (x_new - min_x) / (max_x - min_x)

        # Construct basis matrix for new covariates
        B <- object$B
        B_attr <- attributes(B)
        B_new <- splines::ns(
            x = x_new,
            knots = B_attr$knots,
            Boundary.knots = B_attr$Boundary.knots,
            intercept = B_attr$intercept
        )

        beta <- object$beta
        K <- dim(B)[2]
        L1 <- dim(M_new)[2]
        L2 <- dim(Z_new)[2]
        p <- dim(beta)[2]
        N <- dim(FX_new)[1]

        betas <- array(beta[(L2 + 1):(L2 + K * L1), ], dim = c(K, p, L1))
        betas <- aperm(betas, c(3, 1, 2))  # (L1, K, p)
        Hhat <- matrix(0, N, p)
        for (j in 1:L1) {
            Hhat <- Hhat + M_new[, j] * (B_new %*% betas[j, , ])
        }
        Hhat <- Hhat + Z_new %*% beta[1:L2, ]
    }

    # Add column labels
    d <- dim(FX_new)[2]
    colnames(Hhat) <- paste0("eta", seq_len(choose(d, 2)))
    if (type == "link") {
        return(Hhat)
    }

    # Matrix of estimated correlation coefficients
    Rhat <- t(apply(Hhat, 1, function(v) {
        copula::P2p(vec2cor(v))
    }))

    # Ensure consistent shape of Rhat across all dimensions
    if (d == 2L) {
        Rhat <- t(Rhat)
    }

    # Add column labels
    colnames(Rhat) <- apply(utils::combn(seq_len(d), 2), 2, function(x) {
        paste0("rho", paste(x, collapse = '_'))
    })

    return(Rhat)
}

###############################################################################
