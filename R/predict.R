###############################################################################

#' Predict method for GAM Gaussian copula models
#'
#' @description Predict coefficients for new data from a previously fit GAM
#' Gaussian copula model.
#'
#' @param object A fitted object of class \code{gamGaussianCopula}.
#' @param design A \code{data.frame} specifying the new design matrix. Must
#' include the same columns as in the design matrix in the original
#' \link[gamgc]{fit_gamgc} call.
#' @param type The type of prediction to return The default option
#' \code{"response"} returns a \code{matrix} of correlation coefficients, with
#' rows corresponding to rows in \code{design} and columns corresponding to
#' pairwise coefficients in row-major order. \code{"link"} returns a
#' \code{matrix} of coefficients on the link scale. \code{"response_matrix"}
#' returns a list of pairwise correlation matrices, one per row of
#' \code{design}. \code{"tau"} returns a \code{matrix} of Kendall's tau
#' coefficients with rows corresponding to rows in \code{design} and columns
#' corresponding to pairwise coefficients in row-major order.
#' \code{"tau_matrix"} returns a list of pairwise Kendall's tau matrices, one
#' per row of \code{design}.
#' @param ... Additional arguments.
#'
#' @returns A matrix or list of matrices.
#'
#' @examples
#' \dontrun{
#' library(gamgc)
#' library(copula)
#' set.seed(0)
#'
#' N <- 500
#' # Smooth covariate
#' t <- runif(N)
#' # Design matrix
#' design <- data.frame(t = t)
#'
#' # Covariate-dependent correlation function
#' rho <- function(t) {
#'     0.7 * cos(4 * pi * t)
#' }
#' # Sample from copula
#' U <- t(sapply(seq(N), function(i) {
#'     rho_i <- rho(t[i])
#'     cop <- normalCopula(param = rho_i, dim = 2, dispstr = "un")
#'     return(rCopula(1L, cop))
#' }))
#'
#' # Fit model
#' gc <- fit_gamgc(U, design, formula = ~s(t))
#'
#' # Predict correlation coefficients for original data
#' pred <- predict(gc, design)
#' }
#'
#' @export
#'
#' @method predict gamGaussianCopula
predict.gamGaussianCopula <- function(object,
                                      design,
                                      type = c("response", "link",
                                               "response_matrix", "tau",
                                               "tau_matrix"),
                                      ...) {
    assert_that(
        methods::is(object, "gamGaussianCopula"),
        is.data.frame(design),
        is.character(type)
    )

    type <- match.arg(type)
    N <- nrow(design)

    # Make sure new design matrix has necessary columns
    all_vars <- c(all.vars(object$lin_formula), all.vars(object$int_formula))
    missing_vars <- setdiff(all_vars, colnames(design))
    if (length(missing_vars) > 0L) {
        stop("The following columns are missing in the design matrix: ",
             paste(missing_vars, collapse = ", "))
    }

    # Pad design matrices in case new design matrices are missing factors
    # present in the original design matrix
    Z_new <- stats::model.matrix(object$lin_formula, design)
    Z_col_miss <- setdiff(object$Z_names, colnames(Z_new))
    if (length(Z_col_miss) > 0L) {
        nzn <- dim(Z_new)[2]
        Z_new <- cbind(Z_new, matrix(0, nrow = N, ncol = length(Z_col_miss)))
        colnames(Z_new)[(nzn + 1):(nzn + length(Z_col_miss))] <- Z_col_miss
        Z_new <- Z_new[, object$Z_names]
    }

    smooth_name <- object$smooth_name
    if (!is.null(smooth_name)) {
        M_new <- stats::model.matrix(object$int_formula, design)
        M_col_miss <- setdiff(object$M_names, colnames(M_new))
        if (length(M_col_miss) > 0L) {
            nmn <- dim(M_new)[2]
            M_new <- cbind(M_new, matrix(0, nrow = N, ncol = length(M_col_miss)))
            colnames(M_new)[(nmn + 1):(nmn + length(M_col_miss))] <- M_col_miss
            M_new <- M_new[, object$M_names]
        }
    }

    # Construct matrix of estimated calibration coefficients
    if (is.null(smooth_name)) {
        Hhat <- Z_new %*% object$beta
    } else {
        # Extract new smooth covariate values
        x_new <- design[[smooth_name]]

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

        betas <- array(beta[(L2 + 1):(L2 + K * L1), ], dim = c(K, p, L1))
        betas <- aperm(betas, c(3, 1, 2))  # (L1, K, p)
        Hhat <- matrix(0, N, p)
        for (j in 1:L1) {
            Hhat <- Hhat + M_new[, j] * (B_new %*% betas[j, , ])
        }
        Hhat <- Hhat + Z_new %*% beta[1:L2, ]
    }

    par2tau <- function(x) {
        2 / pi * asin(x)
    }

    d <- object$dim
    if (type == "link") {
        # Matrix of predicted coefficients on the link scale
        colnames(Hhat) <- paste0("eta", seq_len(choose(d, 2)))
        return(Hhat)
    } else if (type == "response_matrix" || type == "tau_matrix") {
        # List of predicted correlation matrices
        Rhat_list <- lapply(seq_len(dim(Hhat)[1]), function(i) {
            R <- vec2cor(Hhat[i, ])
            rownames(R) <- colnames(R) <- object$colnames
            return(R)
        })
        if (type == "response_matrix") {
            return(Rhat_list)
        } else {
            return(lapply(Rhat_list, par2tau))
        }
    } else {
        # Matrix of predicted correlation coefficients
        Rhat <- t(apply(Hhat, 1, function(v) {
            copula::P2p(vec2cor(v))
        }))

        # Ensure consistent shape of Rhat across all dimensions
        if (d == 2L) {
            Rhat <- t(Rhat)
        }

        if (type == "tau") {
            tau_hat <- par2tau(Rhat)

            # Add column labels
            cnames <- object$colnames
            if (is.null(cnames)) {
                cnames <- seq_len(object$dim)
            }
            combs <- t(utils::combn(seq_len(d), 2))
            colnames(tau_hat) <- apply(combs, 1, function(x) {
                paste("tau", cnames[x[1]], cnames[x[2]], sep = '_')
            })

            return(tau_hat)
        } else {
            # Add column labels
            cnames <- object$colnames
            if (is.null(cnames)) {
                cnames <- seq_len(object$dim)
            }
            combs <- t(utils::combn(seq_len(d), 2))
            colnames(Rhat) <- apply(combs, 1, function(x) {
                paste("rho", cnames[x[1]], cnames[x[2]], sep = '_')
            })

            return(Rhat)
        }
    }
}

###############################################################################
