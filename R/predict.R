###############################################################################

#' Predict method for GAM Gaussian copula models
#'
#' @description Predict coefficients for new data from a previously fit GAM
#' Gaussian copula model.
#'
#' @param object A fitted object of class \code{gamGaussianCopula}.
#' @param design A \code{data.frame} specifying the new design matrix. Must
#' include the same columns as in the design matrix in the original
#' \link[DynCopula]{fit_gamgc} call.
#' @param type The type of prediction required. The default option
#' \code{"response"} returns a \code{matrix} of correlation coefficients, with
#' rows corresponding to rows in \code{design}. The \code{"link"} option returns
#' a \code{matrix} of predictions of the coefficients on the link scale. The
#' \code{"response_matrix"} returns a list of correlation matrices, one per row
#' of \code{design}.
#' @param ... Additional arguments.
#'
#' @returns A matrix of coefficients.
#'
#' @examples
#' \dontrun{
#' library(DynCopula)
#' library(copula)
#' set.seed(0)
#'
#' N <- 1000
#' # Smooth covariate
#' t <- runif(N)
#' # Linear covariates
#' x1 <- sample(c("a", "b"), N, replace = TRUE)
#' x2 <- sample(seq(5), N, replace = TRUE)
#' # Design matrix
#' design <- data.frame(t = t, x1 = x1, x2 = x2)
#'
#' # Covariate-dependent correlation function
#' rho <- function(t, x1, x2) {
#'     0.7 * cos(4 * pi * t) + (x1 == "b") * 0.1 - 0.01 * x2
#' }
#' # Sample from copula
#' U <- t(sapply(seq(N), function(i) {
#'     rho_i <- rho(t[i], x1[i], x2[i])
#'     cop <- normalCopula(param = rho_i, dim = 2, dispstr = "un")
#'     return(rCopula(1L, cop))
#' }))
#'
#' # Fit model
#' gc <- fit_gamgc(U, design, formula = ~x1 + x2 + s(t))
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
                                      type = c("response", "link", "response_matrix"),
                                      ...) {
    assert_that(
        methods::is(object, "gamGaussianCopula"),
        is.data.frame(design),
        is.character(type)
    )

    type <- match.arg(type)
    N <- object$nobs

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

    d <- object$dim
    if (type == "link") {
        # Matrix of predicted coefficients on the link scale
        colnames(Hhat) <- paste0("eta", seq_len(choose(d, 2)))
        return(Hhat)
    } else if (type == "response_matrix") {
        # List of predicted correlation matrices
        Rhat <- lapply(seq_len(dim(Hhat)[1]), function(i) {
            R <- vec2cor(Hhat[i, ])
            rownames(R) <- colnames(R) <- object$colnames
            return(R)
        })
        return(Rhat)
    } else {
        # Matrix of predicted correlation coefficients
        Rhat <- t(apply(Hhat, 1, function(v) {
            copula::P2p(vec2cor(v))
        }))

        # Ensure consistent shape of Rhat across all dimensions
        if (d == 2L) {
            Rhat <- t(Rhat)
        }

        # Add column labels
        cnames <- object$colnames
        combs <- t(utils::combn(seq_len(d), 2))
        colnames(Rhat) <- apply(combs, 1, function(x) {
            paste(cnames[x[1]], cnames[x[2]], sep = '_')
        })

        return(Rhat)
    }
}

###############################################################################
