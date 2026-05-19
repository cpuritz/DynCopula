###############################################################################

#' Parse a formula into linear and smooth parts
#'
#' @description Parse a formula into linear and smooth parts. Used internally.
#'
#' @param formula A \code{formula} object.
#'
#' @returns A list with the following components:
#' \itemize{
#'   \item linear Formula for linear covariates.
#'   \item smooth_linear_int Formula for interactions between linear covariates
#'   and the smooth covariate.
#'   \item smooth_name The name of the smooth covariate.
#' }
.parse_formula <- function(formula) {
    terms <- stats::terms(formula)
    vars <- attributes(terms)$variable

    # Check if the formula has no covariates (i.e. ~1)
    if (length(vars) == 1L) {
        return(list(
            linear = ~1,
            smooth_linear_int = NULL,
            smooth_name = NULL
        ))
    }

    varnames <- as.character(vars)[2:length(vars)]
    factors <- attributes(terms)$factors

    # Identify smoothing terms
    s_terms <- unique(varnames[grepl("s(*)", varnames)])
    is_smooth <- (rownames(factors) == s_terms)
    if (length(s_terms) == 0L) {
        int_formula <- NULL
        lin_factors <- factors
        s_name <- NULL
    } else if (length(s_terms) > 1L) {
        stop("You have specified multiple smooth covariates. Only one is ",
             "allowed.")
    } else {
        # Name of the smooth covariate
        s_name <- sub("^s\\((.*)\\)$", "\\1", s_terms)

        # Factors that represent interaction terms
        int_factors <- factors[, grep(":", colnames(factors)), drop = FALSE]
        if (dim(int_factors)[2] == 0L) {
            # There are no interaction terms at all
            int_formula <- ~1
        } else {
            # Interactions that involve the smooth covariate
            which_smooth <- (int_factors[is_smooth, ] > 0)
            s_int <- int_factors[!is_smooth, which_smooth, drop = FALSE]
            if (dim(s_int)[2] == 0L) {
                # There are interaction terms, but none involving the smooth
                # covariate
                int_formula <- ~1
            } else {
                # There are interaction terms involving the smooth covariate
                sd_int <- apply(s_int, 2, function(x) {
                    paste(rownames(s_int)[x > 0], collapse = ":")
                })
                int_terms <- paste("~", paste(sd_int, collapse = "+"))
                int_formula <- stats::as.formula(int_terms)
            }
        }

        # Terms for linear covariates
        lin_factors <- factors[!is_smooth, factors[s_terms, ] == 0, drop = FALSE]
    }

    # Formula for linear covariates
    if (dim(lin_factors)[2] > 0L) {
        lin_terms <- paste("~", paste(colnames(lin_factors), collapse = "+"))
        lin_formula <- stats::as.formula(lin_terms)
    } else {
        lin_formula <- ~1
    }

    return(list(
        linear = lin_formula,
        smooth_linear_int = int_formula,
        smooth_name = s_name
    ))
}

###############################################################################
