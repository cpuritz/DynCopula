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

    varnames <- as.character(vars)[-1]
    factors <- attributes(terms)$factors

    # Identify smoothing term (if any)
    s_term <- unique(varnames[grepl("^s\\(.*\\)$", varnames)])
    is_smooth <- (rownames(factors) == s_term)
    if (length(s_term) == 0L) {
        int_formula <- NULL
        s_name <- NULL
        # Linear covariates
        lin_factors <- colnames(factors)
    } else if (length(s_term) > 1L) {
        stop("You have specified multiple smooth covariates. Only one is ",
             "allowed.")
    } else {
        # Extract smooth function arguments
        cl_args <- as.list(str2lang(s_term))[-1]
        # First argument is the smooth covariate name
        s_name <- as.character(cl_args[[1]])
        # Interaction formula
        if ("by" %in% names(cl_args)) {
            int_formula <- stats::as.formula(paste("~", deparse(cl_args$by)))
        } else {
            int_formula <- ~1
        }
        # Linear covariates
        lin_factors <- colnames(factors)[colnames(factors) != s_term]
    }

    # Formula for linear covariates
    if (length(lin_factors) > 0L) {
        lin_terms <- paste("~", paste(lin_factors, collapse = "+"))
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
