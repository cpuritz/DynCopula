#' Parse a formula into discrete and smooth parts
#'
#' @description Parse a formula into discrete and smooth parts. Used internally.
#'
#' @param formula A \code{formula} object.
#'
#' @returns A list with the following components:
#' \itemize{
#'   \item disc Formula for discrete covariates.
#'   \item int Formula for interactions with smooth covariate.
#' }
.parse_formula <- function(formula) {
    terms <- stats::terms(formula)
    vars <- attributes(terms)$variable
    varnames <- as.character(vars)[2:length(vars)]
    factors <- attributes(terms)$factors

    # Identify smoothing terms
    s_terms <- unique(varnames[grepl("s(*)", varnames)])
    is_smooth <- (rownames(factors) == s_terms)
    if (length(s_terms) == 0L) {
        int_formula <- NULL
        disc_factors <- factors
    } else if (length(s_terms) > 1L) {
        stop("You have specified multiple smooth covariates. Only one is ",
             "allowed.")
    } else {
        s_vars <- sub("^s\\((.*)\\)$", "\\1", s_terms)
        if (s_vars %in% varnames) {
            stop("The smooth covariate cannot be included as a free term.")
        }

        # Factors that represent interaction terms
        int_factors <- factors[, grep(":", colnames(factors)), drop = FALSE]
        if (dim(int_factors)[2] == 0L) {
            # There are no interaction terms at all
            int_formula <- ~1
        } else {
            # Interactions that involve the smooth covariate
            which_smooth <- int_factors[is_smooth, ] > 0
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

        # Terms for discrete covariates
        disc_factors <- factors[!is_smooth, factors[s_terms, ] == 0, drop = FALSE]
    }

    # Formula for discrete covariates
    disc_terms <- paste("~", paste(colnames(disc_factors), collapse = "+"))
    disc_formula <- stats::as.formula(disc_terms)

    return(list(disc = disc_formula, int = int_formula))
}
