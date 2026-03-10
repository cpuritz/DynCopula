###############################################################################

#' Fit a dynamic Gaussian copula to scRNA-seq data
#'
#' @description Fit a dynamic Gaussian copula to scRNA-seq data using local
#' likelihood.
#'
#' @param sce A \code{SingleCellExperiment}.
#' @param formula Formula for covariates. Variable names should be column
#' metadata names.
#' @param lambda Vector of penalty parameters to test. Default is
#' \code{10^(seq(-5, 5, length.out = 7))}.
#' @param K Dimension of the spline basis matrix. Default is \code{30}.
#' @param nfold Number of folds for cross-validation. Default is \code{5}.
#' @param cores Number of cores to use. Default is \code{1}.
#' @param control A \code{list} of control parameters for optimization.
#'
#' @details The formula can include a single smooth covariate and any number of
#' linear covariates. Specify the smooth covariate by \code{s(t)}, replacing
#' \code{t} with the name of the smooth covariate. If interaction terms between
#' linear covariates and the smooth covariate are included, a baseline
#' (intercept) smooth function will be included in the model.
#'
#' If a smooth covariate is included in the formula, it is modeled
#' using penalized splines with a roughness penalty. The linear covariates are
#' penalized by an L2 penalty. The penalty parameters are selected via
#' cross-validation. A separate penalty parameter is used to for each term
#' which includes the smooth covariate. This includes a baseline smooth
#' function as well as any interactions between the smooth covariate and
#' linear covariates.
#'
#' Optimization is performed using L-BFGS. The \code{control} argument
#' is a list that supplies control parameters for optimization. The following
#' parameters can be supplied:
#' \itemize{
#'   \item \code{max_itr} Maximum number of iterations. Default is \code{100}.
#'   \item \code{history_size} History size. Default is \code{30}.
#'   \item \code{tolerance_grad} Termination tolerance for gradient. Default is
#'   \code{1e-7}.
#'   \item \code{tolerance_change} Termination tolerance for log-likelihood.
#'   Default is \code{1e-9}.
#'   \item \code{precision}: Floating point precision for calculations. Either
#'   \code{"float32"} or \code{"float64"}. Default is \code{"float64"}.
#'   \item \code{boundary}: Boundary extension for spline knot boundaries.
#'   Default is \code{0.05}.
#' }
#' Any control parameters not specified are replaced by their default values.
#'
#' @returns The same \code{SingleCellExperiment} as was passed as input, but
#' with the metadata entry \code{copula_fit} updated to include results of the
#' model fitting.
#'
#' @export
fit_gamgc_sce <- function(sce,
                          formula,
                          lambda = 10^(seq(-5, 5, length.out = 7)),
                          K = 30,
                          nfold = 5,
                          cores = 1,
                          control = list()) {
    # Check arguments not validated by fit_gamgc
    assert_that(methods::is(sce, "SingleCellExperiment"))

    if (!"copula_fit" %in% names(metadata(sce))) {
        stop("You must fit margins first using 'fit_margins'.")
    }

    copula_fit <- metadata(sce)$copula_fit
    assert_that("FX" %in% names(copula_fit))

    # Identify all covariates in the formula
    vars <- attributes(stats::terms(formula))$variable
    varnames <- as.character(vars)[2:length(vars)]
    s_terms <- unique(varnames[grepl("s(*)", varnames)])
    s_name <- sub("^s\\((.*)\\)$", "\\1", s_terms)
    covs <- unique(c(varnames[varnames != s_terms], s_name))

    # Design matrix
    design <- SummarizedExperiment::colData(sce)[, covs, drop = FALSE]
    design <- as.data.frame(design)

    # Estimate copula parameters
    res <- fit_gamgc(
        FX = copula_fit$FX,
        design = design,
        formula = formula,
        lambda = lambda,
        K = K,
        nfold = nfold,
        cores = cores,
        control = control
    )

    # Save results in metadata
    metadata(sce)$copula_fit <- c(metadata(sce)$copula_fit, res)

    return(sce)
}

###############################################################################
