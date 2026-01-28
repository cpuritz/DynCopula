###############################################################################

#' Fit a dynamic Gaussian copula to scRNA-seq data
#'
#' @description Fit a dynamic Gaussian copula to scRNA-seq data using local
#' likelihood.
#'
#' @param sce A \code{SingleCellExperiment}.
#' @param bandwidths Vector of global kernel bandwidths to test. Default is
#' \code{seq(0.01, 0.20, 0.01)}.
#' @param variable Whether to use a variable bandwidth. Default is \code{TRUE}.
#' @param alpha Powers to test for the variable bandwidth function. Ignored if
#' \code{variable = FALSE}. Default is \code{seq(0, 1, 0.1)}.
#' @param beta Scale factors to test for the variable bandwidth function.
#' Ignored if \code{variable = FALSE}. Default is \code{seq(0.75, 1.25, 0.1)}.
#' @param control A \code{list} of control parameters for optimization.
#' @param ncv Number of pseudotime values to use for LOOCV. If \code{NULL} (the
#' default), all pseudotime values are used.
#' @param degree Degree of local polynomial approximation. Default is \code{0}.
#'
#' @details This function can only be run after
#' \link[DynCopula]{generate_metacells} has been run.
#'
#' Optimization is performed using L-BFGS. The \code{control} argument
#' is a list that supplies control parameters for optimization. The following
#' parameters can be supplied:
#' \itemize{
#'   \item \code{max_outer}: Maximum number of outer iterations. Default is
#'   \code{1}.
#'   \item \code{max_itr}: Maximum number of inner iterations. Default is
#'   \code{100}.
#'   \item \code{history_size}: History size. Default is \code{30}.
#'   \item \code{tolerance_grad}: Termination tolerance for gradient. Default is
#'   \code{1e-7}.
#'   \item \code{tolerance_change}: Termination tolerance for log-likelihood.
#'   Default is \code{1e-9}.
#' }
#' Any control parameters not specified are assigned their default values.
#'
#' Kernel bandwidth selection is performed using leave-one-out cross-validation.
#' If \code{variable = FALSE}, a single global bandwidth is selected from
#' \code{bandwidths}. The optimal bandwidth is the one that maximizes the
#' cross-validated likelihood criterion. The model returned uses this optimal
#' bandwidth.
#'
#' If \code{variable = TRUE}, a variable bandwidth is selected. A global pilot
#' bandwidth is first selected as discussed above. LOOCV is then used to select
#' parameters for the variable bandwidth function from \code{alpha} and
#' \code{beta}. The variable bandwidth function is defined as
#' \deqn{
#' h(x;\alpha,\beta)=\beta h_{0}\big(\hat{f}_{x}(x)/G\big)^{-\alpha}
#' }
#' where \eqn{h_{0}} is the global pilot bandwidth, \eqn{\hat{f}_{x}} is a
#' kernel density estimator for the density of pseudotimes, and \eqn{G} is the
#' geometric mean of \eqn{\hat{f}_{x}(x)}. As long as \eqn{\alpha=0} and
#' \eqn{\beta=1} are included in the parameter sets, the global bandwidth model
#' is included, and thus the variable approach will perform no worse (in terms
#' of cross-validated likelihood) than the global bandwidth model.
#'
#' The argument \code{ncv} specifies the number of covariate values to use for
#' LOOCV. Full LOOCV corresponds to \code{ncv = length(x)}. If
#' \code{ncv < length(x)}, then only a subset of the covariates are used to
#' reduce the run time. The indices of the chosen covariates are equally spaced
#' along \code{seq_along(x)}. This is only an estimate to full LOOCV and may be
#' quite inaccurate for \code{ncv << length(x)}.
#'
#' @returns The same \code{SingleCellExperiment} as was passed as input, but
#' with the metadata entry \code{dyn_corr} updated to include the following
#' elements:
#' \itemize{
#'   \item \code{rho}: Matrix of estimated copula parameters.
#'   \item \code{bandwidth}: The kernel bandwidth used to estimate copula
#'   parameters.
#' }
#'
#' @export
fit_dyn_corr <- function(sce,
                         bandwidths = seq(0.01, 0.20, 0.01),
                         variable = TRUE,
                         alpha = seq(0, 1, 0.1),
                         beta = seq(0.75, 1.25, 0.1),
                         ncv = NULL,
                         degree = 0,
                         control = list()) {
    assert_that(
        methods::is(sce, "SingleCellExperiment"),
        "dyn_corr" %in% names(metadata(sce)),
        "metacell_sce" %in% names(metadata(sce)$dyn_corr),
        is.null(ncv) || (is.numeric(ncv) && ncv > 1)
    )

    sce_mc <- metadata(sce)$dyn_corr$metacell_sce
    dyn_corr <- metadata(sce_mc)$dyn_corr
    assert_that("margins" %in% names(dyn_corr))

    pseudotimes <- sce_mc[[dyn_corr$time_col]]

    # Construct jittered pseudo-observations
    FX <- dyn_corr$FXm + (dyn_corr$FX - dyn_corr$FXm) * dyn_corr$V

    if (is.null(ncv)) {
        ncv <- length(pseudotimes)
    } else {
        ncv <- as.integer(ncv)
        if (ncv > length(pseudotimes)) {
            ncv <- length(pseudotimes)
            message("Only ", ncv, " unique pseudotimes are available. ",
                    "Using ncv = ", ncv, ".")
        }
    }

    res <- fit_local_gaussian(
        FX = FX,
        x = pseudotimes,
        bandwidths = bandwidths,
        variable = variable,
        alpha = alpha,
        beta = beta,
        ncv = ncv,
        degree = degree,
        control = control,
        cores = dyn_corr$cores
    )

    # Convert numeric labels to gene names
    gene_names <- sapply(colnames(res$rho), function(x) {
        x2 <- unlist(strsplit(x, split = "rho"))[2]
        ix <- as.numeric(unlist(strsplit(x2, split = '_')))
        return(paste(dyn_corr$features[ix], collapse = '_'))
    })
    colnames(res$rho) <- gene_names

    # Save results in metadata
    metadata(sce)$dyn_corr$eta <- res$eta
    metadata(sce)$dyn_corr$rho <- res$rho
    metadata(sce)$dyn_corr$bandwidths <- res$bandwidths

    return(sce)
}

###############################################################################
