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
.pseudo_obs <- function(X, cores = 1L) {
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

#' Fit margins
#'
#' @description Fit parametric models for each gene. Multiple families and
#' formulas can be passed, in which case the best fitting model is chosen across
#' all combinations.
#'
#' @param sce A \code{SingleCellExperiment}.
#' @param family Family to use. Either \code{"NBI"} (negative binomial) or
#' \code{"ZINBI"} (zero-inflated negative binomial).
#' @param mu_formula A vector of strings specifying formulas for the mean
#' parameter.
#' @param sigma_formula A vector of strings specifying formulas for the
#' dispersion parameter.
#'
#' @returns The same \code{SingleCellExperiment} as was passed as input, but
#' modified to include a named metadata entry \code{margins} which is a list
#' of marginal models.
#'
#' @details Formulas should be written as a function of the name of the
#' pseudotime column (what was passed as the argument \code{tcol} to
#' \link[DynCopula]{setup}). Only the righthand side of the formula should be
#' passed. For example (if the pseudotime column is named \code{"pseudotime"}):
#' \itemize{
#'   \item \code{mu_formula = "1"}: mean is constant.
#'   \item \code{mu_formula = "pseudotime"}: mean varies linearly with
#'   pseudotime.
#'   \item \code{mu_formula = "gamlss::pb(pseudotime)"}: mean varies through a
#'   P-spline with pseudotime.
#' }
#' For ZINBI models, the zero probability parameter is assumed to be constant.
#'
#' Model selection is done using AIC.
#'
#' @export
fit_margins <- function(sce,
                        family = c("NBI", "ZINBI"),
                        mu_formula,
                        sigma_formula) {
    assertthat::assert_that(
        methods::is(sce, "SingleCellExperiment"),
        "dyn_corr_info" %in% names(metadata(sce)),
        is.character(mu_formula),
        is.character(sigma_formula)
    )
    family <- match.arg(family, several.ok = TRUE)

    message("Fitting marginal distributions")

    info <- metadata(sce)$dyn_corr_info
    X <- SummarizedExperiment::assay(sce, info$assay)
    X <- Matrix::t(X[info$features, ])
    pseudotimes <- SummarizedExperiment::colData(sce)[[info$tcol]]

    cl <- parallel::makeCluster(info$cores)
    future::plan(future::cluster, workers = cl)
    on.exit({ future::plan(future::sequential); parallel::stopCluster(cl) },
            add = TRUE)

    combs <- expand.grid(
        family = family,
        mu_formula = mu_formula,
        sigma_formula = sigma_formula,
        stringsAsFactors = FALSE,
        KEEP.OUT.ATTRS = FALSE
    )

    progressr::with_progress({
        d <- dim(X)[2]
        pbar <- progressr::progressor(along = seq_len(d))

        X_slices <- apply(X, 2, c, simplify = FALSE)
        margins <- future.apply::future_lapply(
            X = X_slices,
            FUN = function(x) {
                ddata <- data.frame(x, pseudotimes)
                names(ddata) <- c("x", info$tcol)
                models <- lapply(seq_len(dim(combs)[1]), function(i) {
                    fmu <- paste("x ~", combs[i, "mu_formula"])
                    fsigma <- paste("~", combs[i, "sigma_formula"])
                    fname <- combs[i, "family"]
                    fam <- utils::getFromNamespace(fname, "gamlss.dist")
                    fit <- gamlss::gamlss(
                        formula = stats::formula(fmu),
                        sigma.formula = stats::formula(fsigma),
                        nu.formula = stats::formula("~ 1"),
                        data = ddata,
                        family = fam,
                        control = gamlss::gamlss.control(trace = FALSE)
                    )
                    fit$call$family <- as.name(fname)
                    return(fit)
                })
                best_ix <- which.min(sapply(models, '[[', "aic"))
                pbar()
                return(models[[best_ix]])
            },
            future.globals = c("combs", "pbar", "pseudotimes"),
            future.seed = TRUE,
            future.packages = c("gamlss.dist")
        )
    })

    names(margins) <- colnames(X)
    metadata(sce)$margins <- margins
    return(sce)
}

###############################################################################
