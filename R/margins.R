###############################################################################

#' Fit marginal distributions to each gene
#'
#' @description Fit parametric models to transcript counts for each gene.
#'
#' @param sce A \code{SingleCellExperiment}.
#' @param family A character vector specifying families to use. Options are
#' \code{"NBI"} (negative binomial) and \code{"ZINBI"}
#' (zero-inflated negative binomial).
#' @param mu_formula A character vector specifying formulas for the mean
#' parameter.
#' @param sigma_formula A character vector specifying formulas for the
#' dispersion parameter.
#' @param save Whether to save all information about margins. Default is
#' \code{FALSE}.
#'
#' @returns The same \code{SingleCellExperiment} as was passed as input, but
#' with the metadata entry \code{dyn_corr} updated to include the following
#' entries:
#' \itemize{
#'   \item \code{margins} List of marginal models.
#'   \item \code{FX} Matrix of pseudo-observations.
#'   \item \code{FXm} Matrix of left limits of pseudo-observations.
#'   \item \code{V} Jittering matrix.
#' }
#'
#' @details Multiple families and formulas can be passed, in which case the best
#' fitting model is chosen across all combinations by minimizing AIC.
#'
#' Formulas should be written as a function of the name of the
#' pseudotime column (what was passed as the argument \code{time_col} to
#' \link[DynCopula]{setup}). Only the righthand side of the formula should be
#' passed. For example, if the pseudotime column is named \code{"pseudotime"}:
#' \itemize{
#'   \item \code{mu_formula = "1"}: mean is constant.
#'   \item \code{mu_formula = "pseudotime"}: mean varies linearly with
#'   pseudotime.
#'   \item \code{mu_formula = "gamlss::pb(pseudotime)"}: mean varies through a
#'   P-spline with pseudotime. Any spline functions need to be explicitly
#'   scoped.
#' }
#' For ZINBI models, the zero probability parameter is assumed to be constant.
#'
#' If \code{save = TRUE}, all information returned by \link[gamlss]{gamlss} is
#' retained. This will greatly increase the size of the returned
#' \code{SingleCellExperiment}, and thus \code{save = FALSE} is the default. See
#' \link[gamlss]{gamlss.control} for more details.
#'
#' @export
fit_margins <- function(sce,
                        family = c("NBI", "ZINBI"),
                        mu_formula,
                        sigma_formula,
                        save = FALSE) {
    family <- match.arg(family, several.ok = TRUE)
    sce <- .fit_margins(
        sce = sce,
        family = family,
        mu_formula = mu_formula,
        sigma_formula = sigma_formula,
        save = save
    )
    return(sce)
}

###############################################################################

#' Fit marginal distributions to each gene for metacells
#'
#' @description Fit parametric models to metacell transcript counts for each
#' gene.
#'
#' @param sce A \code{SingleCellExperiment}.
#' @param family A character vector specifying families to use. Options are
#' \code{"NBI"} (negative binomial) and \code{"ZINBI"}
#' (zero-inflated negative binomial).
#' @param mu_formula A character vector specifying formulas for the mean
#' parameter.
#' @param sigma_formula A character vector specifying formulas for the
#' dispersion parameter.
#' @param save Whether to save all information about margins. Default is
#' \code{FALSE}.
#'
#' @returns The same \code{SingleCellExperiment} as was passed as input, but
#' with the embedded metacell \code{SingleCellExperiment} updated.
#'
#' @details Multiple families and formulas can be passed, in which case the best
#' fitting model is chosen across all combinations by minimizing AIC.
#'
#' Formulas should be written as a function of the name of the
#' pseudotime column (what was passed as the argument \code{time_col} to
#' \link[DynCopula]{setup}). Only the righthand side of the formula should be
#' passed. For example, if the pseudotime column is named \code{"pseudotime"}:
#' \itemize{
#'   \item \code{mu_formula = "1"}: mean is constant.
#'   \item \code{mu_formula = "pseudotime"}: mean varies linearly with
#'   pseudotime.
#'   \item \code{mu_formula = "gamlss::pb(pseudotime)"}: mean varies through a
#'   P-spline with pseudotime. Any spline functions need to be explicitly
#'   scoped.
#' }
#' For ZINBI models, the zero probability parameter is assumed to be constant.
#'
#' If \code{save = TRUE}, all information returned by \link[gamlss]{gamlss} is
#' retained. This will greatly increase the size of the returned
#' \code{SingleCellExperiment}, and thus \code{save = FALSE} is the default. See
#' \link[gamlss]{gamlss.control} for more details.
#'
#' @export
fit_metacell_margins <- function(sce,
                                 family = c("NBI", "ZINBI"),
                                 mu_formula,
                                 sigma_formula,
                                 save = FALSE) {
    assert_that("metacell_sce" %in% names(metadata(sce)$dyn_corr))
    family <- match.arg(family, several.ok = TRUE)

    mc_sce <- metadata(sce)$dyn_corr$metacell_sce
    mc_sce <- .fit_margins(
        sce = mc_sce,
        family = family,
        mu_formula = mu_formula,
        sigma_formula = sigma_formula,
        save = save
    )

    # Save the metacell sce in the original sce
    metadata(sce)$dyn_corr$metacell_sce <- mc_sce
    return(sce)
}

###############################################################################

#' Fit marginal distributions to each gene
#'
#' @description Fit parametric models to transcript counts for each gene.
#' Internal function.
#'
#' @param sce A \code{SingleCellExperiment}.
#' @param family A character vector specifying families to use. Options are
#' \code{"NBI"} (negative binomial) and \code{"ZINBI"}
#' (zero-inflated negative binomial).
#' @param mu_formula A character vector specifying formulas for the mean
#' parameter.
#' @param sigma_formula A character vector specifying formulas for the
#' dispersion parameter.
#' @param save Whether to save all information about margins. Default is
#' \code{FALSE}.
#'
#' @returns The same \code{SingleCellExperiment} as was passed as input, but
#' with the metadata entry \code{dyn_corr} updated to include the following
#' entries:
#' \itemize{
#'   \item \code{margins} List of marginal models.
#'   \item \code{FX} Matrix of pseudo-observations.
#'   \item \code{FXm} Matrix of left limits of pseudo-observations.
#'   \item \code{V} Jittering matrix.
#' }
.fit_margins <- function(sce,
                         family,
                         mu_formula,
                         sigma_formula,
                         save) {
    assert_that(
        methods::is(sce, "SingleCellExperiment"),
        "dyn_corr" %in% names(metadata(sce)),
        is.character(mu_formula),
        is.character(sigma_formula),
        is.logical(save)
    )

    # Set up futures plan
    cores <- metadata(sce)$dyn_corr$cores
    cl <- parallel::makeCluster(cores)
    future::plan(future::cluster, workers = cl)
    on.exit({ future::plan(future::sequential); parallel::stopCluster(cl) },
            add = TRUE)

    # Save information on margins to metadata
    metadata(sce)$dyn_corr$family <- family
    metadata(sce)$dyn_corr$mu_formula <- mu_formula
    metadata(sce)$dyn_corr$sigma_formula <- sigma_formula

    dyn_corr <- metadata(sce)$dyn_corr
    X <- SummarizedExperiment::assay(sce, dyn_corr$assay)
    X <- Matrix::t(X[dyn_corr$features, ])
    pseudotimes <- sce[[dyn_corr$time_col]]

    # All combination of families/formulas to consider
    combs <- expand.grid(
        family = family,
        mu_formula = mu_formula,
        sigma_formula = sigma_formula,
        stringsAsFactors = FALSE,
        KEEP.OUT.ATTRS = FALSE
    )

    progressr::with_progress({
        pbar <- progressr::progressor(along = seq_len(dim(X)[2]))
        margins <- future.apply::future_apply(
            X = X,
            MARGIN = 2,
            FUN = function(x) {
                ddata <- data.frame(x, pseudotimes)
                names(ddata) <- c("x", dyn_corr$time_col)

                models <- lapply(seq_len(dim(combs)[1]), function(i) {
                    # Create formulas for each parameter. The formulas for mu
                    # and sigma are specified by the user. The formula for nu
                    # is constant. If the family is NBI, then the nu formula is
                    # ignored by gamlss.
                    mu_formula <- combs[i, "mu_formula"]
                    mu_formula <- stats::formula(paste("x ~", mu_formula))
                    sigma_formula <- combs[i, "sigma_formula"]
                    sigma_formula <- stats::formula(paste("~", sigma_formula))
                    nu_formula <- stats::formula("~ 1")

                    # Load the gamlss family function
                    fname <- combs[i, "family"]
                    fam <- getExportedValue("gamlss.dist", fname)

                    # Fit the distribution
                    fit <- gamlss::gamlss(
                        formula = mu_formula,
                        sigma.formula = sigma_formula,
                        nu.formula = nu_formula,
                        data = ddata,
                        family = fam,
                        control = gamlss::gamlss.control(trace = FALSE),
                        save = save
                    )
                    fit$call$family <- as.name(fname)
                    return(fit)
                })

                # Select the model which minimizes AIC
                best_ix <- which.min(sapply(models, '[[', "aic"))
                pbar()
                return(models[[best_ix]])
            },
            future.globals = TRUE,
            future.seed = TRUE,
            future.packages = c("gamlss.dist")
        )
        margins <- future::value(margins)
    })

    names(margins) <- colnames(X)
    metadata(sce)$dyn_corr$margins <- margins

    # Compute jittered pseudo-observations
    pobs <- .pseudo_obs(sce)
    FX <- pobs$FX
    FXm <- pobs$FXm
    V <- matrix(stats::runif(prod(dim(FX))), nrow = nrow(FX))
    FXj <- FXm + (FX - FXm) * V

    # Save the pseudo-observations and jittering matrix separately
    metadata(sce)$dyn_corr$FX <- FX
    metadata(sce)$dyn_corr$FXm <- FXm
    metadata(sce)$dyn_corr$V <- V

    return(sce)
}

###############################################################################

#' Compute pseudo-observations
#'
#' @description Compute pseudo-observations and left-limits of
#' pseudo-observations for count-valued data.
#'
#' @param sce A \code{SingleCellExperiment}.
#'
#' @returns A list containing:
#' \itemize{
#'    \item \code{FX} Matrix of pseudo-observations.
#'    \item \code{FXm} Matrix of left limits of pseudo-observations.
#' }
.pseudo_obs <- function(sce) {
    dyn_corr <- metadata(sce)$dyn_corr
    counts <- SummarizedExperiment::assay(sce, dyn_corr$assay)
    counts <- Matrix::t(counts[dyn_corr$features, ])
    pseudotimes <- SummarizedExperiment::colData(sce)[[dyn_corr$time_col]]

    pobs <- lapply(seq_along(dyn_corr$features), function(i) {
        mdat <- data.frame(counts[, i], pseudotimes)
        names(mdat) <- c("x", dyn_corr$time_col)
        mfit <- dyn_corr$margins[[i]]

        par <- gamlss::predictAll(
            object = mfit,
            type = "response",
            data = mdat
        )
        par$y <- NULL

        # Load the correct quantile function
        pfun <- getExportedValue("gamlss.dist", paste0("p", mfit$family[1]))
        pX <- function(x) {
            do.call(pfun, c(list(q = x), par))
        }
        FX <- pX(counts[, i])
        FXm <- pX(counts[, i] - 1)

        # Push values away from boundaries of unit cube
        eps <- 1e-12
        FX[FX > 1 - eps] <- 1 - eps
        FX[FX < eps] <- eps
        FXm[FXm > 1 - eps] <- 1 - eps
        FXm[FXm < eps] <- eps

        return(list(FX = FX, FXm = FXm))
    })

    FX <- do.call(cbind, lapply(pobs, '[[', "FX"))
    FXm <- do.call(cbind, lapply(pobs, '[[', "FXm"))
    return(list(FX = FX, FXm = FXm))
}

###############################################################################
