###############################################################################

#' Compute pseudo-observations
#'
#' @description Compute pseudo-observations and left-limits of
#' pseudo-observations for count-valued data.
#'
#' @param sce A \code{SingleCellExperiment}.
#'
#' @returns A list containing
#' \itemize{
#'    \item \code{FX} Pseudo-observations.
#'    \item \code{FXm} Left limits of pseudo-observations.
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

        pfun <- getExportedValue("gamlss.dist", paste0("p", mfit$family[1]))
        pX <- function(x) {
            do.call(pfun, c(list(q = x), par))
        }
        FX <- pX(counts[, i])
        FXm <- pX(counts[, i] - 1)

        # Push away from boundaries of unit cube
        eps <- 1e-10
        FX[FX == 1] <- 1 - eps
        FX[FX == 0] <- eps
        FXm[FXm == 1] <- 1 - eps
        FXm[FXm == 0] <- eps

        return(list(FX = FX, FXm = FXm))
    })
    return(list(
        FX = do.call(cbind, lapply(pobs, '[[', "FX")),
        FXm = do.call(cbind, lapply(pobs, '[[', "FXm"))
    ))
}

###############################################################################

#' Fit margins
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
#'
#' @returns The same \code{SingleCellExperiment} as was passed as input, but
#' with the metadata entry \code{dyn_corr} updated to include the following
#' entries:
#' \itemize{
#'   \item \code{margins} List of marginal models.
#'   \item \code{FX} Pseudo-observations.
#'   \item \code{FXm} Left limits of pseudo-observations.
#'   \item \code{V} Jittering matrix.
#' }
#'
#' @details Multiple families and formulas can be passed, in which case the best
#' fitting model is chosen across all combinations using AIC.
#'
#' Formulas should be written as a function of the name of the
#' pseudotime column (what was passed as the argument \code{time_col} to
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
#' @export
fit_margins <- function(sce,
                        family = c("NBI", "ZINBI"),
                        mu_formula,
                        sigma_formula) {
    assert_that(
        methods::is(sce, "SingleCellExperiment"),
        "dyn_corr" %in% names(metadata(sce)),
        is.character(mu_formula),
        is.character(sigma_formula)
    )
    family <- match.arg(family, several.ok = TRUE)

    cores <- metadata(sce)$dyn_corr$cores
    cl <- parallel::makeCluster(cores)
    future::plan(future::cluster, workers = cl)
    on.exit({ future::plan(future::sequential); parallel::stopCluster(cl) },
            add = TRUE)

    message("Fitting marginal distributions")
    metadata(sce)$dyn_corr$family <- family
    metadata(sce)$dyn_corr$mu_formula <- mu_formula
    metadata(sce)$dyn_corr$sigma_formula <- sigma_formula

    dyn_corr <- metadata(sce)$dyn_corr
    X <- SummarizedExperiment::assay(sce, dyn_corr$assay)
    X <- Matrix::t(X[dyn_corr$features, ])
    pseudotimes <- sce[[dyn_corr$time_col]]

    combs <- expand.grid(
        family = family,
        mu_formula = mu_formula,
        sigma_formula = sigma_formula,
        stringsAsFactors = FALSE,
        KEEP.OUT.ATTRS = FALSE
    )

    progressr::with_progress({
        pbar <- progressr::progressor(along = seq_len(dim(X)[2]))

        X_slices <- apply(X, 2, c, simplify = FALSE)
        margins <- future.apply::future_lapply(
            X = X_slices,
            FUN = function(x) {
                ddata <- data.frame(x, pseudotimes)
                names(ddata) <- c("x", dyn_corr$time_col)
                models <- lapply(seq_len(dim(combs)[1]), function(i) {
                    fmu <- paste("x ~", combs[i, "mu_formula"])
                    fsigma <- paste("~", combs[i, "sigma_formula"])
                    fname <- combs[i, "family"]
                    fam <- getExportedValue("gamlss.dist", fname)
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
    metadata(sce)$dyn_corr$margins <- margins

    pobs <- .pseudo_obs(sce)
    FX <- pobs$FX
    FXm <- pobs$FXm
    V <- matrix(stats::runif(prod(dim(FX))), nrow = nrow(FX))
    FXj <- FXm + (FX - FXm) * V

    metadata(sce)$dyn_corr$FX <- FX
    metadata(sce)$dyn_corr$FXm <- FXm
    metadata(sce)$dyn_corr$V <- V

    return(sce)
}

###############################################################################

#' Fit metacell margins
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
#'
#' @returns The same \code{SingleCellExperiment} as was passed as input, but
#' with the metacell \code{SingleCellExperiment} updated.
#'
#' @details Multiple families and formulas can be passed, in which case the best
#' fitting model is chosen across all combinations using AIC.
#'
#' Formulas should be written as a function of the name of the
#' pseudotime column (what was passed as the argument \code{time_col} to
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
#' @export
fit_metacell_margins <- function(sce,
                                 family = c("NBI", "ZINBI"),
                                 mu_formula,
                                 sigma_formula) {
    assert_that(
        methods::is(sce, "SingleCellExperiment"),
        "dyn_corr" %in% names(metadata(sce)),
        is.character(mu_formula),
        is.character(sigma_formula),
        "metacell_sce" %in% names(metadata(sce)$dyn_corr)
    )
    family <- match.arg(family, several.ok = TRUE)

    cores <- metadata(sce)$dyn_corr$cores
    cl <- parallel::makeCluster(cores)
    future::plan(future::cluster, workers = cl)
    on.exit({ future::plan(future::sequential); parallel::stopCluster(cl) },
            add = TRUE)

    message("Fitting marginal distributions for metacells")
    mc_sce <- metadata(sce)$dyn_corr$metacell_sce

    metadata(mc_sce)$dyn_corr$family <- family
    metadata(mc_sce)$dyn_corr$mu_formula <- mu_formula
    metadata(mc_sce)$dyn_corr$sigma_formula <- sigma_formula

    dyn_corr <- metadata(mc_sce)$dyn_corr
    X <- SummarizedExperiment::assay(mc_sce, dyn_corr$assay)
    X <- Matrix::t(X[dyn_corr$features, ])
    pseudotimes <- mc_sce[[dyn_corr$time_col]]

    combs <- expand.grid(
        family = family,
        mu_formula = mu_formula,
        sigma_formula = sigma_formula,
        stringsAsFactors = FALSE,
        KEEP.OUT.ATTRS = FALSE
    )

    progressr::with_progress({
        pbar <- progressr::progressor(along = seq_len(dim(X)[2]))

        X_slices <- apply(X, 2, c, simplify = FALSE)
        margins <- future.apply::future_lapply(
            X = X_slices,
            FUN = function(x) {
                ddata <- data.frame(x, pseudotimes)
                names(ddata) <- c("x", dyn_corr$time_col)
                models <- lapply(seq_len(dim(combs)[1]), function(i) {
                    fmu <- paste("x ~", combs[i, "mu_formula"])
                    fsigma <- paste("~", combs[i, "sigma_formula"])
                    fname <- combs[i, "family"]
                    fam <- getExportedValue("gamlss.dist", fname)
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
    metadata(mc_sce)$dyn_corr$margins <- margins

    pobs <- .pseudo_obs(mc_sce)
    FX <- pobs$FX
    FXm <- pobs$FXm
    V <- matrix(stats::runif(prod(dim(FX))), nrow = nrow(FX))
    FXj <- FXm + (FX - FXm) * V

    metadata(mc_sce)$dyn_corr$FX <- FX
    metadata(mc_sce)$dyn_corr$FXm <- FXm
    metadata(mc_sce)$dyn_corr$V <- V

    metadata(sce)$dyn_corr$metacell_sce <- mc_sce
    return(sce)
}

###############################################################################
