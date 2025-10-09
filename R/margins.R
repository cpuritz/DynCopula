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
#'    \item \code{FX} Pseudo-observations
#'    \item \code{FXm} Left limits of pseudo-observations
#' }
.pseudo_obs <- function(sce) {
    info <- metadata(sce)$dyn_corr_info
    counts <- SummarizedExperiment::assay(sce, info$assay)
    counts <- Matrix::t(counts[info$features, ])
    pseudotimes <- SummarizedExperiment::colData(sce)[[info$time_col]]

    pobs <- lapply(seq_along(info$features), function(i) {
        mdat <- data.frame(counts[, i], pseudotimes)
        names(mdat) <- c("x", info$time_col)
        mfit <- metadata(sce)$margins[[i]]

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
        return(list(FX = pX(counts[, i]), FXm = pX(counts[, i] - 1)))
    })
    return(list(
        FX = do.call(cbind, lapply(pobs, '[[', "FX")),
        FXm = do.call(cbind, lapply(pobs, '[[', "FXm"))
    ))
}

###############################################################################

#' Fit margins
#'
#' @description Fit parametric models for each gene. Multiple families and
#' formulas can be passed, in which case the best fitting model is chosen across
#' all combinations.
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
#' modified to include a named metadata entry \code{margins} which is a list
#' of marginal models.
#'
#' @details Formulas should be written as a function of the name of the
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
#' Model selection is done using AIC.
#'
#' @export
fit_margins <- function(sce,
                        family = c("NBI", "ZINBI"),
                        mu_formula,
                        sigma_formula) {
    assert_that(
        methods::is(sce, "SingleCellExperiment"),
        "dyn_corr_info" %in% names(metadata(sce)),
        is.character(mu_formula),
        is.character(sigma_formula)
    )
    family <- match.arg(family, several.ok = TRUE)

    metadata(sce)$dyn_corr_info$family <- family
    metadata(sce)$dyn_corr_info$mu_formula <- mu_formula
    metadata(sce)$dyn_corr_info$sigma_formula <- sigma_formula

    message("Fitting marginal distributions")

    info <- metadata(sce)$dyn_corr_info
    X <- SummarizedExperiment::assay(sce, info$assay)
    X <- Matrix::t(X[info$features, ])
    pseudotimes <- SummarizedExperiment::colData(sce)[[info$time_col]]

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
                names(ddata) <- c("x", info$time_col)
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
    metadata(sce)$margins <- margins
    return(sce)
}

###############################################################################
