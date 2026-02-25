###############################################################################

#' Fit marginal distributions to each gene
#'
#' @description Fit marginal distributions to transcript counts for each gene.
#' If the assay to model is \code{"counts"}, parametric models are fit. If the
#' assay to model is \code{"logcounts"}, nonparametric models are fit.
#'
#' @param sce A \code{SingleCellExperiment}.
#' @param quick Whether to perform cross-validated bandwidth selection (which is
#' expensive but more accurate) or use a quick plug-in estimator. Only used if
#' \code{assay = "logcounts"}. Default is \code{TRUE}.
#' @param family A character vector specifying families to use. Only used
#' if \code{assay = "counts"}. Options are \code{"NBI"} (negative binomial) and
#' \code{"ZINBI"} (zero-inflated negative binomial).
#' @param mu_formula A character vector specifying formulas for the mean
#' parameter. Only used if \code{assay = "counts"}.
#' @param sigma_formula A character vector specifying formulas for the
#' dispersion parameter. Only used if \code{assay = "counts"}.
#'
#' @returns The same \code{SingleCellExperiment} as was passed as input, but
#' with the metadata entry \code{dyn_corr} updated to include necessary
#' information about the margins.
#'
#' @details If \code{assay = "logcounts"}, kernel-smoothed hurdle models are
#' fit. The only modeling decision to make is whether to perform cross-validated
#' bandwidth selection for the margins. This approach is more accurate than
#' using the plug-in bandwidth selector, but is much more expensive. See the
#' \code{bwmethod} argument of \link[np]{npcdistbw} for details.
#'
#' If \code{assay = "counts"}, parametric models are fit. Multiple families and
#' formulas can be passed, in which case the best fitting model is chosen across
#' all combinations by minimizing AIC. Formulas should be written as a function
#' of the name of the pseudotime column (what was passed as the argument
#' \code{time_col} to \link[DynCopula]{setup}). Only the righthand side of the
#' formula should be passed. For example, if the pseudotime column is named
#' \code{"pseudotime"}:
#' \itemize{
#'   \item \code{mu_formula = "1"}: mean is constant.
#'   \item \code{mu_formula = "pseudotime"}: mean varies linearly with
#'   pseudotime.
#'   \item \code{mu_formula = "gamlss::pb(pseudotime)"}: mean varies through a
#'   P-spline with pseudotime. Any spline functions need to be explicitly
#'   scoped (e.g., \code{gamlss::pb} instead of \code{pb}).
#' }
#' For ZINBI models, the zero probability parameter is assumed to be constant.
#'
#' @export
fit_margins <- function(sce,
                        quick = TRUE,
                        family = c("NBI", "ZINBI"),
                        mu_formula,
                        sigma_formula) {
    assay <- metadata(sce)$dyn_corr$assay
    if (assay == "logcounts") {
        assert_that(is.logical(quick))
        sce <- .fit_hurdle_margins(
            sce = sce,
            quick = quick
        )
    } else {
        family <- match.arg(family, several.ok = TRUE)
        # We need margin information to sample new cells, so 'save' is set to
        # TRUE.
        sce <- .fit_count_margins(
            sce = sce,
            family = family,
            mu_formula = mu_formula,
            sigma_formula = sigma_formula,
            save = TRUE
        )
    }
    return(sce)
}

###############################################################################

#' Fit marginal distributions to log-counts of each gene
#'
#' @description Fit kernel-smoothed hurdle models to log-counts of each gene.
#' Internal function.
#'
#' @param sce A \code{SingleCellExperiment}.
#' @param quick Whether bandwidth selection should be quick or not.
#'
#' @returns The same \code{SingleCellExperiment} as was passed as input, but
#' with the metadata entry \code{dyn_corr} updated to include \code{FX}, the
#' matrix of pseudo-observations.
.fit_hurdle_margins <- function(sce, quick) {
    dyn_corr <- metadata(sce)$dyn_corr
    X <- SummarizedExperiment::assay(sce, "logcounts")
    X <- Matrix::t(X[dyn_corr$features, ])
    time_col <- dyn_corr$time_col
    pseudotimes <- sce[[time_col]]

    # Set up futures plan
    cores <- metadata(sce)$dyn_corr$cores
    run_parallel <- (cores > 1L)
    if (run_parallel) {
        cl <- parallel::makeCluster(cores)
        future::plan(future::cluster, workers = cl)
        on.exit({
            future::plan(future::sequential)
            parallel::stopCluster(cl)
        }, add = TRUE)
    }
    apply_fun <- ifelse(run_parallel, future.apply::future_apply, apply)

    if (quick) {
        bwmethod <- "normal-reference"
    } else {
        bwmethod <- "cv.ls"
    }

    ngenes <- dim(X)[2]
    progressr::with_progress({
        pbar <- progressr::progressor(along = seq_len(ngenes))
        apply_args <- list(
            X = X,
            MARGIN = 2,
            FUN = function(x) {
                is_pos <- (x > 0)
                ind0 <- as.numeric(!is_pos)

                # Estimate P(Y = 0 | t)
                sink <- utils::capture.output(
                    # No option for a plug-in bandwidth estimator, we have to use
                    # a data-driven approach. This is expensive, but much less expensive
                    # than cross-validation for the bivariate np::npcdistbw.
                    bw0 <- np::npregbw(
                        xdat = pseudotimes,
                        ydat = ind0,
                        bwmethod = "cv.ls",
                        nmulti = 1L
                    )
                )
                fit0 <- np::npreg(bws = bw0)
                p0_hat <- stats::predict(fit0)

                # Estimate P(Y <= y | Y > 0, t)
                t_pos <- pseudotimes[is_pos]
                x_pos <- x[is_pos]
                # We use the plug-in bandwidth selector instead of cross-validation. The
                # latter will of course be much more accurate, but it is way too
                # expensive for a large number of genes.
                sink <- utils::capture.output(
                    bwF <- np::npcdistbw(
                        xdat = t_pos,
                        ydat = x_pos,
                        bwmethod = bwmethod
                    )
                )
                fitF <- np::npcdist(bws = bwF)
                Fpos <- stats::predict(fitF)

                # Hurdle CDF
                FX <- p0_hat
                FX[is_pos] <- p0_hat[is_pos] + (1 - p0_hat[is_pos]) * Fpos

                # Jittering at atom at 0
                V <- stats::runif(length(p0_hat), min = 0, max = p0_hat)
                FX[!is_pos] <- V[!is_pos]

                # Push values away from boundaries of unit cube
                eps <- 1e-12
                FX[FX > 1 - eps] <- 1 - eps
                FX[FX < eps] <- eps

                pbar()

                return(FX)
            }
        )
        if (run_parallel) {
            apply_args <- c(apply_args, list(
                future.globals = c("pseudotimes", "bwmethod"),
                future.seed = TRUE,
                future.packages = c("np")
            ))
        }
        FX <- do.call(apply_fun, apply_args)
    })

    colnames(FX) <- colnames(X)
    metadata(sce)$dyn_corr$FX <- FX
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
#' \code{"NBI"} and \code{"ZINBI"}.
#' @param mu_formula A character vector specifying formulas for the mean
#' parameter.
#' @param sigma_formula A character vector specifying formulas for the
#' dispersion parameter.
#' @param save Whether to save all information about margins.
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
.fit_count_margins <- function(sce,
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
    run_parallel <- (cores > 1L)
    if (run_parallel) {
        cl <- parallel::makeCluster(cores)
        future::plan(future::cluster, workers = cl)
        on.exit({
            future::plan(future::sequential)
            parallel::stopCluster(cl)
        }, add = TRUE)
    }

    # Save information on margins to metadata
    metadata(sce)$dyn_corr$family <- family
    metadata(sce)$dyn_corr$mu_formula <- mu_formula
    metadata(sce)$dyn_corr$sigma_formula <- sigma_formula

    dyn_corr <- metadata(sce)$dyn_corr
    X <- SummarizedExperiment::assay(sce, dyn_corr$assay)
    X <- Matrix::t(X[dyn_corr$features, ])
    time_col <- dyn_corr$time_col
    pseudotimes <- sce[[time_col]]

    # All combination of families/formulas to consider
    combs <- expand.grid(
        family = family,
        mu_formula = mu_formula,
        sigma_formula = sigma_formula,
        stringsAsFactors = FALSE,
        KEEP.OUT.ATTRS = FALSE
    )

    ###########################################################################
    # This purpose of this line is solely to avoid triggering a check note about
    # gamlss.dist being imported but not used, since
    # getExportedValue("gamlss.dist", ) is not recognized as using gamlss.dist.
    GAMLSS_DIST_NOTE <- gamlss.dist::dNBI(1)
    ###########################################################################

    apply_fun <- ifelse(run_parallel, future.apply::future_apply, apply)
    apply_args <- list(X = X, MARGIN = 2)

    if (run_parallel) {
        apply_args <- c(apply_args, list(
            future.globals = c("pseudotimes", "time_col", "combs"),
            future.seed = TRUE,
            future.packages = c("gamlss", "gamlss.dist")
        ))
    }

    ngenes <- dim(X)[2]
    progressr::with_progress({
        pbar <- progressr::progressor(along = seq_len(ngenes))
        apply_args <- c(apply_args, list(
            FUN = function(x) {
                ddata <- data.frame(x, pseudotimes)
                names(ddata) <- c("x", time_col)

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
                mfit <- models[[best_ix]]

                # Get model parameters
                par <- gamlss::predictAll(
                    object = mfit,
                    type = "response",
                    data = ddata
                )
                par$y <- NULL

                # Load the correct distribution function
                pname <- paste0("p", mfit$family[1])
                pfun <- getExportedValue("gamlss.dist", pname)
                # Pseudo-observations
                FX <- do.call(pfun, c(list(q = x), par))
                # Left limits of pseudo-observations
                FXm <- do.call(pfun, c(list(q = x - 1), par))

                # Push values away from boundaries of unit cube
                eps <- 1e-12
                FX[FX > 1 - eps] <- 1 - eps
                FX[FX < eps] <- eps
                FXm[FXm > 1 - eps] <- 1 - eps
                FXm[FXm < eps] <- eps

                pbar()
                return(list(model = mfit, FX = FX, FXm = FXm))
            }
        ))
        model_fits <- do.call(apply_fun, apply_args)
    })

    # Save the marginal models
    margins <- lapply(model_fits, '[[', "model")
    names(margins) <- colnames(X)
    metadata(sce)$dyn_corr$margins <- margins

    # Record whether full marginal models were saved
    metadata(sce)$dyn_corr$full_margins <- save

    # Save the pseudo-observations and jittering matrix
    FX <- do.call(cbind, lapply(model_fits, '[[', "FX"))
    FXm <- do.call(cbind, lapply(model_fits, '[[', "FXm"))
    V <- matrix(stats::runif(prod(dim(FX))), nrow = nrow(FX))
    colnames(FX) <- colnames(FXm) <- colnames(V) <- colnames(X)
    metadata(sce)$dyn_corr$FX <- FX
    metadata(sce)$dyn_corr$FXm <- FXm
    metadata(sce)$dyn_corr$V <- V

    return(sce)
}

###############################################################################
