###############################################################################

#' Fit marginal distributions to each gene
#'
#' @description Fit GAMLSS models to transcript counts for each gene.
#'
#' @param sce A \code{SingleCellExperiment}.
#' @param features Which genes to model. Default is all genes.
#' @param family The parametric family to use. Options are \code{"NBI"} and
#' \code{"ZINBI"}.
#' @param mu_formula Formula for mean vector.
#' @param sigma_formula Formula for dispersion parameter.
#' @param nu_formula Formula for zero proportion parameter.
#' @param save Whether to save all model information. Default is \code{TRUE}.
#' @param cores Number of cores to use. Default is \code{1}.
#' @param cl_type Type of cluster for parallel computations. If \code{NULL},
#' the value of \code{snow::getClusterOption("type")} is used. See
#' \link[parallel]{makeCluster} for details.
#'
#' @details Variable names in the formulas should be column metadata names.
#'
#' The argument \code{save} should be set to \code{TRUE} if you plan on sampling
#' cells from the fitted model. If not, setting it to \code{FALSE} will save a
#' large amount of memory. If it is set to \code{FALSE} and you later decide you
#' want to sample cells, you will need to refit the margins with
#' \code{save = TRUE}.
#'
#' @returns The same \code{SingleCellExperiment} as was passed as input, but
#' modified to include a named metadata entry \code{copula_fit} that is used to
#' store all information related to this analysis.
#'
#' @export
fit_margins <- function(sce,
                        features = rownames(sce),
                        family = c("NBI", "ZINBI"),
                        mu_formula,
                        sigma_formula,
                        nu_formula,
                        save = TRUE,
                        cores = 1,
                        cl_type = NULL) {
    assert_that(
        methods::is(sce, "SingleCellExperiment"),
        is.character(features),
        all(features %in% rownames(sce)),
        methods::is(mu_formula, "formula"),
        methods::is(sigma_formula, "formula"),
        is.numeric(cores) && cores >= 1,
        is.logical(save)
    )
    family <- match.arg(family)
    cores <- as.integer(cores)

    # Only ZINBI needs a nu_formula
    if (family == "ZINBI") {
        assert_that(methods::is(nu_formula, "formula"))
    } else {
        nu_formula <- NA
    }

    ###########################################################################
    # This purpose of this line is solely to avoid triggering a check note about
    # gamlss.dist being imported but not used, since
    # getExportedValue("gamlss.dist", ) is not recognized as using gamlss.dist.
    GAMLSS_DIST_NOTE <- gamlss.dist::dNBI(1)
    ###########################################################################

    # Set up futures plan
    run_parallel <- (cores > 1L)
    if (run_parallel) {
        if (is.null(cl_type)) {
            cl_type <- snow::getClusterOption("type")
        }
        cl <- parallel::makeCluster(cores, type = cl_type)
        future::plan(future::cluster, workers = cl)
        on.exit({
            future::plan(future::sequential)
            parallel::stopCluster(cl)
        }, add = TRUE)
    }

    # Extract counts matrix
    X <- SummarizedExperiment::assay(sce, "counts")
    X <- Matrix::t(X[features, ])

    # Design matrix
    design <- as.data.frame(SummarizedExperiment::colData(sce))
    if ("countsforgene" %in% colnames(design)) {
        stop("'countsforgene' cannot be the name of a column in colData.")
    }

    apply_fun <- ifelse(run_parallel, future.apply::future_apply, apply)
    apply_args <- list(X = X, MARGIN = 2)

    if (run_parallel) {
        apply_args <- c(apply_args, list(
            future.globals = c("family", "mu_formula", "sigma_formula",
                               "nu_formula", "save"),
            future.seed = TRUE,
            future.packages = c("gamlss", "gamlss.dist")
        ))
    }

    ngenes <- dim(X)[2]
    progressr::with_progress({
        pbar <- progressr::progressor(along = seq_len(ngenes))
        apply_args <- c(apply_args, list(
            FUN = function(x) {
                ddata <- cbind(design, data.frame(countsforgene = x))

                # The mu formula requires a response variable
                fmu <- stats::update(mu_formula, paste("countsforgene", "~ ."))

                # Load the correct gamlss family object
                gamlss_family <- getExportedValue("gamlss.dist", family)

                # Fit the model
                mfit <- gamlss::gamlss(
                    formula = fmu,
                    sigma.formula = sigma_formula,
                    nu.formula = nu_formula,
                    data = ddata,
                    family = gamlss_family,
                    control = gamlss::gamlss.control(trace = FALSE)
                )
                mfit$call$family <- as.name(family)

                # Get model parameters
                par <- gamlss::predictAll(
                    object = mfit,
                    type = "response",
                    data = ddata
                )
                par$y <- NULL

                # Load the correct distribution function
                pfun <- getExportedValue("gamlss.dist", paste0("p", family))
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

    margins <- lapply(model_fits, '[[', "model")
    names(margins) <- colnames(X)

    # Save information about margins to metadata
    metadata(sce)$copula_fit <- list(
        features = features,
        family = family,
        mu_formula = mu_formula,
        sigma_formula = sigma_formula,
        nu_formula = nu_formula,
        save = save
    )

    # The models are very large and not worth saving unless we need them again
    if (save) {
        metadata(sce)$copula_fit$margins <- margins
    } else {
        metadata(sce)$copula_fit$margins <- NULL
    }

    # Pseudo-observations
    FX <- do.call(cbind, lapply(model_fits, '[[', "FX"))
    # Left-limits of pseudo-observations
    FXm <- do.call(cbind, lapply(model_fits, '[[', "FXm"))
    # Jittering matrix
    V <- matrix(stats::runif(prod(dim(FX))), nrow = nrow(FX))

    # Jittered pseudo-observations
    FX <- FXm + (FX - FXm) * V
    colnames(FX) <- colnames(X)
    metadata(sce)$copula_fit$FX <- FX

    return(sce)
}

###############################################################################
