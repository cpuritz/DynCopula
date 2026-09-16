# Test argument checks for fit_dynamic_gaussian
test_that("fit_gamgc arguments", {
    d <- 3
    N <- 100
    FX <- matrix(rep(0, d * N), ncol = d)
    t <- seq(0, 1, length.out = N)
    x1 <- rep(c("a", "b"), N / 2)
    design <- data.frame(t = t, x1 = x1)

    # Test FX
    expect_error(fit_gamgc(NULL, design, formula = ~1))
    expect_error(fit_gamgc(seq(N), design, formula = ~1))
    expect_error(fit_gamgc(as.data.frame(FX), design, formula = ~1))
    expect_error(fit_gamgc(FX[1:(N - 1), ], design, formula = ~1))

    # Test design
    expect_error(fit_gamgc(FX, NULL, formula = ~1))
    expect_error(fit_gamgc(FX, design[, 1], formula = ~1))
    expect_error(fit_gamgc(FX, design[1:(N - 1), ], formula = ~1))

    ## Test formula
    # Missing formula
    expect_error(fit_gamgc(FX, design))
    expect_error(fit_gamgc(FX, design, formula = NULL))
    # Formula passed as string
    expect_error(fit_gamgc(FX, design, formula = "~1"))
    # Formula includes variable not in design matrix
    expect_error(fit_gamgc(FX, design, formula = ~notavar))
    expect_error(fit_gamgc(FX, design, formula = ~t + notavar))
    # Multiple smooth variables
    expect_error(fit_gamgc(FX, design, formula = ~s(t) + s(x1)))
    # Trying to use exact mgcv style
    #expect_error(fit_gamgc(FX, design, formula = ~s(t, k = 30)))

    # Test lambda
    bad_lambda <- list(NULL, NA, NaN, Inf, 0, -1, c(1, 2, 0), c(-1, 2, 1), "a")
    for (x in bad_lambda) {
        expect_error(fit_gamgc(FX, design, ~1, lambda = x))
    }

    # Test K
    bad_K <- list(NULL, NA, NaN, Inf, "a", 0, 1, 2, 2.99999, "3", "4L")
    for (x in bad_K) {
        expect_error(fit_gamgc(FX, design, ~1, K = x))
    }

    # Test nfold
    bad_nfold <- list(NULL, NA, NaN, Inf, "a", 0, 1, 1.99999, "3", "4L")
    for (x in bad_nfold) {
        expect_error(fit_gamgc(FX, design, ~1, nfold = x))
    }

    # Test cores
    bad_cores <- list(NULL, NA, NaN, Inf, "a", 0, 0.99999, "2", "3L")
    for (x in bad_cores) {
        expect_error(fit_gamgc(FX, design, ~1, cores = x))
    }

    ## Test control ##
    # Control must be a list, not a vector
    expect_error(fit_gamgc(FX, design, ~1, control = c()))
    expect_error(fit_gamgc(FX, design, ~1, control = c(history_size = 1)))

    bad_max_itr <- list(-1, NA, NULL, 0, "a", 0.99, "2", "3L")
    for (x in bad_max_itr) {
        expect_error(
            fit_gamgc(FX, design, ~1, control = list(max_itr = x))
        )
    }

    bad_history_size <- list(-1, NA, NULL, 0, "a", 0.99, "2", "3L")
    for (x in bad_history_size) {
        expect_error(
            fit_gamgc(FX, design, ~1, control = list(history_size = x))
        )
    }

    bad_tol_grad <- list(-1, NA, NULL, -0.001, 0, "a", "0.1")
    for (x in bad_tol_grad) {
        expect_error(
            fit_gamgc(FX, design, ~1, control = list(tolerance_grad = x))
        )
    }

    bad_tol_change <- list(-1, NA, NULL, -0.001, 0, "a", "0.1")
    for (x in bad_tol_change) {
        expect_error(
            fit_gamgc(FX, design, ~1, control = list(tolerance_change = x))
        )
    }

    bad_precision <- list(1, NA, NULL, "float23", "32", 32)
    for (x in bad_precision) {
        expect_error(
            fit_gamgc(FX, design, ~1, control = list(precision = x))
        )
    }

    bad_boundary <- list(-1, -0.1, -0.001, -1e-6, "a", "0", "0.1", NULL, NA)
    for (x in bad_boundary) {
        expect_error(
            fit_gamgc(FX, design, ~1, control = list(boundary = x))
        )
    }
})
