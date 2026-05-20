# Test argument checks for fit_dynamic_gaussian
test_that("fit_gamgc arguments", {
    d <- 3
    N <- 100
    FX <- matrix(rep(0, d * N), ncol = d)
    t <- seq(0, 1, length.out = N)
    x1 <- rep(c("a", "b"), N / 2)
    design <- data.frame(t = t, x1 = x1)

    expect_error(fit_gamgc(FX, design, formula, lambda, K, nfold, cores, cl_type, control))

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
    expect_error(fit_gamgc(FX, design, formula = ~s(t, k = 30)))

    # Test lambda
    expect_error(fit_gamgc(FX, design, ~1, lambda = NULL))
    expect_error(fit_gamgc(FX, design, ~1, lambda = 0))
    expect_error(fit_gamgc(FX, design, ~1, lambda = -1))
    expect_error(fit_gamgc(FX, design, ~1, lambda = c(1, 2, 0)))
    expect_error(fit_gamgc(FX, design, ~1, lambda = c(-1, 2, 1)))

    # Test K
    expect_error(fit_gamgc(FX, design, ~1, K = NULL))
    expect_error(fit_gamgc(FX, design, ~1, K = NA))
    expect_error(fit_gamgc(FX, design, ~1, K = NaN))
    expect_error(fit_gamgc(FX, design, ~1, K = Inf))
    expect_error(fit_gamgc(FX, design, ~1, K = "a"))
    expect_error(fit_gamgc(FX, design, ~1, K = 0))
    expect_error(fit_gamgc(FX, design, ~1, K = 1))
    expect_error(fit_gamgc(FX, design, ~1, K = 2))
    expect_error(fit_gamgc(FX, design, ~1, K = 2.99999))

    # Test nfold
    expect_error(fit_gamgc(FX, design, ~1, nfold = NULL))
    expect_error(fit_gamgc(FX, design, ~1, nfold = NA))
    expect_error(fit_gamgc(FX, design, ~1, nfold = NaN))
    expect_error(fit_gamgc(FX, design, ~1, nfold = Inf))
    expect_error(fit_gamgc(FX, design, ~1, nfold = "a"))
    expect_error(fit_gamgc(FX, design, ~1, nfold = 0))
    expect_error(fit_gamgc(FX, design, ~1, nfold = 1))
    expect_error(fit_gamgc(FX, design, ~1, nfold = 1.99999))

    # Test cores
    expect_error(fit_gamgc(FX, design, ~1, cores = NULL))
    expect_error(fit_gamgc(FX, design, ~1, cores = NA))
    expect_error(fit_gamgc(FX, design, ~1, cores = NaN))
    expect_error(fit_gamgc(FX, design, ~1, cores = Inf))
    expect_error(fit_gamgc(FX, design, ~1, cores = "a"))
    expect_error(fit_gamgc(FX, design, ~1, cores = 0))
    expect_error(fit_gamgc(FX, design, ~1, cores = 0.99999))

    # Test control
    expect_error(fit_gamgc(FX, design, ~1, control = c()))
    expect_error(fit_gamgc(FX, design, ~1, control = c(history_size = 1)))
    expect_error(fit_gamgc(FX, design, ~1, control = list(max_itr = -1)))
    expect_error(fit_gamgc(FX, design, ~1, control = list(max_itr = NA)))
    expect_error(fit_gamgc(FX, design, ~1, control = list(max_itr = NULL)))
    expect_error(fit_gamgc(FX, design, ~1, control = list(max_itr = 0)))
    expect_error(fit_gamgc(FX, design, ~1, control = list(max_itr = "a")))
    expect_error(fit_gamgc(FX, design, ~1, control = list(max_itr = 0.99)))

    expect_error(fit_gamgc(FX, design, ~1, control = list(history_size = -1)))
    expect_error(fit_gamgc(FX, design, ~1, control = list(history_size = NA)))
    expect_error(fit_gamgc(FX, design, ~1, control = list(history_size = NULL)))
    expect_error(fit_gamgc(FX, design, ~1, control = list(history_size = 0)))
    expect_error(fit_gamgc(FX, design, ~1, control = list(history_size = "a")))
    expect_error(fit_gamgc(FX, design, ~1, control = list(history_size = 0.99)))

    expect_error(fit_gamgc(FX, design, ~1, control = list(tolerance_grad = -1)))
    expect_error(fit_gamgc(FX, design, ~1, control = list(tolerance_grad = NA)))
    expect_error(fit_gamgc(FX, design, ~1, control = list(tolerance_grad = NULL)))
    expect_error(fit_gamgc(FX, design, ~1, control = list(tolerance_grad = -0.001)))
    expect_error(fit_gamgc(FX, design, ~1, control = list(tolerance_grad = 0)))
    expect_error(fit_gamgc(FX, design, ~1, control = list(tolerance_grad = "a")))

    expect_error(fit_gamgc(FX, design, ~1, control = list(tolerance_change = -1)))
    expect_error(fit_gamgc(FX, design, ~1, control = list(tolerance_change = NA)))
    expect_error(fit_gamgc(FX, design, ~1, control = list(tolerance_change = NULL)))
    expect_error(fit_gamgc(FX, design, ~1, control = list(tolerance_change = -0.001)))
    expect_error(fit_gamgc(FX, design, ~1, control = list(tolerance_change = 0)))
    expect_error(fit_gamgc(FX, design, ~1, control = list(tolerance_change = "a")))

    expect_error(fit_gamgc(FX, design, ~1, control = list(precision = 1)))
    expect_error(fit_gamgc(FX, design, ~1, control = list(precision = NA)))
    expect_error(fit_gamgc(FX, design, ~1, control = list(precision = NULL)))
    expect_error(fit_gamgc(FX, design, ~1, control = list(precision = "float23")))
    expect_error(fit_gamgc(FX, design, ~1, control = list(precision = "32")))
    expect_error(fit_gamgc(FX, design, ~1, control = list(precision = 32)))

    expect_error(fit_gamgc(FX, design, ~1, control = list(boundary = -1)))
    expect_error(fit_gamgc(FX, design, ~1, control = list(boundary = -0.001)))
    expect_error(fit_gamgc(FX, design, ~1, control = list(boundary = "a")))
    expect_error(fit_gamgc(FX, design, ~1, control = list(boundary = NULL)))
    expect_error(fit_gamgc(FX, design, ~1, control = list(boundary = NA)))
})
