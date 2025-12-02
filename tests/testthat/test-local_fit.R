# Test argument checks for fit_dynamic_gaussian
test_that("fit_dynamic_gaussian arguments", {
    d <- 3
    N <- 100
    FX <- matrix(rep(0, d * N), ncol = d)
    x <- seq(0, 1, length.out = N)
    x0 <- 0.5
    h <- 0.1

    # Test FX
    expect_error(fit_dynamic_gaussian(FX[1:2, ], x, x0, h))
    expect_error(fit_dynamic_gaussian(FX[1:(N - 1), ], x, x0, h))
    expect_error(fit_dynamic_gaussian(FX[, 1, drop = FALSE], x, x0, h))
    expect_error(fit_dynamic_gaussian(1, x, x0, h))
    expect_error(fit_dynamic_gaussian("a", x, x0, h))
    expect_error(fit_dynamic_gaussian(rep("a", N), x, x0, h))
    expect_error(fit_dynamic_gaussian(as.data.frame(FX), x, x0, h))

    # Test x
    expect_error(fit_dynamic_gaussian(FX, 0, x0, h))
    expect_error(fit_dynamic_gaussian(FX, x[1:2], x0, h))
    expect_error(fit_dynamic_gaussian(FX, x[1:(N - 1)], x0, h))
    expect_error(fit_dynamic_gaussian(FX, "a", x0, h))
    expect_error(fit_dynamic_gaussian(FX, rep("a", N), x0, h))
    expect_error(fit_dynamic_gaussian(FX, c(x[1], x[1:99]), x0, h))
    expect_error(fit_dynamic_gaussian(FX, c(x[1:99], -0.1), x0, h))

    # Test x0
    expect_error(fit_dynamic_gaussian(FX, x, "a", h))
    expect_error(fit_dynamic_gaussian(FX, x, c("a", "a"), h))

    # Test h
    expect_error(fit_dynamic_gaussian(FX, x, x0, 0))
    expect_error(fit_dynamic_gaussian(FX, x, x0, 1))
    expect_error(fit_dynamic_gaussian(FX, x, x0, -0.01))
    expect_error(fit_dynamic_gaussian(FX, x, x0, 1.01))
    expect_error(fit_dynamic_gaussian(FX, x, x0, "a"))

    # Test cores
    expect_error(fit_dynamic_gaussian(FX, x, x0, h, cores = 0L))
    expect_error(fit_dynamic_gaussian(FX, x, x0, h, cores = 0.0))
    expect_error(fit_dynamic_gaussian(FX, x, x0, h, cores = -1L))
    expect_error(fit_dynamic_gaussian(FX, x, x0, h, cores = -0.1))
    expect_error(fit_dynamic_gaussian(FX, x, x0, h, cores = -1.1))
    expect_error(fit_dynamic_gaussian(FX, x, x0, h, cores = 0.9999))
    expect_error(fit_dynamic_gaussian(FX, x, x0, h, cores = "a"))

    # Test degree
    expect_error(fit_dynamic_gaussian(FX, x, x0, h, degree = -1L))
    expect_error(fit_dynamic_gaussian(FX, x, x0, h, degree = -0.001))
    expect_error(fit_dynamic_gaussian(FX, x, x0, h, degree = "a"))

    # Test control
    expect_error(fit_dynamic_gaussian(FX, x, x0, h, control = c()))
    expect_error(fit_dynamic_gaussian(FX, x, x0, h, control = c(max_itr = 10L)))

    # Test max_outer
    expect_error(fit_dynamic_gaussian(FX, x, x0, h,
                                      control = list(max_outer = 0L)))
    expect_error(fit_dynamic_gaussian(FX, x, x0, h,
                                      control = list(max_outer = 0.99)))
    expect_error(fit_dynamic_gaussian(FX, x, x0, h,
                                      control = list(max_outer = -1L)))
    expect_error(fit_dynamic_gaussian(FX, x, x0, h,
                                      control = list(max_outer = "a")))

    # Test max_itr
    expect_error(fit_dynamic_gaussian(FX, x, x0, h,
                                      control = list(max_itr = 0L)))
    expect_error(fit_dynamic_gaussian(FX, x, x0, h,
                                      control = list(max_itr = 0.99)))
    expect_error(fit_dynamic_gaussian(FX, x, x0, h,
                                      control = list(max_itr = -1L)))
    expect_error(fit_dynamic_gaussian(FX, x, x0, h,
                                      control = list(max_itr = "a")))

    # Test history_size
    expect_error(fit_dynamic_gaussian(FX, x, x0, h,
                                      control = list(history_size = 0L)))
    expect_error(fit_dynamic_gaussian(FX, x, x0, h,
                                      control = list(history_size = 0.99)))
    expect_error(fit_dynamic_gaussian(FX, x, x0, h,
                                      control = list(history_size = -1L)))
    expect_error(fit_dynamic_gaussian(FX, x, x0, h,
                                      control = list(history_size = "a")))

    # Test tolerance_grad
    expect_error(fit_dynamic_gaussian(FX, x, x0, h,
                                      control = list(tolerance_grad = 0)))
    expect_error(fit_dynamic_gaussian(FX, x, x0, h,
                                      control = list(tolerance_grad = -0.1)))
    expect_error(fit_dynamic_gaussian(FX, x, x0, h,
                                      control = list(tolerance_grad = "a")))

    # Test tolerance_change
    expect_error(fit_dynamic_gaussian(FX, x, x0, h,
                                      control = list(tolerance_change = 0)))
    expect_error(fit_dynamic_gaussian(FX, x, x0, h,
                                      control = list(tolerance_change = -0.1)))
    expect_error(fit_dynamic_gaussian(FX, x, x0, h,
                                      control = list(tolerance_change = "a")))
})
