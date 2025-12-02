# Test that cor2vec is the inverse of vec2cor in two dimensions
test_that("vec2cor 2D", {
    etas <- seq(-10, 10, 0.1)
    for (x in etas) {
        expect_equal(cor2vec(vec2cor(x)), x)
    }
})

# Test that vec2cor is the inverse of cor2vec in two dimensions
test_that("cor2vec 2D", {
    rhos <- seq(-0.99, 0.99, 1e-2)
    for (x in rhos) {
        R <- diag(2) * (1 - x) + x
        expect_equal(vec2cor(cor2vec(R)), R)
    }
})

# Test that cor2vec is the inverse of vec2cor in three dimensions
test_that("vec2cor 3D", {
    etas <- as.matrix(do.call(expand.grid, rep(list(seq(-5, 5, 0.5)), 3)))
    colnames(etas) <- NULL
    for (i in seq_len(dim(etas)[1])) {
        x <- etas[i, ]
        expect_equal(cor2vec(vec2cor(x)), x)
    }
})

# Test that vec2cor is the inverse of cor2vec in three dimensions
test_that("cor2vec 3D", {
    rhos <- as.matrix(do.call(expand.grid, rep(list(seq(-0.9, 0.9, 0.1)), 3)))
    colnames(rhos) <- NULL
    for (i in seq_len(dim(rhos)[1])) {
        R <- copula::p2P(rhos[i, ])
        R <- Matrix::nearPD(x = R, corr = TRUE, base.matrix = TRUE)$mat
        expect_equal(vec2cor(cor2vec(R)), R)
    }
})

# Test that cor2vec is the inverse of vec2cor in four dimensions
test_that("vec2cor 4D", {
    etas <- as.matrix(do.call(expand.grid, rep(list(c(-2, 0, 2)), 6)))
    colnames(etas) <- NULL
    for (i in seq_len(dim(etas)[1])) {
        x <- etas[i, ]
        expect_equal(cor2vec(vec2cor(x)), x)
    }
})

# Test that vec2cor is the inverse of cor2vec in four dimensions
test_that("cor2vec 4D", {
    rhos <- as.matrix(do.call(expand.grid, rep(list(c(-0.5, 0, 0.5)), 6)))
    colnames(rhos) <- NULL
    for (i in seq_len(dim(rhos)[1])) {
        R <- copula::p2P(rhos[i, ])
        R <- Matrix::nearPD(x = R, corr = TRUE, base.matrix = TRUE)$mat
        expect_equal(vec2cor(cor2vec(R)), R)
    }
})

# Test that cor2vec is the inverse of vec2cor in five dimensions
test_that("vec2cor 5D", {
    etas <- seq(-5, 5, 1e-2)
    for (x in etas) {
        eta <- rep(x, choose(5, 2))
        expect_equal(cor2vec(vec2cor(eta)), eta)
    }
})

# Test that vec2cor is the inverse of cor2vec in five dimensions
test_that("cor2vec 5D", {
    rhos <- seq(-0.25, 0.25, 5e-4)
    for (r in rhos) {
        R <- copula::p2P(rep(r, choose(5, 2)))
        expect_equal(vec2cor(cor2vec(R)), R)
    }
})
