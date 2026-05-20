# Test formula parsing
test_that("parse_formula", {
    compare_formulas <- function(f, r1, r2, t) {
        res <- .parse_formula(f)
        f1 <- res$linear
        f2 <- res$smooth_linear_int

        get_terms <- function(f) {
            if (is.null(f)) {
                return(NULL)
            }
            return(sort(attr(terms(f), "term.labels")))
        }

        testthat::expect_true(identical(get_terms(f1), get_terms(r1)))
        testthat::expect_true(identical(get_terms(f2), get_terms(r2)))
        testthat::expect_identical(t, res$smooth_name)
    }

    # Only linear covariates
    compare_formulas(~1, ~1, NULL, NULL)
    compare_formulas(~x, ~x, NULL, NULL)
    compare_formulas(~x*y, ~x+y+x:y, NULL, NULL)
    compare_formulas(~x+y+z+t, ~x+y+z+t, NULL, NULL)

    # Smooth covariate with no interactions
    compare_formulas(~s(t), ~1, ~1, "t")
    compare_formulas(~s(t) + x, ~x, ~1, "t")
    compare_formulas(~x + s(t), ~x, ~1, "t")
    compare_formulas(~s(time), ~1, ~1, "time")
    compare_formulas(~x+y+z+t+s(t), ~x+y+z+t, ~1, "t")
    compare_formulas(~x*y+z+t+s(t), ~x+y+x:y+z+t, ~1, "t")
    compare_formulas(~x*y*t+s(t), ~x*y*t, ~1, "t")
    compare_formulas(~x*y*t+s(x), ~x*y*t, ~1, "x")

    # Interactions with smooth covariate
    compare_formulas(~x:s(t) + s(t), ~1, ~x, "t")
    compare_formulas(~x:s(t) + y:s(t) + s(t), ~1, ~x + y, "t")
    compare_formulas(~x:s(t) + y:s(t), ~1, ~x + y, "t")

    # Specify multiple smooth covariates (not allowed)
    testthat::expect_error(.parse_formula(~s(t1) + s(t2)))
    testthat::expect_error(.parse_formula(~s(x) * s(y)))
    testthat::expect_error(.parse_formula(~s(xy) + s(x)))
})
