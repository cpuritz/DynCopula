# # Test argument checks for fit_dyn_corr
# test_that("fit_dyn_corr arguments", {
#     sce <- SingleCellExperiment::SingleCellExperiment()
#     metadata(sce)$dyn_corr <- list(metacell_sce = list())
#     bandwidth <- 0.1
#     control <- list()
#     ncv <- NULL
#     metadata(sce)$dyn_corr <- list()
#
#     # Test sce
#     expect_error(fit_dyn_corr("sce", bandwidth))
#     expect_error(fit_dyn_corr(5, bandwidth))
#     expect_error(fit_dyn_corr(NULL, bandwidth))
#
#     # Test bandwidth
#     expect_error(fit_dyn_corr(sce, -0.1))
#     expect_error(fit_dyn_corr(sce, 0))
#     expect_error(fit_dyn_corr(sce, "a"))
#
#     # Test control
#     expect_error(fit_dyn_corr(sce, bandwidth, control = c()))
#     expect_error(fit_dyn_corr(sce, bandwidth, control = NULL))
#
#     # Test ncv
#     expect_error(fit_dyn_corr(sce, bandwidth, ncv = -1))
#     expect_error(fit_dyn_corr(sce, bandwidth, ncv = 0))
#     expect_error(fit_dyn_corr(sce, bandwidth, ncv = 1))
#     expect_error(fit_dyn_corr(sce, bandwidth, ncv = "a"))
#     expect_error(fit_dyn_corr(sce, bandwidth, ncv = NA))
#
#     # Test sce metadata
#     metadata(sce) <- list()
#     expect_error(fit_dyn_corr(sce, bandwidth))
#     metadata(sce)$dyn_corr <- list()
#     expect_error(fit_dyn_corr(sce, bandwidth))
# })
