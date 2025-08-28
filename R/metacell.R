###############################################################################

#' Bin a vector
#'
#' @description Construct a sequence of bins for a vector with a maximum bin
#' width and a desired number of points per bin.
#'
#' @param x Vector of points.
#' @param N Target number of points per bin.
#' @param max_width Maximum allowable bin width as a fraction of the range of
#' \code{x}.
#'
#' @details \code{max_width} is a strict upper bound on the width of each bin.
#' Bins may end up having more or less than \code{N} points.
#'
.bin_vector <- function(x, N, max_width) {
    N <- as.integer(N)
    n <- length(x)

    # Standardize to [0, 1]
    x <- sort(x)
    min_x <- x[1]
    max_x <- x[length(x)]
    x <- (x - min_x) / (max_x - min_x)

    cuts <- c(0)
    bins <- list()
    i <- 1L

    while (i <= n) {
        # Number of points left to allocate
        remaining <- n - i + 1L
        # Expected number of bins left to create
        bins_left <- max(1L, ceiling(remaining / N))

        # Target number of points to allocate to the bin
        target <- as.integer(remaining / bins_left)
        lo <- max(1L, as.integer(floor(0.5 * N)))
        hi <- as.integer(ceiling(1.5 * N))
        target <- max(lo, min(hi, target))
        # Leave at least one point per future bin
        target <- min(target, remaining - (bins_left - 1L))

        # Index that satisfies the width constraint but allocates as close to
        # 'target' points to the bin
        ix_by_count <- min(n, i + target - 1L)
        ix_by_width <- utils::tail(seq_along(x)[(x - x[i]) <= max_width], 1)
        end_ix <- max(min(ix_by_count, ix_by_width), i)

        # Keep enough points for future bins
        min_needed <- (bins_left - 1L) * lo
        while ((n - end_ix) < min_needed && end_ix > i) {
            end_ix <- end_ix - 1L
            if (x[end_ix] - x[i] > max_width) {
                break
            }
        }

        # Shrink bin if it is wider than max_width
        while (x[end_ix] - x[i] > max_width && end_ix > i) {
            end_ix <- end_ix - 1L
        }

        bins[[length(bins) + 1L]] <- c(i, end_ix)
        cuts <- c(cuts, ifelse(end_ix == n, 1, x[end_ix]))

        if (end_ix == n) {
            break
        } else {
            i <- end_ix + 1L
        }
    }

    # Ensure last cut-point is exactly 1
    cuts[length(cuts)] <- 1

    # Rescale cut points to original length scale
    cuts <- (max_x - min_x) * cuts + min_x

    # Counts per bin
    counts <- sapply(bins, function(b) { b[2] - b[1] + 1 })

    return(list(cut_points = cuts, counts = counts))
}

###############################################################################

#' Generate metacells
#'
#' @description Generate metacells by binning cells based on their pseudotimes.
#'
#' @param sce A SingleCellExperiment.
#' @param col The name of the \code{colData} column containing pseudotimes.
#' @param N Desired number of cells per metacell.
#' @param max_width Maximum difference in pseudotimes allowed for cells assigned
#' to the same metacell. Expressed as a fraction of the total range of
#' pseudotimes.
#' @param agg Method for assigning pseudotimes to metacells. Default is
#' \code{"mean"}.
#'
#' @details Cells are first binned by their pseudotimes. \code{N} sets the
#' target number of cells per bin, but the number of cells per bin may be
#' smaller or larger. \code{max_width} sets an upper bound on the width of each
#' bin. This bound may be violated if there are a large number of ties in
#' pseudotimes. Each bin forms a single metacell. The raw transcript counts for
#' cells in a single bin are summed to form the metacell's counts.
#'
#' @export
generate_metacells <- function(sce,
                               col,
                               N,
                               max_width,
                               agg = c("mean", "min", "max")) {
    assertthat::assert_that(
        methods::is(sce, "SingleCellExperiment"),
        is.character(col),
        col %in% colnames(SingleCellExperiment::colData(sce)),
        is.numeric(N) && N > 0,
        is.numeric(max_width) && 0 < max_width && max_width <= 1
    )
    N <- as.integer(N)

    agg <- match.arg(agg)
    agg <- methods::getFunction(agg, where = getNamespace("base"))

    # Get cut points for binning
    times <- sce[[col]]
    bins <- .bin_vector(unique(times), N, max_width)
    cut_pts <- bins$cut_points

    # Bin pseudotimes
    intervals <- cut(
        x = times,
        breaks = bins$cut_points,
        include.lowest = TRUE,
        right = TRUE,
        labels = FALSE
    )

    # Build metacell count matrix
    uint <- seq_len(max(intervals))
    sce_counts <- SingleCellExperiment::counts(sce)
    mc_counts <- sapply(uint, function(i) {
        Matrix::rowSums(sce_counts[, intervals == i, drop = FALSE])
    })
    mc_counts <- methods::as(mc_counts, "dgCMatrix")

    # Construct SingleCellExperiment object for metacells
    mc_sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = mc_counts)
    )
    colnames(mc_sce) <- paste("metacell", uint, sep = '_')

    return(list(
        sce_metacell = mc_sce,
        metacell_assignment = stats::setNames(intervals, colnames(sce))
    ))
}

###############################################################################
