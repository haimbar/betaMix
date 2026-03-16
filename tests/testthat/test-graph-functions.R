library(Matrix)

# Helper: build a dense symmetric 0/1 Matrix with given density
make_adj_matrix <- function(n, density = 0.15, seed = 1) {
  set.seed(seed)
  A <- matrix(0L, n, n)
  upper <- which(upper.tri(A))
  A[upper[sample(length(upper), round(density * length(upper)))]] <- 1L
  A <- A + t(A)
  Matrix(A)
}

# ── clusteringCoef ────────────────────────────────────────────────────────────

test_that("clusteringCoef returns 1 for every node in a complete graph", {
  A <- Matrix(1L, 4L, 4L)
  diag(A) <- 0L
  cc <- clusteringCoef(A)
  expect_equal(length(cc), 4L)
  expect_true(all(abs(cc - 1) < 1e-10))
})

test_that("clusteringCoef returns 0 for all nodes of a star graph", {
  # Node 1 = center, nodes 2-5 = leaves; no leaf-to-leaf edges
  A <- Matrix(0L, 5L, 5L)
  A[1, 2:5] <- 1L
  A[2:5, 1] <- 1L
  cc <- clusteringCoef(A)
  expect_equal(length(cc), 5L)
  expect_true(all(cc == 0))
})

test_that("clusteringCoef values are in [0, 1]", {
  A <- make_adj_matrix(30, density = 0.2)
  cc <- clusteringCoef(A)
  expect_true(all(cc >= 0))
  expect_true(all(cc <= 1))
})

test_that("clusteringCoef returns 0 for an isolated node", {
  A <- Matrix(0L, 4L, 4L)
  A[1, 2] <- A[2, 1] <- 1L   # edge only between 1 and 2
  cc <- clusteringCoef(A)
  # nodes 3 and 4 are isolated
  expect_equal(cc[3], 0)
  expect_equal(cc[4], 0)
})

test_that("clusteringCoef returns a vector of length nrow(A)", {
  A <- make_adj_matrix(20)
  cc <- clusteringCoef(A)
  expect_equal(length(cc), 20L)
})

# ── graphComponents ────────────────────────────────────────────────────────────

test_that("graphComponents returns a data frame with required columns", {
  A <- Matrix(1L, 6L, 6L)
  diag(A) <- 0L
  result <- graphComponents(A, minCtr = 1)
  expected_cols <- c("labels", "degree", "cc", "ctr", "clustNo",
                     "iscenter", "intEdges", "extEdges", "distCenter")
  expect_true(all(expected_cols %in% colnames(result)))
  expect_equal(nrow(result), 6L)
})

test_that("graphComponents requires a Matrix (not a base matrix)", {
  M <- matrix(1L, 5L, 5L)
  expect_error(graphComponents(M))
})

test_that("graphComponents clustNo is non-negative integer", {
  A <- make_adj_matrix(20, density = 0.3)
  result <- graphComponents(A, minCtr = 2)
  expect_true(all(result$clustNo >= 0))
})

test_that("graphComponents degree column matches Matrix::rowSums", {
  A <- make_adj_matrix(15, density = 0.2)
  result <- graphComponents(A)
  expect_equal(result$degree, as.numeric(Matrix::rowSums(A)))
})

test_that("graphComponents cc column matches clusteringCoef", {
  A <- make_adj_matrix(12, density = 0.25)
  cc_ref <- clusteringCoef(abs(A))
  result  <- graphComponents(A)
  expect_equal(result$cc, cc_ref)
})

test_that("graphComponents type=0 centrality equals degree", {
  A <- make_adj_matrix(10, density = 0.3)
  result <- graphComponents(A, type = 0)
  expect_equal(result$ctr, result$degree)
})

# ── sphericalCaps ─────────────────────────────────────────────────────────────

test_that("sphericalCaps returns NULL for a graph with no edges", {
  A <- Matrix(0L, 5L, 5L)
  expect_null(sphericalCaps(A))
})

test_that("sphericalCaps requires a Matrix (not a base matrix)", {
  A <- matrix(0L, 5L, 5L)
  expect_error(sphericalCaps(A))
})

test_that("sphericalCaps returns a data frame with required columns", {
  A <- Matrix(1L, 5L, 5L)
  diag(A) <- 0L
  result <- sphericalCaps(A)
  expect_true(is.data.frame(result))
  expect_true(all(c("node", "capNum", "isCtr", "deg", "cc") %in% colnames(result)))
})

test_that("sphericalCaps capNum values are positive integers", {
  A <- make_adj_matrix(15, density = 0.3)
  result <- sphericalCaps(A)
  if (!is.null(result)) {
    expect_true(all(result$capNum >= 1L))
  }
})

test_that("sphericalCaps isCtr is 0 or 1", {
  A <- Matrix(1L, 6L, 6L)
  diag(A) <- 0L
  result <- sphericalCaps(A)
  expect_true(all(result$isCtr %in% c(0, 1)))
})

# ── summarizeClusters ─────────────────────────────────────────────────────────

test_that("summarizeClusters returns a matrix with 12 columns", {
  A <- make_adj_matrix(20, density = 0.3, seed = 7)
  ci <- graphComponents(A, minCtr = 2)
  if (max(ci$clustNo) > 0) {
    result <- summarizeClusters(ci)
    expect_true(is.matrix(result))
    expect_equal(ncol(result), 12L)
  }
})

test_that("summarizeClusters returns NULL when no clusters", {
  A <- Matrix(0L, 5L, 5L)
  ci <- graphComponents(A, minCtr = 1)
  result <- summarizeClusters(ci)
  expect_null(result)
})

# ── collapsedGraph ────────────────────────────────────────────────────────────

test_that("collapsedGraph returns a symmetric Matrix", {
  A <- make_adj_matrix(20, density = 0.3, seed = 11)
  ci <- graphComponents(A, minCtr = 2)
  if (max(ci$clustNo) > 0) {
    cg <- collapsedGraph(A, ci)
    expect_true(is(cg, "Matrix"))
    cg_dense <- as.matrix(cg)
    expect_equal(cg_dense, t(cg_dense))
  }
})
