library(Matrix)

# Helper: base matrix with normal random data
make_base_matrix <- function(n = 60, p = 12, seed = 777) {
  set.seed(seed)
  matrix(rnorm(n * p), nrow = n, ncol = p)
}

# ── (a) All-zeros column ──────────────────────────────────────────────────────

test_that("betaMix warns about zero-variance column (all zeros)", {
  M <- make_base_matrix()
  M[, 3] <- 0
  expect_warning(
    betaMix(M, msg = FALSE),
    regexp = "zero-variance.*column|column.*zero-variance",
    ignore.case = TRUE
  )
})

test_that("betaMix with all-zeros column warns about column index", {
  M <- make_base_matrix()
  M[, 3] <- 0
  expect_warning(
    betaMix(M, msg = FALSE),
    regexp = "\\b3\\b"
  )
})

test_that("betaMix with all-zeros column produces no NaN in m0", {
  M <- make_base_matrix()
  M[, 3] <- 0
  res <- suppressWarnings(betaMix(M, msg = FALSE))
  expect_false(any(is.nan(res$m0)))
})

test_that("betaMix with all-zeros column produces finite ppthr", {
  M <- make_base_matrix()
  M[, 3] <- 0
  res <- suppressWarnings(betaMix(M, msg = FALSE))
  expect_true(is.finite(res$ppthr))
})

test_that("betaMix with all-zeros column: constant column has no edges", {
  M <- make_base_matrix()
  M[, 3] <- 0
  res <- suppressWarnings(betaMix(M, msg = FALSE))
  # ppthr=1 forces all non-zero correlations to be edges; column 3 was zeroed
  adj <- getAdjMat(res, ppthr = 1)
  expect_true(all(adj[3, ] == 0))
  expect_true(all(adj[, 3] == 0))
})

# ── (b) All-ones column ───────────────────────────────────────────────────────

test_that("betaMix warns about zero-variance column (all ones)", {
  M <- make_base_matrix()
  M[, 5] <- 1
  expect_warning(
    betaMix(M, msg = FALSE),
    regexp = "zero-variance.*column|column.*zero-variance",
    ignore.case = TRUE
  )
})

test_that("betaMix with all-ones column warns about column index", {
  M <- make_base_matrix()
  M[, 5] <- 1
  expect_warning(
    betaMix(M, msg = FALSE),
    regexp = "\\b5\\b"
  )
})

test_that("betaMix with all-ones column produces no NaN in m0", {
  M <- make_base_matrix()
  M[, 5] <- 1
  res <- suppressWarnings(betaMix(M, msg = FALSE))
  expect_false(any(is.nan(res$m0)))
})

test_that("betaMix with all-ones column: constant column has no edges", {
  M <- make_base_matrix()
  M[, 5] <- 1
  res <- suppressWarnings(betaMix(M, msg = FALSE))
  adj <- getAdjMat(res, ppthr = 1)
  expect_true(all(adj[5, ] == 0))
  expect_true(all(adj[, 5] == 0))
})

# ── (c) Perfectly correlated columns ─────────────────────────────────────────

test_that("betaMix runs without error with perfectly correlated columns", {
  M <- make_base_matrix()
  M[, 7] <- M[, 2]  # identical columns -> cor = 1, z_j clamped to calcAcc
  expect_no_error(suppressWarnings(betaMix(M, msg = FALSE)))
})

test_that("betaMix m0 is NaN-free with perfectly correlated columns", {
  M <- make_base_matrix()
  M[, 7] <- M[, 2]
  res <- suppressWarnings(betaMix(M, msg = FALSE))
  expect_false(any(is.nan(res$m0)))
})

test_that("betaMix m0 values are in [0, 1] with perfectly correlated columns", {
  M <- make_base_matrix()
  M[, 7] <- M[, 2]
  res <- suppressWarnings(betaMix(M, msg = FALSE))
  expect_true(all(res$m0 >= 0 & res$m0 <= 1))
})

test_that("betaMix detects perfectly correlated pair as an edge", {
  M <- make_base_matrix()
  M[, 7] <- M[, 2]
  res <- suppressWarnings(betaMix(M, msg = FALSE, maxalpha = 0.05, ppr = 0.2))
  adj <- getAdjMat(res)
  expect_true(as.logical(adj[2, 7]) || as.logical(adj[7, 2]))
})

# ── Combined: all three degenerate cases together ─────────────────────────────

test_that("betaMix handles all three degenerate cases simultaneously", {
  M <- make_base_matrix(n = 80, p = 15)
  M[, 1]  <- 0        # all zeros
  M[, 2]  <- 1        # all ones
  M[, 10] <- M[, 5]   # perfectly correlated pair

  # warns about the two constant columns
  expect_warning(
    res <- betaMix(M, msg = FALSE),
    regexp = "zero-variance",
    ignore.case = TRUE
  )
  expect_false(any(is.nan(res$m0)))
  expect_true(is.finite(res$ppthr))

  adj <- getAdjMat(res, ppthr = 1)
  # constant columns have no edges
  expect_true(all(adj[1, ] == 0))
  expect_true(all(adj[2, ] == 0))
  # perfectly correlated pair is connected
  expect_true(as.logical(adj[5, 10]) || as.logical(adj[10, 5]))
})
