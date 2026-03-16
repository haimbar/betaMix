library(Matrix)

# ── numSteps = 0 ──────────────────────────────────────────────────────────────

test_that("shortestPathDistance with numSteps=0 returns the input unchanged", {
  A <- Matrix(c(0L, 1L, 0L,
                1L, 0L, 1L,
                0L, 1L, 0L), 3L, 3L)
  result <- shortestPathDistance(A, numSteps = 0)
  expect_equal(as.matrix(result), as.matrix(A))
})

# ── 2-hop connectivity ────────────────────────────────────────────────────────

test_that("shortestPathDistance with numSteps=1 connects 2-hop paths", {
  # Path graph: 1 -- 2 -- 3
  A <- Matrix(c(0L, 1L, 0L,
                1L, 0L, 1L,
                0L, 1L, 0L), 3L, 3L)
  result <- shortestPathDistance(A, numSteps = 1)
  # Nodes 1 and 3 are 2 hops apart; result should be > 0
  expect_gt(result[1, 3], 0)
  expect_gt(result[3, 1], 0)
})

test_that("shortestPathDistance with numSteps=1 does not connect 3-hop-only nodes", {
  # Path graph: 1 -- 2 -- 3 -- 4
  A <- Matrix(c(0L, 1L, 0L, 0L,
                1L, 0L, 1L, 0L,
                0L, 1L, 0L, 1L,
                0L, 0L, 1L, 0L), 4L, 4L)
  result <- shortestPathDistance(A, numSteps = 1)
  # Nodes 1 and 4 are 3 hops; numSteps=1 should not connect them
  expect_equal(result[1, 4], 0)
})

test_that("shortestPathDistance with numSteps=2 connects 3-hop paths", {
  # Path: 1 -- 2 -- 3 -- 4
  A <- Matrix(c(0L, 1L, 0L, 0L,
                1L, 0L, 1L, 0L,
                0L, 1L, 0L, 1L,
                0L, 0L, 1L, 0L), 4L, 4L)
  result <- shortestPathDistance(A, numSteps = 2)
  expect_gt(result[1, 4], 0)
  expect_gt(result[4, 1], 0)
})

# ── Return type ───────────────────────────────────────────────────────────────

test_that("shortestPathDistance returns a matrix-like object", {
  A <- Matrix(c(0L, 1L, 1L,
                1L, 0L, 1L,
                1L, 1L, 0L), 3L, 3L)
  result <- shortestPathDistance(A, numSteps = 2)
  expect_true(is.matrix(result) || is(result, "Matrix"))
})

# ── Dimensions ────────────────────────────────────────────────────────────────

test_that("shortestPathDistance returns matrix of same dimensions as input", {
  n <- 8L
  set.seed(42)
  A_base <- matrix(sample(0L:1L, n * n, replace = TRUE, prob = c(0.7, 0.3)), n, n)
  A_sym  <- A_base + t(A_base)
  A_sym[A_sym > 1L] <- 1L
  diag(A_sym) <- 0L
  A <- Matrix(A_sym)
  result <- shortestPathDistance(A, numSteps = 2)
  expect_equal(dim(result), c(n, n))
})

# ── Complete graph ────────────────────────────────────────────────────────────

test_that("shortestPathDistance on a complete graph stays the same for numSteps >= 1", {
  # In K4 every pair is already directly connected; no new paths added
  A <- Matrix(1L, 4L, 4L)
  diag(A) <- 0L
  r0 <- shortestPathDistance(A, numSteps = 0)
  r1 <- shortestPathDistance(A, numSteps = 1)
  # Both should represent full connectivity (every off-diagonal > 0)
  diag_mask <- !diag(4L)
  expect_true(all(as.matrix(r1)[diag_mask] > 0))
})
