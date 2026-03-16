library(Matrix)

# ── Return structure ─────────────────────────────────────────────────────────

test_that("betaMix returns a list with all expected fields", {
  set.seed(123)
  M <- matrix(rnorm(50 * 20), nrow = 50, ncol = 20)
  res <- betaMix(M, msg = FALSE)

  expect_type(res, "list")
  expect_named(
    res,
    c("angleMat", "z_j", "m0", "p0", "N", "ahat", "bhat",
      "etahat", "bmax", "ppthr", "nodes", "edges", "cnt"),
    ignore.order = TRUE
  )
})

test_that("betaMix reports correct N and nodes", {
  set.seed(1)
  M <- matrix(rnorm(60 * 25), nrow = 60, ncol = 25)
  res <- betaMix(M, msg = FALSE)

  expect_equal(res$N, 60L)
  expect_equal(res$nodes, 25L)
})

# ── Parameter bounds ─────────────────────────────────────────────────────────

test_that("betaMix p0 is in [0, 1]", {
  set.seed(2)
  M <- matrix(rnorm(50 * 15), nrow = 50, ncol = 15)
  res <- betaMix(M, msg = FALSE)
  expect_gte(res$p0, 0)
  expect_lte(res$p0, 1)
})

test_that("betaMix ahat and bhat are positive", {
  set.seed(3)
  M <- matrix(rnorm(50 * 15), nrow = 50, ncol = 15)
  res <- betaMix(M, msg = FALSE)
  expect_gt(res$ahat, 0)
  expect_gt(res$bhat, 0)
})

test_that("betaMix edges is non-negative and at most choose(P,2)", {
  set.seed(4)
  M <- matrix(rnorm(50 * 20), nrow = 50, ncol = 20)
  res <- betaMix(M, msg = FALSE)
  expect_gte(res$edges, 0)
  expect_lte(res$edges, choose(20L, 2L))
})

test_that("betaMix ppthr is in (0, 1)", {
  set.seed(5)
  M <- matrix(rnorm(50 * 15), nrow = 50, ncol = 15)
  res <- betaMix(M, msg = FALSE)
  expect_gt(res$ppthr, 0)
  expect_lt(res$ppthr, 1)
})

test_that("betaMix posterior probabilities m0 are in [0, 1]", {
  set.seed(6)
  M <- matrix(rnorm(50 * 20), nrow = 50, ncol = 20)
  res <- betaMix(M, msg = FALSE)
  expect_true(all(res$m0 >= 0))
  expect_true(all(res$m0 <= 1))
})

# ── z_j range ────────────────────────────────────────────────────────────────

test_that("betaMix z_j values are strictly in (0, 1)", {
  set.seed(7)
  M <- matrix(rnorm(80 * 30), nrow = 80, ncol = 30)
  res <- betaMix(M, msg = FALSE)
  expect_true(all(res$z_j > 0))
  expect_true(all(res$z_j < 1))
})

# ── ind=FALSE ────────────────────────────────────────────────────────────────

test_that("betaMix runs without error when ind=FALSE", {
  set.seed(8)
  M <- matrix(rnorm(50 * 12), nrow = 50, ncol = 12)
  expect_no_error(betaMix(M, ind = FALSE, msg = FALSE))
})

test_that("betaMix ind=FALSE etahat differs from (N-1)/2 (or converges reasonably)", {
  set.seed(9)
  M <- matrix(rnorm(50 * 12), nrow = 50, ncol = 12)
  res <- betaMix(M, ind = FALSE, msg = FALSE)
  # etahat must be > 0 and <= (N-1)/2
  expect_gt(res$etahat, 0)
  expect_lte(res$etahat, (50 - 1) / 2)
})

# ── Hub structure ─────────────────────────────────────────────────────────────

test_that("betaMix detects edges in a dataset with hub structure", {
  set.seed(42)
  n <- 100
  p <- 20
  hub <- rnorm(n)
  M <- matrix(rnorm(n * p), nrow = n, ncol = p)
  # Strong positive correlation among first 5 columns
  M[, 1:5] <- M[, 1:5] + 3 * hub
  res <- betaMix(M, msg = FALSE, maxalpha = 0.05, ppr = 0.2)
  expect_gt(res$edges, 0)
})

# ── bmax parameter ───────────────────────────────────────────────────────────

test_that("betaMix respects custom bmax", {
  set.seed(10)
  M <- matrix(rnorm(50 * 15), nrow = 50, ncol = 15)
  res <- betaMix(M, bmax = 0.95, msg = FALSE)
  expect_equal(res$bmax, 0.95)
})

# ── angleMat dimensions ──────────────────────────────────────────────────────

test_that("betaMix angleMat is a P x P matrix", {
  set.seed(11)
  M <- matrix(rnorm(50 * 18), nrow = 50, ncol = 18)
  res <- betaMix(M, msg = FALSE)
  expect_equal(dim(res$angleMat), c(18L, 18L))
})
