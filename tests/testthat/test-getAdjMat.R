library(Matrix)

# Shared fixture: small betaMix result
local_res <- local({
  set.seed(101)
  M <- matrix(rnorm(60 * 18), nrow = 60, ncol = 18)
  betaMix(M, msg = FALSE)
})

# ── Return type ───────────────────────────────────────────────────────────────

test_that("getAdjMat returns a matrix-like object", {
  adj <- getAdjMat(local_res)
  expect_true(is.matrix(adj) || is(adj, "Matrix"))
})

test_that("getAdjMat has correct dimensions", {
  adj <- getAdjMat(local_res)
  expect_equal(nrow(adj), local_res$nodes)
  expect_equal(ncol(adj), local_res$nodes)
})

# ── Diagonal ─────────────────────────────────────────────────────────────────

test_that("getAdjMat unsigned has zero/FALSE diagonal", {
  adj <- getAdjMat(local_res)
  expect_true(all(!diag(adj)))
})

test_that("getAdjMat signed has zero/FALSE diagonal", {
  adj <- getAdjMat(local_res, signed = TRUE)
  expect_true(all(!diag(adj)))
})

# ── Symmetry ─────────────────────────────────────────────────────────────────

test_that("getAdjMat unsigned is symmetric", {
  adj <- getAdjMat(local_res)
  adj_m <- as.matrix(adj)
  expect_equal(adj_m, t(adj_m))
})

# ── Signed values ─────────────────────────────────────────────────────────────

test_that("getAdjMat signed contains only -1, 0, or 1 (or logicals)", {
  adj <- getAdjMat(local_res, signed = TRUE)
  vals <- unique(as.vector(adj))
  expect_true(all(vals %in% c(-1L, 0L, 1L, FALSE, TRUE)))
})

# ── Custom ppthr ──────────────────────────────────────────────────────────────

test_that("getAdjMat with ppthr=1 returns all-edges matrix (except diagonal)", {
  adj_all <- getAdjMat(local_res, ppthr = 1)
  # All off-diagonal entries should be edges
  diag_mask <- !diag(nrow(adj_all))
  expect_true(all(adj_all[diag_mask]))
})

test_that("getAdjMat with ppthr=0 returns no-edges matrix", {
  adj_none <- getAdjMat(local_res, ppthr = 0)
  expect_true(all(!adj_none))
})

# ── Node subset ───────────────────────────────────────────────────────────────

test_that("getAdjMat nodes subset returns matrix with at least those nodes", {
  nodes_sel <- 1:3
  adj_sub <- getAdjMat(local_res, nodes = nodes_sel)
  expect_gte(nrow(adj_sub), length(nodes_sel))
})
