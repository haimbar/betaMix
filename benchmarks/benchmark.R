## betaMix performance benchmarks
## Run from the package root:  Rscript benchmarks/benchmark.R

# ── Load package ─────────────────────────────────────────────────────────────
if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".", quiet = TRUE)
} else {
  library(betaMix)
}
library(Matrix)

# ── Benchmark helper ──────────────────────────────────────────────────────────
# Accepts a zero-argument function; repeats 'times' calls and reports
# median/min/max elapsed time per call.
bench <- function(fun, times = 10, label = "") {
  gc(verbose = FALSE)
  # One warm-up call (not counted)
  fun()
  elapsed <- numeric(times)
  for (i in seq_len(times)) {
    t0 <- proc.time()[["elapsed"]]
    fun()
    elapsed[i] <- proc.time()[["elapsed"]] - t0
  }
  # For very fast ops proc.time() resolution is ~10ms; use a bulk run instead.
  if (max(elapsed) == 0) {
    bulk <- 200L
    t0   <- proc.time()[["elapsed"]]
    for (i in seq_len(bulk)) fun()
    elapsed[] <- (proc.time()[["elapsed"]] - t0) / bulk
  }
  cat(sprintf("  %-52s  median %8.4fs  [min %8.4fs  max %8.4fs]\n",
              label,
              median(elapsed), min(elapsed), max(elapsed)))
  invisible(elapsed)
}

sep <- function(title)
  cat(sprintf("\n── %s %s\n", title, strrep("─", max(0, 62 - nchar(title)))))

# ── Data preparation ──────────────────────────────────────────────────────────
sep("Loading datasets")
data(DrySeeds, package = "betaMix")    # 50 rows x  68 cols
data(SIM,      package = "betaMix")    # 200 rows x 1000 cols

set.seed(42)
SIM_small <- SIM[, 1:100]             # 200 rows x 100 cols

cat(sprintf("  DrySeeds  : %d rows x %d cols  (%d pairs)\n",
            nrow(DrySeeds), ncol(DrySeeds), choose(ncol(DrySeeds), 2)))
cat(sprintf("  SIM_small : %d rows x %d cols  (%d pairs)\n",
            nrow(SIM_small), ncol(SIM_small), choose(ncol(SIM_small), 2)))
cat(sprintf("  SIM (full): %d rows x %d cols  (%d pairs, subsampled)\n",
            nrow(SIM), ncol(SIM), choose(ncol(SIM), 2)))

# Helper: random symmetric 0/1 sparse adjacency Matrix with given density
make_adj <- function(n, density = 0.10, seed = 1) {
  set.seed(seed)
  A <- matrix(0L, n, n)
  ui <- which(upper.tri(A))
  A[sample(ui, round(density * length(ui)))] <- 1L
  A <- A + t(A)
  diag(A) <- 0L
  Matrix(A)
}

A50  <- make_adj(50)
A100 <- make_adj(100)
A200 <- make_adj(200)

# ── 1. betaMix ────────────────────────────────────────────────────────────────
sep("1. betaMix() -- EM model fitting")
bench(function() betaMix(DrySeeds,  msg = FALSE, ind = TRUE),
      times = 5,  label = "DrySeeds  (68 nodes, 50 samples,   2278 pairs)")
bench(function() betaMix(SIM_small, msg = FALSE, ind = TRUE),
      times = 3,  label = "SIM_small (100 nodes, 200 samples, 4950 pairs)")
bench(function() betaMix(SIM, msg = FALSE, maxalpha = 1e-6, ppr = 0.01,
                          ind = TRUE),
      times = 3,  label = "SIM full  (1000 nodes, 200 samples, subsampled)")

# ── Pre-compute results for downstream benchmarks ────────────────────────────
cat("\n  [pre-computing fitted models...]\n")
res_dry <- betaMix(DrySeeds,  msg = FALSE)
res_sim <- betaMix(SIM, msg = FALSE, maxalpha = 1e-6, ppr = 0.01)

adj_dry <- getAdjMat(res_dry)
adj_sim <- getAdjMat(res_sim)

# ── 2. getAdjMat ──────────────────────────────────────────────────────────────
sep("2. getAdjMat() -- adjacency matrix extraction")
bench(function() getAdjMat(res_dry),                 times = 50, label = "DrySeeds unsigned")
bench(function() getAdjMat(res_dry, signed = TRUE),  times = 50, label = "DrySeeds signed")
bench(function() getAdjMat(res_sim),                 times = 20, label = "SIM unsigned")
bench(function() getAdjMat(res_sim, signed = TRUE),  times = 20, label = "SIM signed")

# ── 3. clusteringCoef ─────────────────────────────────────────────────────────
sep("3. clusteringCoef() -- clustering coefficients")
bench(function() clusteringCoef(A50),     times = 20, label = "random 50-node graph")
bench(function() clusteringCoef(A100),    times = 10, label = "random 100-node graph")
bench(function() clusteringCoef(A200),    times = 10, label = "random 200-node graph")
bench(function() clusteringCoef(adj_dry), times = 20, label = "DrySeeds adjacency (68 nodes)")
bench(function() clusteringCoef(adj_sim), times = 10, label = "SIM adjacency (1000 nodes)")

# ── 4. graphComponents ────────────────────────────────────────────────────────
sep("4. graphComponents() -- cluster detection")
bench(function() graphComponents(A50,    minCtr = 3), times = 10, label = "random 50-node graph")
bench(function() graphComponents(A100,   minCtr = 3), times = 10, label = "random 100-node graph")
bench(function() graphComponents(A200,   minCtr = 3), times =  5, label = "random 200-node graph")
bench(function() graphComponents(adj_dry),             times = 10, label = "DrySeeds adjacency")
bench(function() graphComponents(adj_sim),             times =  3, label = "SIM adjacency (1000 nodes)")

# ── 5. shortestPathDistance ───────────────────────────────────────────────────
sep("5. shortestPathDistance() -- path distances")
bench(function() shortestPathDistance(A50,  numSteps = 1), times = 20, label = "n=50  steps=1")
bench(function() shortestPathDistance(A50,  numSteps = 2), times = 20, label = "n=50  steps=2")
bench(function() shortestPathDistance(A50,  numSteps = 3), times = 20, label = "n=50  steps=3")
bench(function() shortestPathDistance(A100, numSteps = 1), times = 10, label = "n=100 steps=1")
bench(function() shortestPathDistance(A100, numSteps = 2), times = 10, label = "n=100 steps=2")
bench(function() shortestPathDistance(A200, numSteps = 1), times =  5, label = "n=200 steps=1")

# ── 6. sphericalCaps ──────────────────────────────────────────────────────────
sep("6. sphericalCaps() -- hub detection")
bench(function() sphericalCaps(adj_dry), times = 20, label = "DrySeeds adjacency (68 nodes)")
bench(function() sphericalCaps(adj_sim), times =  5, label = "SIM adjacency (1000 nodes)")

# ── 7. summarizeClusters ──────────────────────────────────────────────────────
sep("7. summarizeClusters() -- cluster summary statistics")
ci_dry <- graphComponents(adj_dry)
ci_sim <- graphComponents(adj_sim)
bench(function() capture.output(summarizeClusters(ci_dry)), times = 50, label = "DrySeeds clusters")
bench(function() capture.output(summarizeClusters(ci_sim)), times = 20, label = "SIM clusters")

sep("Benchmark complete")
cat(sprintf("  R version: %s\n",   R.version$version.string))
cat(sprintf("  Platform : %s\n\n", R.version$platform))
