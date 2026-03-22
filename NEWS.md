# betaMix 0.3.2

## Bug fixes

- `shortestPathDistance()`: fixed two correctness bugs in the previous
  implementation.  (1) Off-by-one: k-hop paths were labelled with distance
  k−1 instead of k, so 2-hop paths and direct edges were indistinguishable.
  (2) Diagonal: A² has nonzero diagonals (every node has a 2-step path back
  to itself), causing self-distances to be set to 1 instead of 0.

## Performance

- `shortestPathDistance()`: rewrote the matrix-power loop to work with binary
  (logical) sparse matrices throughout.  The previous implementation
  accumulated integer path counts, causing rapid fill-in — a P=1000 network
  at 1 % density goes from 1 % dense (A¹) to 60 % dense (A³) to 100 % dense
  (A⁴).  The new frontier-BFS approach casts to binary after each step,
  keeping the matrices sparse and giving a 2–6× speed-up depending on P and
  numSteps.

---

# betaMix 0.3.1

## Bug fixes

- Fixed `jacmle()`: all four entries of the 2×2 Jacobian matrix had the wrong
  sign, causing `nleqslv` to fail with `termcd = 3` (maximum iterations
  exceeded) on every EM M-step. As a result, `ahat` and `bhat` never moved
  from their initial values (8, 3) in any dataset.
- M-step now only updates `ahat`/`bhat` when `nleqslv` reports successful
  convergence (`termcd` ∈ {1, 2}). Non-converged (`termcd = 3`) and
  step-back (`termcd = 5`) results are silently ignored, preventing the EM
  oscillation cycle that appeared for SixHourImbibed once the Jacobian sign
  was corrected.
- Vignette images and example output regenerated with the corrected estimates.

---

# betaMix 0.3.0

## New features

- `assessFit()`: new exported function producing a 2×2 goodness-of-fit
  diagnostic plot and a printed assessment for a fitted betaMix model.
  - Panel 1: fitted mixture histogram (`plotFittedBetaMix`).
  - Panel 2: Q-Q plot of null-assigned observations (m0 ≥ 0.5) vs
    Beta(η̂, 0.5); points are blue when KS D ≤ 0.05, red when D > 0.05.
  - Panel 3: Q-Q plot of non-null-assigned observations (m0 < 0.5, scaled by
    bmax) vs Beta(â, b̂).
  - Panel 4: kernel density of non-null z_j overlaid with the fitted Beta
    density; orange dotted lines mark extra modes when multimodal.
  - Printed assessment covers KS D for both components, mode count, heavy
    left-tail flag, no-signal flag, and an overall verdict ("Good fit",
    "Acceptable fit", or "Deviations detected"). Contextual hints suggest
    `ind = FALSE` on null deviation and warn of multimodality on non-null
    deviation.
  - Uses KS D > 0.05 as the practical threshold (raw p-values are unreliable
    for large numbers of pairs).
- `plotROC()`: new exported function that computes and plots a model-based ROC
  curve by sweeping the edge-detection threshold τ from 0 to 1.  FPR and TPR
  are derived directly from the fitted Beta distributions — no ground-truth
  labels are required.  Marks the current threshold (orange) and the
  Youden-index optimum (red); returns AUC and `tau_youden` for downstream use
  via `getAdjMat(res, ppthr = roc$tau_youden)`.
- `betaMixReport()`: the single fitted-mixture plot is replaced by a "Model
  Fit Assessment" section (`assessFit()`, 7 × 7 in) followed by a "ROC Curve
  and Threshold Selection" section (`plotROC()`, 4.5 × 4.5 in).
- `betaMixReport()`: added `betaMixReport()` function (introduced in 0.2.11)
  for generating self-contained PDF reports.
- Vignette: added "Goodness-of-fit assessment" and "Model-based ROC curve"
  sections documenting the new functions.

---

# betaMix 0.2.14

## Documentation

- Vignette: added "Automated PDF Report" section with a `betaMixReport()`
  invocation example and a description of its key arguments.

---

# betaMix 0.2.13

## Bug fixes

- `betaMixReport()`: cluster summary table headers no longer spill across
  columns in the PDF. Long column names (`degreeMedian`, `pctInClstMedian`,
  etc.) are replaced by short abbreviations (`Deg. Med.`, `Pct. Med.`, …) and
  `kable_styling(latex_options = "scale_down")` is applied so the table always
  fits within the page width.

---

# betaMix 0.2.12

# betaMix 0.2.11

## New features

- `betaMixReport()`: new exported function that runs the full betaMix analysis
  pipeline on a data matrix and renders a self-contained PDF report.  The
  report includes:
  - Dataset overview (dimensions, variable names, summary statistics table for
    P ≤ 30, or a distributional summary for larger P).
  - Estimated model parameters table (â, b̂, η̂, p₀, τ).
  - Explanation of the posterior probability threshold τ and the equivalent
    minimum |r| required to declare an edge.
  - `plotFittedBetaMix()`: histogram of z_ij with fitted null/non-null
    components and the orange significance region.
  - `plotDegCC()`: scatter of node degree vs. degree × clustering coefficient.
  - `plotBitmapCC()`: adjacency bitmap ordered by cluster membership.
  - `plotCluster()`: network plot for a representative cluster (the largest
    cluster in the 5–30 node range; falls back to the globally largest cluster
    when no cluster satisfies this criterion).
  - A `rep_cluster` argument lets the user override the automatic selection.
  - Output is written to `output_dir` (default `"reports/"`); the directory is
    created automatically if absent.
- Added `inst/report_template.Rmd` (parametrised Rmd used by `betaMixReport`).
- Added `rmarkdown` to `Suggests` in `DESCRIPTION`.
- Added `reports/`, `betaMix_*.pdf`, and related intermediate-file patterns to
  `.gitignore`.

---

# betaMix 0.2.10

## Cross-platform portability fixes in C++ (`writeLargeTable.cpp`)

- `calcCorr()`: changed `fopen()` modes from `"r"`/`"w"` to `"rb"`/`"wb"`.
  The file format uses fixed-width records addressed by arithmetic byte offsets
  via `fseek`. Text mode on Windows expands `\n` to `\r\n`, corrupting the
  offset arithmetic and producing wrong correlations or a crash.
- `calcCorr()`: replaced `fseek(fp, i * recSize, SEEK_SET)` (both calls) with
  a new `pkg_fseek()` helper that uses `_fseeki64` on Windows and `fseeko`
  (with `off_t`) elsewhere. The old `int × int` multiplication overflows for
  `P × recSize > 2 × 10⁹`; `pkg_fseek` computes the offset as `long long`.
- Added `src/Makevars.win` (mirrors `src/Makevars`) so Windows builds receive
  explicit `-D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE` preprocessor flags.
- Added `-D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE` to `src/Makevars`
  (no-op on 64-bit Linux/macOS; ensures `off_t` is 64-bit on 32-bit POSIX).

## Vignette

- Corrected `ind` default description: was "FALSE (the default)", now correctly
  states "TRUE (the default)".
- Fixed copy-paste error: `adjMat2 <- getAdjMat(res1)` → `getAdjMat(res2)`.
- Added section "Using Spearman's rank correlation" with a worked example.
- Minor wording and formatting improvements.

---

# betaMix 0.2.9

## New features

- `betaMix()` gains a `method` argument (`"pearson"` (default) or
  `"spearman"`). Spearman's rank correlation is computed by passing
  `method = "spearman"` to `cor()` in the in-memory path, and by ranking
  each column in-place before standardisation in the C++ (`calcCorr()`)
  file-backed path. Average ranks are used for ties. The chosen method is
  stored in the returned list as `$method`.
- `calcCorr()` gains a `spearman` argument (default `FALSE`). When `TRUE`,
  columns are ranked via a new internal `rankVec()` helper (O(n log n),
  average-rank tie handling) before Pearson standardisation, giving exact
  Spearman correlations.

## Tests

- Added 6 tests for the `method` parameter: default is `"pearson"`, Spearman
  is recorded correctly, valid list structure, `angleMat` differs between
  methods, Spearman detects monotone non-linear relationships that Pearson
  may miss, and invalid method values are rejected.

---

# betaMix 0.2.8

## Performance / memory

- `sphericalCaps()`: replaced incremental `rbind()` inside the while-loop with
  a pre-allocated list of chunks and a single `do.call(rbind, …)` at the end,
  reducing peak memory from O(n²) to O(n) allocations.
- `summarizeClusters()`: replaced per-cluster data-frame subset
  (`clustersInfo[which(…), ]`) with a logical index vector, eliminating one
  temporary frame copy per cluster.
- `collapsedGraph()`: pre-computed all cluster-membership index vectors once
  before the nested loops (`cluster_idx <- lapply(…)`); also hoisted the
  invariant `ni <- length(notInCluster)` out of the outer loop. Removes
  O(k²) repeated full-scan `which()` calls.
- `plotCluster()`: moved `edgecol` allocation outside the per-node loop and
  reset only the modified entries after each iteration, cutting per-iteration
  allocation from O(n) to O(1).
- `shortestPathDistance()`: removed the redundant third matrix copy in
  `An <- Ap <- minDist <- AdjMat`; `An` is unconditionally overwritten on
  the first line of the loop body.

## Housekeeping

- `benchmarks/benchmark.R`: removed stale comment and redundant
  `Matrix(getAdjMat(…) * 1L)` double-wrapping; `getAdjMat()` already returns
  a `Matrix` object since 0.2.7.

---

# betaMix 0.2.7

## Bug fixes

- Fixed `getAdjMat()` in-memory path: now always returns a `Matrix` object
  (was returning a base R matrix), consistent with the SQLite path and with
  the `is(A, "Matrix")` requirement in `graphComponents()` and
  `sphericalCaps()`.
- Fixed `summarizeClusters()`: `intEdges / degree` no longer produces a
  transient `NaN`/`Inf` vector before the zero-degree correction; replaced
  with a single `ifelse()`.
- Fixed `plotCluster()`: division by `max(degree)` and by per-node `degree`
  are now guarded — emits a warning and returns `NULL` if all cluster nodes
  have degree 0.
- Fixed EM E-step (three sites): `p0f0_n / (p0f0_n + p1f1_n)` now guards
  against a zero denominator (both beta densities underflow to 0) by treating
  affected observations as null (`m0 = 1`) and issuing a warning.

## Improvements

- `betaMix()` now warns when constant (zero-variance) columns are detected
  before calling `cor()`, identifying the offending column indices; their
  correlations are set to 0 as before.

## Defensive fixes in C++ (`writeLargeTable.cpp`)

- `stdvec()`: added an early-return guard with a warning when standard
  deviation is zero, preventing division by zero (safe by call-site contract,
  but now self-documenting).
- `corr()`: added guard against `n <= 1`, returning 0 with a warning instead
  of dividing by `n - 1 = 0`.

## Tests

- Added `tests/testthat/test-degenerate-inputs.R` (21 tests) covering
  all-zeros columns, all-ones columns, perfectly correlated column pairs,
  and a combined scenario; verifies correct warning messages, NaN-free `m0`,
  finite `ppthr`, and correct edge detection.
- Fixed `test-getAdjMat.R` symmetry test: `t(adj)` on a `Matrix` object
  dispatches to `t.default` via S3; fixed by converting to base matrix first.

---

# betaMix 0.2.6

## Bug fixes

- Fixed `betaMix()` with `ind=FALSE`: when `uniroot()` failed to find the root
  for `etahat`, the error object was stored back into `etahat`, causing
  `dbeta()` to throw "Non-numeric argument to mathematical function" on the next
  EM iteration. The function now retains the previous `etahat` value and emits
  a warning message instead.
- Fixed two `&` (bitwise) operators used instead of `&&` (logical short-circuit)
  in `if` and `while` conditions; the old form caused an extra EM iteration at
  convergence and could evaluate side-effects unnecessarily.
- Replaced `class(x) == "try-error"` checks with `inherits(x, "try-error")` for
  R 4.x compatibility.
- Fixed silent `stopifnot()` no-op in `sphericalCaps()` and `graphComponents()`:
  `grep("Matrix", class(A)) > 0` returns `logical(0)` for non-Matrix input
  (vacuously TRUE), replaced with `is(A, "Matrix")`.
- Fixed inconsistent boundary condition in `betaMix()`: initial
  `inNonNullSupport` used `<=` while all subsequent iterations used `<`.
- Fixed empty-vector crash in `sphericalCaps()`: added `length(orddeg) > 0`
  guard before `max(deg[orddeg])`.
- Fixed `clusteringCoef()`: neighbour lookup used `== 1` instead of `!= 0`,
  causing signed-edge graphs to miss negative-correlation neighbours.
- Fixed colour-opacity recycling bug in `plotCluster()`: opacity vector (cluster
  length) was applied before subsetting `nodecol` to cluster nodes, producing
  arbitrary per-node colours.
- Fixed edge-colour index-space bug in `plotCluster()`: `which(tmpA[i,nbrs]==-1)`
  returns positions within the subvector, not column indices; fixed to
  `nbrs[which(...)]`.
- Fixed `plotDegCC()`: was accessing non-existent field `betamixobj$AdjMat`
  (returns NULL silently); corrected to `getAdjMat(betamixobj)`.

## Performance

- Second-round optimisations across R and C++:
  - `angleMat` now stores the raw correlation matrix `corM` instead of
    `acos(corM)`, eliminating O(P²) `acos()` calls in `betaMix()` and the
    matching `sin()` calls in `getAdjMat()`; the identity
    `sin²(arccos(r)) = 1 − r²` is used in both places.
  - `getAdjMat()` signed-edge detection changed from
    `acos(r) > π/2` to the equivalent `r < 0`, avoiding a further O(E)
    `acos()` evaluation.
  - `etafun` (called by `uniroot`) now accepts precomputed scalar sums
    `sum_m0` and `sum_m0·log(z_j)` instead of recomputing them on every
    function evaluation (~10–50 evaluations per EM step for the `ind=FALSE`
    path).
  - `MLEfun` and `jacmle` (called by `nleqslv`) now accept precomputed
    scalar sums `sm_val`, `swt_logz`, `swt_log1mz`; each is now O(1)
    scalar arithmetic (~7 evaluations per EM step, all saved).
  - `dbeta(z_j[nns], etahat, 0.5)` precomputed once before the `while`
    loop when `ind=TRUE` (constant null density); refreshed only when
    `etahat` changes in the `ind=FALSE` path.
  - `p0new` recomputed as O(|nns|) rather than O(n), exploiting the fact
    that `m0[-nns] = 1` always.
  - Initial `p0` estimate uses `findInterval()` for O(log n) binary search
    on sorted `z_j` instead of an O(n) comparison.
  - `calcCorr.cpp`: `pow(sin(acos(r)),2)` replaced by `1.0 − r*r`;
    `stdvec()` changed to pass `NumericVector&` by reference (eliminates a
    copy per column); `corr()` uses a scalar accumulator loop instead of
    `sum(x*y)` (eliminates a temporary O(n) vector per pair); VLAs
    `char pname[]` / `char val[]` replaced with `std::vector<char>`.
  - `calcCorr.cpp` bug fix: output-file null check used `finptr` instead
    of `foutptr`, so a failed `fopen()` for the output file was never
    detected.

- `betaMix()` is now 1.9–2.7× faster on typical datasets after profiling
  identified and eliminated four major bottlenecks in the EM loop:
  - `intersect(which(), which())` in `MLEfun` (the dominant cost, ~44% of
    runtime) replaced by `which(& &)` — **47× faster** per call.
  - Dead code removed from `jacmle`: two filter/scale operations on `z_j0`
    were performed but never used (the Jacobian depends only on `sum(m-1)`
    and trigamma values).
  - `inNonNullSupport` was recomputed inside every EM iteration despite
    `z_j` and `bmax` being constant; it is now precomputed once with
    `findInterval()` (O(log n) vs O(n) per iteration).
  - `MLEfun`/`jacmle` re-filtered and re-scaled `z_j` on every `nleqslv`
    function evaluation (~7 per EM step); the filtered/scaled slice and its
    log values are now precomputed before the loop.
  - E-step `dbeta()` is now called only on the nonnull-support slice (z_j <
    bmax); entries ≥ bmax have posterior m0 = 1 by definition and no longer
    require a density evaluation (~33% fewer `dbeta` calls for SIM N=200).
- `sin^2(arccos(r))` for the z_j statistic replaced by the algebraically
  equivalent `1 - r^2`, avoiding one `sin()` call per lower-triangle pair.
- Minor: `corM[which(is.na(corM))]` simplified to `corM[is.na(corM)]`;
  `length(which(z_j > t))` replaced by `sum(z_j > t)`.

## Infrastructure

- Renamed `Data/` directory to `data/` (lowercase) so that R's standard
  `LazyData` mechanism correctly installs and lazy-loads the package datasets
  (`DrySeeds`, `SixHourImbibed`, `SIM`).
- Added `tests/testthat/` test suite (68 tests covering `betaMix()`,
  `getAdjMat()`, `clusteringCoef()`, `graphComponents()`, `sphericalCaps()`,
  `summarizeClusters()`, `collapsedGraph()`, and `shortestPathDistance()`).
- Added `benchmarks/benchmark.R` performance benchmark script.
- Added `testthat (>= 3.0.0)` to `Suggests` in `DESCRIPTION`.

## Removed dead code

- Removed unused assignments `cccsig` (`plotFittedBetaMix()`), `inCluster`
  (`collapsedGraph()`), and `degs` (`shortestPathDistance()`).

## Documentation

- Corrected `@return` field names in `betaMix()`: `m_0` → `m0`, `p_0` → `p0`,
  `P` → `nodes` (matching actual list names).
- Corrected `@param subsamplesize` threshold: "greater than 20000" →
  "at least 20000".
- `getAdjMat()`: clarified that the function handles both in-memory and
  SQLite-backed correlation data, and that `res` is always required.
- `shortSummary()`: expanded description to list all printed values.
- `sphericalCaps()`: fixed typo "apherical" → "spherical".
- `graphComponents()`: expanded `@param type` with the formula used.
- `plotDegCC()`: corrected y-axis description (degree × CC, not a proportion).
- `plotCluster()`: fixed `@param edgecols` description (said "edgecol").
- `shortestPathDistance()`: removed stray "expMat" reference; fixed typo
  "shortset" → "shortest".

---

# betaMix 0.2.3

- `bmax` is now a user-adjustable parameter with default 0.999.

# betaMix 0.2.2

- Fixed bug in determination of the nonnull support.
- Improved estimation of `etahat` for dependent samples.
- Added Jacobian for faster MLE of the nonnull beta parameters.

# betaMix 0.2.1

- Added `sphericalCaps()` for detecting overlapping hub structures.
- Fixed `getAdjMat()` for single-node neighbourhood queries.
- Faster adjacency matrix construction for node subsets.

# betaMix 0.1.1

- Initial CRAN-ready release.
