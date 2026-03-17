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
