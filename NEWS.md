# betaMix 0.2.5

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
