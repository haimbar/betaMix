# betaMix 0.2.4

## Bug fixes

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
