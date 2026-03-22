# betaMix
The betaMix package detects edges in co-expression networks from gene expression (or other multivariate) data. Given an N × P expression matrix it identifies which pairs among the P predictors are significantly correlated — the edges of the underlying graphical model.

**Statistical model.** The method is built on the geometric result of Frankl and Maehara (1990): in high-dimensional space, two random (unassociated) vectors are nearly perpendicular with high probability. The pairwise correlations equal the cosines of the angles between column vectors, and z_j = sin²(θ) is computed for each of the P(P-1)/2 pairs. Under the null (no association), z_j ~ Beta((N-1)/2, 1/2); correlated pairs follow Beta(a, b) on a rescaled support [0, b_max]. A two-component beta mixture is fitted by the EM algorithm, simultaneously estimating a, b, and the null proportion p0.

**One-step efficiency.** betaMix fits a single model to all P(P-1)/2 pairwise statistics simultaneously. Many competing methods estimate the network by running P separate penalized regressions, each designating one variable as the response and the remaining P-1 as predictors. betaMix avoids this entirely: one EM run produces the fitted mixture, and a single threshold determines all edges at once.

**False positive control.** betaMix provides two complementary safeguards enforced simultaneously:
- **`maxalpha`** (frequentist): the probability that any truly unassociated pair is declared an edge is bounded by `maxalpha` (default 1e-4). For large P this should be set much smaller (e.g. 1e-6) to account for the large number of tests.
- **`ppr`** (Bayesian): only pairs whose posterior probability of belonging to the null component is below `ppr` (default 0.05) are declared edges.

The edge-detection threshold is set to satisfy whichever criterion is more stringent.

The input to the program is a normalized expression matrix, with predictors/genes (nodes) in the columns, and samples in the rows.

With a large number of predictors, P, the estimation may be slow, so it is recommended to set the parameter subsamplesize to something smaller than choose(P,2). The minimum allowed by the program is 20,000. Using anything smaller will cause betaMix to fit the model to all choose(P,2) pairs.

If the N samples can be assumed to be independent, set the parameter ind to TRUE (the default). If it is set to FALSE, the null set follows a Beta((nu-1)/2, 1/2) distribution and nu (the effective sample size) is estimated from the data in the EM algorithm.

The betaMix package depends on the Matrix package, to allow for efficient storage and computation of large adjacency matrices. If P is very large, the correlation matrix may require an unattainable amount of memory. For such cases, the betaMix package has a procedure to store the correlation data in a SQLite database, thus never loading more than two columns into memory. For this purpose, the package depends on DBI and RSQLite. For the EM algorithm, betaMix requires the nleqslv package.
