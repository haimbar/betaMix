# betaMix
The betaMix package is used to find edges in gene networks using co-expression data. The betaMix method is built on the insightful results of Frankl and Maehara (1990) who showed that two random vectors are approximately perpendicular with high probability, if the dimension of the space is sufficiently large. The pair-wise correlations between pairs of predictors are equal to the cosine of the angles between the pairs. From these angles we compute $z_j=\sin^2(\theta)$ and fit a mixture of two beta distributions to the $z_j$'s. For pairs of random vectors (the null set) the distribution of $z_j$ is Beta((N-1)/2, 1/2), where N is the sample size. The nonnull set is assumed to follow a Beta(a,b) distribution, and using the EM algorithm we estimate a,b, and the proportion, p0, of the null set of pairs. The betaMix function determines a threshold which will control the error rate given by the user. Any $z_j$ below that threshold corresponds to a significantly correlated pair of predictors (an edge in the graphical model.)

The input to the program is a normalized expression matrix, with predictors/genes (nodes) in the columns, and samples in the rows.

With a large number of predictors, P, the estimation may be slow, so it is recommended to set the parameter subsamplesize to something smaller than choose(P,2). The minimum allowed by the program is 20,000. Using anything smaller will cause betaMix to fit the model to all choose(P,2) pairs. 

If the N samples can be assumed to be independent, set the parameter ind to TRUE. If it is set to FALSE (the default), the null set follows a Beta((nu-1)/2, 1/2) distribution and nu (the effective sample size) is estimated from the data in the EM algorithm.

The betaMix package depends on the 'Matrix' package, to allow for efficient storage and computation of large co-occurrence matrices. If P is very large, the correlation matrix may require an unattainable amount of memory. For such cases, the betaMix package has a procedure to store the correlation data in a SQLite database, thus never loading more than two columns into memory. For this purpose, the package depends on DBI and RSQLite. For the EM algorithm, betaMix requires the nleqslv package.
