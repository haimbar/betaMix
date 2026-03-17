#' @useDynLib betaMix
#' @importFrom Rcpp sourceCpp evalCpp
#' @export calcCorr
NULL


#' Fit a two-component beta mixture model to the matrix of all pairwise correlations.
#'
#' From the pairwise correlations, the function calculates the statistics z_j=sin^2(arccos(cor(y_i,y_j))) and fits the two-component model using the EM algorithm.
#' @param M A matrix with N rows (samples) and P columns (variables).
#' @param dbname The sqlite database, if one is used to store the pairwise correlation data instead of using the cor function and storing the cor(M) matrix in memory (for situations in which P is very large).
#' @param tol The convergence threshold for the EM algorithm (default=1e-4, but taken to be the maximum of the user's input and 1/(P(P-1)/2)).
#' @param calcAcc The calculation accuracy threshold (to avoid values greater than 1 when calling asin.) Default=1e-9.
#' @param maxalpha The probability of Type I error (default=1e-4). For a large P, use a much smaller value.
#' @param ppr The null posterior probability threshold (default=0.05).
#' @param mxcnt The maximum number of EM iterations (default=200).
#' @param ahat The initial value for the first parameter of the nonnull beta distribution (default=8).
#' @param bhat The initial value for the second parameter of the nonnull beta distribution (default=3).
#' @param bmax The RHS of the support of the non-null component (default=0.999)
#' @param subsamplesize If at least 20000 and the number of pairs exceeds this
#'   value, take a random sample of this size to fit the model. Otherwise, use
#'   all the data (default=50000).
#' @param seed The random seed to use if selecting a subset with the subsamplesize parameter (default=912469).
#' @param ind Whether the N samples should be assumed to be independent (default=TRUE).
#' @param msg Whether to print intermediate output messages (default=TRUE).
#' @param method Correlation method: \code{"pearson"} (default) or \code{"spearman"}.
#'   Spearman's rank correlation is more robust to outliers and non-linear
#'   monotone relationships. In the SQLite path this controls \code{calcCorr()};
#'   in the in-memory path it is passed directly to \code{cor()}.
#' @return A list with the following:
#' \describe{
#' \item{angleMat}{A PxP matrix with angles between pairs of vectors. If the correlation data is stored in SQLite, then the returned value is database name.}
#' \item{z_j}{The statistics z_j=sin^2(angles).}
#' \item{m0}{The posterior null probabilities.}
#' \item{p0}{The estimated probability of the null component.}
#' \item{ahat}{The estimated first parameter of the nonnull beta component.}
#' \item{bhat}{The estimated second parameter of the nonnull beta component.}
#' \item{N}{The sample size.}
#' \item{etahat}{If the samples are not assumed to be independent, this corresponds to the effective sample size, ESS=2*etahat+1}
#' \item{bmax}{The user-defined right-hand side of the support of the non-null component.}
#' \item{ppthr}{The estimated posterior probability threshold, under which all the z_j correspond to nonnull edges.}
#' \item{nodes}{The number of nodes.}
#' \item{edges}{The number of edges found.}
#' \item{cnt}{The number of EM iterations.}
#' \item{method}{The correlation method used (\code{"pearson"} or \code{"spearman"}).}
#' }
#' @seealso \code{\link{betaMix-package}} for an overview and index of all functions and datasets in this package.
#' @export
#' @import DBI stats nleqslv
#' @importFrom stats cor dbeta pbeta qbeta optimize
#' @examples
#' \dontrun{
#' data(SIM, package = "betaMix") # variables correspond to columns, samples to rows
#' res <- betaMix(betaMix::SIM, maxalpha = 1e-6, ppr = 0.01, subsamplesize = 30000, ind = TRUE)
#' plotFittedBetaMix(res)
#' }
betaMix <- function(M, dbname = NULL, tol = 1e-4, calcAcc = 1e-9, maxalpha = 1e-4,
                    ppr = 0.05, mxcnt = 200, ahat = 8, bhat = 3, bmax = 0.999,
                    subsamplesize = 50000, seed = 912469, ind = TRUE, msg = TRUE,
                    method = c("pearson", "spearman")) {
  method <- match.arg(method)
  if (!is.numeric(bmax) || length(bmax) != 1 || bmax <= 0 || bmax >= 1)
    stop("betaMix: 'bmax' must be a single number in (0, 1); got bmax = ", bmax)
  if (!is.numeric(calcAcc) || length(calcAcc) != 1 || calcAcc <= 0 || calcAcc >= 0.5)
    stop("betaMix: 'calcAcc' must be a single number in (0, 0.5); got calcAcc = ", calcAcc)
  if (msg) cat("Generating the z_ij statistics...\n")
  if (!is.null(dbname)) {
    if (!file.exists(dbname)) {
      message(paste("SQLite database", dbname, "not found!\n"))
      return(NULL)
    }
    con <- dbConnect(RSQLite::SQLite(), dbname)
    res <- dbSendQuery(con, "SELECT * FROM metadata")
    metadata <- dbFetch(res)
    dbClearResult(res)
    P <- metadata$p
    N <- metadata$n
    etahat <- (N - 1) / 2
    if (subsamplesize < 20000) {
      res <- dbSendQuery(con, sprintf("SELECT * FROM correlations"))
    } else {
      res <- dbSendQuery(con, sprintf("SELECT * FROM correlations ORDER BY random() LIMIT %d", subsamplesize))
    }
    subtable <- dbFetch(res)
    dbClearResult(res)
    z_j <- sort(pmin(1 - calcAcc, pmax(calcAcc, subtable$zij)))
    angleMat <- dbname
  } else {
    N <- nrow(M)
    P <- ncol(M)
    if (N < 2) stop("betaMix: M must have at least 2 rows (samples); got N = ", N)
    if (P < 2) stop("betaMix: M must have at least 2 columns (variables); got P = ", P)
    etahat <- (N - 1) / 2
    const_cols <- which(apply(M, 2, sd) == 0)
    if (length(const_cols) > 0) {
      warning("betaMix: ", length(const_cols), " constant (zero-variance) column(s) detected at index ",
              paste(const_cols, collapse = ", "), "; their correlations will be set to 0")
    }
    corM <- cor(M, method = method)
    if (any(is.na(corM))) {
      corM[is.na(corM)] <- 0
    }
    # Store corM directly (not acos(corM)): getAdjMat uses (1 - corM^2) < ppthr
    # instead of sin(acos(r))^2 — same value, eliminates P^2 trig calls.
    angleMat <- corM
    # sin^2(arccos(r)) == 1 - r^2; avoids sin() on the lower triangle
    z_j <- pmin(1 - calcAcc, pmax(calcAcc, 1 - corM[lower.tri(corM)]^2))
    z_jall <- c()
    if ((length(z_j) > subsamplesize) && (subsamplesize >= 20000)) {
      z_jall <- sort(z_j)
      set.seed(seed)
      z_j <- sort(z_j[sample(length(z_j), subsamplesize)])
    } else {
      z_j <- sort(z_j) # sort so findInterval() works below
    }
  }
  tol <- max(tol, 1 / length(z_j))
  # O(log n) initial p0: z_j is sorted so findInterval gives a binary search
  thr0 <- qbeta(0.1, etahat, 0.5)
  p0 <- min((length(z_j) - findInterval(thr0, z_j)) / (0.9 * length(z_j)), 1)
  if (msg) cat("Fitting the model...\n")
  # Precompute nonnull support once: z_j is sorted and bmax is constant
  nns <- seq_len(findInterval(bmax, z_j))
  z_j_nns <- z_j[nns] / bmax
  # Interior mask and log values for MLEfun — also constant throughout EM
  inc <- which(z_j_nns > 1e-6 & z_j_nns < 1 - 1e-6)
  log_z <- log(z_j_nns[inc])
  log1mz <- log(1 - z_j_nns[inc])
  # Precompute null density once; stays constant when ind=TRUE (etahat fixed).
  # When ind=FALSE it is refreshed after each successful eta update below.
  # m0 is initialised here and only m0[nns] is ever updated — m0[-nns] == 1.
  m0 <- rep(1, length(z_j))
  if (length(nns) > 0) {
    dbeta0_nns <- dbeta(z_j[nns], etahat, 0.5)
    p0f0_n <- p0 * dbeta0_nns
    p1f1_n <- (1 - p0) * dbeta(z_j_nns, ahat, bhat)
    .denom <- p0f0_n + p1f1_n
    if (any(.denom == 0)) warning("betaMix E-step: zero total density for ", sum(.denom == 0), " observation(s); treating as null")
    m0[nns] <- ifelse(.denom > 0, p0f0_n / .denom, 1)
  }
  p0new <- p0 - 10 * tol
  cnt <- 0
  while (abs(p0 - p0new) > tol && (cnt <- cnt + 1) < mxcnt) {
    p0 <- p0new
    if (!ind) {
      # Precompute O(n) sums once; uniroot calls etafun ~10-50x per solve
      sum_m0 <- sum(m0)
      sum_m0_logz <- sum(m0 * log(z_j))
      etahat_new <- try(uniroot(etafun, c(1, (N - 1) / 2),
        sum_m0 = sum_m0, sum_m0_logz = sum_m0_logz,
        lower = 1, upper = (N - 1) / 2
      )$root, silent = TRUE)
      if (inherits(etahat_new, "try-error")) {
        message("betaMix error when estimating eta:", etahat_new, "\n")
      } else {
        etahat <- etahat_new
        # Refresh null density after successful eta update
        if (length(nns) > 0) {
          dbeta0_nns <- dbeta(z_j[nns], etahat, 0.5)
        }
      }
    }
    # M-step: precompute O(|nns|) quantities once per EM step so that
    # nleqslv's ~7 evaluations of MLEfun/jacmle are O(1) scalar arithmetic.
    m_nns_cur <- m0[nns]
    sm_val <- sum(1 - m_nns_cur)
    if (!any(is.na(m_nns_cur)) && sm_val >= 1e-10) {
      wt <- 1 - m_nns_cur[inc]
      swt_logz <- sum(wt * log_z)
      swt_log1mz <- sum(wt * log1mz)
      ests_full <- try(nleqslv(c(ahat, bhat), MLEfun,
        jac = jacmle,
        sm_val = sm_val, swt_logz = swt_logz,
        swt_log1mz = swt_log1mz
      ), silent = TRUE)
      if (inherits(ests_full, "try-error")) {
        message("betaMix error when estimating a and b:", ests_full, "\n")
      } else if (ests_full$termcd %in% c(1L, 2L)) {
        ahat <- min(max(ests_full$x[1], 0.5), 1000)
        bhat <- min(max(ests_full$x[2], 0.5), 1000)
      }
    }
    # E-step: reuse precomputed dbeta0_nns; m0[-nns] stays 1 (no NA possible)
    if (length(nns) > 0) {
      p0f0_n <- p0 * dbeta0_nns
      p1f1_n <- (1 - p0) * dbeta(z_j_nns, ahat, bhat)
      .denom <- p0f0_n + p1f1_n
      if (any(.denom == 0)) warning("betaMix E-step: zero total density for ", sum(.denom == 0), " observation(s); treating as null")
      m0[nns] <- ifelse(.denom > 0, p0f0_n / .denom, 1)
    }
    # O(|nns|) mean: m0[-nns] == 1 always, no NA after pmax/pmin
    p0new <- (sum(m0[nns]) + length(z_j) - length(nns)) / length(z_j)
  }
  if (!is.null(dbname)) {
    ppthr <- qbeta(maxalpha, etahat, 0.5)
    critPP <- which(m0 < ppr)
    if (length(critPP) > 0) {
      ppthr <- max(z_j[critPP])
    }
    p0 <- mean(m0)
    res <- dbSendQuery(con, sprintf("SELECT * FROM correlations WHERE zij < %f", ppthr))
    selected <- dbFetch(res)
    edges <- nrow(selected)
    dbClearResult(res)
    dbDisconnect(con)
  } else {
    if (length(z_jall) > subsamplesize) {
      z_j <- z_jall
    }
    # Final E-step on the full z_j (recompute nns since z_j may have expanded)
    nns_f <- seq_len(findInterval(bmax, z_j))
    m0 <- rep(1, length(z_j))
    if (length(nns_f) > 0) {
      z_j_nns_f <- z_j[nns_f] / bmax
      p0f0_n <- p0 * dbeta(z_j[nns_f], etahat, 0.5)
      p1f1_n <- (1 - p0) * dbeta(z_j_nns_f, ahat, bhat)
      .denom <- p0f0_n + p1f1_n
      if (any(.denom == 0)) warning("betaMix final E-step: zero total density for ", sum(.denom == 0), " observation(s); treating as null")
      m0[nns_f] <- ifelse(.denom > 0, p0f0_n / .denom, 1)
    }
    p0 <- mean(m0)
    ppthr <- qbeta(maxalpha, etahat, 0.5)
    critPP <- which(m0 < ppr)
    if (length(critPP) > 0) {
      ppthr <- max(z_j[critPP])
    }
    nonnull <- which(z_j <= ppthr)
    edges <- length(nonnull)
  }
  if (msg) cat("Done.\n")
  list(
    angleMat = angleMat, z_j = z_j, m0 = m0, p0 = p0, N = N, ahat = ahat, bhat = bhat,
    etahat = etahat, bmax = bmax, ppthr = ppthr, nodes = P, edges = edges, cnt = cnt,
    method = method
  )
}


#' Get the adjacency matrix from a fitted betaMix model.
#'
#' Works for both in-memory correlation matrices and SQLite-backed storage.
#' SQLite storage should be used if P is large and calculating and storing the
#' full correlation matrix would require too much memory.
#' @param res The object returned by betaMix(). Always required: used to obtain
#'   the posterior probability threshold and node count regardless of whether
#'   the correlation data is in memory or in a database.
#' @param dbname The sqlite database with the pairwise correlations, if the data is not stored in memory (for very large P). Default=NULL.
#' @param ppthr A threshold to select edges. Default=NULL, in which case the returned value from betaMix() is used (available in res).
#' @param signed If TRUE, the returned matrix will contain -1 for negatively correlated pairs. Otherwise, all correlated pairs will have 1 in the returned matrix (Default).
#' @param nodes An optional parameter which allows to create a sparse adjacency matrix containing only the neighbors of the selected nodes (default=NULL).
#' @importFrom Matrix Matrix rowSums colSums diag
#' @seealso \code{\link{betaMix-package}} for an overview and index of all functions and datasets in this package.
#' @export
#' @examples
#' \dontrun{
#' res <- betaMix(betaMix::SIM, maxalpha = 1e-5, ppr = 0.01, subsamplesize = 30000, ind = TRUE)
#' adjMat <- getAdjMat(res)
#' image(adjMat[1:80, 1:80])
#' # See online documentation for usage when SQLite is used.
#' }
getAdjMat <- function(res, dbname = NULL, ppthr = NULL, signed = FALSE, nodes = NULL) {
  if (!is.null(nodes)) {
    if (!is.numeric(nodes) || any(nodes < 1) || any(nodes > res$nodes))
      stop("getAdjMat: 'nodes' must be integers in [1, ", res$nodes, "]")
  }
  if (is.null(ppthr)) {
    ppthr <- res$ppthr
  }
  if (is.null(dbname)) {
    if (!is.null(nodes)) {
      if (length(nodes) == 1) {
        Atmp <- Matrix(res$angleMat[nodes, ], nrow = 1)
      } else {
        Atmp <- Matrix(res$angleMat[nodes, ])
      }
      nbrs <- which(colSums((1 - Atmp^2) < ppthr) > 0)
      selected <- sort(union(nodes, nbrs))
      angMat <- res$angleMat[selected, selected]
    } else {
      angMat <- res$angleMat
    }
    # angleMat now holds corM; sin(acos(r))^2 == 1-r^2, acos(r)>pi/2 == r<0
    A <- (1 - angMat^2) < ppthr
    if (signed) {
      negedges <- intersect(which(A > 0), which(angMat < 0))
      if (length(negedges) > 0) {
        A[negedges] <- -1
      }
    }
    diag(A) <- FALSE
    return(Matrix::Matrix(A))
  }
  if (!file.exists(dbname)) {
    message(paste("SQLite database", dbname, "not found!\n"))
    return(NULL)
  }
  con <- dbConnect(RSQLite::SQLite(), dbname)
  res <- dbSendQuery(con, "SELECT * FROM metadata")
  metadata <- dbFetch(res)
  dbClearResult(res)
  P <- metadata$p
  res <- dbSendQuery(con, "select name from predictors")
  pnames <- unlist(dbFetch(res))
  dbClearResult(res)
  A <- Matrix::Matrix(FALSE, P, P)
  colnames(A) <- pnames
  if (is.null(nodes)) {
    res <- dbSendQuery(con, sprintf("SELECT * FROM correlations WHERE zij < %f", ppthr))
  } else {
    sset <- paste(nodes, collapse = ",")
    res <- dbSendQuery(con, sprintf("SELECT * FROM correlations WHERE (node1 IN (%s) OR node2 IN (%s)) AND zij < %f", sset, sset, ppthr))
  }
  selected <- dbFetch(res)
  dbClearResult(res)
  dbDisconnect(con)
  if (nrow(selected) == 0) {
    message("Zero rows selected.\n")
    return(NULL)
  }
  selectedCells <- P * (selected$node1 - 1) + selected$node2
  A[selectedCells] <- TRUE
  if (signed) {
    A[selectedCells[which(selected$corr < 0)]] <- -1
  }
  return(A + Matrix::t(A))
}


# the Jacobian, to speed up the MLE calculation of a and b in the non-null
# component.
# sm_val: sum(1 - m0[nns]) — precomputed before nleqslv so it is not
# recomputed on each of nleqslv's ~7 function evaluations per EM step.
# '...' absorbs the swt_logz / swt_log1mz scalars forwarded by nleqslv.
jacmle <- function(par, sm_val, ...) {
  a <- par[1]
  b <- par[2]
  tab <- trigamma(a + b) # computed once; used three times below
  # J[i,j] = d(MLEfun[i])/d(par[j]):
  #   J[1,1] = sm_val*(trigamma(a) - trigamma(a+b))  > 0
  #   J[1,2] = J[2,1] = -sm_val*trigamma(a+b)        < 0
  #   J[2,2] = sm_val*(trigamma(b) - trigamma(a+b))  > 0
  sm_val * matrix(c(
    trigamma(a) - tab, -tab,
    -tab, trigamma(b) - tab
  ), 2, 2)
}


# The maximum likelihood estimation of a,b of the nonnull component.
# nu (eta = (nu-1)/2) is estimated separately.
# All O(|nns|) sums are precomputed once before each nleqslv call and passed
# as scalars, so MLEfun is pure scalar arithmetic (~7x fewer vector ops/step).
MLEfun <- function(par, sm_val, swt_logz, swt_log1mz) {
  a <- par[1]
  b <- par[2]
  dab <- digamma(a + b) # computed once; used twice below
  c(
    sm_val * (digamma(a) - dab) - swt_logz,
    sm_val * (digamma(b) - dab) - swt_log1mz
  )
}


# The maximum likelihood estimation of the null parameter (if samples are dependent).
# sum_m0 and sum_m0_logz are O(n) sums that are constant within a uniroot solve;
# precomputing them before uniroot avoids recomputing them on each evaluation.
etafun <- function(eta, sum_m0, sum_m0_logz) {
  -(digamma(eta) - digamma(eta + 0.5)) * sum_m0 + sum_m0_logz
}


#' Model-based ROC curve for edge-detection threshold selection.
#'
#' Computes and plots a model-based ROC curve by sweeping the edge-detection
#' threshold \eqn{\tau} from 0 to 1.  At each \eqn{\tau}:
#' \itemize{
#'   \item \strong{FPR}(\eqn{\tau}) = P(z < \eqn{\tau} | null)
#'         = \code{pbeta(tau, etahat, 0.5)}
#'   \item \strong{TPR}(\eqn{\tau}) = P(z < \eqn{\tau} | non-null)
#'         = \code{pbeta(min(tau, bmax) / bmax, ahat, bhat)}
#' }
#' No ground-truth labels are needed; both rates are derived entirely from the
#' fitted model distributions.  The curve is therefore exact under the
#' assumed model.
#'
#' Two reference thresholds are marked:
#' \itemize{
#'   \item \strong{Orange} — the current \code{ppthr} from \code{\link{betaMix}}.
#'   \item \strong{Red} — the Youden-index optimum,
#'     \eqn{\tau^* = \arg\max_\tau [\text{TPR}(\tau) - \text{FPR}(\tau)]},
#'     which maximises the sum of sensitivity and specificity.
#' }
#' A grey vertical dotted line marks the FPR at \eqn{\tau = b_{\max}}; for
#' \eqn{\tau > b_{\max}} the TPR is 1 (all non-null pairs detected) and the
#' curve runs horizontally to (1, 1).
#'
#' @param betaMixObj The list returned by \code{\link{betaMix}}.
#' @return Invisibly, a named list:
#' \describe{
#'   \item{tau}{Threshold grid (length 2001, from 0 to 1).}
#'   \item{fpr}{False positive rate at each \eqn{\tau}.}
#'   \item{tpr}{True positive rate at each \eqn{\tau}.}
#'   \item{auc}{Area under the ROC curve (trapezoidal rule).}
#'   \item{tau_youden}{Threshold that maximises TPR \eqn{-} FPR.}
#'   \item{tau_current}{Current threshold (\code{ppthr}) from the fitted model.}
#' }
#' @seealso \code{\link{betaMix-package}} for an overview and index of all functions and datasets in this package.
#' @export
#' @import stats graphics
#' @examples
#' \dontrun{
#' data(SIM, package = "betaMix")
#' res <- betaMix(betaMix::SIM, maxalpha = 1e-5, ppr = 0.01,
#'                subsamplesize = 30000, ind = TRUE)
#' roc <- plotROC(res)
#' roc$tau_youden     # Youden-optimal threshold
#' # Apply the new threshold:
#' adjMat <- getAdjMat(res, ppthr = roc$tau_youden)
#' }
plotROC <- function(betaMixObj) {
  with(betaMixObj, {

    tau_grid <- seq(0, 1, length.out = 2001)

    # FPR from null; TPR from non-null (clamped at bmax)
    fpr <- pbeta(tau_grid, etahat, 0.5)
    tpr <- pbeta(pmin(tau_grid, bmax) / bmax, ahat, bhat)

    # AUC via trapezoidal rule
    auc <- sum(diff(fpr) * (tpr[-1L] + tpr[-length(tpr)]) / 2)

    # Youden index: argmax(TPR - FPR)
    j_idx      <- which.max(tpr - fpr)
    tau_youden <- tau_grid[j_idx]
    fpr_youden <- fpr[j_idx]
    tpr_youden <- tpr[j_idx]

    # Current threshold
    fpr_cur <- pbeta(ppthr, etahat, 0.5)
    tpr_cur <- pbeta(min(ppthr, bmax) / bmax, ahat, bhat)

    plot(fpr, tpr, type = "l", lwd = 2, col = "steelblue",
         xlim = c(0, 1), ylim = c(0, 1),
         xlab = "False Positive Rate  P(z < tau | null)",
         ylab = "True Positive Rate  P(z < tau | non-null)",
         main = sprintf("Model-based ROC  (AUC = %.4f)", auc))
    abline(0, 1, lty = 2, col = "grey60")
    # mark where non-null support ends (tau = bmax)
    abline(v = pbeta(bmax, etahat, 0.5), lty = 3, col = "grey70")

    points(fpr_cur,    tpr_cur,    pch = 21, bg = "orange", cex = 1.8)
    points(fpr_youden, tpr_youden, pch = 21, bg = "tomato", cex = 1.8)
    legend("bottomright",
           c(sprintf("Current  tau = %.3g  (FPR=%.3f, TPR=%.3f)",
                     ppthr, fpr_cur, tpr_cur),
             sprintf("Youden   tau = %.3g  (FPR=%.3f, TPR=%.3f)",
                     tau_youden, fpr_youden, tpr_youden)),
           pch = 21, pt.bg = c("orange", "tomato"), pt.cex = 1.4,
           cex = 0.75, bty = "n")

    cat("=== plotROC: Model-based ROC ===\n")
    cat(sprintf("  AUC = %.4f\n", auc))
    cat(sprintf("  Current  tau = %.4g  FPR = %.4f  TPR = %.4f  edges = %d\n",
                ppthr, fpr_cur, tpr_cur, edges))
    cat(sprintf("  Youden   tau = %.4g  FPR = %.4f  TPR = %.4f\n",
                tau_youden, fpr_youden, tpr_youden))
    cat("  To apply: getAdjMat(res, ppthr = roc$tau_youden)\n")
    cat("================================\n")

    invisible(list(
      tau         = tau_grid,
      fpr         = fpr,
      tpr         = tpr,
      auc         = auc,
      tau_youden  = tau_youden,
      tau_current = ppthr
    ))
  })
}


#' Plot the kernel density of non-null observations with fitted Beta overlay.
#'
#' Plots the empirical kernel density of non-null-assigned \eqn{z_j} values
#' overlaid with the fitted \eqn{\text{Beta}(\hat{a}, \hat{b})} density.
#' When more than one mode is detected, orange dotted vertical lines mark each
#' additional peak.
#'
#' @param betaMixObj The list returned by \code{\link{betaMix}}.
#' @return Invisibly, a named list:
#' \describe{
#'   \item{n_modes}{Integer number of detected density peaks.}
#'   \item{peaks_x}{x-positions (on the \eqn{z_j} scale) of the detected peaks.}
#' }
#' @seealso \code{\link{betaMix-package}} for an overview and index of all functions and datasets in this package.
#' @export
#' @import graphics stats
#' @examples
#' \dontrun{
#' data(SIM, package = "betaMix")
#' res <- betaMix(betaMix::SIM, maxalpha = 1e-5, ppr = 0.01,
#'                subsamplesize = 30000, ind = TRUE)
#' plotNonNullDensity(res)
#' }
plotNonNullDensity <- function(betaMixObj) {
  with(betaMixObj, {
    nonnull_idx <- which(m0 < 0.5)
    n_nonnull   <- length(nonnull_idx)

    if (n_nonnull < 10) {
      plot.new()
      text(0.5, 0.5, "Insufficient non-null\nobservations", cex = 1.1)
      title(main = "Non-null density  (no data)")
      return(invisible(list(n_modes = NA_integer_, peaks_x = numeric(0))))
    }

    dens_nn        <- density(z_j[nonnull_idx])
    peak_threshold <- 0.10 * max(dens_nn$y)
    dy             <- dens_nn$y
    n_d            <- length(dy)
    is_peak        <- c(FALSE,
                        dy[seq(2, n_d - 1)] > dy[seq(1, n_d - 2)] &
                        dy[seq(2, n_d - 1)] > dy[seq(3, n_d)],
                        FALSE)
    peaks_above_thr <- which(is_peak & dy >= peak_threshold)
    n_modes         <- length(peaks_above_thr)

    xgrid     <- seq(1e-4, bmax - 1e-4, length.out = 512)
    beta_dens <- dbeta(xgrid / bmax, ahat, bhat) / bmax   # Jacobian correction
    ylim_p    <- c(0, max(max(dens_nn$y), max(beta_dens)) * 1.05)

    plot(dens_nn,
         ylim = ylim_p,
         main = sprintf("Non-null density  (%d mode%s)",
                        n_modes, if (n_modes == 1L) "" else "s"),
         xlab = expression(z[j]),
         ylab = "Density",
         col = "steelblue", lwd = 2)
    lines(xgrid, beta_dens, col = "firebrick", lwd = 2, lty = 2)
    if (length(peaks_above_thr) > 1L)
      abline(v = dens_nn$x[peaks_above_thr], col = "orange", lty = 3, lwd = 1.5)
    legend("topright",
           legend = c("Kernel density",
                      sprintf("Beta(%.2f, %.2f)/bmax", ahat, bhat)),
           col    = c("steelblue", "firebrick"),
           lty    = c(1, 2), lwd = 2, cex = 0.75, bty = "n")

    invisible(list(n_modes = n_modes,
                   peaks_x = dens_nn$x[peaks_above_thr]))
  })
}


#' Plot the histogram of the z_j and the fitted mixture distribution.
#'
#' @param betaMixObj An object returned from betaMix()
#' @param yLim The maximum value on the y-axis (default=5)
#' @seealso \code{\link{betaMix-package}} for an overview and index of all functions and datasets in this package.
#' @export
#' @examples
#' \dontrun{
#' data(SIM, package = "betaMix")
#' res <- betaMix(betaMix::SIM, maxalpha = 1e-5, ppr = 0.01, subsamplesize = 30000, ind = TRUE)
#' plotFittedBetaMix(res)
#' }
plotFittedBetaMix <- function(betaMixObj, yLim = 5) {
  with(betaMixObj, {
    ccc <- seq(0.001, 0.999, length = 1000)
    hist(z_j,
      freq = FALSE, breaks = 300, border = "grey", main = "", ylim = c(0, yLim),
      xlim = c(0, 1), xlab = expression(sin^{
        2
      } ~ (theta))
    )
    lines(ccc, p0 * dbeta(ccc, etahat, 0.5), col = 5, lwd = 4, lty = 2)
    lines(ccc, (1 - p0) * dbeta(ccc / bmax, ahat, bhat), lwd = 1, col = 2)
    lines(ccc, (1 - p0) * dbeta(ccc / bmax, ahat, bhat) +
      p0 * dbeta(ccc, etahat, 0.5), col = 4, lwd = 2)
    rect(0, 0, ppthr, yLim, col = "#FF7F5020", border = "orange")
  })
}

#' Assess the goodness of fit of a betaMix model.
#'
#' Produces a 2x2 diagnostic plot and a printed verbal assessment of how well
#' the fitted two-component beta mixture describes the data.
#'
#' The four panels are:
#' \enumerate{
#'   \item The standard \code{\link{plotFittedBetaMix}} histogram overlay.
#'   \item Q-Q plot for the null component: observations with posterior null
#'     probability \code{m0 >= qq_purity} (default 0.9), compared against
#'     \code{Beta(etahat, 0.5)} \emph{conditional on} the sub-range spanned
#'     by that pure-null subset.  Restricting to this high-purity upper tail
#'     removes contamination from the overlap region with the non-null
#'     component.  Falls back to all hard-assigned observations
#'     (\code{m0 >= 0.5}) if fewer than 10 pass the threshold.
#'   \item Q-Q plot for the non-null component: observations with posterior
#'     non-null probability \code{m0 <= 1 - qq_purity} (default
#'     \eqn{\le 0.1}), scaled by \code{bmax}, compared against
#'     \code{Beta(ahat, bhat)} \emph{conditional on} the sub-range of the
#'     pure-non-null subset.  Restricting to this high-purity lower tail
#'     removes contamination from the right-side overlap with the null.
#'     Falls back to all \code{m0 < 0.5} observations if fewer than 10
#'     pass.
#'   \item Model-based ROC curve (see \code{\link{plotROC}}): FPR and TPR swept
#'     over all thresholds using the fitted Beta distributions.  The orange
#'     point marks the current threshold; the red point marks the Youden-index
#'     optimum.
#' }
#'
#' Q-Q points are blue when the Kolmogorov-Smirnov D statistic is
#' \eqn{\le 0.05} and red when \eqn{D > 0.05}.  Raw p-values are not used
#' as the decision threshold because with tens of thousands of observations
#' even \eqn{D \approx 0.01} yields \eqn{p \approx 0}.
#'
#' Multimodality is assessed by counting local maxima of the kernel density
#' that reach at least 10\% of the global maximum.
#'
#' @param betaMixObj The list returned by \code{\link{betaMix}}.
#' @param yLim The maximum y-axis value for panel 1 (passed to
#'   \code{\link{plotFittedBetaMix}}).  Default \code{5}.
#' @param qq_purity Posterior probability threshold for the high-purity Q-Q
#'   subsets.  Only observations with \code{m0 >= qq_purity} are used for
#'   the null Q-Q panel, and only those with \code{m0 <= 1 - qq_purity} for
#'   the non-null panel.  The reference distribution in each case is the
#'   fitted Beta \emph{conditional on} the sub-range spanned by the pure
#'   subset, excluding the overlap region.  Must be in (0.5, 1).  Default
#'   \code{0.9}.  Falls back to hard assignment if fewer than 10
#'   observations pass.
#' @return Invisibly, a named list:
#' \describe{
#'   \item{ks_null}{Result of \code{ks.test} for the null component
#'     (\code{NULL} if fewer than 2 null-assigned observations).}
#'   \item{ks_nonnull}{Result of \code{ks.test} for the non-null component
#'     (\code{NULL} if fewer than 2 non-null-assigned observations after
#'     boundary filtering).}
#'   \item{ks_null_pure}{Conditional KS test for pure-null observations
#'     (\code{m0 >= qq_purity}) vs \code{Beta(etahat, 0.5)} restricted to
#'     their observed range.  \code{NULL} if fewer than 2 such observations.}
#'   \item{ks_nonnull_pure}{Conditional KS test for pure-non-null observations
#'     (\code{m0 <= 1 - qq_purity}, scaled by \code{bmax}) vs
#'     \code{Beta(ahat, bhat)} restricted to their observed range.
#'     \code{NULL} if fewer than 2 such observations.}
#'   \item{n_null}{Number of hard-assigned null observations.}
#'   \item{n_nonnull}{Number of hard-assigned non-null observations.}
#'   \item{n_null_pure}{Number of pure-null observations (\code{m0 >= qq_purity}).}
#'   \item{n_nonnull_pure}{Number of pure-non-null observations
#'     (\code{m0 <= 1 - qq_purity}).}
#'   \item{n_modes}{Number of detected density peaks in the non-null component;
#'     \code{NA} if fewer than 10 non-null observations are available.}
#'   \item{assessment}{Overall character verdict: \code{"Good fit"},
#'     \code{"Acceptable fit"}, or \code{"Deviations detected"}.}
#' }
#' @seealso \code{\link{betaMix-package}} for an overview and index of all functions and datasets in this package.
#' @export
#' @import stats graphics
#' @examples
#' \dontrun{
#' data(SIM, package = "betaMix")
#' res <- betaMix(betaMix::SIM, maxalpha = 1e-5, ppr = 0.01,
#'                subsamplesize = 30000, ind = TRUE)
#' fit_diag <- assessFit(res)
#' fit_diag$assessment
#' }
assessFit <- function(betaMixObj, yLim = 5, qq_purity = 0.9) {
  with(betaMixObj, {

    ## ── Hard assignment ───────────────────────────────────────────────────────
    null_idx    <- which(m0 >= 0.5)
    nonnull_idx <- which(m0 <  0.5)
    n_null      <- length(null_idx)
    n_nonnull   <- length(nonnull_idx)

    ## ── KS test: null component ───────────────────────────────────────────────
    ks_null <- NULL
    if (n_null >= 2) {
      ks_null <- ks.test(z_j[null_idx], "pbeta", etahat, 0.5)
    } else {
      warning("assessFit: fewer than 2 null-assigned observations; null KS test skipped.")
    }

    ## ── KS test: non-null component ───────────────────────────────────────────
    z_nn_raw <- if (n_nonnull > 0) z_j[nonnull_idx] / bmax else numeric(0)
    z_nn     <- z_nn_raw[z_nn_raw > 1e-6 & z_nn_raw < 1 - 1e-6]
    ks_nonnull <- NULL
    if (length(z_nn) >= 2) {
      ks_nonnull <- ks.test(z_nn, "pbeta", ahat, bhat)
    } else {
      warning("assessFit: fewer than 2 non-null observations after boundary filtering; ",
              "non-null KS test skipped.")
    }

    ## ── High-purity subsets for conditional Q-Q plots ─────────────────────────
    # Pure-null: observations firmly in the null (upper) tail, away from the
    # non-null overlap.  Pure-non-null: observations firmly in the non-null
    # (lower) tail, away from the null overlap.
    null_pure_idx    <- which(m0 >= qq_purity)
    nonnull_pure_idx <- which(m0 <= 1 - qq_purity)
    n_null_pure      <- length(null_pure_idx)

    z_nn_pure_raw  <- if (length(nonnull_pure_idx) > 0)
                        z_j[nonnull_pure_idx] / bmax else numeric(0)
    z_nn_pure      <- z_nn_pure_raw[z_nn_pure_raw > 1e-6 & z_nn_pure_raw < 1 - 1e-6]
    n_nonnull_pure <- length(z_nn_pure)

    # KS tests against the conditional (truncated) reference distribution.
    ks_null_pure    <- NULL
    null_D_pure     <- NA_real_
    ks_nonnull_pure <- NULL
    nonnull_D_pure  <- NA_real_

    if (n_null_pure >= 2) {
      z_np  <- sort(z_j[null_pure_idx])
      p_low <- pbeta(z_np[1L], etahat, 0.5)
      if (p_low < 1 - .Machine$double.eps) {
        cond_null_cdf <- function(z)
          (pbeta(z, etahat, 0.5) - p_low) / (1 - p_low)
        ks_null_pure <- ks.test(z_np, cond_null_cdf)
        null_D_pure  <- unname(ks_null_pure$statistic)
      }
    }

    if (n_nonnull_pure >= 2) {
      z_np   <- sort(z_nn_pure)
      p_high <- pbeta(z_np[length(z_np)], ahat, bhat)
      if (p_high > .Machine$double.eps) {
        cond_nonnull_cdf <- function(z) pbeta(z, ahat, bhat) / p_high
        ks_nonnull_pure  <- ks.test(z_np, cond_nonnull_cdf)
        nonnull_D_pure   <- unname(ks_nonnull_pure$statistic)
      }
    }

    ## ── Multimodality in non-null ─────────────────────────────────────────────
    dens_nn         <- NULL
    n_modes         <- NA_integer_
    peaks_above_thr <- integer(0)
    if (n_nonnull >= 10) {
      dens_nn        <- density(z_j[nonnull_idx])
      peak_threshold <- 0.10 * max(dens_nn$y)
      dy             <- dens_nn$y
      n_d            <- length(dy)
      is_peak        <- c(FALSE,
                          dy[seq(2, n_d - 1)] > dy[seq(1, n_d - 2)] &
                          dy[seq(2, n_d - 1)] > dy[seq(3, n_d)],
                          FALSE)
      peaks_above_thr <- which(is_peak & dy >= peak_threshold)
      n_modes         <- length(peaks_above_thr)
    }

    ## ── Heavy left tail ───────────────────────────────────────────────────────
    left_tail_frac  <- if (n_nonnull > 0) mean(z_j[nonnull_idx] < 0.05) else 0
    heavy_left_tail <- left_tail_frac > 0.20

    ## ── Deviation flags ───────────────────────────────────────────────────────
    null_D    <- if (!is.null(ks_null))    unname(ks_null$statistic)    else NA_real_
    nonnull_D <- if (!is.null(ks_nonnull)) unname(ks_nonnull$statistic) else NA_real_
    null_dev    <- !is.na(null_D)    && null_D    > 0.05
    nonnull_dev <- !is.na(nonnull_D) && nonnull_D > 0.05
    multimodal  <- !is.na(n_modes) && n_modes > 1
    no_signal   <- edges == 0L
    extreme_p0  <- p0 < 0.01

    n_issues   <- sum(null_dev, nonnull_dev, multimodal, heavy_left_tail, no_signal, extreme_p0)
    assessment <- if (n_issues == 0) "Good fit" else
                  if (n_issues == 1) "Acceptable fit" else
                  "Deviations detected"

    ## ── 2x2 plot panel ────────────────────────────────────────────────────────
    op <- par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1))
    on.exit(par(op), add = TRUE)

    ## Panel 1: standard fitted mixture
    plotFittedBetaMix(betaMixObj, yLim = yLim)
    title(main = "Fitted mixture", line = 0.5)

    ## Panel 2: Q-Q plot, null component
    # Prefer the pure-null (m0 >= qq_purity) conditional Q-Q; the reference
    # distribution is Beta(etahat, 0.5) truncated to [z_min, 1] so the plot
    # covers only the part of the null distribution free of non-null overlap.
    use_null_pure <- n_null_pure >= 10
    if (use_null_pure || n_null >= 2) {
      if (use_null_pure) {
        z_null_sorted <- sort(z_j[null_pure_idx])
        n_qq_null     <- n_null_pure
        # Conditional quantiles: spread ppoints uniformly over [p_low, 1]
        p_low_qq      <- min(pbeta(z_null_sorted[1L], etahat, 0.5),
                             1 - .Machine$double.eps * 100)
        q_theo_null   <- qbeta(p_low_qq + ppoints(n_qq_null) * (1 - p_low_qq),
                               etahat, 0.5)
        d_null_show   <- null_D_pure
        qq_note_null  <- sprintf("m0 >= %.2f  (n = %d of %d null-assigned)",
                                 qq_purity, n_null_pure, n_null)
      } else {
        z_null_sorted <- sort(z_j[null_idx])
        n_qq_null     <- n_null
        q_theo_null   <- qbeta(ppoints(n_qq_null), etahat, 0.5)
        d_null_show   <- null_D
        qq_note_null  <- sprintf("m0 >= 0.50 (fallback, n = %d)", n_null)
      }
      col_null <- if (!is.na(d_null_show) && d_null_show > 0.05) "tomato" else "steelblue"
      plot(q_theo_null, z_null_sorted,
           pch = 16, cex = 0.4, col = col_null,
           xlab = sprintf("Theoretical quantiles  Beta(%.2f, 0.5)", etahat),
           ylab = "Sample quantiles (null)",
           main = sprintf("Q-Q: Null  (D = %s)",
                          if (is.na(d_null_show)) "NA"
                          else sprintf("%.4f", d_null_show)))
      mtext(qq_note_null, side = 3, line = 0.25, cex = 0.62, col = "grey35")
      abline(0, 1, lwd = 1.5, lty = 2, col = "grey40")
    } else {
      plot.new()
      text(0.5, 0.5, "Insufficient null\nobservations", cex = 1.1)
      title(main = "Q-Q: Null  (no data)")
    }

    ## Panel 3: Q-Q plot, non-null component (scaled by bmax)
    # Prefer the pure-non-null (m0 <= 1 - qq_purity) conditional Q-Q; the
    # reference distribution is Beta(ahat, bhat) truncated to [0, z_max] so
    # the plot covers only the lower tail free of null overlap.
    use_nonnull_pure <- n_nonnull_pure >= 10
    if (use_nonnull_pure || length(z_nn) >= 2) {
      if (use_nonnull_pure) {
        z_nn_sorted  <- sort(z_nn_pure)
        n_qq_nn      <- n_nonnull_pure
        # Conditional quantiles: spread ppoints uniformly over [0, p_high]
        p_high_qq    <- max(pbeta(z_nn_sorted[n_qq_nn], ahat, bhat),
                            .Machine$double.eps * 100)
        q_theo_nn    <- qbeta(ppoints(n_qq_nn) * p_high_qq, ahat, bhat)
        d_nn_show    <- nonnull_D_pure
        qq_note_nn   <- sprintf("m0 <= %.2f  (n = %d of %d non-null-assigned)",
                                1 - qq_purity, n_nonnull_pure, length(z_nn))
      } else {
        z_nn_sorted  <- sort(z_nn)
        n_qq_nn      <- length(z_nn)
        q_theo_nn    <- qbeta(ppoints(n_qq_nn), ahat, bhat)
        d_nn_show    <- nonnull_D
        qq_note_nn   <- sprintf("m0 < 0.50 (fallback, n = %d)", n_qq_nn)
      }
      col_nn <- if (!is.na(d_nn_show) && d_nn_show > 0.05) "tomato" else "steelblue"
      plot(q_theo_nn, z_nn_sorted,
           pch = 16, cex = 0.4, col = col_nn,
           xlab = sprintf("Theoretical quantiles  Beta(%.2f, %.2f)", ahat, bhat),
           ylab = sprintf("Sample quantiles (non-null, z_j/%.3f)", bmax),
           main = sprintf("Q-Q: Non-null  (D = %s)",
                          if (is.na(d_nn_show)) "NA"
                          else sprintf("%.4f", d_nn_show)))
      mtext(qq_note_nn, side = 3, line = 0.25, cex = 0.62, col = "grey35")
      abline(0, 1, lwd = 1.5, lty = 2, col = "grey40")
    } else {
      plot.new()
      text(0.5, 0.5, "Insufficient non-null\nobservations", cex = 1.1)
      title(main = "Q-Q: Non-null  (no data)")
    }

    ## Panel 4: model-based ROC curve (plot only, no console output)
    tau_grid   <- seq(0, 1, length.out = 2001)
    fpr_roc    <- pbeta(tau_grid, etahat, 0.5)
    tpr_roc    <- pbeta(pmin(tau_grid, bmax) / bmax, ahat, bhat)
    auc_roc    <- sum(diff(fpr_roc) * (tpr_roc[-1L] + tpr_roc[-length(tpr_roc)]) / 2)
    j_idx      <- which.max(tpr_roc - fpr_roc)
    tau_youden <- tau_grid[j_idx]
    fpr_youden <- fpr_roc[j_idx]
    tpr_youden <- tpr_roc[j_idx]
    fpr_cur    <- pbeta(ppthr, etahat, 0.5)
    tpr_cur    <- pbeta(min(ppthr, bmax) / bmax, ahat, bhat)

    plot(fpr_roc, tpr_roc, type = "l", lwd = 2, col = "steelblue",
         xlim = c(0, 1), ylim = c(0, 1),
         xlab = "False Positive Rate",
         ylab = "True Positive Rate",
         main = sprintf("ROC  (AUC = %.4f)", auc_roc))
    abline(0, 1, lty = 2, col = "grey60")
    abline(v = pbeta(bmax, etahat, 0.5), lty = 3, col = "grey70")
    points(fpr_cur,    tpr_cur,    pch = 21, bg = "orange", cex = 1.8)
    points(fpr_youden, tpr_youden, pch = 21, bg = "tomato",  cex = 1.8)
    legend("bottomright",
           c(sprintf("Current  tau=%.3g", ppthr),
             sprintf("Youden   tau=%.3g", tau_youden)),
           pch = 21, pt.bg = c("orange", "tomato"), pt.cex = 1.4,
           cex = 0.75, bty = "n")

    ## ── Verbal assessment ─────────────────────────────────────────────────────
    cat("=== assessFit: Goodness-of-Fit Summary ===\n")
    cat(sprintf("  Hard-assigned null: %d,  non-null: %d\n", n_null, n_nonnull))
    cat(sprintf("  Null component   KS D = %s  [%s]\n",
                if (is.na(null_D)) "NA" else sprintf("%.4f", null_D),
                if      (is.na(null_D)) "skipped"
                else if (null_dev)      "NOTABLE - D > 0.05"
                else                    "OK"))
    if (null_dev)
      cat("    Hint: null deviation may indicate dependent samples;",
          "try betaMix(..., ind = FALSE).\n")
    cat(sprintf("  Non-null component   KS D = %s  [%s]\n",
                if (is.na(nonnull_D)) "NA" else sprintf("%.4f", nonnull_D),
                if      (is.na(nonnull_D)) "skipped"
                else if (nonnull_dev)      "NOTABLE - D > 0.05"
                else                       "OK"))
    if (nonnull_dev)
      cat("    Hint: non-null deviation may indicate a multimodal or",
          "heavy-tailed non-null distribution.\n")
    cat(sprintf("  Conditional Q-Q (pure subsets, m0 >= %.2f / m0 <= %.2f):\n",
                qq_purity, 1 - qq_purity))
    cat(sprintf("    Null D = %s  (n = %d),  Non-null D = %s  (n = %d)\n",
                if (is.na(null_D_pure))    "NA" else sprintf("%.4f", null_D_pure),
                n_null_pure,
                if (is.na(nonnull_D_pure)) "NA" else sprintf("%.4f", nonnull_D_pure),
                n_nonnull_pure))
    if (!is.na(n_modes))
      cat(sprintf("  Non-null modes: %d  [%s]\n",
                  n_modes,
                  if (multimodal) "NOTABLE - multimodal non-null component" else "OK"))
    cat(sprintf("  Heavy left tail (z_j < 0.05): %.1f%%  [%s]\n",
                100 * left_tail_frac,
                if (heavy_left_tail) "NOTABLE - excess mass near 0" else "OK"))
    if (no_signal)
      cat("  No edges detected  [WARNING: no non-null signal found]\n")
    cat(sprintf("  p0 = %.4f  [%s]\n", p0,
                if   (p0 < 0.01) "WARNING: ill-conditioned null (nearly all pairs non-null)"
                else              "OK"))
    cat(sprintf("  Overall verdict: >>> %s <<<\n", assessment))
    cat("==========================================\n")

    invisible(list(
      ks_null         = ks_null,
      ks_nonnull      = ks_nonnull,
      ks_null_pure    = ks_null_pure,
      ks_nonnull_pure = ks_nonnull_pure,
      n_null          = n_null,
      n_nonnull       = n_nonnull,
      n_null_pure     = n_null_pure,
      n_nonnull_pure  = n_nonnull_pure,
      n_modes         = n_modes,
      assessment      = assessment
    ))
  })
}


#' Print a short summary of the fitted mixture model.
#'
#' Prints the estimated model parameters (ahat, bhat, etahat), the posterior
#' probability threshold, sample size, node count, maximum possible edges,
#' detected edges, and estimated null proportion p0.
#' @param betamixobj The object (list) returned from the betaMix function.
#' @seealso \code{\link{betaMix-package}} for an overview and index of all functions and datasets in this package.
#' @export
#' @examples
#' \dontrun{
#' data(SIM, package = "betaMix")
#' res <- betaMix(betaMix::SIM, maxalpha = 1e-6, ppr = 0.01, subsamplesize = 30000, ind = TRUE)
#' shortSummary(res)
#' }
shortSummary <- function(betamixobj) {
  with(betamixobj, {
    cat(paste0(
      "Nonnull support = [0,", format(bmax, digits = 3), "]\nahat = ",
      format(ahat, digits = 2), ", bhat = ", format(bhat, digits = 2),
      "\netahat = ", format(etahat, digits = 2),
      "\nPost. Pr. threshold = ", format(ppthr, digits = 2)
    ), "\n")
    cat("Sample size =", N, "\n")
    cat("No. nodes =", prettyNum(nodes, big.mark = ","), "\n")
    cat("Max no. edges =", prettyNum(choose(nodes, 2), big.mark = ","), "\n")
    cat("No. edges detected =", prettyNum(edges, big.mark = ","), "\n")
    cat("p0 =", format(p0, digits = 3), "\n")
  })
}


#' Find spherical caps with more than one node.
#'
#' Takes a PxP adjacency Matrix as input and find spherical caps with more than
#' one node in R^n. Cap centers can only appear in one cap, but other nodes are
#' allowed to appear in multiple caps.
#' @param A An adjacency Matrix(0/1), zeros on the main diagonal.
#' @return A data frame with the following columns:
#' \describe{
#' \item{node}{Node number (row/column in A.)}
#' \item{capNum}{The cap number.}
#' \item{isCtr}{1 if the node is the center of the cap, 0 otherwise.}
#' \item{deg}{Node degree.}
#' \item{cc}{Clustering coefficient.}
#' }
#' @importFrom Matrix Matrix rowSums diag
#' @importFrom methods is
#' @seealso \code{\link{betaMix-package}} for an overview and index of all functions and datasets in this package.
#' @export
#' @examples
#' \dontrun{
#' data(SIM, package = "betaMix")
#' res <- betaMix(betaMix::SIM, maxalpha = 1e-6, ppr = 0.01, subsamplesize = 30000, ind = TRUE)
#' adjMat <- getAdjMat(res)
#' caps <- sphericalCaps(adjMat)
#' head(caps)
#' am <- getAdjMat(res, signed = TRUE, nodes = caps$node[which(caps$capNum == 1)])
#' plotCluster(am, 1, edgecols = c("blue", "red"), labels = TRUE)
#' }
sphericalCaps <- function(A) {
  stopifnot(is(A, "Matrix"))
  A <- abs(A)
  diag(A) <- FALSE
  deg <- Matrix::rowSums(A)
  possibleCtrs <- which(deg > 0)
  if (length(possibleCtrs) == 0) {
    cat("No edges!\n")
    return(NULL)
  }
  CC <- clusteringCoef(A)
  chunks <- list()
  orddeg <- order(deg, decreasing = TRUE)
  capNum <- 1
  while (length(orddeg) > 0 && max(deg[orddeg]) > 0) {
    capCtr <- orddeg[1]
    if (deg[capCtr] == 0) {
      break
    }
    nbrs <- setdiff(which(A[capCtr, ] != 0), capCtr)
    orddeg <- setdiff(orddeg, union(capCtr, nbrs))
    tmpdf <- cbind(
      c(capCtr, nbrs),
      rep(capNum, length(c(capCtr, nbrs))),
      c(1, rep(0, length(nbrs))),
      deg[c(capCtr, nbrs)],
      CC[c(capCtr, nbrs)]
    )
    rownames(tmpdf) <- sprintf("%s_%04d", rownames(tmpdf), tmpdf[, 2])
    chunks[[capNum]] <- tmpdf
    capNum <- capNum + 1
    if (length(orddeg) == 0) {
      break
    }
  }
  retdf <- as.data.frame(do.call(rbind, chunks))
  colnames(retdf) <- c("node", "capNum", "isCtr", "deg", "cc")
  return(retdf)
}


#' Find clusters, and return node characteristics.
#'
#' Take an adjacency Matrix as input and find clusters. For each node, find the degree and clustering coefficient (CC). Then, calculate a centrality measure (type\*CC+1)\*deg. For type=0, it's just the degree. Note that setting type=1 means that we assign a higher value to nodes that not only have many neighbors, but the neighbors are highly interconnected. For example, suppose we have two components with k nodes, one has a star shape, and the other is a complete graph. With type=0 both graphs will get the same value, but with type=1 the complete graph will be picked by the algorithm first. Setting type to a negative value gives CC\*deg as the centrality measure.
#' @param A An adjacency Matrix(0/1).
#' @param minCtr The minimum centrality value to be considered for a cluster center (default=5).
#' @param type Determines how the centrality measure is computed: positive
#'   values use (type*CC+1)*deg; negative values use CC*deg; 0 gives plain
#'   degree.
#' @return A data frame with the following columns:
#' \describe{
#' \item{labels}{Node label (e.g. gene names).}
#' \item{degree}{Node degree.}
#' \item{cc}{Node clustering coefficient.}
#' \item{ctr}{Node centrality measure: (type*CC+1)*deg, or CC*deg if type is negative.}
#' \item{clustNo}{Cluster number.}
#' \item{iscenter}{1 if the node was chosen as the cluster's center, 0 otherwise.}
#' \item{intEdges}{Number of edges from the node to nodes in the same cluster.}
#' \item{extEdges}{Number of edges from the node to nodes NOT in the same cluster.}
#' \item{distCenter}{Standardized Manhattan distance to the central node.}
#' }
#' @importFrom methods is
#' @seealso \code{\link{betaMix-package}} for an overview and index of all functions and datasets in this package.
#' @export
#' @examples
#' \dontrun{
#' data(SIM, package = "betaMix")
#' res <- betaMix(betaMix::SIM, maxalpha = 1e-6, ppr = 0.01, subsamplesize = 30000, ind = TRUE)
#' adjMat <- getAdjMat(res)
#' SimComp <- graphComponents(adjMat)
#' head(SimComp)
#' }
graphComponents <- function(A, minCtr = 5, type = 1) {
  stopifnot(is(A, "Matrix"))
  A <- abs(A)
  Vn <- ncol(A)
  ctrs <- rep(2 * Vn, Vn)
  labels <- 1:Vn
  if (!is.null(rownames(A))) {
    labels <- rownames(A)
  }
  deg <- Matrix::rowSums(A)
  CC <- clusteringCoef(A)
  ctrs <- (type * CC + 1) * deg
  if (type < 0) {
    ctrs <- CC * deg
  }
  clustersInfo <- data.frame(
    labels = labels, degree = deg, cc = CC, ctr = ctrs,
    clustNo = rep(0, Vn), iscenter = rep(0, Vn),
    intEdges = rep(0, Vn), extEdges = rep(0, Vn),
    distCenter = rep(0, Vn)
  )
  clustNo <- 1
  clustered <- which(deg < 1)
  while (length(clustered) < Vn) {
    notInCluster <- setdiff(1:Vn, clustered)
    if (max(ctrs[notInCluster]) < minCtr) {
      return(clustersInfo)
    }
    ctrnode <- notInCluster[which.max(ctrs[notInCluster])]
    # candidate cluster neighbors
    nbrs <- setdiff(sort(c(ctrnode, which(A[ctrnode, ] != 0))), clustered)
    if (length(nbrs) > 1) {
      if (length(nbrs) > minCtr) {
        clustersInfo$iscenter[ctrnode] <- 1
        clustersInfo$clustNo[union(ctrnode, nbrs)] <- clustNo
        clustersInfo$intEdges[nbrs] <- Matrix::rowSums(A[nbrs, nbrs])
        if (length(nbrs) < ncol(A)) {
          clustersInfo$extEdges[nbrs] <- Matrix::rowSums(A[nbrs, -nbrs])
        } else {
          clustersInfo$extEdges[nbrs] <- 0
        }
        for (i in seq_along(nbrs)) {
          clustersInfo$distCenter[nbrs[i]] <- mean(xor(A[ctrnode, ], A[nbrs[i], ]))
        }
        clustNo <- clustNo + 1
      } else {
        nbrs <- c()
      }
    } else {
      nbrs <- c()
    }
    clustered <- union(clustered, c(nbrs, ctrnode))
  }
  return(clustersInfo)
}


#' Show cluster characteristics.
#'
#' Takes an object obtained from graphComponents and prints and returns summary statistics.
#' @param clustersInfo Obtained from graphComponents.
#' @return A matrix with cluster number, number of nodes, and fivenum summaries for the degrees of nodes in the cluster, and the percentage of edges that are within the cluster.
#' @seealso \code{\link{betaMix-package}} for an overview and index of all functions and datasets in this package.
#' @export
#' @examples
#' \dontrun{
#' data(SIM, package = "betaMix")
#' res <- betaMix(betaMix::SIM, maxalpha = 1e-6, ppr = 0.01, subsamplesize = 30000, ind = TRUE)
#' adjMat <- getAdjMat(res)
#' SimComp <- graphComponents(adjMat)
#' (summtab <- summarizeClusters(SimComp))
#' }
summarizeClusters <- function(clustersInfo) {
  cat("Num of nodes:", nrow(clustersInfo), "\n")
  cat("Num of edges:", sum(clustersInfo$degree) / 2, "\n")
  cat("Num of clusters:", max(clustersInfo$clustNo), "\n")
  cat("Num of unclustered nodes:", length(which(clustersInfo$clustNo == 0)), "\n")
  percentInCluster <- ifelse(clustersInfo$degree > 0, clustersInfo$intEdges / clustersInfo$degree, 0)
  if (max(clustersInfo$clustNo) == 0) {
    return(NULL)
  }
  tab <- matrix(0, nrow = max(clustersInfo$clustNo), ncol = 12)
  for (cnum in 1:max(clustersInfo$clustNo)) {
    idx <- clustersInfo$clustNo == cnum
    tab[cnum, ] <- c(
      cnum, sum(idx), fivenum(clustersInfo$degree[idx]),
      fivenum(percentInCluster[idx])
    )
  }
  colnames(tab) <- c(
    "Cluster", "Nodes", "degreeMin", "degreeQ25", "degreeMedian",
    "degreeQ75", "degreeMax", "pctInClstMin", "pctInClstQ25",
    "pctInClstMedian", "pctInClstQ75", "pctInClstMax"
  )
  tab
}


#' Reduce an adjacency matrix to include only central nodes and nodes not assigned to any cluster.
#'
#' Return an adjacency matrix after collapsing clusters into their central nodes.
#' @param A An adjacency Matrix.
#' @param clustersInfo Obtained from graphComponents
#' @return A weighted adjacency matrix between clusters and unclustered nodes.
#' @seealso \code{\link{betaMix-package}} for an overview and index of all functions and datasets in this package.
#' @export
#' @examples
#' \dontrun{
#' data(DrySeeds)
#' res <- betaMix(DrySeeds, maxalpha = 1e-4, ppr = 0.01, ind = TRUE)
#' adjMat <- getAdjMat(res)
#' rownames(adjMat) <- colnames(DrySeeds)
#' SimComp <- graphComponents(adjMat)
#' Adj1 <- collapsedGraph(adjMat, SimComp) > 0
#' # with the igraph package, we can do:
#' # plot(graph.adjacency(Adj1, mode="undirected"), vertex.label=rownames(Adj1),
#' # vertex.label.cex=0.7, vertex.size=0, vertex.frame.color="white",
#' # vertex.label.color='blue', edge.color='grey80',  asp=1)
#' }
collapsedGraph <- function(A, clustersInfo) {
  collDim <- length(which(clustersInfo$clustNo == 0)) + max(clustersInfo$clustNo)
  collA <- Matrix::Matrix(0, ncol = collDim, nrow = collDim)
  notInCluster <- which(clustersInfo$clustNo == 0)
  if (length(notInCluster) > 0) {
    ni <- length(notInCluster)
    collA[seq_len(ni), seq_len(ni)] <- A[notInCluster, notInCluster] > 0
  }
  if (length(rownames(A)) != nrow(A)) {
    rownames(A) <- seq_len(nrow(A))
  }
  nc_max <- max(clustersInfo$clustNo)
  rownames(collA) <- c(
    rownames(A)[notInCluster],
    paste0("CLS", seq_len(nc_max))
  )
  cluster_idx <- lapply(seq_len(nc_max), function(k) which(clustersInfo$clustNo == k))
  ni <- length(notInCluster)
  for (i in seq_len(nc_max)) {
    Ci <- cluster_idx[[i]]
    if (ni > 0) {
      Atmp <- matrix(A[notInCluster, Ci], nrow = ni, ncol = length(Ci))
      collA[i + ni, seq_len(ni)] <- Matrix::rowSums(Atmp)
    }
    if (i < nc_max) {
      for (j in (i + 1):nc_max) {
        collA[i + ni, j + ni] <- sum(A[Ci, cluster_idx[[j]]])
      }
    }
  }
  collA + Matrix::t(collA)
}


#' Calculate the clustering coefficient of each node.
#'
#' @param A an adjacency Matrix (0/1).
#' @return A vector with the clustering coefficient of each node.
#' @seealso \code{\link{betaMix-package}} for an overview and index of all functions and datasets in this package.
#' @export
#' @examples
#' \dontrun{
#' data(SIM, package = "betaMix")
#' res <- betaMix(betaMix::SIM, maxalpha = 1e-6, ppr = 0.01, subsamplesize = 30000, ind = TRUE)
#' adjMat <- getAdjMat(res)
#' clusteringCoef(adjMat)
#' }
#'
clusteringCoef <- function(A) {
  rsum <- Matrix::rowSums(A)
  cc <- rep(0, nrow(A))
  for (i in seq_len(nrow(A))) {
    if (rsum[i] <= 1) {
      cc[i] <- 0
    } else {
      nbrs <- which(A[i, ] != 0)
      At <- A[nbrs, nbrs]
      cc[i] <- 0.5 * sum(At) / choose(rsum[i], 2)
    }
  }
  cc
}


#' Plot the degree of nodes versus the degree times the clustering coefficient.
#'
#' The x-axis represents node degree (number of neighbors). The y-axis
#' represents degree multiplied by the clustering coefficient (CC*degree), a
#' centrality-like measure that weights high-degree nodes by how interconnected
#' their neighborhood is.
#' @param betamixobj The object (list) returned by betaMix
#' @param clusterInfo obtained from graphComponents. If not provided by the user, it will be computed on the fly.
#' @param highlightNodes A vector of node-numbers which will be shown in red. Default is NULL.
#' @seealso \code{\link{betaMix-package}} for an overview and index of all functions and datasets in this package.
#' @export
#' @import stats graphics
#' @examples
#' \dontrun{
#' data(SIM, package = "betaMix")
#' res <- betaMix(betaMix::SIM, maxalpha = 1e-6, ppr = 0.01, subsamplesize = 30000, ind = TRUE)
#' adjMat <- getAdjMat(res)
#' SimComp <- graphComponents(adjMat)
#' plotDegCC(res, SimComp)
#' }
plotDegCC <- function(betamixobj, clusterInfo = NULL, highlightNodes = NULL) {
  if (is.null(clusterInfo)) {
    clusterInfo <- graphComponents(getAdjMat(betamixobj))
  }
  cc0 <- clusterInfo$cc
  deg0 <- clusterInfo$degree
  if (length(deg0) == 0 || max(deg0) == 0) {
    warning("plotDegCC: no edges in graph — nothing to plot")
    return(invisible(NULL))
  }
  plot(deg0, deg0 * cc0,
    axes = FALSE, xlim = c(0, max(deg0)),
    ylim = c(0, 1.1 * max(deg0 * cc0)), main = "",
    xlab = bquote("degree"), ylab = bquote("CC*degree"),
    col = "thistle", pch = 24, cex = 0.5
  )
  axis(1)
  axis(2)
  grid()
  abline(0, 1, col = "seagreen1", lwd = 2)
  if (!is.null(highlightNodes)) {
    points(deg0[highlightNodes], (deg0 * cc0)[highlightNodes], col = 2, pch = 24, cex = 0.5)
  }
}


#' Edge-indicator bitmap plot.
#'
#' Plot a bitmap in which a black dot corresponds to a pair of highly correlated genes (an edge in the graph).
#' The default is to show the nodes according to their order in the input.
#' By setting orderByCluster=TRUE, it is possible to cluster nodes and show
#' them in increasing degree order (from left to right).
#' @param AdjMat An adjacency Matrix (0/1).
#' @param clusterInfo obtained from graphComponents. If not provided by the user, it will be computed on the fly.
#' @param orderByCluster If FALSE, show the bitmap in the original node order.
#'   If TRUE, show nodes by clusters, sorted by distance from the cluster center.
#' @param showMinDegree Non-negative integer indicating the minimum degree of nodes that should be displayed. Default=0 (all nodes).
#' @seealso \code{\link{betaMix-package}} for an overview and index of all functions and datasets in this package.
#' @export
#' @examples
#' \dontrun{
#' data(SIM, package = "betaMix")
#' res <- betaMix(betaMix::SIM, maxalpha = 1e-6, ppr = 0.01, subsamplesize = 30000, ind = TRUE)
#' adjMat <- getAdjMat(res)
#' SimComp <- graphComponents(adjMat)
#' plotBitmapCC(adjMat)
#' plotBitmapCC(adjMat, SimComp, orderByCluster = TRUE)
#' }
plotBitmapCC <- function(AdjMat, clusterInfo = NULL, orderByCluster = FALSE, showMinDegree = 0) {
  if (!is.null(clusterInfo)) {
    orderByCluster <- TRUE
  }
  if (orderByCluster) {
    if (is.null(clusterInfo)) {
      clusterInfo <- graphComponents(AdjMat)
    }
    nodeOrder <- order(clusterInfo$clustNo, clusterInfo$distCenter)
    AdjMat <- AdjMat[nodeOrder, nodeOrder]
  }
  showNodes <- which(Matrix::rowSums(AdjMat) >= showMinDegree)
  Matrix::image(AdjMat[showNodes, showNodes])
}


#' Plot cluster network
#'
#' Plot a cluster network with all the nodes and edges - the central node is marked by a black circle. The radius of each point corresponds to its degree. The opacity corresponds to the percentage of edges from the node that is in the cluster (the darker it is, the larger the percentage of edges is within the cluster.) The distance from the center corresponds to the relative dissimilarity with the central node. This is computed as the number of neighbors the node and the central node do not have in common.
#' @param AdjMat An adjacency Matrix (0/1).
#' @param clustNo The chosen cluster.
#' @param clusterInfo Obtained from graphComponents.
#' @param labels If set to TRUE, show node names (default=FALSE).
#' @param nodecol The color(s) of the nodes. Can be a single value or a vector of length equal to the number of rows in AdjMat
#' @param labelsize Text size of node labels.
#' @param figtitle The title of the plot (default=NULL).
#' @param edgecols The colors to be used for edges. Default="grey88". If one
#'   value is given, all edges will be drawn using this color. If edgecols
#'   contains two valid colors, the first is used for positive correlations and
#'   the second for negative ones.
#' @seealso \code{\link{betaMix-package}} for an overview and index of all functions and datasets in this package.
#' @export
#' @importFrom grDevices col2rgb colours rgb
#' @examples
#' \dontrun{
#' data(SIM, package = "betaMix")
#' res <- betaMix(betaMix::SIM, maxalpha = 1e-5, ppr = 0.01, subsamplesize = 30000, ind = TRUE)
#' adjMat <- getAdjMat(res, signed = TRUE)
#' SimComp <- graphComponents(adjMat)
#' plotCluster(adjMat, 2, SimComp)
#' plotCluster(adjMat, 2, SimComp, edgecols = c("blue", "orange"), labels = TRUE)
#' }
plotCluster <- function(AdjMat, clustNo, clusterInfo = NULL, labels = FALSE, nodecol = "blue", labelsize = 1, figtitle = NULL, edgecols = "grey88") {
  if (is.null(clusterInfo)) {
    clusterInfo <- graphComponents(AdjMat, minCtr = 2, type = 0)
  }
  if (length(nodecol) < nrow(AdjMat)) {
    nodecol <- rep(nodecol[1], length = nrow(AdjMat))
  }
  ids <- which(clusterInfo$clustNo == clustNo)
  if (length(ids) > 0) {
    tmpA <- AdjMat[ids, ids]
    tmpclusterInfo <- clusterInfo[ids, ]
    rads <- round(10 * tmpclusterInfo$distCenter /
      max(0.1, max(tmpclusterInfo$distCenter)))
    thetas <- rep(0, length(rads))
    intvls <- findInterval(rads, seq(1, 10))
    for (intvl in unique(sort(intvls))) {
      pts <- which(intvls == intvl)
      thetas[pts] <- 3 * intvl * pi / max(0.01, max(intvls)) + seq(0, 1.9 * pi, length = length(pts))
    }
    max_deg <- max(tmpclusterInfo$degree)
    if (max_deg == 0) {
      warning("plotCluster: all nodes in cluster ", clustNo, " have degree 0; cannot render plot")
      return(invisible(NULL))
    }
    sizes <- pmax(0.3, tmpclusterInfo$degree / max_deg)
    opacity <- ifelse(tmpclusterInfo$degree > 0, 0.25 + tmpclusterInfo$intEdges / tmpclusterInfo$degree, 0.25)
    opacity <- opacity / max(opacity)
    nodecol <- rgb(t(col2rgb(nodecol[ids]) / 255), alpha = opacity)
    plot(rads * cos(thetas), rads * sin(thetas),
      cex = sizes * 3, pch = 19, axes = FALSE,
      xlab = "", ylab = "", col = nodecol, main = figtitle,
      ylim = c(min(rads * sin(thetas)), 1.1 * max(rads * sin(thetas)))
    )
    edgecol <- rep(edgecols[1], ncol(tmpA))
    has_neg_col <- length(edgecols) >= 2 && edgecols[2] %in% colours()
    for (i in seq_len(ncol(tmpA))) {
      nbrs <- setdiff(which(abs(tmpA[i, ]) == 1), seq_len(i))
      if (length(nbrs) > 0) {
        neg_nbrs <- if (has_neg_col) nbrs[which(tmpA[i, nbrs] == -1)] else integer(0)
        edgecol[neg_nbrs] <- edgecols[2]
        for (j in nbrs) {
          lines(c(rads[i] * cos(thetas[i]), rads[j] * cos(thetas[j])),
            c(rads[i] * sin(thetas[i]), rads[j] * sin(thetas[j])),
            col = edgecol[j], lwd = 0.5
          )
        }
        edgecol[neg_nbrs] <- edgecols[1]
      }
    }
    points(rads * cos(thetas), rads * sin(thetas), cex = sizes * 3, pch = 19, col = nodecol)
    if (labels) {
      text(rads * cos(thetas), rads * sin(thetas), tmpclusterInfo$labels, pos = 3, cex = labelsize)
    }
    ctr <- which(tmpclusterInfo$iscenter == 1)
    points(rads[ctr] * cos(thetas[ctr]), rads[ctr] * sin(thetas[ctr]),
      pch = 21,
      cex = sizes[ctr] * 3, col = "black", lwd = 2
    )
  } else {
    cat("Invalid cluster number\n")
  }
}


#' Return a Matrix with the shortest path distance between nodes.
#'
#' Returns a matrix whose entry (i,j) is the minimum number of steps from node
#' i to node j, considering paths up to numSteps edges long.
#' @param AdjMat An adjacency Matrix (0/1).
#' @param numSteps The maximum number of edges between pairs of nodes. If numSteps=0, returns the input matrix. numSteps=1 adds neighbors of direct neighbors, etc.
#' @return A Matrix containing the shortest path lengths between nodes i and j
#' @seealso \code{\link{betaMix-package}} for an overview and index of all functions and datasets in this package.
#' @export
#' @examples
#' \dontrun{
#' data(SIM, package = "betaMix")
#' res <- betaMix(betaMix::SIM, maxalpha = 1e-6, ppr = 0.01, subsamplesize = 30000, ind = TRUE)
#' adjMat <- getAdjMat(res)
#' AdjMat <- shortestPathDistance(adjMat, numSteps = 2)
#' Matrix::image((AdjMat > 0)[1:200, 1:200])
#' adjMat1 <- AdjMat > 0
#' SimComp <- graphComponents(adjMat1)
#' head(summarizeClusters(SimComp))
#' }
shortestPathDistance <- function(AdjMat, numSteps = 0) {
  if (!is.numeric(numSteps) || length(numSteps) != 1 || numSteps != round(numSteps) || numSteps < 0)
    stop("shortestPathDistance: 'numSteps' must be a non-negative integer; got ", numSteps)
  if (numSteps == 0) {
    return(AdjMat)
  }
  Ap <- minDist <- AdjMat
  for (i in 1:numSteps) {
    An <- Ap %*% AdjMat
    if (sum((An | Ap) - (An & Ap)) == 0) {
      break
    }
    minDist[(An > 0) & (Ap == 0) & (minDist == 0)] <- i
    Ap <- An
  }
  rownames(minDist) <- colnames(minDist) <- rownames(AdjMat)
  minDist
}



#' Metabolite Expression data for the dry seed group.
#'
#' DrySeeds is a matrix with normalized metabolite data containing 68 named metabolites, and 50 samples.
#'
#' @docType data
#' @keywords datasets
#' @name DrySeeds
#' @format A matrix with 50 rows and 68 columns
#' @references \url{https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST000561&StudyType=MS&ResultType=1}
NULL

#' Metabolite Expression data for the 6-hour imbibed seed group.
#'
#' SixHourImbibed is a matrix with normalized metabolite data containing 68 named metabolites, and 50 samples.
#'
#' @docType data
#' @keywords datasets
#' @name SixHourImbibed
#' @format A matrix with 50 rows and 68 columns
#' @references \url{https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST000561&StudyType=MS&ResultType=1}
NULL


#' Simulated gene Expression data using the huge package
#'
#' SIM is a simulated dataset with a hub structure, consisting of 1000 nodes and 50 hubs.
#'
#' @docType data
#' @keywords datasets
#' @name SIM
#' @usage data(SIM,package = "betaMix")
#' @format A 1000 by 200 matrix, representing 50 hubs
NULL
