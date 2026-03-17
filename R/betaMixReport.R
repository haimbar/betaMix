#' Generate a betaMix analysis report
#'
#' Runs the full betaMix analysis pipeline on a data matrix and renders a
#' self-contained PDF report that includes dataset overview, summary
#' statistics, fitted mixture model parameters, all standard diagnostic plots,
#' and a network plot for a representative cluster.
#'
#' @param M Numeric matrix (N × P) of expression/measurement data (samples in
#'   rows, variables in columns).
#' @param dataset_name Character string used in the report title and output
#'   filename.  Defaults to the name of the object passed as \code{M}.
#' @param output_dir Directory where the PDF will be written.  Created
#'   automatically if it does not exist.  Default \code{"reports"}.
#' @param ind Logical; passed to \code{\link{betaMix}}.  Default \code{TRUE}.
#' @param maxalpha Numeric; passed to \code{\link{betaMix}}.  Default
#'   \code{1e-4}.
#' @param ppr Numeric; passed to \code{\link{betaMix}}.  Default \code{0.05}.
#' @param method Correlation method: \code{"pearson"} (default) or
#'   \code{"spearman"}.
#' @param subsamplesize Integer; passed to \code{\link{betaMix}}.  Default
#'   \code{50000}.
#' @param signed Logical; whether to use signed edges (positive/negative
#'   correlations distinguished by colour).  Default \code{FALSE}.
#' @param nodecol Node colour for the cluster network plot.  Default
#'   \code{"steelblue"}.
#' @param edgecols Edge colour(s) for the cluster network plot.  A single
#'   colour for all edges, or two colours for positive/negative correlations
#'   (requires \code{signed = TRUE}).  Default \code{"grey55"}.
#' @param rep_cluster Integer; which cluster to display in the network plot.
#'   If \code{NULL} (default), the function picks the largest cluster whose
#'   size falls in the range [5, 30]; if no cluster satisfies this, the
#'   largest cluster overall is used.
#' @return Invisibly, the path to the generated PDF file.
#' @export
#' @examples
#' \dontrun{
#' data(SixHourImbibed, package = "betaMix")
#' betaMixReport(SixHourImbibed, dataset_name = "SixHourImbibed",
#'               output_dir = "reports", ind = TRUE)
#' }
betaMixReport <- function(M,
                           dataset_name  = deparse(substitute(M)),
                           output_dir    = "reports",
                           ind           = TRUE,
                           maxalpha      = 1e-4,
                           ppr           = 0.05,
                           method        = "pearson",
                           subsamplesize = 50000,
                           signed        = FALSE,
                           nodecol       = "steelblue",
                           edgecols      = "grey55",
                           rep_cluster   = NULL) {

  if (!requireNamespace("rmarkdown", quietly = TRUE))
    stop("Package 'rmarkdown' is required. Install with install.packages('rmarkdown').")

  if (!dir.exists(output_dir))
    dir.create(output_dir, recursive = TRUE)

  # ── Run betaMix pipeline ────────────────────────────────────────────────────
  message("Running betaMix on '", dataset_name, "'...")
  res <- betaMix(M, ind = ind, maxalpha = maxalpha, ppr = ppr,
                 method = method, subsamplesize = subsamplesize, msg = FALSE)

  message("Building adjacency matrix and graph components...")
  adjMat    <- getAdjMat(res, signed = signed)
  graphComp <- graphComponents(adjMat, minCtr = 2, type = 1)

  # ── Pick representative cluster ─────────────────────────────────────────────
  # Prefer clusters with 5–30 nodes (best for visualisation); fall back to
  # the largest cluster if none fall in that range.
  if (is.null(rep_cluster)) {
    clust_nos    <- graphComp$clustNo[graphComp$clustNo > 0]
    cluster_sizes <- if (length(clust_nos) > 0) table(clust_nos) else integer(0)
    if (length(cluster_sizes) > 0) {
      sweet <- cluster_sizes[cluster_sizes >= 5 & cluster_sizes <= 30]
      if (length(sweet) > 0) {
        rep_cluster <- as.integer(names(which.max(sweet)))
      } else {
        rep_cluster <- as.integer(names(which.max(cluster_sizes)))
      }
    } else {
      rep_cluster <- NA_integer_
    }
  }

  # ── Per-column summary statistics (computed once, passed to template) ───────
  col_stats <- data.frame(
    Min    = apply(M, 2, min,             na.rm = TRUE),
    Q1     = apply(M, 2, quantile, 0.25, na.rm = TRUE),
    Median = apply(M, 2, median,         na.rm = TRUE),
    Mean   = apply(M, 2, mean,           na.rm = TRUE),
    Q3     = apply(M, 2, quantile, 0.75, na.rm = TRUE),
    Max    = apply(M, 2, max,            na.rm = TRUE),
    SD     = apply(M, 2, sd,             na.rm = TRUE),
    row.names = colnames(M)
  )

  # ── Locate template and render ───────────────────────────────────────────────
  template <- system.file("report_template.Rmd", package = "betaMix")
  if (!nzchar(template))
    stop("Report template not found. Re-install betaMix from source.")

  safe_name <- gsub("[^A-Za-z0-9_-]", "_", dataset_name)
  out_file  <- file.path(normalizePath(output_dir, mustWork = FALSE),
                         paste0("betaMix_", safe_name, ".pdf"))

  message("Rendering report to: ", out_file)
  rmarkdown::render(
    template,
    output_file = out_file,
    params = list(
      dataset_name = dataset_name,
      M            = M,
      res          = res,
      adjMat       = adjMat,
      graphComp    = graphComp,
      col_stats    = col_stats,
      rep_cluster  = rep_cluster,
      signed       = signed,
      nodecol      = nodecol,
      edgecols     = edgecols
    ),
    quiet = TRUE
  )

  message("Done.")
  invisible(out_file)
}
