#' Cluster SNPs Based on Genomic Distance and P-value
#'
#' This function clusters SNPs by chromosome and genomic distance. Within each cluster,
#' the SNP with the smallest p-value is marked as the "top" SNP.
#'
#' @param data A data.frame containing SNP information.
#' @param chr_col Character. Name of the chromosome column. Default is "chr".
#' @param pos_col Character. Name of the position column. Default is "pos".
#' @param rsid_col Character. Name of the SNP ID column. Default is "rsid".
#' @param pvalue_col Character. Name of the p-value column. Default is "pvalue".
#' @param gene_col Character. Name of the gene column (optional). Default is "gene".
#' @param pvalue_threshold Numeric. Maximum p-value for SNPs to be considered. Default is 0.05.
#' @param distance_threshold Numeric. Maximum distance (bp) between SNPs to belong to the same cluster. Default is 1e6.
#'
#' @return A data.frame identical to \code{data}, with two new columns:
#' \code{cluster} (cluster assignment) and \code{top} (indicates the top SNP in each cluster).
#'
#' @import igraph
#' @export
#'
#' @examples
#' \dontrun{
#' clustered_snps <- cluster_snps(data = snp_data)
#' }
cluster_snps <- function(data, chr_col = "chr", pos_col = "pos", rsid_col = "rsid",
                         pvalue_col = "pvalue", gene_col = "gene",
                         pvalue_threshold = 0.05, distance_threshold = 1e6) {
  # Ensure required package is installed
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required. Install it via install.packages('igraph').")
  }
  # Copy original data
  result_data <- data
  result_data$cluster <- NA_character_
  result_data$top <- NA_character_
  # Filter SNPs by p-value threshold
  filtered_data <- data[data[[pvalue_col]] <= pvalue_threshold, ]
  # If no SNPs pass threshold, return original data
  if (nrow(filtered_data) == 0) return(result_data)
  # Process SNPs chromosome by chromosome
  for (chr in unique(filtered_data[[chr_col]])) {
    chr_snps <- filtered_data[filtered_data[[chr_col]] == chr, ]
    # Handle single SNP case
    if (nrow(chr_snps) == 1) {
      cluster_id <- paste0("chr", chr, "_cluster1")
      idx <- which(result_data[[rsid_col]] == chr_snps[[rsid_col]])
      result_data$cluster[idx] <- cluster_id
      result_data$top[idx] <- "top"
      next
    }
    # Compute distance matrix
    positions <- chr_snps[[pos_col]]
    dist_matrix <- abs(outer(positions, positions, "-"))
    # Build adjacency matrix based on distance threshold
    adj_matrix <- ifelse(dist_matrix <= distance_threshold, 1, 0)
    diag(adj_matrix) <- 0  # Remove self-loops
    # Construct graph and identify clusters
    g <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
    clusters <- igraph::components(g)
    # Assign cluster IDs and identify top SNPs
    for (i in seq_len(clusters$no)) {
      cluster_indices <- which(clusters$membership == i)
      top_idx <- cluster_indices[which.min(chr_snps[[pvalue_col]][cluster_indices])]
      cluster_id <- paste0("chr", chr, "_cluster", i)
      for (idx in cluster_indices) {
        rsid <- chr_snps[[rsid_col]][idx]
        result_idx <- which(result_data[[rsid_col]] == rsid)
        result_data$cluster[result_idx] <- cluster_id
        if (idx == top_idx) result_data$top[result_idx] <- "top"
      }
    }
  }
  return(result_data)
}
