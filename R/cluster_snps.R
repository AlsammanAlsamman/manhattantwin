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


############## Other Versions ################
#' Cluster SNPs Based on Genomic Distance and P-value (Optimized)
#'
#' This function clusters SNPs by chromosome and genomic distance using an efficient
#' algorithm that avoids large distance matrices.
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
#' @export
#'
#' @examples
#' \dontrun{
#' clustered_snps <- cluster_snps_fast(data = snp_data)
#' }
cluster_snps_fast <- function(data, chr_col = "chr", pos_col = "pos", rsid_col = "rsid",
                              pvalue_col = "pvalue", gene_col = "gene",
                              pvalue_threshold = 0.05, distance_threshold = 1e6) {

  # Create result dataframe
  result_data <- data
  result_data$cluster <- NA_character_
  result_data$top <- NA_character_

  # Filter SNPs by p-value threshold
  filtered_data <- data[data[[pvalue_col]] <= pvalue_threshold, ]

  # If no SNPs pass threshold, return original data
  if (nrow(filtered_data) == 0) {
    return(result_data)
  }

  # Process each chromosome separately
  chromosomes <- unique(filtered_data[[chr_col]])

  for (chr in chromosomes) {
    # Get SNPs for this chromosome and sort by position
    chr_idx <- which(filtered_data[[chr_col]] == chr)
    if (length(chr_idx) == 0) next

    chr_data <- filtered_data[chr_idx, ]
    chr_data <- chr_data[order(chr_data[[pos_col]]), ]

    # Get original indices for mapping back
    original_indices <- which(data[[chr_col]] == chr & data[[pvalue_col]] <= pvalue_threshold)
    original_indices <- original_indices[order(chr_data[[pos_col]])]

    # Efficient clustering using single pass algorithm
    clusters <- list()
    current_cluster <- 1
    cluster_start <- 1

    if (nrow(chr_data) == 1) {
      # Single SNP case
      cluster_id <- paste0("chr", chr, "_cluster1")
      result_data$cluster[original_indices[1]] <- cluster_id
      result_data$top[original_indices[1]] <- "top"
      next
    }

    # Initialize first cluster
    clusters[[current_cluster]] <- 1

    # Cluster SNPs based on distance threshold
    for (i in 2:nrow(chr_data)) {
      current_pos <- chr_data[[pos_col]][i]
      prev_pos <- chr_data[[pos_col]][i-1]

      if (current_pos - chr_data[[pos_col]][cluster_start] <= distance_threshold) {
        # SNP belongs to current cluster
        clusters[[current_cluster]] <- c(clusters[[current_cluster]], i)
      } else {
        # Start new cluster
        current_cluster <- current_cluster + 1
        clusters[[current_cluster]] <- i
        cluster_start <- i
      }
    }

    # Assign clusters and identify top SNPs
    for (cluster_num in seq_along(clusters)) {
      cluster_indices <- clusters[[cluster_num]]

      if (length(cluster_indices) == 0) next

      # Find SNP with smallest p-value in cluster
      cluster_pvalues <- chr_data[[pvalue_col]][cluster_indices]
      top_idx_in_cluster <- cluster_indices[which.min(cluster_pvalues)]

      cluster_id <- paste0("chr", chr, "_cluster", cluster_num)

      # Map back to original indices
      original_cluster_indices <- original_indices[cluster_indices]
      original_top_idx <- original_indices[top_idx_in_cluster]

      # Assign cluster and top SNP
      result_data$cluster[original_cluster_indices] <- cluster_id
      result_data$top[original_top_idx] <- "top"
    }
  }

  return(result_data)
}


########## Other Versions ##########
#' Cluster SNPs Based on Genomic Distance and P-value (Data.Table Version)
#'
#' This function clusters SNPs by chromosome and genomic distance using an efficient
#' algorithm that avoids large distance matrices.
#'
#' @param data A data.frame containing SNP information.
#' @param chr_col Character. Name of the chromosome column. Default is "chr".
#' @param pos_col Character. Name of the position column. Default is "pos".
#' @param rsid_col Character. Name of the SNP ID column. Default is "rsid".
#' @param pvalue_col Character. Name of the p-value column. Default is "pvalue".
#' @param gene_col Character. Name of the gene column (optional). Default is "gene".
#' @param pvalue_threshold Numeric. Maximum p-value for SNPs to be considered. Default is 0.05.
#' @param distance_threshold Numeric. Maximum distance (bp) between SNPs to belong to the same cluster. Default is 1e6.
#' @return A data.table with cluster assignments.
#'
#' @import data.table
#' @export
cluster_snps_dt <- function(data, chr_col = "chr", pos_col = "pos", rsid_col = "rsid",
                            pvalue_col = "pvalue", gene_col = "gene",
                            pvalue_threshold = 0.05, distance_threshold = 1e6) {

  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required for this function.")
  }

  # Convert to data.table
  dt <- data.table::as.data.table(data)

  # Add cluster and top columns
  dt[, cluster := NA_character_]
  dt[, top := NA_character_]

  # Filter by p-value
  filtered_dt <- dt[get(pvalue_col) <= pvalue_threshold]

  if (nrow(filtered_dt) == 0) {
    return(dt)
  }

  # Process by chromosome
  chromosomes <- unique(filtered_dt[[chr_col]])

  for (chr in chromosomes) {
    chr_dt <- filtered_dt[get(chr_col) == chr]
    chr_dt <- chr_dt[order(get(pos_col))]

    if (nrow(chr_dt) == 0) next

    # Efficient clustering
    positions <- chr_dt[[pos_col]]
    cluster_id <- 1
    cluster_starts <- 1

    if (length(positions) > 1) {
      cluster_assignments <- rep(1, length(positions))

      for (i in 2:length(positions)) {
        if (positions[i] - positions[cluster_starts] > distance_threshold) {
          cluster_id <- cluster_id + 1
          cluster_starts <- i
        }
        cluster_assignments[i] <- cluster_id
      }
    } else {
      cluster_assignments <- 1
    }

    # Assign clusters and find top SNPs
    chr_dt[, temp_cluster := cluster_assignments]

    # Find top SNP in each cluster (minimum p-value)
    top_snps <- chr_dt[, .SD[which.min(get(pvalue_col))], by = temp_cluster]

    # Create cluster IDs
    chr_dt[, cluster := paste0("chr", chr, "_cluster", temp_cluster)]

    # Mark top SNPs
    chr_dt[get(rsid_col) %in% top_snps[[rsid_col]], top := "top"]

    # Update main data.table
    dt[get(rsid_col) %in% chr_dt[[rsid_col]],
       c("cluster", "top") := .(chr_dt$cluster[match(get(rsid_col), chr_dt[[rsid_col]])],
                                chr_dt$top[match(get(rsid_col), chr_dt[[rsid_col]])])]
  }

  return(dt)
}
