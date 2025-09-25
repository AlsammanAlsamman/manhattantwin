#' Annotate GWAS Data with Genomic Loci and Identify Top SNPs
#'
#' This function annotates GWAS data with genomic loci information and identifies the top SNP for each locus.
#' It can use either pre-defined top SNP positions from the loci data or dynamically determine the top SNP
#' by finding the variant with the smallest p-value within each locus. All column names are user-specifiable.
#'
#' @param gwas_data A data.table or data.frame containing GWAS results.
#' @param loci_data A data.table or data.frame defining genomic loci.
#' @param loci_name A prefix for the new column names (e.g., "fuma8").
#' @param p_col Optional character string specifying the p-value column in `gwas_data` to use for
#'   dynamically identifying the top SNP if not pre-defined in `loci_data`.
#' @param gwas_chr_col Column name for chromosome in `gwas_data` (default: "CHR").
#' @param gwas_pos_col Column name for base pair position in `gwas_data` (default: "BP").
#' @param loci_chr_col Column name for chromosome in `loci_data` (default: "chr").
#' @param loci_start_col Column name for start position in `loci_data` (default: "start").
#' @param loci_end_col Column name for end position in `loci_data` (default: "end").
#' @param loci_name_col Column name for locus name in `loci_data` (default: "locusname").
#' @param loci_pos_col Column name for pre-defined top SNP position in `loci_data` (default: "pos").
#'
#' @return The modified gwas_data data.table with two new columns:
#'   \itemize{
#'     \item A locus annotation column named "loci.{loci_name}" indicating which locus each SNP belongs to
#'     \item A top SNP flag column named "{loci_name}_topsnp" marking the top SNP in each locus
#'   }
#'
#' @examples
#' # Create example GWAS data with custom column names
#' gwas_data <- data.table::data.table(
#'   chromosome = c("1", "1", "1", "2", "2", "3"),
#'   position = c(110000, 120000, 510000, 760000, 770000, 220000),
#'   p_value = c(0.001, 0.0001, 0.01, 0.00001, 0.00005, 0.0002),
#'   SNP = c("rs111", "rs112", "rs113", "rs114", "rs115", "rs116")
#' )
#' 
#' # Create example loci data with custom column names
#' loci_data <- data.table::data.table(
#'   chrom = c("1", "1", "2", "3"),
#'   begin = c(100000, 500000, 750000, 200000),
#'   stop = c(150000, 550000, 800000, 250000),
#'   locus_id = c("LOC1", "LOC2", "LOC3", "LOC4"),
#'   top_position = c(120000, 520000, 770000, 220000)  # Pre-defined top SNP positions
#' )
#' 
#' # Annotate using custom column names
#' annotated_data <- annotate_gwas_with_loci(
#'   gwas_data = gwas_data,
#'   loci_data = loci_data,
#'   loci_name = "custom",
#'   p_col = "p_value",
#'   gwas_chr_col = "chromosome",
#'   gwas_pos_col = "position",
#'   loci_chr_col = "chrom",
#'   loci_start_col = "begin",
#'   loci_end_col = "stop",
#'   loci_name_col = "locus_id",
#'   loci_pos_col = "top_position"
#' )
#' print(annotated_data)
#'
#' @export
#' @importFrom data.table data.table setDT setnames copy
annotate_gwas_with_loci <- function(gwas_data, loci_data, loci_name, p_col = NULL,
                                   gwas_chr_col = "CHR", gwas_pos_col = "BP",
                                   loci_chr_col = "chr", loci_start_col = "start",
                                   loci_end_col = "end", loci_name_col = "locusname",
                                   loci_pos_col = "pos") {
  
  # --- 1. Input Validation ---
  if (!is.data.table(gwas_data)) data.table::setDT(gwas_data)
  if (!is.data.table(loci_data)) data.table::setDT(loci_data)
  
  if (!is.character(loci_name) || length(loci_name) != 1 || loci_name == "") {
    stop("'loci_name' must be a single, non-empty string.")
  }
  
  # Validate column names in gwas_data
  required_gwas_cols <- c(gwas_chr_col, gwas_pos_col)
  missing_gwas_cols <- setdiff(required_gwas_cols, names(gwas_data))
  if (length(missing_gwas_cols) > 0) {
    stop(sprintf("Required columns not found in gwas_data: %s", paste(missing_gwas_cols, collapse = ", ")))
  }
  
  # Validate column names in loci_data
  required_loci_cols <- c(loci_chr_col, loci_start_col, loci_end_col, loci_name_col)
  missing_loci_cols <- setdiff(required_loci_cols, names(loci_data))
  if (length(missing_loci_cols) > 0) {
    stop(sprintf("Required columns not found in loci_data: %s", paste(missing_loci_cols, collapse = ", ")))
  }
  
  if (!is.null(p_col)) {
    if (!is.character(p_col) || length(p_col) != 1 || p_col == "") {
      stop("'p_col' must be a single, non-empty string.")
    }
    if (!p_col %in% names(gwas_data)) {
      stop(sprintf("The specified p-value column '%s' was not found in gwas_data.", p_col))
    }
  }
  
  # --- 2. Standardize Column Types ---
  message("Standardizing input tables for locus annotation...")
  
  # Harmonize chromosome types
  gwas_data[, (gwas_chr_col) := as.character(get(gwas_chr_col))]
  loci_data[, (loci_chr_col) := as.character(get(loci_chr_col))]
  
  # --- 3. Annotate Locus Regions ---
  new_col_name <- paste0("loci.", loci_name)
  message(sprintf("Annotating with loci. New column will be named '%s'.", new_col_name))
  
  # Create a copy of loci_data with standardized column names
  loci_copy <- data.table::copy(loci_data)
  data.table::setnames(loci_copy, 
                      old = c(loci_chr_col, loci_start_col, loci_end_col, loci_name_col),
                      new = c("TEMP_CHR", "TEMP_START", "TEMP_END", "TEMP_LOCUS_NAME"))
  
  # Define join conditions
  join_on_conditions <- c(
    paste0("TEMP_CHR==", gwas_chr_col),
    paste0("TEMP_START<=", gwas_pos_col),
    paste0("TEMP_END>=", gwas_pos_col)
  )
  
  # Perform the join to annotate loci
  gwas_data[, (new_col_name) := loci_copy[.SD, on = join_on_conditions, x.TEMP_LOCUS_NAME, mult = "first"]]
  
  annotated_rows <- sum(!is.na(gwas_data[[new_col_name]]))
  message(sprintf("Locus region annotation complete. %d variants were assigned to a locus.", annotated_rows))
  
  # --- 4. Identify and Flag Top SNP ---
  topsnp_col_name <- paste0(loci_name, "_topsnp")
  gwas_data[, (topsnp_col_name) := NA_character_]  # Initialize column
  
  # Method 1: Use pre-defined top SNP positions if they exist in loci_data
  if (loci_pos_col %in% names(loci_data)) {
    message(sprintf("Found '%s' column in loci data. Flagging pre-defined top SNPs in new column '%s'.", 
                   loci_pos_col, topsnp_col_name))
    
    # Create a lookup table of CHR-BP for top SNPs
    topsnp_lookup <- unique(loci_data[, .(chr = get(loci_chr_col), pos = get(loci_pos_col))])
    
    # Perform a join to flag the matching SNPs in gwas_data
    join_condition <- setNames(c("chr", "pos"), c(gwas_chr_col, gwas_pos_col))
    gwas_data[topsnp_lookup, on = join_condition, (topsnp_col_name) := "top"]
    
    # Method 2: Dynamically find top SNP using provided p-value column
  } else if (!is.null(p_col)) {
    message(sprintf("Dynamically identifying top SNP for each locus using p-value column '%s'.", p_col))
    
    # Subset to only the annotated SNPs that have a p-value
    annotated_snps <- gwas_data[!is.na(get(new_col_name)) & !is.na(get(p_col))]
    
    if (nrow(annotated_snps) > 0) {
      # For each locus, find the row (SNP) with the minimum p-value
      top_snps_per_locus <- annotated_snps[, .SD[which.min(get(p_col))], by = new_col_name]
      top_snps_lookup <- top_snps_per_locus[, .(CHR = get(gwas_chr_col), BP = get(gwas_pos_col))]
      
      # Join this lookup table back to the main gwas_data to flag the top SNPs
      join_condition <- setNames(c("CHR", "BP"), c(gwas_chr_col, gwas_pos_col))
      gwas_data[top_snps_lookup, on = join_condition, (topsnp_col_name) := "top"]
    } else {
      message("  -> Note: No SNPs with both locus annotation and p-value found. Skipping dynamic top SNP search.")
    }
    
  } else {
    message(sprintf("No '%s' column in loci_data and no 'p_col' provided. Skipping top SNP annotation.", loci_pos_col))
  }
  
  topsnp_annotated_rows <- sum(gwas_data[[topsnp_col_name]] == "top", na.rm = TRUE)
  if (topsnp_annotated_rows > 0) {
    message(sprintf("Top SNP annotation complete. %d variants were flagged as 'top'.", topsnp_annotated_rows))
  }
  
  # --- 5. Return the modified data.table ---
  message("All annotations are complete.")
  return(gwas_data)
}