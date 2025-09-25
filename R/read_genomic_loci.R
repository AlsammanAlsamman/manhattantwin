#' Read Genomic Loci Definition File
#'
#' This function reads a file containing genomic loci definitions from different sources
#' (e.g., simple format or FUMA) and standardizes the output into a consistent format.
#' It efficiently handles large files using data.table for fast reading.
#'
#' @param file_path Character string. The path to the loci file.
#' @param source Character string. The format of the file. Options:
#'   \itemize{
#'     \item "default": A simple 4-column file (chr, start, end, locusname)
#'     \item "fuma": A 'GenomicRiskLoci.txt' file from FUMA. This will create
#'         a 'locusname' column (e.g., 'fuma_1') and add 'pos' and 'topsnp'
#'         columns from the 'pos' and 'rsID' columns in the file.
#'   }
#'   Default: "default"
#'
#' @return A data.table containing the loci information with columns:
#'   \itemize{
#'     \item chr: Chromosome number
#'     \item start: Start position of the locus
#'     \item end: End position of the locus
#'     \item locusname: Name/identifier of the locus
#'     \item pos: (For FUMA source) Position of the top SNP
#'     \item topsnp: (For FUMA source) rsID of the top SNP
#'   }
#'
#'
#' @examples
#' # Example for default format
#' default_loci_data <- data.table::data.table(
#'   chr = c("1", "1", "2", "3", "5"),
#'   start = c(100000, 500000, 750000, 200000, 300000),
#'   end = c(150000, 550000, 800000, 250000, 350000),
#'   locusname = c("LOC1", "LOC2", "LOC3", "LOC4", "LOC5")
#' )
#' 
#' default_file <- tempfile(fileext = ".txt")
#' data.table::fwrite(default_loci_data, file = default_file, sep = "\t", quote = FALSE)
#' default_result <- read_genomic_loci(default_file, source = "default")
#' print(default_result)
#' 
#' # Example for FUMA format
#' fuma_loci_data <- data.table::data.table(
#'   GenomicLocus = c(1, 2, 3, 4, 5),
#'   rsID = c("rs12345", "rs67890", "rs11111", "rs22222", "rs33333"),
#'   pos = c(120000, 520000, 770000, 220000, 320000),
#'   chr = c("1", "1", "2", "3", "5"),
#'   start = c(100000, 500000, 750000, 200000, 300000),
#'   end = c(150000, 550000, 800000, 250000, 350000)
#' )
#' 
#' fuma_file <- tempfile(fileext = ".txt")
#' data.table::fwrite(fuma_loci_data, file = fuma_file, sep = "\t", quote = FALSE)
#' fuma_result <- read_genomic_loci(fuma_file, source = "fuma")
#' print(fuma_result)
#'
#' @export
#' @importFrom data.table fread
read_genomic_loci <- function(file_path, source = "default") {
  
  # Validate file existence
  if (!file.exists(file_path)) {
    stop("File not found at the specified path: ", file_path)
  }
  
  # Validate source parameter
  if (!source %in% c("default", "fuma")) {
    stop("Invalid 'source' specified. Please use 'default' or 'fuma'.")
  }
  
  message("Reading loci file: ", file_path, " (format: ", source, ")")
  
  if (source == "default") {
    # Read simple 4-column format
    loci_data <- data.table::fread(
      file = file_path,
      col.names = c("chr", "start", "end", "locusname")
    )
    
  } else if (source == "fuma") {
    # Read FUMA format
    fuma_data <- data.table::fread(file = file_path)
    
    # Validate required columns for FUMA format
    required_cols <- c("GenomicLocus", "rsID", "pos", "chr", "start", "end")
    missing_cols <- setdiff(required_cols, names(fuma_data))
    
    if (length(missing_cols) > 0) {
      stop(paste("Input file specified as 'fuma' is missing required column(s):",
                 paste(missing_cols, collapse = ", ")))
    }
    
    # Process FUMA data
    loci_data <- fuma_data[, .(
      chr = chr,
      pos = pos,
      start = start,
      end = end,
      locusname = paste0("fuma_", GenomicLocus),
      topsnp = rsID
    )]
  }
  
  # Validate that we have the expected columns
  expected_cols <- c("chr", "start", "end", "locusname")
  if (source == "fuma") {
    expected_cols <- c(expected_cols, "pos", "topsnp")
  }
  
  missing_final_cols <- setdiff(expected_cols, names(loci_data))
  if (length(missing_final_cols) > 0) {
    stop("Internal error: Missing expected columns in output: ",
         paste(missing_final_cols, collapse = ", "))
  }
  
  message("Successfully read ", nrow(loci_data), " loci.")
  
  return(loci_data)
}