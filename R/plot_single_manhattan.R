#' Publication-Quality Manhattan Plot with Advanced Features
#'
#' (Version 17: Implemented dynamic Y-axis breaks for better visualization)
#'
#' @import data.table
#' @import ggplot2
#' @import ggrepel
#' @import scales
#' @import stringr
#' @import ragg
#' @import ggnewscale
#' @import extrafont
#'
#' @param gwas_data A data.frame or data.table containing GWAS summary statistics.
#' @param chr_col The name of the column for chromosome information. Default is "chr".
#' @param bp_col The name of the column for base pair position. Default is "pos".
#' @param p_col The name of the column for p-values.
#' @param gene_col The name of the column for gene/SNP labels. Set to NULL to disable labeling.
#' @param group_col The column used to group SNPs for labeling. Defaults to gene_col.
#' @param loci_to_label Optional vector of locus/group names to label.
#' @param genes_to_label Optional vector of gene names to label.
#' @param custom_gene_colors Optional named vector of custom colors for specific gene labels.
#' @param n_cases Number of cases in the GWAS.
#' @param n_controls Number of controls in the GWAS.
#' @param plot_title_prefix Prefix for the plot title.
#' @param lambda Optional genomic inflation factor. If NULL, it is calculated from the full dataset.
#' @param total_snps_in_study Optional total number of SNPs in the study for title display.
#' @param add_date_to_title Logical; if TRUE, appends the current date to the plot title.
#' @param plot_pval_threshold Maximum p-value to include in the plot (default 1).
#' @param sig_lines Named vector of significance lines (e.g., c("red" = 5e-8)).
#' @param label_threshold_colors Named vector of p-value thresholds for label colors.
#' @param y_axis_squish_threshold Y-axis value above which to compress the scale.
#' @param gene_label_size Size of gene/SNP label text.
#' @param label_orientation The orientation of gene labels. Can be "vertical" or "horizontal" (default).
#' @param label_alpha The transparency of the gene label background (0=transparent, 1=opaque). Default is 0.7.
#' @param output_folder Folder to save output plots.
#' @param file_name_prefix Prefix for output file names (PNG and PDF).
#' @param font_family Font family to use for the plot.
#' @param output_width The width of the saved plot in inches. Default is 18.
#' @param output_height The height of the saved plot in inches. Default is 8.
#' @export
#' @return Invisibly returns NULL. Saves a PNG and a PDF file of the plot.
#'
plot_single_manhattan <- function(gwas_data,
                                  chr_col = "chr",
                                  bp_col = "pos",
                                  p_col= "pvalue",
                                  gene_col = NULL,
                                  group_col = NULL,
                                  loci_to_label = NULL,
                                  genes_to_label = NULL,
                                  custom_gene_colors = NULL,
                                  n_cases,
                                  n_controls,
                                  plot_title_prefix,
                                  lambda = NULL,
                                  total_snps_in_study = NULL,
                                  add_date_to_title = FALSE,
                                  plot_pval_threshold = 1,
                                  sig_lines = c("red" = 5e-8, "blue" = 1e-5),
                                  label_threshold_colors = c("red" = 5e-8, "orange" = 1e-6, "darkblue" = 1e-5),
                                  y_axis_squish_threshold = 30,
                                  gene_label_size = 3.5,
                                  label_orientation = "horizontal",
                                  label_alpha = 0.7,
                                  output_folder = "Single_Manhattan_Plots",
                                  file_name_prefix,
                                  font_family = "Arial Unicode MS",
                                  output_width = 18,
                                  output_height = 8) {
  
  # --- 1. Input Validation ---
  message(paste("--- Generating Manhattan plot for:", plot_title_prefix, "---"))
  required_cols <- c(chr_col, bp_col, p_col)
  labeling_enabled <- !is.null(gene_col)
  
  if (labeling_enabled) {
    if (is.null(group_col)) group_col <- gene_col
    required_cols <- c(required_cols, gene_col, group_col)
  }
  
  if (!all(required_cols %in% names(gwas_data))) {
    missing <- required_cols[!required_cols %in% names(gwas_data)]
    stop(paste("Input data is missing required columns:", paste(missing, collapse=", ")))
  }
  
  # --- 2. Calculate Lambda from FULL dataset ---
  if (is.null(lambda)) {
    message("Lambda not provided, calculating from full dataset...")
    p_values_for_lambda <- gwas_data[[p_col]]
    p_values_for_lambda <- p_values_for_lambda[!is.na(p_values_for_lambda) & p_values_for_lambda > 0 & p_values_for_lambda <= 1]
    chisq <- stats::qchisq(1 - p_values_for_lambda, 1)
    lambda_val <- round(median(chisq, na.rm = TRUE) / stats::qchisq(0.5, 1), 4)
    message(paste("Calculated Lambda GC:", lambda_val))
  } else {
    lambda_val <- lambda
    message(paste("Using provided Lambda GC:", lambda_val))
  }
  
  # --- 3. Filter Data for Plotting and Prepare Columns ---
  plot_data <- data.table::as.data.table(gwas_data)
  plot_data <- plot_data[!is.na(get(p_col)) & get(p_col) > 0 & get(p_col) <= plot_pval_threshold]
  
  if (nrow(plot_data) == 0) {
    message("No SNPs remained after p-value filtering. Skipping plot generation.")
    return(invisible(NULL))
  }
  
  message(sprintf("Filtered to %d SNPs with p <= %.2e for plotting.", nrow(plot_data), plot_pval_threshold))
  plot_data[, `:=` (pvalue = as.numeric(get(p_col)), CHR = as.character(get(chr_col)), BP = as.numeric(get(bp_col)))]
  plot_data[, log_p := -log10(pvalue)]
  
  # --- 4. Prepare Cumulative Plotting Data ---
  plot_data[, CHR_num := as.numeric(factor(CHR, levels = unique(stringr::str_sort(unique(CHR), numeric = TRUE))))]
  plot_data <- plot_data[order(CHR_num, BP)]
  
  cumulative_summary <- plot_data[, .(max_bp = max(BP)), by = CHR_num][order(CHR_num)]
  cumulative_summary$bp_add <- c(0, cumsum(as.numeric(cumulative_summary$max_bp[-nrow(cumulative_summary)])))
  
  plot_data <- merge(plot_data, cumulative_summary[, .(CHR_num, bp_add)], by = "CHR_num")
  plot_data[, bp_cum := BP + bp_add]
  axis_data <- plot_data[, .(center = (max(bp_cum) + min(bp_cum)) / 2), by = CHR]
  
  # --- 5. Prepare Gene Labels ---
  snps_to_label <- data.table::data.table()
  if (labeling_enabled) {
    message("Gene labeling is enabled.")
    snps_for_labeling <- plot_data[!is.na(get(group_col))]
    
    if (!is.null(loci_to_label)) snps_for_labeling <- snps_for_labeling[get(group_col) %in% loci_to_label]
    else if (!is.null(genes_to_label)) snps_for_labeling <- snps_for_labeling[get(gene_col) %in% genes_to_label]
    else snps_for_labeling <- snps_for_labeling[pvalue < max(label_threshold_colors)]
    
    if (nrow(snps_for_labeling) > 0) {
      snps_to_label <- snps_for_labeling[, .SD[which.min(pvalue)], by = group_col]
      sorted_thresholds <- sort(label_threshold_colors, decreasing = TRUE)
      snps_to_label[, label_color := names(sorted_thresholds)[1]]
      if (length(sorted_thresholds) > 1) {
        for(i in 2:length(sorted_thresholds)) snps_to_label[pvalue <= sorted_thresholds[i], label_color := names(sorted_thresholds)[i]]
      }
      if (!is.null(custom_gene_colors)) {
        custom_color_dt <- data.table::data.table(gene_name = names(custom_gene_colors), new_color = as.character(custom_gene_colors))
        join_on_condition <- setNames("gene_name", gene_col)
        snps_to_label[custom_color_dt, on = join_on_condition, label_color := i.new_color]
      }
    }
  } else {
    message("Gene labeling is disabled (gene_col is NULL).")
  }
  
  # --- 6. Custom Y-Axis Transformation and Dynamic Breaks ---
  compression_factor <- 0.1
  squish_forward <- function(x) ifelse(x <= y_axis_squish_threshold, x, y_axis_squish_threshold + (x - y_axis_squish_threshold) * compression_factor)
  squish_inverse <- function(x) ifelse(x <= y_axis_squish_threshold, x, (x - y_axis_squish_threshold) / compression_factor + y_axis_squish_threshold)
  squish_trans <- scales::trans_new("squish", squish_forward, squish_inverse)
  
  # **FIXED**: Dynamically generate Y-axis breaks based on data range
  max_log_p <- max(plot_data$log_p, na.rm = TRUE)
  linear_breaks <- seq(0, min(y_axis_squish_threshold, max_log_p), by = 10)
  
  y_breaks <- if (max_log_p > y_axis_squish_threshold) {
    squished_breaks <- pretty(c(y_axis_squish_threshold, max_log_p), n = 3)
    unique(sort(c(linear_breaks, squished_breaks)))
  } else {
    linear_breaks
  }
  
  # --- 7. Generate Plot Title ---
  snps_title_label <- "Total SNPs"
  snps_for_title <- if(is.null(total_snps_in_study)) format(nrow(gwas_data), big.mark=",") else {
    snps_title_label <- "Total SNPs in Study"
    format(total_snps_in_study, big.mark=",")
  }
  title_string <- sprintf("%s\n(%s: %s | \U03BB = %s | Cases: %s, Controls: %s)",
                          plot_title_prefix, snps_title_label, snps_for_title, lambda_val,
                          format(n_cases, big.mark=","), format(n_controls, big.mark=","))
  if (add_date_to_title) title_string <- paste0(title_string, " | ", Sys.Date())
  
  # --- 8. Generate the Plot ---
  rainbow_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#F781BF", "#A65628", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3")
  manhattan_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = bp_cum, y = log_p)) +
    ggplot2::geom_point(ggplot2::aes(color = as.factor(CHR_num)), alpha = 0.8, size = 1.5) +
    ggplot2::scale_color_manual(values = rep(rainbow_palette, 2), guide = "none") +
    ggplot2::geom_hline(yintercept = -log10(sig_lines), color = names(sig_lines), linetype = "dashed", linewidth = 0.8)
  
  if (labeling_enabled && nrow(snps_to_label) > 0) {
    repel_layer <- NULL
    if (label_orientation == "vertical") {
      repel_layer <- ggrepel::geom_label_repel(
        data = snps_to_label, ggplot2::aes(label = get(gene_col), color = label_color),
        size = gene_label_size, fontface = "bold.italic", nudge_y = 15, angle = 90, hjust = 0,
        direction = "y", segment.color = 'grey50', segment.linetype = "dashed",
        min.segment.length = 0, max.overlaps = Inf, box.padding = 0.8, label.size = 0.25,
        fill = ggplot2::alpha("white", label_alpha)
      )
    } else {
      repel_layer <- ggrepel::geom_label_repel(
        data = snps_to_label, ggplot2::aes(label = get(gene_col), color = label_color),
        size = gene_label_size, fontface = "bold.italic", nudge_y = 15, angle = 0,
        direction = "x", segment.color = 'grey50', segment.linetype = "dashed",
        min.segment.length = 0, max.overlaps = Inf, box.padding = 0.5, label.size = 0.25,
        fill = ggplot2::alpha("white", label_alpha)
      )
    }
    manhattan_plot <- manhattan_plot + ggnewscale::new_scale_color() + repel_layer + ggplot2::scale_color_identity(guide = "none")
  }
  
  manhattan_plot <- manhattan_plot +
    ggplot2::scale_x_continuous(label = axis_data$CHR, breaks = axis_data$center) +
    # Use the new dynamic breaks vector
    ggplot2::scale_y_continuous(trans = squish_trans, expand = ggplot2::expansion(mult = c(0, 0.1)),
                                breaks = y_breaks) +
    ggplot2::labs(x = "Chromosome", y = expression(-log[10](italic(P))), title = title_string) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      text = ggplot2::element_text(family = font_family),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold", color = "black"),
      panel.grid.major.x = ggplot2::element_blank(), panel.grid.minor.x = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(color = "black"),
      axis.title = ggplot2::element_text(color = "black", face = "bold")
    )
  
  # --- 9. Save the Plot ---
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  file_path_png <- file.path(output_folder, paste0(file_name_prefix, "_pub.png"))
  file_path_pdf <- file.path(output_folder, paste0(file_name_prefix, "_pub.pdf"))
  
  tryCatch({
    grDevices::png(filename = file_path_png, width = output_width, height = output_height, units = "in", res = 300, type = "cairo", bg = "white")
    print(manhattan_plot)
    grDevices::dev.off()
    grDevices::cairo_pdf(filename = file_path_pdf, width = output_width, height = output_height)
    print(manhattan_plot)
    grDevices::dev.off()
  }, error = function(e) {
    message("Cairo devices failed, falling back to standard devices.")
    ggplot2::ggsave(file_path_png, plot = manhattan_plot, device = ragg::agg_png, width = output_width, height = output_height, units = "in", res = 300, bg = "white")
    ggplot2::ggsave(file_path_pdf, plot = manhattan_plot, device = "pdf", width = output_width, height = output_height, units = "in")
  })
  
  message(paste("Successfully saved plots to:", file_path_png, "and", file_path_pdf))
  return(invisible(NULL))
}
