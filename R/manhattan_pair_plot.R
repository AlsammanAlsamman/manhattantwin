#' Helper function to generate a formatted plot title
#'
#' @param plot_title_prefix The prefix for the plot title (e.g., study name).
#' @param n_cases The number of cases in the study.
#' @param n_controls The number of controls in the study.
#' @param lambda The genomic inflation factor (lambda GC).
#' @param total_snps_in_study The total number of SNPs in the study.
#' @param add_date_to_title A logical flag to append the current date to the title.
#' @param gwas_data The GWAS dataset.
#' @param font_family The font family to use for the title.
#'
#' @return A formatted string for the plot title.
construct_plot_title <- function(plot_title_prefix, n_cases, n_controls, lambda,
                                 total_snps_in_study, add_date_to_title, gwas_data,
                                 font_family = "Arial Unicode MS") {
  
  if (is.null(total_snps_in_study)) {
    snps_for_title <- format(nrow(gwas_data), big.mark = ",", scientific = FALSE)
    snps_title_label <- "Total SNPs"
  } else {
    snps_for_title <- format(total_snps_in_study, big.mark = ",", scientific = FALSE)
    snps_title_label <- "Total SNPs in Study"
  }
  
  title_string <- sprintf("%s\n(%s: %s | \U03BB = %s | Cases: %s, Controls: %s)",
                          plot_title_prefix, snps_title_label, snps_for_title,
                          lambda, format(n_cases, big.mark=","), format(n_controls, big.mark=","))
  
  if (add_date_to_title) {
    title_string <- paste0(title_string, " | ", Sys.Date())
  }
  return(title_string)
}

#' Helper function to preprocess GWAS data for Manhattan plotting
#'
#' @param gwas_data The GWAS dataset (a data.table or data.frame).
#' @param chr_col The name of the chromosome column.
#' @param bp_col The name of the base pair position column.
#' @param p_col The name of the p-value column.
#' @param gene_col The name of the gene column.
#' @param group_col The name of the grouping column for labeling.
#' @param plot_pval_threshold The p-value threshold for filtering SNPs to plot.
#' @param all_chrs A vector of all chromosome names to ensure consistent ordering.
#'
#' @return A list containing the processed plot data (`plot_data`) and axis annotation data (`axis_data`).
prepare_gwas_data_for_plotting <- function(gwas_data, chr_col, bp_col, p_col, gene_col,
                                           group_col, plot_pval_threshold, all_chrs) {
  
  plot_data <- data.table::as.data.table(gwas_data)
  plot_data <- plot_data[!is.na(get(p_col)) & get(p_col) > 0 & get(p_col) <= plot_pval_threshold]
  
  if (nrow(plot_data) == 0) {
    return(list(plot_data = data.table(), axis_data = data.table()))
  }
  
  plot_data[, `:=` (pvalue = as.numeric(get(p_col)),
                    CHR = as.character(get(chr_col)),
                    BP = as.numeric(get(bp_col)))]
  plot_data[, log_p := -log10(pvalue)]
  
  plot_data[, CHR_num := as.numeric(factor(CHR, levels = unique(stringr::str_sort(all_chrs, numeric = TRUE))))]
  plot_data <- plot_data[order(CHR_num, BP)]
  
  cumulative_summary <- plot_data[, .(max_bp = max(BP)), by = CHR_num][order(CHR_num)]
  cumulative_summary$bp_add <- c(0, cumsum(as.numeric(cumulative_summary$max_bp[-nrow(cumulative_summary)])))
  
  plot_data <- merge(plot_data, cumulative_summary[, .(CHR_num, bp_add)], by = "CHR_num")
  plot_data[, bp_cum := BP + bp_add]
  
  axis_data <- plot_data[, .(center = (max(bp_cum) + min(bp_cum)) / 2), by = CHR]
  
  return(list(plot_data = plot_data, axis_data = axis_data))
}

#' Helper function to generate a single Manhattan plot (normal or inverted)
#'
#' (Internal function for create_mirrored_manhattan_plot)
#'
generate_manhattan_subplot <- function(plot_data, axis_data, plot_title, sig_lines,
                                       label_threshold_colors, y_axis_squish_threshold,
                                       gene_label_size, gene_col, group_col,
                                       loci_to_label, genes_to_label, custom_gene_colors,
                                       inverted, font_family = "Arial Unicode MS",
                                       label_orientation = "horizontal",
                                       label_alpha = 0.7) {
  
  snps_to_label <- data.table::data.table()
  if (!is.null(gene_col) && gene_col %in% names(plot_data) && nrow(plot_data) > 0) {
    if (is.null(group_col)) group_col <- gene_col
    if (group_col %in% names(plot_data)) {
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
          data.table::setnames(custom_color_dt, "gene_name", gene_col)
          snps_to_label[custom_color_dt, on = gene_col, label_color := i.new_color]
        }
        snps_to_label[, label_text := get(gene_col)]
      }
    }
  }
  
  # --- **FIXED**: Define squish transformation WITH inverse function ---
  compression_factor <- 0.1
  squish_forward <- function(x) ifelse(x <= y_axis_squish_threshold, x, y_axis_squish_threshold + (x - y_axis_squish_threshold) * compression_factor)
  squish_inverse <- function(x) ifelse(x <= y_axis_squish_threshold, x, (x - y_axis_squish_threshold) / compression_factor + y_axis_squish_threshold)
  squish_trans <- scales::trans_new("squish", squish_forward, squish_inverse) # <-- THE FIX IS HERE
  
  x_axis_expansion <- ggplot2::expansion(mult = c(0.015, 0.015))
  custom_chromosome_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78")
  
  max_log_p <- if(nrow(plot_data) > 0) max(plot_data$log_p, na.rm = TRUE) else 0
  linear_breaks <- seq(0, min(y_axis_squish_threshold, max_log_p), by = 10)
  y_breaks <- if (max_log_p > y_axis_squish_threshold) {
    squished_breaks <- pretty(c(y_axis_squish_threshold, max_log_p), n = 3)
    unique(sort(c(linear_breaks, squished_breaks)))
  } else {
    linear_breaks
  }
  
  if (inverted) {
    inverted_data <- data.table::copy(plot_data)
    inverted_data[, transformed_log_p := -squish_forward(log_p)]
    transformed_breaks <- -squish_forward(y_breaks)
    
    p <- ggplot2::ggplot(inverted_data) +
      ggplot2::geom_point(ggplot2::aes(x = bp_cum, y = transformed_log_p, color = factor(CHR_num)), alpha = 0.8, size = 1.5) +
      ggplot2::scale_color_manual(values = rep(custom_chromosome_colors, 2), guide = "none") +
      ggplot2::geom_hline(yintercept = -squish_forward(-log10(sig_lines)), color = names(sig_lines), linetype = "dashed", linewidth = 0.8) +
      ggplot2::scale_x_continuous(label = axis_data$CHR, breaks = axis_data$center, expand = x_axis_expansion) +
      ggplot2::scale_y_continuous(breaks = transformed_breaks, labels = as.character(y_breaks), position = "left") +
      ggplot2::labs(x = "Chromosome", y = expression(-log[10](italic(P))), caption = plot_title) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        text = ggplot2::element_text(family = font_family), panel.grid.major.x = ggplot2::element_blank(), panel.grid.minor.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(color = "black"), axis.title.x = ggplot2::element_text(color = "black", face = "bold"),
        plot.caption = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold", color = "black", margin = ggplot2::margin(t = 10)),
        axis.text.y = ggplot2::element_text(angle = 90, hjust = 0.5, color = "black"),
        axis.title.y = ggplot2::element_text(angle = 90, vjust = 0.5, color = "black", face = "bold"),
        plot.margin = ggplot2::margin(t = -15, r = 5.5, b = 5.5, l = 5.5, unit = "pt")
      )
  } else {
    p <- ggplot2::ggplot(plot_data) +
      ggplot2::geom_point(ggplot2::aes(x = bp_cum, y = log_p, color = factor(CHR_num)), alpha = 0.8, size = 1.5) +
      ggplot2::scale_color_manual(values = rep(custom_chromosome_colors, 2), guide = "none") +
      ggplot2::geom_hline(yintercept = -log10(sig_lines), color = names(sig_lines), linetype = "dashed", linewidth = 0.8) +
      ggplot2::scale_x_continuous(breaks = axis_data$center, expand = x_axis_expansion) +
      ggplot2::scale_y_continuous(trans = squish_trans, breaks = y_breaks, expand = ggplot2::expansion(mult = c(0, 0.1))) +
      ggplot2::labs(x = NULL, y = expression(-log[10](italic(P))), title = plot_title) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        text = ggplot2::element_text(family = font_family), plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold", color = "black"),
        panel.grid.major.x = ggplot2::element_blank(), panel.grid.minor.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(color = "black"), axis.title.y = ggplot2::element_text(color = "black", face = "bold"),
        axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(t = 5.5, r = 5.5, b = -15, l = 5.5, unit = "pt")
      )
  }
  
  if (nrow(snps_to_label) > 0) {
    label_data <- if(inverted) data.table::copy(snps_to_label)[, transformed_log_p := -squish_forward(log_p)] else snps_to_label
    y_aes <- if(inverted) "transformed_log_p" else "log_p"
    nudge_val <- if(inverted) -15 else 15
    repel_layer <- NULL
    
    if (label_orientation == "vertical") {
      repel_layer <- ggrepel::geom_label_repel(
        data = label_data, ggplot2::aes_string(x = "bp_cum", y = y_aes, label = "label_text", color = "label_color"),
        size = gene_label_size, fontface = "bold.italic", nudge_y = nudge_val, angle = 90, hjust = 0.5,
        direction = "y", segment.color = 'grey50', segment.linetype = "dashed", min.segment.length = 0,
        max.overlaps = Inf, box.padding = 0.8, label.size = 0.25, fill = ggplot2::alpha("white", label_alpha)
      )
    } else {
      repel_layer <- ggrepel::geom_label_repel(
        data = label_data, ggplot2::aes_string(x = "bp_cum", y = y_aes, label = "label_text", color = "label_color"),
        size = gene_label_size, fontface = "bold.italic", nudge_y = nudge_val, angle = 0,
        direction = "x", segment.color = 'grey50', segment.linetype = "dashed", min.segment.length = 0,
        max.overlaps = Inf, box.padding = 0.5, label.size = 0.25, fill = ggplot2::alpha("white", label_alpha)
      )
    }
    p <- p + ggnewscale::new_scale_color() + repel_layer + ggplot2::scale_color_identity(guide = "none")
  }
  return(p)
}


#' Create a Pair of Inverted Manhattan Plots with Advanced Features
#'
#' (Version 18: Synced features with single plot function)
#'
#' @param gwas_data1 The first GWAS dataset.
#' @param gwas_data2 The second GWAS dataset.
#' @param chr_col Column name for chromosome.
#' @param bp_col Column name for base pair position.
#' @param p_col Column name for p-value.
#' @param rsid_col Optional: Column name for SNP identifiers (e.g., RSIDs) for filtering.
#' @param rsids_to_remove Optional: A vector of SNP IDs to remove from both datasets before analysis.
#' @param gene_col Column name for gene labels.
#' @param group_col Column for grouping labels. Defaults to gene_col.
#' @param ... other parameters
#' @param n_cases1 Number of cases in dataset 1.
#' @param n_controls1 Number of controls in dataset 1.
#' @param n_cases2 Number of cases in dataset 2.
#' @param n_controls2 Number of controls in dataset 2.
#' @param plot_title1 Title for the top plot.
#' @param plot_title2 Title for the bottom plot.
#' @param lambda1 Optional lambda for dataset 1.
#' @param lambda2 Optional lambda for dataset 2.
#' @param label_orientation Orientation of gene labels ("vertical" or "horizontal").
#' @param label_alpha Transparency of gene label backgrounds (0-1).
#' @param output_width Width of the output plot in inches.
#' @param output_height Height of the output plot in inches.
#' @export
create_mirrored_manhattan_plot <- function(gwas_data1, gwas_data2,
                                           chr_col = "chr",
                                           bp_col = "pos",
                                           p_col = "pvalue",
                                           rsid_col = NULL,
                                           rsids_to_remove = NULL,
                                           gene_col = "gene",
                                           group_col = NULL,
                                           loci_to_label = NULL,
                                           genes_to_label = NULL,
                                           custom_gene_colors = NULL,
                                           n_cases1, n_controls1,
                                           n_cases2, n_controls2,
                                           plot_title1 = "Dataset 1",
                                           plot_title2 = "Dataset 2",
                                           lambda1 = NULL, lambda2 = NULL,
                                           total_snps_in_study1 = NULL, total_snps_in_study2 = NULL,
                                           add_date_to_title = FALSE,
                                           plot_pval_threshold = 1,
                                           sig_lines = c("red" = 5e-8, "blue" = 1e-5),
                                           label_threshold_colors = c("red" = 5e-8, "orange" = 1e-6, "darkblue" = 1e-5),
                                           y_axis_squish_threshold = 30,
                                           gene_label_size = 3.5,
                                           label_orientation = "horizontal",
                                           label_alpha = 0.7,
                                           output_folder = "Inverted_Manhattan_Plots",
                                           file_name_prefix = "inverted_manhattan",
                                           font_family = "Arial Unicode MS",
                                           output_width = 18,
                                           output_height = 12) {
  
  data1 <- data.table::as.data.table(gwas_data1)
  data2 <- data.table::as.data.table(gwas_data2)
  
  if (!is.null(rsid_col) && !is.null(rsids_to_remove)) {
    if (rsid_col %in% names(data1) && rsid_col %in% names(data2)) {
      initial_count1 <- nrow(data1)
      data1 <- data1[!(get(rsid_col) %in% rsids_to_remove)]
      message(sprintf("Dataset 1: Removed %d SNPs based on the provided rsid list.", initial_count1 - nrow(data1)))
      
      initial_count2 <- nrow(data2)
      data2 <- data2[!(get(rsid_col) %in% rsids_to_remove)]
      message(sprintf("Dataset 2: Removed %d SNPs based on the provided rsid list.", initial_count2 - nrow(data2)))
    } else {
      warning(paste0("`rsid_col` '", rsid_col, "' not found in one or both datasets. Skipping removal."))
    }
  }
  
  all_chrs <- sort(unique(c(data1[[chr_col]], data2[[chr_col]])))
  
  if (is.null(lambda1)) {
    p1 <- data1[[p_col]][!is.na(data1[[p_col]]) & data1[[p_col]] > 0 & data1[[p_col]] <= 1]
    chisq1 <- stats::qchisq(1 - p1, 1)
    lambda1 <- round(median(chisq1, na.rm = TRUE) / stats::qchisq(0.5, 1), 4)
    message(paste("Calculated Lambda GC for Dataset 1:", lambda1))
  }
  if (is.null(lambda2)) {
    p2 <- data2[[p_col]][!is.na(data2[[p_col]]) & data2[[p_col]] > 0 & data2[[p_col]] <= 1]
    chisq2 <- stats::qchisq(1 - p2, 1)
    lambda2 <- round(median(chisq2, na.rm = TRUE) / stats::qchisq(0.5, 1), 4)
    message(paste("Calculated Lambda GC for Dataset 2:", lambda2))
  }
  
  processed_data1 <- prepare_gwas_data_for_plotting(
    gwas_data = data1, chr_col = chr_col, bp_col = bp_col, p_col = p_col, gene_col = gene_col,
    group_col = group_col, plot_pval_threshold = plot_pval_threshold, all_chrs = all_chrs
  )
  processed_data2 <- prepare_gwas_data_for_plotting(
    gwas_data = data2, chr_col = chr_col, bp_col = bp_col, p_col = p_col, gene_col = gene_col,
    group_col = group_col, plot_pval_threshold = plot_pval_threshold, all_chrs = all_chrs
  )
  
  title1 <- construct_plot_title(
    plot_title_prefix = plot_title1, n_cases = n_cases1, n_controls = n_controls1, lambda = lambda1,
    total_snps_in_study = total_snps_in_study1, add_date_to_title = add_date_to_title, gwas_data = data1, font_family = font_family
  )
  title2 <- construct_plot_title(
    plot_title_prefix = plot_title2, n_cases = n_cases2, n_controls = n_controls2, lambda = lambda2,
    total_snps_in_study = total_snps_in_study2, add_date_to_title = add_date_to_title, gwas_data = data2, font_family = font_family
  )
  
  top_plot <- generate_manhattan_subplot(
    plot_data = processed_data1$plot_data, axis_data = processed_data1$axis_data, plot_title = title1,
    sig_lines = sig_lines, label_threshold_colors = label_threshold_colors, y_axis_squish_threshold = y_axis_squish_threshold,
    gene_label_size = gene_label_size, gene_col = gene_col, group_col = group_col, loci_to_label = loci_to_label,
    genes_to_label = genes_to_label, custom_gene_colors = custom_gene_colors, inverted = FALSE, font_family = font_family,
    label_orientation = label_orientation, label_alpha = label_alpha
  )
  bottom_plot <- generate_manhattan_subplot(
    plot_data = processed_data2$plot_data, axis_data = processed_data2$axis_data, plot_title = title2,
    sig_lines = sig_lines, label_threshold_colors = label_threshold_colors, y_axis_squish_threshold = y_axis_squish_threshold,
    gene_label_size = gene_label_size, gene_col = gene_col, group_col = group_col, loci_to_label = loci_to_label,
    genes_to_label = genes_to_label, custom_gene_colors = custom_gene_colors, inverted = TRUE, font_family = font_family,
    label_orientation = label_orientation, label_alpha = label_alpha
  )
  
  combined_plot <- top_plot / bottom_plot + patchwork::plot_layout(heights = c(1, 1))
  
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  file_path_png <- file.path(output_folder, paste0(file_name_prefix, ".png"))
  file_path_pdf <- file.path(output_folder, paste0(file_name_prefix, ".pdf"))
  
  tryCatch({
    grDevices::png(filename = file_path_png, width = output_width, height = output_height, units = "in", res = 300, type = "cairo", bg = "white")
    print(combined_plot)
    grDevices::dev.off()
    grDevices::cairo_pdf(filename = file_path_pdf, width = output_width, height = output_height)
    print(combined_plot)
    grDevices::dev.off()
  }, error = function(e) {
    message("Cairo devices failed, falling back to standard ggsave devices.")
    ggplot2::ggsave(file_path_png, plot = combined_plot, device = ragg::agg_png, width = output_width, height = output_height, units = "in", res = 300, bg = "white")
    ggplot2::ggsave(file_path_pdf, plot = combined_plot, device = "pdf", width = output_width, height = output_height, units = "in")
  })
  
  message(paste("Successfully saved plots to:", file_path_png, "and", file_path_pdf))
  return(invisible(combined_plot))
}

