
library(manhattantwin)
data("gwasdataseta")


gwasdataseta <- manhattantwin::cluster_snps(gwasdataseta,chr_col = "chr", pos_col = "pos", pvalue_col = "pvalue", rsid_col = "rsid", pvalue_threshold = 5e-8, distance_threshold = 250000)

gwasdatasetb <- manhattantwin::cluster_snps(gwasdataseta,chr_col = "chr", pos_col = "pos", pvalue_col = "pvalue", rsid_col = "rsid", pvalue_threshold = 5e-5, distance_threshold = 250000)

manhattantwin::plot_single_manhattan(
  gwasdataseta,
  plot_title_prefix = "Example Dataset A",
  p_col = "pvalue",
  n_cases = 1000,
  n_controls = 1000,
  file_name_prefix = "example_a",group_col = "cluster",
  gene_col = "gene",
  output_folder="single_plots",
  y_axis_squish_threshold = 10
)




create_mirrored_manhattan_plot(
  gwasdataseta,
  gwasdatasetb,
  chr_col = "chr",
  bp_col = "pos",
  p_col = "pvalue",
  n_cases1 = 1000,
  n_controls1 = 1000,
  n_cases2 = 1000,
  n_controls2 = 1000,
  file_name_prefix = "example_a_vs_b",
  group_col = "cluster",
  gene_col = "gene",
  label_threshold_colors = c("red" = 5e-8, "orange" = 1e-6, "darkblue" = 1e-5),
  output_folder="Inverted_Manhattan_Plots/",
  y_axis_squish_threshold = 10,

)



################################ New Code #####################################


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

  # Filter data based on p-value and ensure no missing values
  plot_data <- gwas_data[!is.na(get(p_col)) & get(p_col) > 0 & get(p_col) <= plot_pval_threshold]

  # Standardize column names and calculate -log10(p-value)
  plot_data[, `:=` (pvalue = as.numeric(get(p_col)),
                    CHR = as.character(get(chr_col)),
                    BP = as.numeric(get(bp_col)))]
  plot_data[, log_p := -log10(pvalue)]

  # Sort chromosomes numerically and create a numeric factor for plotting
  plot_data[, CHR_num := as.numeric(factor(CHR, levels = unique(stringr::str_sort(unique(CHR), numeric = TRUE))))]
  plot_data <- plot_data[order(CHR_num, BP)]

  # Calculate cumulative base pair positions for a continuous x-axis
  cumulative_summary <- plot_data[, .(max_bp = max(BP)), by = CHR_num]
  cumulative_summary <- cumulative_summary[order(CHR_num)]
  cumulative_summary$bp_add <- c(0, cumsum(as.numeric(cumulative_summary$max_bp[-nrow(cumulative_summary)])))

  plot_data <- merge(plot_data, cumulative_summary[, .(CHR_num, bp_add)], by = "CHR_num")
  plot_data[, bp_cum := BP + bp_add]

  # Calculate the center position for each chromosome label on the x-axis
  axis_data <- plot_data[, .(center = (max(bp_cum) + min(bp_cum)) / 2), by = CHR]

  return(list(plot_data = plot_data, axis_data = axis_data))
}

#' Helper function to generate a single Manhattan plot (normal or inverted)
#'
#' @param plot_data Processed GWAS data from `prepare_gwas_data_for_plotting`.
#' @param axis_data Axis data for plotting chromosome labels.
#' @param plot_title The main title for the plot.
#' @param sig_lines A named vector defining significance lines and their colors.
#' @param label_threshold_colors A named vector for p-value thresholds and corresponding label colors.
#' @param y_axis_squish_threshold The -log10(p) value at which to start compressing the y-axis.
#' @param gene_label_size The font size for gene labels.
#' @param gene_col The name of the gene column.
#' @param group_col The name of the grouping column for labeling.
#' @param loci_to_label A vector of specific loci to label.
#' @param genes_to_label A vector of specific genes to label.
#' @param custom_gene_colors A named vector for custom colors for specific genes.
#' @param inverted A logical flag to create an inverted (downward-facing) plot.
#' @param font_family The font family to use for plot text.
#'
#' @return A ggplot object representing a single Manhattan plot.
generate_manhattan_subplot <- function(plot_data, axis_data, plot_title, sig_lines,
                                       label_threshold_colors, y_axis_squish_threshold,
                                       gene_label_size, gene_col, group_col,
                                       loci_to_label, genes_to_label, custom_gene_colors,
                                       inverted, font_family = "Arial Unicode MS") {

  # --- Prepare gene labels based on significance or specified lists ---
  snps_to_label <- data.table::data.table()
  if (!is.null(gene_col) && gene_col %in% names(plot_data)) {
    if (is.null(group_col)) group_col <- gene_col
    if (group_col %in% names(plot_data)) {
      # Filter for SNPs that can be labeled
      snps_for_labeling <- plot_data[!is.na(get(group_col))]

      # Prioritize labeling based on user-provided lists or significance thresholds
      if (!is.null(loci_to_label)) {
        snps_for_labeling <- snps_for_labeling[get(group_col) %in% loci_to_label]
      } else if (!is.null(genes_to_label)) {
        snps_for_labeling <- snps_for_labeling[get(gene_col) %in% genes_to_label]
      } else {
        min_label_thresh <- max(label_threshold_colors, na.rm = TRUE)
        snps_for_labeling <- snps_for_labeling[pvalue < min_label_thresh]
      }

      # Select the most significant SNP per group to label
      if (nrow(snps_for_labeling) > 0) {
        snps_to_label <- snps_for_labeling[, .SD[which.min(pvalue)], by = group_col]
      }

      # Assign colors to labels based on significance or custom settings
      if (nrow(snps_to_label) > 0) {
        sorted_thresholds <- sort(label_threshold_colors, decreasing = TRUE)
        snps_to_label[, label_color := names(sorted_thresholds)[1]]
        if (length(sorted_thresholds) > 1) {
          for(i in 2:length(sorted_thresholds)) {
            snps_to_label[pvalue <= sorted_thresholds[i], label_color := names(sorted_thresholds)[i]]
          }
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

  # --- Define y-axis transformation for compressing large values ---
  compression_factor <- 0.1
  squish_forward <- function(x) ifelse(x <= y_axis_squish_threshold, x, y_axis_squish_threshold + (x - y_axis_squish_threshold) * compression_factor)
  squish_inverse <- function(x) ifelse(x <= y_axis_squish_threshold, x, (x - y_axis_squish_threshold) / compression_factor + y_axis_squish_threshold)
  squish_trans <- scales::trans_new("squish", squish_forward, squish_inverse)
  x_axis_expansion <- ggplot2::expansion(mult = c(0.015, 0.015))

  # Define a custom color palette for chromosomes
  custom_chromosome_colors <- c(
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
    "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78"
  )

  # --- Build the plot (inverted or normal) ---
  if (inverted) {
    # Prepare data for inverted plot
    inverted_data <- data.table::copy(plot_data)
    inverted_data[, transformed_log_p := -squish_forward(log_p)]
    desired_breaks <- c(0, 10, 20, 30, 60, 90, 120)
    transformed_breaks <- -squish_forward(desired_breaks)

    p <- ggplot2::ggplot(inverted_data) +
      ggplot2::geom_point(ggplot2::aes(x = bp_cum, y = transformed_log_p, color = factor(CHR_num)), alpha = 0.8, size = 1.5) +
      ggplot2::scale_color_manual(values = rep(custom_chromosome_colors, 2), guide = "none") +
      ggplot2::geom_hline(yintercept = -squish_forward(-log10(sig_lines)), color = names(sig_lines), linetype = "dashed", linewidth = 0.8) +
      ggplot2::scale_x_continuous(label = axis_data$CHR, breaks = axis_data$center, expand = x_axis_expansion) +
      ggplot2::scale_y_continuous(breaks = transformed_breaks, labels = as.character(desired_breaks), position = "left") +
      ggplot2::labs(x = "Chromosome", y = expression(-log[10](italic(P))), caption = plot_title) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        text = ggplot2::element_text(family = font_family),
        panel.grid.major.x = ggplot2::element_blank(), panel.grid.minor.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(color = "black"),
        axis.title.x = ggplot2::element_text(color = "black", face = "bold"),
        plot.caption = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold", color = "black", margin = ggplot2::margin(t = 10)),
        axis.text.y = ggplot2::element_text(angle = 90, hjust = 0.5, color = "black"),
        axis.title.y = ggplot2::element_text(angle = 90, vjust = 0.5, color = "black", face = "bold"),
        plot.margin = ggplot2::margin(t = -15, r = 5.5, b = 5.5, l = 5.5, unit = "pt")
      )

    if (nrow(snps_to_label) > 0) {
      inverted_snps_to_label <- data.table::copy(snps_to_label)
      inverted_snps_to_label[, transformed_log_p := -squish_forward(log_p)]
      p <- p +
        ggnewscale::new_scale_color() +
        # **MODIFIED:** Use geom_label_repel for boxed labels
        ggrepel::geom_label_repel(
          data = inverted_snps_to_label,
          ggplot2::aes(x = bp_cum, y = transformed_log_p, label = label_text, color = label_color),
          size = gene_label_size,
          fontface = "bold.italic",
          nudge_y = -10,
          direction = "x",
          angle = 0, # Horizontal labels are more readable in boxes
          segment.color = 'grey50',
          segment.linetype = "dashed",
          min.segment.length = 0,
          max.overlaps = Inf,
          box.padding = 0.5,
          label.size = 0.25,
          fill = ggplot2::alpha("white", 0.7)
        ) + ggplot2::scale_color_identity(guide = "none")
    }

  } else {
    # Create the normal (top) plot
    p <- ggplot2::ggplot(plot_data) +
      ggplot2::geom_point(ggplot2::aes(x = bp_cum, y = log_p, color = factor(CHR_num)), alpha = 0.8, size = 1.5) +
      ggplot2::scale_color_manual(values = rep(custom_chromosome_colors, 2), guide = "none") +
      ggplot2::geom_hline(yintercept = -log10(sig_lines), color = names(sig_lines), linetype = "dashed", linewidth = 0.8) +
      ggplot2::scale_x_continuous(breaks = axis_data$center, expand = x_axis_expansion) +
      ggplot2::scale_y_continuous(trans = squish_trans, breaks = c(0, 10, 20, 30, 60, 90, 120), expand = ggplot2::expansion(mult = c(0, 0.1))) +
      ggplot2::labs(x = NULL, y = expression(-log[10](italic(P))), title = plot_title) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        text = ggplot2::element_text(family = font_family),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold", color = "black"),
        panel.grid.major.x = ggplot2::element_blank(), panel.grid.minor.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(color = "black"),
        axis.title.y = ggplot2::element_text(color = "black", face = "bold"),
        axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(t = 5.5, r = 5.5, b = -15, l = 5.5, unit = "pt")
      )

    if (nrow(snps_to_label) > 0) {
      p <- p +
        ggnewscale::new_scale_color() +
        # **MODIFIED:** Use geom_label_repel for boxed labels
        ggrepel::geom_label_repel(
          data = snps_to_label,
          ggplot2::aes(x = bp_cum, y = log_p, label = label_text, color = label_color),
          size = gene_label_size,
          fontface = "bold.italic",
          nudge_y = 15,
          direction = "x",
          angle = 0, # Horizontal labels are more readable in boxes
          segment.color = 'grey50',
          segment.linetype = "dashed",
          min.segment.length = 0,
          max.overlaps = Inf,
          box.padding = 0.5,
          label.size = 0.25,
          fill = ggplot2::alpha("white", 0.7)
        ) + ggplot2::scale_color_identity(guide = "none")
    }
  }
  return(p)
}


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

  # Construct the title string with study metadata
  title_string <- sprintf("%s\n(%s: %s | \U03BB = %s | Cases: %s, Controls: %s)",
                          plot_title_prefix, snps_title_label, snps_for_title,
                          lambda, format(n_cases, big.mark=","), format(n_controls, big.mark=","))

  if (add_date_to_title) {
    title_string <- paste0(title_string, " | ", Sys.Date())
  }
  return(title_string)
}


#' Create a Pair of Inverted Manhattan Plots from Two GWAS Datasets
#'
#' Generates and saves a mirrored Manhattan plot comparison between two GWAS datasets.
#'
#' @import data.table
#' @import ggplot2
#' @import ggrepel
#' @import scales
#' @import stringr
#' @import ragg
#' @import ggnewscale
#' @import patchwork
#' @import extrafont
#' @param gwas_data1 The first GWAS dataset (a data.frame or data.table).
#' @param gwas_data2 The second GWAS dataset (a data.frame or data.table).
#' @param chr_col The name of the chromosome column (default: "CHR").
#' @param bp_col The name of the base pair position column (default: "BP").
#' @param p_col The name of the p-value column (default: "P").
#' @param gene_col The name of the gene column used for labeling (default: "GENE").
#' @param group_col The name of a grouping column to select one label per group (e.g., locus name). Defaults to `gene_col`.
#' @param loci_to_label An optional vector of specific locus names to be labeled.
#' @param genes_to_label An optional vector of specific gene names to be labeled.
#' @param custom_gene_colors An optional named vector to apply custom colors to specific gene labels.
#' @param n_cases1 The number of cases in the first dataset.
#' @param n_controls1 The number of controls in the first dataset.
#' @param n_cases2 The number of cases in the second dataset.
#' @param n_controls2 The number of controls in the second dataset.
#' @param plot_title1 A title for the first (top) plot (default: "Dataset 1").
#' @param plot_title2 A title for the second (bottom) plot (default: "Dataset 2").
#' @param lambda1 The genomic inflation factor for the first dataset (auto-calculated if NULL).
#' @param lambda2 The genomic inflation factor for the second dataset (auto-calculated if NULL).
#' @param total_snps_in_study1 The total number of SNPs in the first study (uses SNPs in data if NULL).
#' @param total_snps_in_study2 The total number of SNPs in the second study (uses SNPs in data if NULL).
#' @param add_date_to_title A logical flag to add the current date to the plot titles (default: FALSE).
#' @param plot_pval_threshold The maximum p-value to include in the plot (default: 1).
#' @param sig_lines A named vector for significance lines and their colors (default: `c("red" = 5e-8, "blue" = 1e-5)`).
#' @param label_threshold_colors A named vector for p-value thresholds that trigger labeling and the color of the label (default: `c("red" = 5e-8, "orange" = 1e-6, "darkblue" = 1e-5)`).
#' @param y_axis_squish_threshold The -log10(p) value at which to begin compressing the y-axis (default: 30).
#' @param gene_label_size The font size for gene labels (default: 3.5).
#' @param output_folder The folder where the output plots will be saved (default: "Inverted_Manhattan_Plots").
#' @param file_name_prefix The file name prefix for the saved PNG and PDF files (default: "inverted_manhattan").
#' @param font_family The font family to use for all plot text, chosen for Unicode support (default: "Arial Unicode MS").
#' @export
#' @return Invisibly returns the combined plot object created by `patchwork`.
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' data("gwasdataseta")
#' data1 <- gwasdataseta
#' data2 <- data.table::copy(gwasdataseta)
#' create_mirrored_manhattan_plot(
#'   gwas_data1 = data1,
#'   gwas_data2 = data2,
#'   n_cases1 = 10000, n_controls1 = 10000,
#'   n_cases2 = 8000, n_controls2 = 8000,
#'   plot_title1 = "Study 1", plot_title2 = "Study 2"
#' )
#' }
create_mirrored_manhattan_plot <- function(gwas_data1, gwas_data2,
                                           chr_col = "chr",
                                           bp_col = "pos",
                                           p_col = "pvalue",
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
                                           output_folder = "Inverted_Manhattan_Plots",
                                           file_name_prefix = "inverted_manhattan",
                                           font_family = "Arial Unicode MS") {

  # --- Initial Setup ---
  data1 <- data.table::as.data.table(gwas_data1)
  data2 <- data.table::as.data.table(gwas_data2)

  # Use a common set of chromosomes from both datasets to align plots
  all_chrs <- sort(unique(c(data1[[chr_col]], data2[[chr_col]])))

  # --- Process Datasets ---
  processed_data1 <- prepare_gwas_data_for_plotting(
    gwas_data = data1, chr_col = chr_col, bp_col = bp_col, p_col = p_col, gene_col = gene_col,
    group_col = group_col, plot_pval_threshold = plot_pval_threshold, all_chrs = all_chrs
  )
  processed_data2 <- prepare_gwas_data_for_plotting(
    gwas_data = data2, chr_col = chr_col, bp_col = bp_col, p_col = p_col, gene_col = gene_col,
    group_col = group_col, plot_pval_threshold = plot_pval_threshold, all_chrs = all_chrs
  )

  # --- Calculate Lambda GC if not provided ---
  if (is.null(lambda1)) {
    chisq1 <- stats::qchisq(1 - processed_data1$plot_data$pvalue, 1)
    lambda1 <- round(median(chisq1, na.rm = TRUE) / stats::qchisq(0.5, 1), 4)
  }
  if (is.null(lambda2)) {
    chisq2 <- stats::qchisq(1 - processed_data2$plot_data$pvalue, 1)
    lambda2 <- round(median(chisq2, na.rm = TRUE) / stats::qchisq(0.5, 1), 4)
  }

  # --- Generate Plot Titles ---
  title1 <- construct_plot_title(
    plot_title_prefix = plot_title1, n_cases = n_cases1, n_controls = n_controls1, lambda = lambda1,
    total_snps_in_study = total_snps_in_study1, add_date_to_title = add_date_to_title, gwas_data = data1, font_family = font_family
  )
  title2 <- construct_plot_title(
    plot_title_prefix = plot_title2, n_cases = n_cases2, n_controls = n_controls2, lambda = lambda2,
    total_snps_in_study = total_snps_in_study2, add_date_to_title = add_date_to_title, gwas_data = data2, font_family = font_family
  )

  # --- Create Subplots ---
  top_plot <- generate_manhattan_subplot(
    plot_data = processed_data1$plot_data, axis_data = processed_data1$axis_data, plot_title = title1,
    sig_lines = sig_lines, label_threshold_colors = label_threshold_colors, y_axis_squish_threshold = y_axis_squish_threshold,
    gene_label_size = gene_label_size, gene_col = gene_col, group_col = group_col, loci_to_label = loci_to_label,
    genes_to_label = genes_to_label, custom_gene_colors = custom_gene_colors, inverted = FALSE, font_family = font_family
  )
  bottom_plot <- generate_manhattan_subplot(
    plot_data = processed_data2$plot_data, axis_data = processed_data2$axis_data, plot_title = title2,
    sig_lines = sig_lines, label_threshold_colors = label_threshold_colors, y_axis_squish_threshold = y_axis_squish_threshold,
    gene_label_size = gene_label_size, gene_col = gene_col, group_col = group_col, loci_to_label = loci_to_label,
    genes_to_label = genes_to_label, custom_gene_colors = custom_gene_colors, inverted = TRUE, font_family = font_family
  )

  # --- Combine and Save Plots ---
  combined_plot <- top_plot / bottom_plot + patchwork::plot_layout(heights = c(1, 1))

  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  file_path_png <- file.path(output_folder, paste0(file_name_prefix, ".png"))
  file_path_pdf <- file.path(output_folder, paste0(file_name_prefix, ".pdf"))

  # Attempt to save with Cairo for best compatibility, with fallbacks
  tryCatch({
    grDevices::png(filename = file_path_png, width = 18, height = 12, units = "in", res = 300, type = "cairo", bg = "white")
    print(combined_plot)
    grDevices::dev.off()
    grDevices::cairo_pdf(filename = file_path_pdf, width = 18, height = 12)
    print(combined_plot)
    grDevices::dev.off()
  }, error = function(e) {
    message("Cairo devices failed, falling back to standard ggsave devices.")
    ggplot2::ggsave(file_path_png, plot = combined_plot, device = ragg::agg_png, width = 18, height = 12, units = "in", res = 300, bg = "white")
    ggplot2::ggsave(file_path_pdf, plot = combined_plot, device = "pdf", width = 18, height = 12, units = "in")
  })

  message(paste("Successfully saved plots to:", file_path_png, "and", file_path_pdf))
  return(invisible(combined_plot))
}




################### Vertical (nice Code)

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

  # Filter data based on p-value and ensure no missing values
  plot_data <- gwas_data[!is.na(get(p_col)) & get(p_col) > 0 & get(p_col) <= plot_pval_threshold]

  # Standardize column names and calculate -log10(p-value)
  plot_data[, `:=` (pvalue = as.numeric(get(p_col)),
                    CHR = as.character(get(chr_col)),
                    BP = as.numeric(get(bp_col)))]
  plot_data[, log_p := -log10(pvalue)]

  # Sort chromosomes numerically and create a numeric factor for plotting
  plot_data[, CHR_num := as.numeric(factor(CHR, levels = unique(stringr::str_sort(unique(CHR), numeric = TRUE))))]
  plot_data <- plot_data[order(CHR_num, BP)]

  # Calculate cumulative base pair positions for a continuous x-axis
  cumulative_summary <- plot_data[, .(max_bp = max(BP)), by = CHR_num]
  cumulative_summary <- cumulative_summary[order(CHR_num)]
  cumulative_summary$bp_add <- c(0, cumsum(as.numeric(cumulative_summary$max_bp[-nrow(cumulative_summary)])))

  plot_data <- merge(plot_data, cumulative_summary[, .(CHR_num, bp_add)], by = "CHR_num")
  plot_data[, bp_cum := BP + bp_add]

  # Calculate the center position for each chromosome label on the x-axis
  axis_data <- plot_data[, .(center = (max(bp_cum) + min(bp_cum)) / 2), by = CHR]

  return(list(plot_data = plot_data, axis_data = axis_data))
}

#' Helper function to generate a single Manhattan plot (normal or inverted)
#'
#' @param plot_data Processed GWAS data from `prepare_gwas_data_for_plotting`.
#' @param axis_data Axis data for plotting chromosome labels.
#' @param plot_title The main title for the plot.
#' @param sig_lines A named vector defining significance lines and their colors.
#' @param label_threshold_colors A named vector for p-value thresholds and corresponding label colors.
#' @param y_axis_squish_threshold The -log10(p) value at which to start compressing the y-axis.
#' @param gene_label_size The font size for gene labels.
#' @param gene_col The name of the gene column.
#' @param group_col The name of the grouping column for labeling.
#' @param loci_to_label A vector of specific loci to label.
#' @param genes_to_label A vector of specific genes to label.
#' @param custom_gene_colors A named vector for custom colors for specific genes.
#' @param inverted A logical flag to create an inverted (downward-facing) plot.
#' @param font_family The font family to use for plot text.
#'
#' @return A ggplot object representing a single Manhattan plot.
generate_manhattan_subplot <- function(plot_data, axis_data, plot_title, sig_lines,
                                       label_threshold_colors, y_axis_squish_threshold,
                                       gene_label_size, gene_col, group_col,
                                       loci_to_label, genes_to_label, custom_gene_colors,
                                       inverted, font_family = "Arial Unicode MS") {

  # --- Prepare gene labels ---
  snps_to_label <- data.table::data.table()
  if (!is.null(gene_col) && gene_col %in% names(plot_data)) {
    if (is.null(group_col)) group_col <- gene_col
    if (group_col %in% names(plot_data)) {
      snps_for_labeling <- plot_data[!is.na(get(group_col))]
      if (!is.null(loci_to_label)) {
        snps_for_labeling <- snps_for_labeling[get(group_col) %in% loci_to_label]
      } else if (!is.null(genes_to_label)) {
        snps_for_labeling <- snps_for_labeling[get(gene_col) %in% genes_to_label]
      } else {
        min_label_thresh <- max(label_threshold_colors, na.rm = TRUE)
        snps_for_labeling <- snps_for_labeling[pvalue < min_label_thresh]
      }
      if (nrow(snps_for_labeling) > 0) {
        snps_to_label <- snps_for_labeling[, .SD[which.min(pvalue)], by = group_col]
      }
      if (nrow(snps_to_label) > 0) {
        sorted_thresholds <- sort(label_threshold_colors, decreasing = TRUE)
        snps_to_label[, label_color := names(sorted_thresholds)[1]]
        if (length(sorted_thresholds) > 1) {
          for(i in 2:length(sorted_thresholds)) {
            snps_to_label[pvalue <= sorted_thresholds[i], label_color := names(sorted_thresholds)[i]]
          }
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

  # --- Define plot parameters ---
  compression_factor <- 0.1
  squish_forward <- function(x) ifelse(x <= y_axis_squish_threshold, x, y_axis_squish_threshold + (x - y_axis_squish_threshold) * compression_factor)
  squish_inverse <- function(x) ifelse(x <= y_axis_squish_threshold, x, (x - y_axis_squish_threshold) / compression_factor + y_axis_squish_threshold)
  squish_trans <- scales::trans_new("squish", squish_forward, squish_inverse)
  x_axis_expansion <- ggplot2::expansion(mult = c(0.015, 0.015))
  custom_chromosome_colors <- c(
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
    "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78"
  )

  # --- Build the plot ---
  if (inverted) {
    inverted_data <- data.table::copy(plot_data)
    inverted_data[, transformed_log_p := -squish_forward(log_p)]
    desired_breaks <- c(0, 10, 20, 30, 60, 90, 120)
    transformed_breaks <- -squish_forward(desired_breaks)

    p <- ggplot2::ggplot(inverted_data) +
      ggplot2::geom_point(ggplot2::aes(x = bp_cum, y = transformed_log_p, color = factor(CHR_num)), alpha = 0.8, size = 1.5) +
      ggplot2::scale_color_manual(values = rep(custom_chromosome_colors, 2), guide = "none") +
      ggplot2::geom_hline(yintercept = -squish_forward(-log10(sig_lines)), color = names(sig_lines), linetype = "dashed", linewidth = 0.8) +
      ggplot2::scale_x_continuous(label = axis_data$CHR, breaks = axis_data$center, expand = x_axis_expansion) +
      ggplot2::scale_y_continuous(breaks = transformed_breaks, labels = as.character(desired_breaks), position = "left") +
      ggplot2::labs(x = "Chromosome", y = expression(-log[10](italic(P))), caption = plot_title) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        text = ggplot2::element_text(family = font_family),
        panel.grid.major.x = ggplot2::element_blank(), panel.grid.minor.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(color = "black"),
        axis.title.x = ggplot2::element_text(color = "black", face = "bold"),
        plot.caption = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold", color = "black", margin = ggplot2::margin(t = 10)),
        axis.text.y = ggplot2::element_text(angle = 90, hjust = 0.5, color = "black"),
        axis.title.y = ggplot2::element_text(angle = 90, vjust = 0.5, color = "black", face = "bold"),
        plot.margin = ggplot2::margin(t = -15, r = 5.5, b = 5.5, l = 5.5, unit = "pt")
      )

    if (nrow(snps_to_label) > 0) {
      inverted_snps_to_label <- data.table::copy(snps_to_label)
      inverted_snps_to_label[, transformed_log_p := -squish_forward(log_p)]
      p <- p +
        ggnewscale::new_scale_color() +
        ggrepel::geom_label_repel(
          data = inverted_snps_to_label,
          ggplot2::aes(x = bp_cum, y = transformed_log_p, label = label_text, color = label_color),
          size = gene_label_size,
          fontface = "bold.italic",
          nudge_y = -15,
          angle = 90,
          hjust = 1,
          segment.color = 'grey50',
          segment.linetype = "dashed",
          min.segment.length = 0,
          max.overlaps = Inf,
          box.padding = 0.5,
          label.size = 0.25,
          fill = ggplot2::alpha("white", 0.7)
        ) + ggplot2::scale_color_identity(guide = "none")
    }

  } else {
    p <- ggplot2::ggplot(plot_data) +
      ggplot2::geom_point(ggplot2::aes(x = bp_cum, y = log_p, color = factor(CHR_num)), alpha = 0.8, size = 1.5) +
      ggplot2::scale_color_manual(values = rep(custom_chromosome_colors, 2), guide = "none") +
      ggplot2::geom_hline(yintercept = -log10(sig_lines), color = names(sig_lines), linetype = "dashed", linewidth = 0.8) +
      ggplot2::scale_x_continuous(breaks = axis_data$center, expand = x_axis_expansion) +
      ggplot2::scale_y_continuous(trans = squish_trans, breaks = c(0, 10, 20, 30, 60, 90, 120), expand = ggplot2::expansion(mult = c(0, 0.1))) +
      ggplot2::labs(x = NULL, y = expression(-log[10](italic(P))), title = plot_title) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        text = ggplot2::element_text(family = font_family),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold", color = "black"),
        panel.grid.major.x = ggplot2::element_blank(), panel.grid.minor.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(color = "black"),
        axis.title.y = ggplot2::element_text(color = "black", face = "bold"),
        axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(t = 5.5, r = 5.5, b = -15, l = 5.5, unit = "pt")
      )

    if (nrow(snps_to_label) > 0) {
      p <- p +
        ggnewscale::new_scale_color() +
        ggrepel::geom_label_repel(
          data = snps_to_label,
          ggplot2::aes(x = bp_cum, y = log_p, label = label_text, color = label_color),
          size = gene_label_size,
          fontface = "bold.italic",
          nudge_y = 15,
          angle = 90,
          hjust = 0,
          segment.color = 'grey50',
          segment.linetype = "dashed",
          min.segment.length = 0,
          max.overlaps = Inf,
          box.padding = 0.5,
          label.size = 0.25,
          fill = ggplot2::alpha("white", 0.7)
        ) + ggplot2::scale_color_identity(guide = "none")
    }
  }
  return(p)
}


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

  # Construct the title string with study metadata
  title_string <- sprintf("%s\n(%s: %s | \U03BB = %s | Cases: %s, Controls: %s)",
                          plot_title_prefix, snps_title_label, snps_for_title,
                          lambda, format(n_cases, big.mark=","), format(n_controls, big.mark=","))

  if (add_date_to_title) {
    title_string <- paste0(title_string, " | ", Sys.Date())
  }
  return(title_string)
}


#' Create a Pair of Inverted Manhattan Plots from Two GWAS Datasets
#'
#' Generates and saves a mirrored Manhattan plot comparison between two GWAS datasets.
#'
#' @import data.table
#' @import ggplot2
#' @import ggrepel
#' @import scales
#' @import stringr
#' @import ragg
#' @import ggnewscale
#' @import patchwork
#' @import extrafont
#' @param gwas_data1 The first GWAS dataset (a data.frame or data.table).
#' @param gwas_data2 The second GWAS dataset (a data.frame or data.table).
#' @param chr_col The name of the chromosome column (default: "CHR").
#' @param bp_col The name of the base pair position column (default: "BP").
#' @param p_col The name of the p-value column (default: "P").
#' @param gene_col The name of the gene column used for labeling (default: "GENE").
#' @param group_col The name of a grouping column to select one label per group (e.g., locus name). Defaults to `gene_col`.
#' @param loci_to_label An optional vector of specific locus names to be labeled.
#' @param genes_to_label An optional vector of specific gene names to be labeled.
#' @param custom_gene_colors An optional named vector to apply custom colors to specific gene labels.
#' @param n_cases1 The number of cases in the first dataset.
#' @param n_controls1 The number of controls in the first dataset.
#' @param n_cases2 The number of cases in the second dataset.
#' @param n_controls2 The number of controls in the second dataset.
#' @param plot_title1 A title for the first (top) plot (default: "Dataset 1").
#' @param plot_title2 A title for the second (bottom) plot (default: "Dataset 2").
#' @param lambda1 The genomic inflation factor for the first dataset (auto-calculated if NULL).
#' @param lambda2 The genomic inflation factor for the second dataset (auto-calculated if NULL).
#' @param total_snps_in_study1 The total number of SNPs in the first study (uses SNPs in data if NULL).
#' @param total_snps_in_study2 The total number of SNPs in the second study (uses SNPs in data if NULL).
#' @param add_date_to_title A logical flag to add the current date to the plot titles (default: FALSE).
#' @param plot_pval_threshold The maximum p-value to include in the plot (default: 1).
#' @param sig_lines A named vector for significance lines and their colors (default: `c("red" = 5e-8, "blue" = 1e-5)`).
#' @param label_threshold_colors A named vector for p-value thresholds that trigger labeling and the color of the label (default: `c("red" = 5e-8, "orange" = 1e-6, "darkblue" = 1e-5)`).
#' @param y_axis_squish_threshold The -log10(p) value at which to begin compressing the y-axis (default: 30).
#' @param gene_label_size The font size for gene labels (default: 3.5).
#' @param output_folder The folder where the output plots will be saved (default: "Inverted_Manhattan_Plots").
#' @param file_name_prefix The file name prefix for the saved PNG and PDF files (default: "inverted_manhattan").
#' @param font_family The font family to use for all plot text, chosen for Unicode support (default: "Arial Unicode MS").
#' @export
#' @return Invisibly returns the combined plot object created by `patchwork`.
#'
create_mirrored_manhattan_plot <- function(gwas_data1, gwas_data2,
                                           chr_col = "chr",
                                           bp_col = "pos",
                                           p_col = "pvalue",
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
                                           output_folder = "Inverted_Manhattan_Plots",
                                           file_name_prefix = "inverted_manhattan",
                                           font_family = "Arial Unicode MS") {

  # --- Initial Setup ---
  data1 <- data.table::as.data.table(gwas_data1)
  data2 <- data.table::as.data.table(gwas_data2)

  # Use a common set of chromosomes from both datasets to align plots
  all_chrs <- sort(unique(c(data1[[chr_col]], data2[[chr_col]])))

  # --- Process Datasets ---
  processed_data1 <- prepare_gwas_data_for_plotting(
    gwas_data = data1, chr_col = chr_col, bp_col = bp_col, p_col = p_col, gene_col = gene_col,
    group_col = group_col, plot_pval_threshold = plot_pval_threshold, all_chrs = all_chrs
  )
  processed_data2 <- prepare_gwas_data_for_plotting(
    gwas_data = data2, chr_col = chr_col, bp_col = bp_col, p_col = p_col, gene_col = gene_col,
    group_col = group_col, plot_pval_threshold = plot_pval_threshold, all_chrs = all_chrs
  )

  # --- Calculate Lambda GC if not provided ---
  if (is.null(lambda1)) {
    chisq1 <- stats::qchisq(1 - processed_data1$plot_data$pvalue, 1)
    lambda1 <- round(median(chisq1, na.rm = TRUE) / stats::qchisq(0.5, 1), 4)
  }
  if (is.null(lambda2)) {
    chisq2 <- stats::qchisq(1 - processed_data2$plot_data$pvalue, 1)
    lambda2 <- round(median(chisq2, na.rm = TRUE) / stats::qchisq(0.5, 1), 4)
  }

  # --- Generate Plot Titles ---
  title1 <- construct_plot_title(
    plot_title_prefix = plot_title1, n_cases = n_cases1, n_controls = n_controls1, lambda = lambda1,
    total_snps_in_study = total_snps_in_study1, add_date_to_title = add_date_to_title, gwas_data = data1, font_family = font_family
  )
  title2 <- construct_plot_title(
    plot_title_prefix = plot_title2, n_cases = n_cases2, n_controls = n_controls2, lambda = lambda2,
    total_snps_in_study = total_snps_in_study2, add_date_to_title = add_date_to_title, gwas_data = data2, font_family = font_family
  )

  # --- Create Subplots ---
  top_plot <- generate_manhattan_subplot(
    plot_data = processed_data1$plot_data, axis_data = processed_data1$axis_data, plot_title = title1,
    sig_lines = sig_lines, label_threshold_colors = label_threshold_colors, y_axis_squish_threshold = y_axis_squish_threshold,
    gene_label_size = gene_label_size, gene_col = gene_col, group_col = group_col, loci_to_label = loci_to_label,
    genes_to_label = genes_to_label, custom_gene_colors = custom_gene_colors, inverted = FALSE, font_family = font_family
  )
  bottom_plot <- generate_manhattan_subplot(
    plot_data = processed_data2$plot_data, axis_data = processed_data2$axis_data, plot_title = title2,
    sig_lines = sig_lines, label_threshold_colors = label_threshold_colors, y_axis_squish_threshold = y_axis_squish_threshold,
    gene_label_size = gene_label_size, gene_col = gene_col, group_col = group_col, loci_to_label = loci_to_label,
    genes_to_label = genes_to_label, custom_gene_colors = custom_gene_colors, inverted = TRUE, font_family = font_family
  )

  # --- Combine and Save Plots ---
  combined_plot <- top_plot / bottom_plot + patchwork::plot_layout(heights = c(1, 1))

  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  file_path_png <- file.path(output_folder, paste0(file_name_prefix, ".png"))
  file_path_pdf <- file.path(output_folder, paste0(file_name_prefix, ".pdf"))

  # Attempt to save with Cairo for best compatibility, with fallbacks
  tryCatch({
    grDevices::png(filename = file_path_png, width = 18, height = 12, units = "in", res = 300, type = "cairo", bg = "white")
    print(combined_plot)
    grDevices::dev.off()
    grDevices::cairo_pdf(filename = file_path_pdf, width = 18, height = 12)
    print(combined_plot)
    grDevices::dev.off()
  }, error = function(e) {
    message("Cairo devices failed, falling back to standard ggsave devices.")
    ggplot2::ggsave(file_path_png, plot = combined_plot, device = ragg::agg_png, width = 18, height = 12, units = "in", res = 300, bg = "white")
    ggplot2::ggsave(file_path_pdf, plot = combined_plot, device = "pdf", width = 18, height = 12, units = "in")
  })

  message(paste("Successfully saved plots to:", file_path_png, "and", file_path_pdf))
  return(invisible(combined_plot))
}

######################### NEW Codes ( this nice too but the horzinatiln is related to the dispertion of allels)
#' Helper function to generate a single Manhattan plot (normal or inverted)
#'
#' @param plot_data Processed GWAS data from `prepare_gwas_data_for_plotting`.
#' @param axis_data Axis data for plotting chromosome labels.
#' @param plot_title The main title for the plot.
#' @param sig_lines A named vector defining significance lines and their colors.
#' @param label_threshold_colors A named vector for p-value thresholds and corresponding label colors.
#' @param y_axis_squish_threshold The -log10(p) value at which to start compressing the y-axis.
#' @param gene_label_size The font size for gene labels.
#' @param gene_col The name of the gene column.
#' @param group_col The name of the grouping column for labeling.
#' @param loci_to_label A vector of specific loci to label.
#' @param genes_to_label A vector of specific genes to label.
#' @param custom_gene_colors A named vector for custom colors for specific genes.
#' @param inverted A logical flag to create an inverted (downward-facing) plot.
#' @param font_family The font family to use for plot text.
#' @param label_orientation The orientation of gene labels ("vertical" or "horizontal").
#'
#' @return A ggplot object representing a single Manhattan plot.
generate_manhattan_subplot <- function(plot_data, axis_data, plot_title, sig_lines,
                                       label_threshold_colors, y_axis_squish_threshold,
                                       gene_label_size, gene_col, group_col,
                                       loci_to_label, genes_to_label, custom_gene_colors,
                                       inverted, font_family = "Arial Unicode MS",
                                       label_orientation = "horizontal") { # <-- New parameter

  # --- Prepare gene labels ---
  snps_to_label <- data.table::data.table()
  if (!is.null(gene_col) && gene_col %in% names(plot_data)) {
    if (is.null(group_col)) group_col <- gene_col
    if (group_col %in% names(plot_data)) {
      snps_for_labeling <- plot_data[!is.na(get(group_col))]
      if (!is.null(loci_to_label)) {
        snps_for_labeling <- snps_for_labeling[get(group_col) %in% loci_to_label]
      } else if (!is.null(genes_to_label)) {
        snps_for_labeling <- snps_for_labeling[get(gene_col) %in% genes_to_label]
      } else {
        min_label_thresh <- max(label_threshold_colors, na.rm = TRUE)
        snps_for_labeling <- snps_for_labeling[pvalue < min_label_thresh]
      }
      if (nrow(snps_for_labeling) > 0) {
        snps_to_label <- snps_for_labeling[, .SD[which.min(pvalue)], by = group_col]
      }
      if (nrow(snps_to_label) > 0) {
        sorted_thresholds <- sort(label_threshold_colors, decreasing = TRUE)
        snps_to_label[, label_color := names(sorted_thresholds)[1]]
        if (length(sorted_thresholds) > 1) {
          for(i in 2:length(sorted_thresholds)) {
            snps_to_label[pvalue <= sorted_thresholds[i], label_color := names(sorted_thresholds)[i]]
          }
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

  # --- Define plot parameters ---
  compression_factor <- 0.1
  squish_forward <- function(x) ifelse(x <= y_axis_squish_threshold, x, y_axis_squish_threshold + (x - y_axis_squish_threshold) * compression_factor)
  squish_inverse <- function(x) ifelse(x <= y_axis_squish_threshold, x, (x - y_axis_squish_threshold) / compression_factor + y_axis_squish_threshold)
  squish_trans <- scales::trans_new("squish", squish_forward, squish_inverse)
  x_axis_expansion <- ggplot2::expansion(mult = c(0.015, 0.015))
  custom_chromosome_colors <- c(
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
    "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78"
  )

  # --- Build the plot ---
  if (inverted) {
    inverted_data <- data.table::copy(plot_data)
    inverted_data[, transformed_log_p := -squish_forward(log_p)]
    desired_breaks <- c(0, 10, 20, 30, 60, 90, 120)
    transformed_breaks <- -squish_forward(desired_breaks)

    p <- ggplot2::ggplot(inverted_data) +
      ggplot2::geom_point(ggplot2::aes(x = bp_cum, y = transformed_log_p, color = factor(CHR_num)), alpha = 0.8, size = 1.5) +
      ggplot2::scale_color_manual(values = rep(custom_chromosome_colors, 2), guide = "none") +
      ggplot2::geom_hline(yintercept = -squish_forward(-log10(sig_lines)), color = names(sig_lines), linetype = "dashed", linewidth = 0.8) +
      ggplot2::scale_x_continuous(label = axis_data$CHR, breaks = axis_data$center, expand = x_axis_expansion) +
      ggplot2::scale_y_continuous(breaks = transformed_breaks, labels = as.character(desired_breaks), position = "left") +
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
      ggplot2::scale_y_continuous(trans = squish_trans, breaks = c(0, 10, 20, 30, 60, 90, 120), expand = ggplot2::expansion(mult = c(0, 0.1))) +
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

  # --- Add Gene Labels Based on Orientation Choice ---
  if (nrow(snps_to_label) > 0) {

    label_data <- if(inverted) {
      data.table::copy(snps_to_label)[, transformed_log_p := -squish_forward(log_p)]
    } else {
      snps_to_label
    }

    y_aes <- if(inverted) "transformed_log_p" else "log_p"
    nudge_val <- if(inverted) -15 else 15

    # Define the base repel layer
    repel_layer <- NULL

    if (label_orientation == "vertical") {
      repel_layer <- ggrepel::geom_label_repel(
        data = label_data,
        ggplot2::aes_string(x = "bp_cum", y = y_aes, label = "label_text", color = "label_color"),
        size = gene_label_size,
        fontface = "bold.italic",
        nudge_y = nudge_val,
        angle = 90,
        hjust = 0.5,
        direction = "y",
        segment.color = 'grey50',
        segment.linetype = "dashed",
        min.segment.length = 0,
        max.overlaps = Inf,
        box.padding = 0.8,
        label.size = 0.25,
        fill = ggplot2::alpha("white", 0.7)
      )
    } else { # Default to horizontal
      repel_layer <- ggrepel::geom_label_repel(
        data = label_data,
        ggplot2::aes_string(x = "bp_cum", y = y_aes, label = "label_text", color = "label_color"),
        size = gene_label_size,
        fontface = "bold.italic",
        nudge_y = nudge_val,
        angle = 0,
        direction = "x",
        segment.color = 'grey50',
        segment.linetype = "dashed",
        min.segment.length = 0,
        max.overlaps = Inf,
        box.padding = 0.5,
        label.size = 0.25,
        fill = ggplot2::alpha("white", 0.7)
      )
    }

    p <- p + ggnewscale::new_scale_color() + repel_layer + ggplot2::scale_color_identity(guide = "none")
  }

  return(p)
}


#' Create a Pair of Inverted Manhattan Plots from Two GWAS Datasets
#'
#' Generates and saves a mirrored Manhattan plot comparison between two GWAS datasets.
#'
#' @param gwas_data1 The first GWAS dataset (a data.frame or data.table).
#' @param gwas_data2 The second GWAS dataset (a data.frame or data.table).
#' @param chr_col The name of the chromosome column (default: "CHR").
#' @param bp_col The name of the base pair position column (default: "BP").
#' @param p_col The name of the p-value column (default: "P").
#' @param gene_col The name of the gene column used for labeling (default: "GENE").
#' @param label_orientation The orientation of gene labels. Can be "vertical" or "horizontal" (default).
#' @param ... other parameters
#' @export
create_mirrored_manhattan_plot <- function(gwas_data1, gwas_data2,
                                           chr_col = "chr",
                                           bp_col = "pos",
                                           p_col = "pvalue",
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
                                           label_orientation = "vertical", # <-- New parameter
                                           output_folder = "Inverted_Manhattan_Plots",
                                           file_name_prefix = "inverted_manhattan",
                                           font_family = "Arial Unicode MS") {

  # --- Initial Setup ---
  data1 <- data.table::as.data.table(gwas_data1)
  data2 <- data.table::as.data.table(gwas_data2)
  all_chrs <- sort(unique(c(data1[[chr_col]], data2[[chr_col]])))

  # --- Process Datasets ---
  processed_data1 <- prepare_gwas_data_for_plotting(
    gwas_data = data1, chr_col = chr_col, bp_col = bp_col, p_col = p_col, gene_col = gene_col,
    group_col = group_col, plot_pval_threshold = plot_pval_threshold, all_chrs = all_chrs
  )
  processed_data2 <- prepare_gwas_data_for_plotting(
    gwas_data = data2, chr_col = chr_col, bp_col = bp_col, p_col = p_col, gene_col = gene_col,
    group_col = group_col, plot_pval_threshold = plot_pval_threshold, all_chrs = all_chrs
  )

  # --- Calculate Lambda GC ---
  if (is.null(lambda1)) {
    chisq1 <- stats::qchisq(1 - processed_data1$plot_data$pvalue, 1)
    lambda1 <- round(median(chisq1, na.rm = TRUE) / stats::qchisq(0.5, 1), 4)
  }
  if (is.null(lambda2)) {
    chisq2 <- stats::qchisq(1 - processed_data2$plot_data$pvalue, 1)
    lambda2 <- round(median(chisq2, na.rm = TRUE) / stats::qchisq(0.5, 1), 4)
  }

  # --- Generate Plot Titles ---
  title1 <- construct_plot_title(
    plot_title_prefix = plot_title1, n_cases = n_cases1, n_controls = n_controls1, lambda = lambda1,
    total_snps_in_study = total_snps_in_study1, add_date_to_title = add_date_to_title, gwas_data = data1, font_family = font_family
  )
  title2 <- construct_plot_title(
    plot_title_prefix = plot_title2, n_cases = n_cases2, n_controls = n_controls2, lambda = lambda2,
    total_snps_in_study = total_snps_in_study2, add_date_to_title = add_date_to_title, gwas_data = data2, font_family = font_family
  )

  # --- Create Subplots ---
  top_plot <- generate_manhattan_subplot(
    plot_data = processed_data1$plot_data, axis_data = processed_data1$axis_data, plot_title = title1,
    sig_lines = sig_lines, label_threshold_colors = label_threshold_colors, y_axis_squish_threshold = y_axis_squish_threshold,
    gene_label_size = gene_label_size, gene_col = gene_col, group_col = group_col, loci_to_label = loci_to_label,
    genes_to_label = genes_to_label, custom_gene_colors = custom_gene_colors, inverted = FALSE, font_family = font_family,
    label_orientation = label_orientation # Pass parameter
  )
  bottom_plot <- generate_manhattan_subplot(
    plot_data = processed_data2$plot_data, axis_data = processed_data2$axis_data, plot_title = title2,
    sig_lines = sig_lines, label_threshold_colors = label_threshold_colors, y_axis_squish_threshold = y_axis_squish_threshold,
    gene_label_size = gene_label_size, gene_col = gene_col, group_col = group_col, loci_to_label = loci_to_label,
    genes_to_label = genes_to_label, custom_gene_colors = custom_gene_colors, inverted = TRUE, font_family = font_family,
    label_orientation = label_orientation # Pass parameter
  )

  # --- Combine and Save Plots ---
  combined_plot <- top_plot / bottom_plot + patchwork::plot_layout(heights = c(1, 1))

  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

  file_path_png <- file.path(output_folder, paste0(file_name_prefix, ".png"))
  file_path_pdf <- file.path(output_folder, paste0(file_name_prefix, ".pdf"))

  tryCatch({
    grDevices::png(filename = file_path_png, width = 18, height = 12, units = "in", res = 300, type = "cairo", bg = "white")
    print(combined_plot)
    grDevices::dev.off()
    grDevices::cairo_pdf(filename = file_path_pdf, width = 18, height = 12)
    print(combined_plot)
    grDevices::dev.off()
  }, error = function(e) {
    message("Cairo devices failed, falling back to standard ggsave devices.")
    ggplot2::ggsave(file_path_png, plot = combined_plot, device = ragg::agg_png, width = 18, height = 12, units = "in", res = 300, bg = "white")
    ggplot2::ggsave(file_path_pdf, plot = combined_plot, device = "pdf", width = 18, height = 12, units = "in")
  })

  message(paste("Successfully saved plots to:", file_path_png, "and", file_path_pdf))
  return(invisible(combined_plot))
}
