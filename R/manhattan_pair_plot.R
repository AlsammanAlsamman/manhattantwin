#' Create Inverted Manhattan Plot Pair from Data Objects
#'
#' Creates a pair of Manhattan plots with one inverted and aligned for comparison.
#' Takes two data objects (data.frames or data.tables) as input.
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
#' @param gwas_data1 First GWAS dataset (data.frame or data.table)
#' @param gwas_data2 Second GWAS dataset (data.frame or data.table)
#' @param chr_col Name of chromosome column (default: "CHR")
#' @param bp_col Name of base pair position column (default: "BP")
#' @param p_col Name of p-value column (default: "P")
#' @param gene_col Name of gene column (default: "GENE")
#' @param group_col Name of grouping column (default: NULL)
#' @param loci_to_label Optional vector of locus names to label
#' @param genes_to_label Optional vector of gene names to label
#' @param custom_gene_colors Optional named vector of custom colors for genes
#' @param n_cases1 Number of cases in first dataset
#' @param n_controls1 Number of controls in first dataset
#' @param n_cases2 Number of cases in second dataset
#' @param n_controls2 Number of controls in second dataset
#' @param plot_title1 Title for first plot (default: "Dataset 1")
#' @param plot_title2 Title for second plot (default: "Dataset 2")
#' @param lambda1 Genomic inflation factor for first dataset (default: NULL)
#' @param lambda2 Genomic inflation factor for second dataset (default: NULL)
#' @param total_snps_in_study1 Total SNPs in first study (default: NULL)
#' @param total_snps_in_study2 Total SNPs in second study (default: NULL)
#' @param add_date_to_title Logical to add date to title (default: FALSE)
#' @param plot_pval_threshold Maximum p-value to plot (default: 1)
#' @param sig_lines Named vector of significance lines (default: c("red" = 5e-8, "blue" = 1e-5))
#' @param label_threshold_colors Named vector of label thresholds (default: c("red" = 5e-8, "orange" = 1e-6, "darkblue" = 1e-5))
#' @param y_axis_squish_threshold Y-axis squish threshold (default: 30)
#' @param gene_label_size Size of gene labels (default: 3.5)
#' @param output_folder Output folder (default: "Inverted_Manhattan_Plots")
#' @param file_name_prefix File name prefix (default: "inverted_manhattan")
#' @param font_family Font family for Unicode support (default: "Arial Unicode MS")
#' @export
#' @return Invisibly returns the combined plot object
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' data("gwasdataseta")
#' data1 <- gwasdataseta
#' data2 <- data.table::copy(gwasdataseta)
#' manhattan_pair_plot(
#'   gwas_data1 = data1,
#'   gwas_data2 = data2,
#'   n_cases1 = 10000, n_controls1 = 10000,
#'   n_cases2 = 8000, n_controls2 = 8000,
#'   plot_title1 = "Study 1", plot_title2 = "Study 2"
#' )
#' }
manhattan_pair_plot <- function(gwas_data1, gwas_data2,
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

  # Convert input data to data.table
  data1 <- data.table::as.data.table(gwas_data1)
  data2 <- data.table::as.data.table(gwas_data2)

  # Process both datasets using the same chromosome ordering
  all_chrs <- unique(c(data1[[chr_col]], data2[[chr_col]]))
  all_chrs <- sort(unique(all_chrs))

  # Process dataset 1
  processed_data1 <- preprocess_manhattan_data(
    gwas_data = data1,
    chr_col = chr_col,
    bp_col = bp_col,
    p_col = p_col,
    gene_col = gene_col,
    group_col = group_col,
    plot_pval_threshold = plot_pval_threshold,
    all_chrs = all_chrs
  )

  # Process dataset 2
  processed_data2 <- preprocess_manhattan_data(
    gwas_data = data2,
    chr_col = chr_col,
    bp_col = bp_col,
    p_col = p_col,
    gene_col = gene_col,
    group_col = group_col,
    plot_pval_threshold = plot_pval_threshold,
    all_chrs = all_chrs
  )

  # Calculate lambda values if not provided
  if (is.null(lambda1)) {
    chisq1 <- stats::qchisq(1 - processed_data1$plot_data$pvalue, 1)
    lambda1 <- round(median(chisq1, na.rm = TRUE) / stats::qchisq(0.5, 1), 4)
  }

  if (is.null(lambda2)) {
    chisq2 <- stats::qchisq(1 - processed_data2$plot_data$pvalue, 1)
    lambda2 <- round(median(chisq2, na.rm = TRUE) / stats::qchisq(0.5, 1), 4)
  }

  # Generate plot titles with lambda symbol
  title1 <- generate_plot_title(
    plot_title_prefix = plot_title1,
    n_cases = n_cases1,
    n_controls = n_controls1,
    lambda = lambda1,
    total_snps_in_study = total_snps_in_study1,
    add_date_to_title = add_date_to_title,
    gwas_data = data1,
    font_family = font_family
  )

  title2 <- generate_plot_title(
    plot_title_prefix = plot_title2,
    n_cases = n_cases2,
    n_controls = n_controls2,
    lambda = lambda2,
    total_snps_in_study = total_snps_in_study2,
    add_date_to_title = add_date_to_title,
    gwas_data = data2,
    font_family = font_family
  )

  # Create the top plot (normal orientation)
  top_plot <- create_single_manhattan(
    plot_data = processed_data1$plot_data,
    axis_data = processed_data1$axis_data,
    plot_title = title1,
    sig_lines = sig_lines,
    label_threshold_colors = label_threshold_colors,
    y_axis_squish_threshold = y_axis_squish_threshold,
    gene_label_size = gene_label_size,
    gene_col = gene_col,
    group_col = group_col,
    loci_to_label = loci_to_label,
    genes_to_label = genes_to_label,
    custom_gene_colors = custom_gene_colors,
    inverted = FALSE,
    font_family = font_family
  )

  # Create the bottom plot (inverted orientation)
  bottom_plot <- create_single_manhattan(
    plot_data = processed_data2$plot_data,
    axis_data = processed_data2$axis_data,
    plot_title = title2,
    sig_lines = sig_lines,
    label_threshold_colors = label_threshold_colors,
    y_axis_squish_threshold = y_axis_squish_threshold,
    gene_label_size = gene_label_size,
    gene_col = gene_col,
    group_col = group_col,
    loci_to_label = loci_to_label,
    genes_to_label = genes_to_label,
    custom_gene_colors = custom_gene_colors,
    inverted = TRUE,
    font_family = font_family
  )

  # Combine plots
  combined_plot <- top_plot / bottom_plot +
    patchwork::plot_layout(heights = c(1, 1))

  # Save plots
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }

  file_path_png <- file.path(output_folder, paste0(file_name_prefix, ".png"))
  file_path_pdf <- file.path(output_folder, paste0(file_name_prefix, ".pdf"))

  # Load fonts if available
  if (requireNamespace("extrafont", quietly = TRUE)) {
    tryCatch({
      extrafont::loadfonts()
    }, error = function(e) {
      message("Could not load fonts. Using default fonts.")
    })
  }

  # Save with Cairo devices for better Unicode support
  tryCatch({
    # Save PNG with Cairo
    grDevices::png(filename = file_path_png, width = 18, height = 12, units = "in", res = 300, type = "cairo", bg = "white")
    print(combined_plot)
    grDevices::dev.off()

    # Save PDF with Cairo
    grDevices::cairo_pdf(filename = file_path_pdf, width = 18, height = 12)
    print(combined_plot)
    grDevices::dev.off()
  }, error = function(e) {
    message("Cairo devices failed, falling back to standard devices")
    # Fallback to ragg and standard pdf
    ggplot2::ggsave(file_path_png, plot = combined_plot, device = ragg::agg_png, width = 18, height = 12, units = "in", res = 300, bg = "white")
    ggplot2::ggsave(file_path_pdf, plot = combined_plot, device = "pdf", width = 18, height = 12, units = "in")
  })

  message(paste("Successfully saved plots to:", file_path_png, "and", file_path_pdf))
  return(invisible(combined_plot))
}

#' Helper function to preprocess GWAS data
#'
#' @param gwas_data GWAS data
#' @param chr_col Chromosome column name
#' @param bp_col Base pair column name
#' @param p_col P-value column name
#' @param gene_col Gene column name
#' @param group_col Group column name
#' @param plot_pval_threshold P-value threshold
#' @param all_chrs All chromosomes to include
#'
#' @return List with processed plot data and axis data

preprocess_manhattan_data <- function(gwas_data, chr_col, bp_col, p_col, gene_col,
                                      group_col, plot_pval_threshold, all_chrs) {

  # Filter data
  plot_data <- gwas_data[!is.na(get(p_col)) & get(p_col) > 0 & get(p_col) <= plot_pval_threshold]

  # Process data
  plot_data[, `:=` (pvalue = as.numeric(get(p_col)),
                    CHR = as.character(get(chr_col)),
                    BP = as.numeric(get(bp_col)))]
  plot_data[, log_p := -log10(pvalue)]

  # Prepare chromosome numbering
  plot_data[, CHR_num := as.numeric(factor(CHR, levels = unique(stringr::str_sort(unique(CHR), numeric = TRUE))))]
  plot_data <- plot_data[order(CHR_num, BP)]

  # Calculate cumulative positions
  cumulative_summary <- plot_data[, .(max_bp = max(BP)), by = CHR_num]
  cumulative_summary <- cumulative_summary[order(CHR_num)]
  cumulative_summary$bp_add <- c(0, cumsum(cumulative_summary$max_bp[-nrow(cumulative_summary)]))

  plot_data <- merge(plot_data, cumulative_summary[, .(CHR_num, bp_add)], by = "CHR_num")
  plot_data[, bp_cum := BP + bp_add]
  axis_data <- plot_data[, .(center = (max(bp_cum) + min(bp_cum)) / 2), by = CHR]

  return(list(plot_data = plot_data, axis_data = axis_data))
}

#' Helper function to create a single Manhattan plot
#'
#' @param plot_data Processed GWAS data
#' @param axis_data Axis data for plotting
#' @param plot_title Plot title
#' @param sig_lines Significance lines
#' @param label_threshold_colors Label threshold colors
#' @param y_axis_squish_threshold Y-axis squish threshold
#' @param gene_label_size Gene label size
#' @param gene_col Gene column name
#' @param group_col Group column name
#' @param loci_to_label Loci to label
#' @param genes_to_label Genes to label
#' @param custom_gene_colors Custom gene colors
#' @param inverted Whether to invert the plot
#' @param font_family Font family for Unicode support
#'
#' @return A ggplot object

create_single_manhattan <- function(plot_data, axis_data, plot_title, sig_lines,
                                    label_threshold_colors, y_axis_squish_threshold,
                                    gene_label_size, gene_col, group_col,
                                    loci_to_label, genes_to_label, custom_gene_colors,
                                    inverted, font_family = "Arial Unicode MS") {

  # Prepare gene labels
  snps_to_label <- data.table::data.table()
  if (!is.null(gene_col) && gene_col %in% names(plot_data)) {
    if (is.null(group_col)) group_col <- gene_col
    if (group_col %in% names(plot_data)) {
      snps_for_labeling <- plot_data[!is.na(get(group_col))]
      if (!is.null(loci_to_label)) snps_for_labeling <- snps_for_labeling[get(group_col) %in% loci_to_label]
      else if (!is.null(genes_to_label)) snps_for_labeling <- snps_for_labeling[get(gene_col) %in% genes_to_label]
      else {
        min_label_thresh <- max(label_threshold_colors, na.rm = TRUE)
        snps_for_labeling <- snps_for_labeling[pvalue < min_label_thresh]
      }
      if (nrow(snps_for_labeling) > 0) snps_to_label <- snps_for_labeling[, .SD[which.min(pvalue)], by = group_col]
      if (nrow(snps_to_label) > 0) {
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

  # Define transformation functions
  compression_factor <- 0.1
  squish_forward <- function(x) ifelse(x <= y_axis_squish_threshold, x, y_axis_squish_threshold + (x - y_axis_squish_threshold) * compression_factor)
  squish_inverse <- function(x) ifelse(x <= y_axis_squish_threshold, x, (x - y_axis_squish_threshold) / compression_factor + y_axis_squish_threshold)
  squish_trans <- scales::trans_new("squish", squish_forward, squish_inverse)
  x_axis_expansion <- ggplot2::expansion(mult = c(0.015, 0.015))

  # Define custom chromosome colors
  custom_chromosome_colors <- c(
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
    "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78"
  )

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
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor.x = ggplot2::element_blank(),
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
        ggrepel::geom_text_repel(
          data = inverted_snps_to_label,
          ggplot2::aes(x = bp_cum, y = transformed_log_p, label = label_text, color = label_color),
          size = gene_label_size, fontface = "bold.italic", nudge_y = -10,
          direction = "x", angle = 90, hjust = 0, segment.color = 'grey50',
          segment.linetype = "dashed", min.segment.length = 0, max.overlaps = Inf
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
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(color = "black"),
        axis.title.y = ggplot2::element_text(color = "black", face = "bold"),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(t = 5.5, r = 5.5, b = -15, l = 5.5, unit = "pt")
      )

    if (nrow(snps_to_label) > 0) {
      p <- p +
        ggnewscale::new_scale_color() +
        ggrepel::geom_text_repel(
          data = snps_to_label,
          ggplot2::aes(x = bp_cum, y = log_p, label = label_text, color = label_color),
          size = gene_label_size, fontface = "bold.italic", nudge_y = 15,
          direction = "x", angle = 90, hjust = 0, segment.color = 'grey50',
          segment.linetype = "dashed", min.segment.length = 0, max.overlaps = Inf
        ) + ggplot2::scale_color_identity(guide = "none")
    }
  }

  return(p)
}

#' Helper function to generate plot titles with lambda symbol
#'
#' @param plot_title_prefix Prefix for the plot title
#' @param n_cases Number of cases
#' @param n_controls Number of controls
#' @param lambda Genomic inflation factor
#' @param total_snps_in_study Total SNPs in study
#' @param add_date_to_title Whether to add date to title
#' @param gwas_data GWAS data
#' @param font_family Font family for Unicode support
#'
#' @return Formatted plot title string

generate_plot_title <- function(plot_title_prefix, n_cases, n_controls, lambda,
                                total_snps_in_study, add_date_to_title, gwas_data,
                                font_family = "Arial Unicode MS") {

  if (is.null(total_snps_in_study)) {
    snps_for_title <- format(nrow(gwas_data), big.mark = ",", scientific = FALSE)
    snps_title_label <- "Total SNPs"
  } else {
    snps_for_title <- format(total_snps_in_study, big.mark = ",", scientific = FALSE)
    snps_title_label <- "Total SNPs in Study"
  }

  title_string <- sprintf("%s\n(%s: %s | λ = %s | Cases: %s, Controls: %s)",
                         plot_title_prefix, snps_title_label, snps_for_title,
                         lambda, n_cases, n_controls)

  if (add_date_to_title) {
    title_string <- paste0(title_string, " | ", Sys.Date())
  }

  return(title_string)
}

