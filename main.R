library(manhattantwin)
data("gwasdataseta")

# Define metadata for two example cohorts
cohortA <- list(
  name = "Simulated European Cohort",
  n_cases = 1200,
  n_controls = 2300,
  note = "Example: European ancestry"
)
cohortB <- list(
  name = "Simulated Admixed Cohort",
  n_cases = 900,
  n_controls = 2100,
  note = "Example: Admixed ancestry"
)

# Cluster SNPs for each cohort with different p-value thresholds
gwasA <- manhattantwin::cluster_snps(
  gwasdataseta,
  chr_col = "chr",
  pos_col = "pos",
  pvalue_col = "pvalue",
  rsid_col = "rsid",
  pvalue_threshold = 5e-8,
  distance_threshold = 250000
)
gwasB <- manhattantwin::cluster_snps(
  gwasdataseta,
  chr_col = "chr",
  pos_col = "pos",
  pvalue_col = "pvalue",
  rsid_col = "rsid",
  pvalue_threshold = 1e-6,
  distance_threshold = 250000
)

# Custom gene highlight colors
highlight_genes <- c(
  "Gene1" = "red",
  "Gene141" = "blue",
  "Gene393" = "darkgreen",
  "Gene81" = "purple",
  "Gene505" = "orange",
  "Gene56" = "brown",
  "Gene670" = "cyan",
  "Gene187" = "magenta",
  "Gene219" = "gold",
  "Gene210" = "darkred",
  "Gene196" = "navy",
  "Gene344" = "forestgreen",
  "Gene690" = "darkorange",
  "Gene26" = "deepskyblue"
)

# Plot single Manhattan for cohort A
manhattantwin::plot_single_manhattan(
  gwasA,
  plot_title_prefix = cohortA$name,
  p_col = "pvalue",
  n_cases = cohortA$n_cases,
  n_controls = cohortA$n_controls,
  file_name_prefix = "sim_eur",
  group_col = "cluster",
  gene_col = "gene",
  output_folder = "single_plots",
  y_axis_squish_threshold = 30,
  label_threshold_colors = c("red" = 5e-8, "orange" = 1e-6, "darkblue" = 1e-5),
  custom_gene_colors = highlight_genes
)

# Plot mirrored Manhattan for comparison
create_mirrored_manhattan_plot(
  gwasA,
  gwasB,
  chr_col = "chr",
  bp_col = "pos",
  p_col = "pvalue",
  n_cases1 = cohortA$n_cases,
  n_controls1 = cohortA$n_controls,
  n_cases2 = cohortB$n_cases,
  n_controls2 = cohortB$n_controls,
  file_name_prefix = "sim_eur_vs_admixed",
  group_col = "cluster",
  gene_col = "gene",
  label_threshold_colors = c("red" = 5e-8, "orange" = 1e-6, "darkblue" = 1e-5),
  output_folder = "Inverted_Manhattan_Plots",
  y_axis_squish_threshold = 30,
  plot_title1 = cohortA$name,
  plot_title2 = cohortB$name,
  custom_gene_colors = highlight_genes,label_alpha = 0.7,label_orientation = "vertical",output_width = 15
)



#################### Building Vignette ####################
# All-in-one command (most common workflow)
devtools::document(); devtools::install(); renv::snapshot()
