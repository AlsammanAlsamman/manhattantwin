
library(manhattantwin)
data("gwasdataseta")

gwasdataseta <- manhattantwin::cluster_snps(gwasdataseta,chr_col = "chr", pos_col = "pos", pvalue_col = "pvalue", rsid_col = "rsid", pvalue_threshold = 5e-8, distance_threshold = 250000)

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

gwasdatasetb <- manhattantwin::cluster_snps(gwasdataseta,chr_col = "chr", pos_col = "pos", pvalue_col = "pvalue", rsid_col = "rsid", pvalue_threshold = 5e-5, distance_threshold = 250000)




manhattantwin::manhattan_pair_plot(
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





