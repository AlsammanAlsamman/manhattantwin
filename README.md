# manhattantwin

**manhattantwin** is an R package for generating publication-quality Manhattan plots and inverted Manhattan pair plots for GWAS data. It provides advanced features for SNP clustering, gene labeling, and customizable plot outputs.

## Features


## Features

- **SNP clustering:** Automatically clusters SNPs by genomic distance and p-value, selecting the most significant SNP (lowest p-value) in each cluster and linking it to the associated gene. This helps highlight the best candidate SNP for each locus.
- **Region expansion for highly significant SNPs:** The package can expand and visually emphasize regions associated with extremely significant SNPs (e.g., those with p-values less than 5e-10), making it easier to interpret and present key findings. This is controlled by the `y_axis_squish_threshold` parameter (e.g., `y_axis_squish_threshold = 10`).
- Generate single Manhattan plots with gene/group labeling
- Create inverted Manhattan pair plots for comparative studies
- Publication-ready PNG and PDF outputs
- **Customizable gene label colors:** Use different colors for each gene threshold, making it easy to visually distinguish significance levels on your plots.

## Installation

You can install the latest version of **manhattantwin** directly from GitHub using the `devtools` package:

```r
# Install devtools if you don't have it
install.packages("devtools")

# Install manhattantwin from GitHub
devtools::install_github("AlsammanAlsamman/manhattantwin")
```

## Example Usage

### 1. Single Manhattan Plot

```r
library(manhattantwin)
data("gwasdataseta")

# Cluster SNPs
gwasdataseta <- manhattantwin::cluster_snps(
  gwasdataseta,
  chr_col = "chr",
  pos_col = "pos",
  pvalue_col = "pvalue",
  rsid_col = "rsid",
  pvalue_threshold = 5e-8,
  distance_threshold = 250000
)

 # Generate single Manhattan plot
 manhattantwin::plot_single_manhattan(
   gwasdataseta,
   plot_title_prefix = "Example Dataset A",
   p_col = "pvalue",
   n_cases = 1000,
   n_controls = 1000,
   file_name_prefix = "example_a",
   group_col = "cluster",
   gene_col = "gene",
   output_folder = "single_plots",
   # You can customize label colors for different gene thresholds:
   label_threshold_colors = c("red" = 5e-8, "orange" = 1e-6, "darkblue" = 1e-5),
   y_axis_squish_threshold = 10 # expandin snps less than 10 -logpvalue
 )
```

#### Example Output

![Single Manhattan Plot](single_plots/example_a_pub.png)


---

### 2. Inverted Manhattan Pair Plot

```r
# Cluster a second dataset (for demonstration, using the same data)
gwasdatasetb <- manhattantwin::cluster_snps(
  gwasdataseta,
  chr_col = "chr",
  pos_col = "pos",
  pvalue_col = "pvalue",
  rsid_col = "rsid",
  pvalue_threshold = 5e-8,
  distance_threshold = 250000
)

# Generate inverted Manhattan pair plot
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
  output_folder = "pair_plots",
  label_threshold_colors = c("red" = 5e-8, "orange" = 1e-6, "darkblue" = 1e-5),
  y_axis_squish_threshold = 10 # expandin snps less than 10 -logpvalue
)
```

![Inverted Manhattan Plots](Inverted_Manhattan_Plots/example_a_vs_b.png)

*In the plot above, different colors are used for gene labels based on their significance thresholds (e.g., red for p < 5e-8, orange for p < 1e-6, dark blue for p < 1e-5).*

## Data

The package includes an example GWAS dataset `gwasdataseta` for demonstration purposes.

## License

MIT License. See `LICENSE` file for details.

## Authors

See package DESCRIPTION for author information.
