### BUSCO R plot

### Input: a TSV table with one row per sample, with at least the columns:
# 1. "sample_id" > like A1, A2, A3, ... 
# 2. "species" > like Mus musculus, Homo sapiens, ...
# 3. "BUSCO_M10" > or name this column as you want, but it must contain the BUSCO summary string (e.g., C:47.6%[S:37.4%,D:10.2%],F:20.4%,M:32.0%,n:954)

# Example input:
# sample_id  species  BUSCO_M10
# A1  Anaspides richardsoni  C:80.2%[S:51.3%,D:28.9%],F:9.4%,M:10.4%,n:954
# A2  Anaspides swaini  C:90.0%[S:50.9%,D:39.1%],F:5.9%,M:4.1%,n:954

### Output: a stacked horizontal bar plot of single-copy (S), duplicated (D), fragmented (F) and missing (M) BUSCOs for the chosen lineage column, with species sorted by completeness (or by a user-defined order).

###########################################################################
# load libs
library(ggplot2)

### busco results
# load data
busco_data <- read.table("busco_result.tsv", header = T, sep = "\t")
# header: sample_id,source,fastq_r1,fastq_r2,species,date,location,cave_adaptation,binary,num_seqs,sum_len

# define the column for BUSCO completeness. In some cases, it is useful to run BUSCO with different lineage datasets, so the column name may vary. Here, we assume the column is named "BUSCO_M10" for the M10 lineage dataset, but edit the column name as needed.
busco_col <- "BUSCO_M10"
x <- busco_data[[busco_col]]
busco_data$C <- as.numeric(sub(".*C:([0-9.]+)%.*", "\\1", x))
busco_data$S <- as.numeric(sub(".*S:([0-9.]+)%.*", "\\1", x))
busco_data$D <- as.numeric(sub(".*D:([0-9.]+)%.*", "\\1", x))
busco_data$F <- as.numeric(sub(".*F:([0-9.]+)%.*", "\\1", x))
busco_data$M <- as.numeric(sub(".*M:([0-9.]+)%.*", "\\1", x))

# create a long format data frame for plotting
bars_df <- stack(busco_data[, c("S", "D", "F", "M")])
bars_df$species <- rep(busco_data$species, 4)

# chose your colors for the barplot
# see here for ideas: https://r-graph-gallery.com/color-palette-finder
busco_colors <- c(S = "steelblue4", D = "cadetblue", F = "tan1", M = "lightgrey")

# if you want to specify the order of species in the plot, create a text file with one species per line, for the top to the bottom of the plot. If you don't provide a file, species will be sorted by completeness (S + D), with the most complete on top.
order_file <- "" # write here the name of the file with the species order, or leave empty to sort by completeness
if (file.exists(order_file)) {
  species_order <- trimws(readLines(order_file))
  species_order <- species_order[species_order != ""]
  } else 
  {species_order <- busco_data$species[order(busco_data$C, decreasing = TRUE)]}
bars_df$species <- factor(bars_df$species, levels = rev(species_order))

# plot
bar_totals <- data.frame(species = factor(busco_data$species, levels = rev(species_order)), total   = busco_data$S + busco_data$D + busco_data$F + busco_data$M)
bar_width <- 0.8

p_bars <- ggplot(bars_df, aes(x = values, y = species, fill = ind)) +
  geom_col(position = position_stack(reverse = TRUE), width = bar_width) +
  geom_vline(xintercept = c(25, 50, 75), color = "white", linewidth = 0.5) +
  geom_col(data = bar_totals, aes(x = total, y = species), inherit.aes = FALSE,
          fill = NA, color = "black", width = bar_width, linewidth = 0.4) +
  scale_fill_manual(values = busco_colors) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100), expand = expansion(mult = c(0, 0.01))) +
  labs(title = busco_col, x = "% BUSCOs", y = NULL) +
  theme_classic() +
  theme(
    text = element_text(family = "Helvetica Neue", size = 12),
    axis.text  = element_text(size = 12, color = "black"),
    axis.line  = element_line(color = "black", linewidth = 0.5),
    axis.line.y  = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_blank())

print(p_bars)

# save, edit width and height as you like
plot_width  <- 6
plot_height <- 1 + 0.2 * nrow(busco_data) # grows with the number of species
ggsave(paste0("buscoBars_", busco_col, ".png"), p_bars, width = plot_width, height = plot_height, dpi = 300)
ggsave(paste0("buscoBars_", busco_col, ".svg"), p_bars, width = plot_width, height = plot_height)
