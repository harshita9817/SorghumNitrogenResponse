peak<- read.csv("totalcombined_snps_summary_HN_sorghum.csv")
head(peak)
library(ggplot2)
library(dplyr)
library(viridis)  # For viridis color palette
library(gtools)  # Load the package
library(patchwork)
# Ensure CHROM and chr are ordered correctly
peak <- peak %>%
  mutate(CHROM = factor(CHROM, levels = mixedsort(unique(CHROM))),  # Order chromosomes
         chr = factor(chr, levels = mixedsort(unique(chr))))  # Order gene chromosomes

# Compute cumulative positions for CHROM (X-axis)
chrom_cum <- peak %>%
  group_by(CHROM) %>%
  summarise(max_bp = max(POS), .groups = "drop") %>%
  mutate(bp_add_x = lag(cumsum(max_bp), default = 0)) %>%
  select(CHROM, bp_add_x)

# Compute cumulative positions for chr (Y-axis)
peak_filtered <- peak %>% filter(!is.na(start))

chr_cum <- peak_filtered %>%
  group_by(chr) %>%
  summarise(max_bp = max(start), .groups = "drop") %>%
  mutate(bp_add_y = lag(cumsum(max_bp), default = 0))

# Merge cumulative positions into peak dataset
peak_cum <- peak %>%
  left_join(chrom_cum, by = "CHROM") %>%
  mutate(bp_cum_x = POS + bp_add_x) %>%
  left_join(chr_cum, by = "chr") %>%
  mutate(bp_cum_y = start + bp_add_y)

library(ggplot2)
library(dplyr)

# Compute chromosome midpoints for axis labels
axis_x_set <- peak_cum %>%
  group_by(CHROM) %>%
  summarize(center = mean(bp_cum_x))

axis_y_set <- peak_cum %>%
  group_by(chr) %>%
  summarize(center = mean(bp_cum_y))

# Manually define chromosome colors: Odd → Purple, Even → Green
chrom_colors <- setNames(ifelse(as.numeric(gsub("Chr", "", levels(peak$CHROM))) %% 2 == 1, 
                                "#440154FF",  # Extreme purple for odd
                                "#22A884FF"), # Middle green for even
                         levels(peak$CHROM))

# Plot
ggplot(peak_cum, aes(x = bp_cum_x, y = bp_cum_y, color = CHROM)) + 
  geom_point(alpha = 0.8, size = 3) +
  scale_color_manual(values = chrom_colors) +  # Custom color mapping
  scale_x_continuous(breaks = axis_x_set$center, labels = axis_x_set$CHROM) +  # Add chromosome labels on x-axis
  scale_y_continuous(breaks = axis_y_set$center, labels = axis_y_set$chr) +  # Add chromosome labels on y-axis
  labs(title = "Peak Positions vs. Gene Positions",
       x = "Chromosome and Peak Position",
       y = "Chromosome and Gene Position") +
  theme_minimal(base_size = 14) +  # Cleaner theme
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    panel.background = element_rect(fill = "white", color = NA),  # Plain white background
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
    axis.text.y = element_text(size = 12),  # Adjust y-axis label size
    legend.position = "none"  # Remove legend
  )
# For combining plots

# Compute chromosome midpoints for axis labels
axis_x_set <- peak_cum %>%
  group_by(CHROM) %>%
  summarize(center = mean(bp_cum_x))

axis_y_set <- peak_cum %>%
  group_by(chr) %>%
  summarize(center = mean(bp_cum_y))

# Manually define chromosome colors: Odd → Purple, Even → Green
chrom_colors <- setNames(ifelse(as.numeric(gsub("Chr", "", levels(peak$CHROM))) %% 2 == 1, 
                                "#440154FF",  # Extreme purple for odd
                                "#22A884FF"), # Middle green for even
                         levels(peak$CHROM))

# Scatterplot: Peak Positions vs. Gene Positions
scatter_plot <- ggplot(peak_cum, aes(x = bp_cum_x, y = bp_cum_y, color = CHROM)) + 
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = chrom_colors) +  # Custom color mapping
  scale_x_continuous(breaks = axis_x_set$center, labels = axis_x_set$CHROM) +  # Add chromosome labels on x-axis
  scale_y_continuous(breaks = axis_y_set$center, labels = axis_y_set$chr) +  # Add chromosome labels on y-axis
  labs(title = "Peak Positions vs. Gene Positions in High Nitrogen",
       x = "Chromosome and Peak Position",
       y = "Chromosome and Gene Position") +
  theme_minimal(base_size = 14) +  # Cleaner theme
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    panel.background = element_rect(fill = "white", color = NA),  # Plain white background
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
    axis.text.y = element_text(size = 12),  # Adjust y-axis label size
    legend.position = "none"  # Remove legend
  )

peak_LN <- read.csv("totalcombined_snps_summary_LN_sorghum.csv")

# Repeat same preprocessing steps for LN data
peak_LN <- peak_LN %>%
  mutate(CHROM = factor(CHROM, levels = mixedsort(unique(CHROM))),
         chr = factor(chr, levels = mixedsort(unique(chr))))

chrom_cum_LN <- peak_LN %>%
  group_by(CHROM) %>%
  summarise(max_bp = max(POS), .groups = "drop") %>%
  mutate(bp_add_x = lag(cumsum(max_bp), default = 0)) %>%
  select(CHROM, bp_add_x)

peak_filtered_LN <- peak_LN %>% filter(!is.na(start))

chr_cum_LN <- peak_filtered_LN %>%
  group_by(chr) %>%
  summarise(max_bp = max(start), .groups = "drop") %>%
  mutate(bp_add_y = lag(cumsum(max_bp), default = 0))

peak_cum_LN <- peak_LN %>%
  left_join(chrom_cum_LN, by = "CHROM") %>%
  mutate(bp_cum_x = POS + bp_add_x) %>%
  left_join(chr_cum_LN, by = "chr") %>%
  mutate(bp_cum_y = start + bp_add_y)

axis_x_set_LN <- peak_cum_LN %>%
  group_by(CHROM) %>%
  summarize(center = mean(bp_cum_x))

axis_y_set_LN <- peak_cum_LN %>%
  group_by(chr) %>%
  summarize(center = mean(bp_cum_y))

chrom_colors_LN <- chrom_colors  # Use same color scheme as HN

# Scatterplot for LN
# Adjusted scatter plot for HN
scatter_plot <- ggplot(peak_cum, aes(x = bp_cum_x, y = bp_cum_y, color = CHROM)) + 
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = chrom_colors) +
  scale_x_continuous(breaks = axis_x_set$center, labels = axis_x_set$CHROM) +
  scale_y_continuous(breaks = axis_y_set$center, labels = axis_y_set$chr) +
  labs(title = "High Nitrogen", x = "Peak Position", y = "Gene Position") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 0, size = 9),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 11, face = "bold"),
    legend.position = "none"
  )

# Adjusted scatter plot for LN
scatter_plot_LN <- ggplot(peak_cum_LN, aes(x = bp_cum_x, y = bp_cum_y, color = CHROM)) + 
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = chrom_colors_LN) +
  scale_x_continuous(breaks = axis_x_set_LN$center, labels = axis_x_set_LN$CHROM) +
  scale_y_continuous(breaks = axis_y_set_LN$center, labels = axis_y_set_LN$chr) +
  labs(title = "Low Nitrogen", x = "Peak Position", y = "Gene Position") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 0, size = 9),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 11, face = "bold"),
    legend.position = "none"
  )


# Combined histogram plot
peak_cum_HN_hist <- peak_cum %>% mutate(condition = "High N")
peak_cum_LN_hist <- peak_cum_LN %>% mutate(condition = "Low N")
combined_hist_data <- bind_rows(peak_cum_HN_hist, peak_cum_LN_hist)

histogram_combined <- ggplot(combined_hist_data, aes(x = bp_cum_x, fill = condition)) +
  geom_histogram(binwidth = 1e6, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = c("High N" = "purple", "Low N" = "darkgreen")) +
  scale_x_continuous(breaks = axis_x_set$center, labels = axis_x_set$CHROM) +
  labs(title = "Peak Density Across Genome", x = "Chromosome and Peak Position", y = "Count") +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 0, hjust = 1),
    axis.text.y = element_text(size = 12)
  )

# Combine all three panels using patchwork
library(patchwork)

# Stack vertically: HN → LN → Histogram
final_combined_plot <- scatter_plot / scatter_plot_LN / histogram_combined +
  plot_layout(heights = c(1, 1, 1))  # Histogram slightly shorter

# Show the combined plot
print(final_combined_plot)
