library(dplyr)
library(ggplot2)
library(viridis)

# Step 1: List all files starting with "Z"
file_names <- list.files(pattern = "^Z.*\\.csv$")

# Step 2: Create labels for traits (remove Z, _signals, and .csv)
trait_labels <- gsub("^Z|_signals|\\.csv", "", file_names)

# Step 3: Safe named list for reading and labeling
files <- setNames(file_names, trait_labels)
df_list <- lapply(names(files), function(label) {
  df <- read.csv(files[[label]])
  df$Trait <- label
  df
})

# Step 4: Combine into one dataframe
df_all <- bind_rows(df_list)

# Step 5: Determine nitrogen level and mirror LN support values
df_all <- df_all %>%
  mutate(
    CHROM = as.character(CHROM),
    Nitrogen = ifelse(grepl("^LN", Trait), "LN", "HN"),
    support = ifelse(Nitrogen == "LN", -support, support)
  )

# Step 6: Compute cumulative position for chromosomes
data_cum <- df_all %>%
  group_by(CHROM) %>%
  summarise(max_bp = max(POS), .groups = "drop") %>%
  mutate(bp_add = lag(cumsum(max_bp), default = 0))

gwas_data <- df_all %>%
  left_join(data_cum, by = "CHROM") %>%
  mutate(bp_cum = POS + bp_add)

# Step 7: Get axis labels for chromosomes
axis_set <- gwas_data %>%
  group_by(CHROM) %>%
  summarize(center = mean(bp_cum), .groups = "drop")
annotate("text", x = max(gwas_data$bp_cum) * 0.95, y = 0.55, label = "HN", size = 6, fontface = "bold") +
annotate("text", x = max(gwas_data$bp_cum) * 0.95, y = -0.55, label = "LN", size = 6, fontface = "bold")
geom_vline(xintercept = data_cum$bp_add, color = "grey80", linetype = "dashed", lwd = 0.2) +
  
  ggplot(gwas_data, aes(x = bp_cum, y = support, color = TraitGroup)) +
  geom_point(alpha = 0.8, size = 3) +
  geom_hline(yintercept = 0, color = "black", linetype = "solid") +
  geom_hline(yintercept = c(0.2, -0.2), color = "grey50", linetype = "dashed") +
  annotate("text", x = max(gwas_data$bp_cum) * 0.95, y = 0.55, label = "HN", size = 6, fontface = "bold") +
  annotate("text", x = max(gwas_data$bp_cum) * 0.95, y = -0.55, label = "LN", size = 6, fontface = "bold") +
  scale_color_manual(values = c(
    "Plant Height" = "#E64B35FF",
    "Flowering" = "#3C5488FF",
    "Chlorophyll" = "#00A087FF",
    "Grain Weight" = "#FFB000",
    "Other" = "grey70"
  )) +
  scale_x_continuous(breaks = axis_set$center, labels = axis_set$CHROM) +
  scale_y_continuous(
    breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6),
    labels = abs(c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6))
  ) +
  labs(x = "Chromosome", y = "Support", color = "Trait Group") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1),
    panel.grid = element_blank()
  )
