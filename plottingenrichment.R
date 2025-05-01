# Load required libraries
library(ggplot2)
library(dplyr)

# Read the data
up <- read.csv("goa_output_down.csv", sep = ",", header = TRUE)

# Filter significant GO terms (FDR < 0.1, optional threshold)
up_sig <- up %>%
  filter(p_fdr_bh < 0.05)

# Prepare -log10 FDR values for plotting
up_sig <- up_sig %>%
  mutate(logFDR = -log10(p_fdr_bh))

# Arrange by significance and pick top N terms (optional, e.g., top 20)
top_terms <- up_sig %>%
  arrange(p_fdr_bh) %>%
  slice_head(n = 20)

# Plot
ggplot(top_terms, aes(x = logFDR, y = reorder(name, logFDR), fill = NS)) +
  geom_bar(stat = "identity") +
  labs(
    x = "-log10(FDR-adjusted p-value)",
    y = "GO Term",
    title = "downregulated Genes"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_fill_brewer(palette = "Set1", name = "GO Domain")


