library(ggplot2)
library(dplyr)

# Helper function to classify genes per condition
classify_genes <- function(peak_data) {
  cis_genes <- peak_data %>%
    filter(cis_trans == "cis") %>%
    pull(gene) %>%
    unique()
  
  trans_genes <- peak_data %>%
    filter(cis_trans == "trans") %>%
    pull(gene) %>%
    unique()
  
  both_genes <- intersect(cis_genes, trans_genes)
  cis_only <- setdiff(cis_genes, both_genes)
  trans_only <- setdiff(trans_genes, both_genes)
  
  data.frame(
    Category = c("Cis", "Trans", "Both"),
    Count = c(length(cis_only), length(trans_only), length(both_genes))
  )
}

# Apply to both conditions
hn_counts <- classify_genes(peak_HN) %>%
  mutate(Condition = "High N")

ln_counts <- classify_genes(peak_LN) %>%
  mutate(Condition = "Low N")

# Combine
all_counts <- bind_rows(hn_counts, ln_counts)

# Plot
ggplot(all_counts, aes(x = Condition, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(aes(label = Count), 
            position = position_dodge(width = 0.75), 
            vjust = -0.4, size = 5) +
  scale_fill_manual(values = c("Cis" = "#1b9e77", "Trans" = "#d95f02", "Both" = "#7570b3")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Unique Genes with Cis, Trans, or Both Associations",
    subtitle = "Grouped by Condition (HN vs LN)",
    x = "Condition",
    y = "Number of Unique Genes"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.position = "top"
  )



