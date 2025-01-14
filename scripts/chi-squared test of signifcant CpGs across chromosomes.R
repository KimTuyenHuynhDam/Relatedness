# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readxl)
library(ggrepel)


# Load the data
data <- read.xlsx("./DNAm_BW/nested_anova_final_hu_pman_2.1.100_anno.xlsx", sheet = "Sheet 1") %>%

 filter(seqnames %in% c(as.character(1:22), "X")) 


# Calculate the overall ratio of significant CpGs (FDR < 0.1)
total_CpGs <- nrow(data)  # Total CpGs
significant_CpGs <- data %>%
  filter(FDR < 0.1) %>%
  nrow()  # Total significant CpGs

overall_ratio <- significant_CpGs / total_CpGs  # Ratio of significant CpGs




#  Calculate chromosome-specific ratios, expected values, and test deviation direction
prepared_data <- data %>%
  distinct() %>%
  group_by(seqnames) %>% 
  summarize(
    total_CpGs = n(),
    significant_CpGs = sum(FDR < 0.1, na.rm = TRUE),
    non_significant_CpGs = total_CpGs - significant_CpGs,
    ratio = significant_CpGs / total_CpGs,  # Actual ratio for this chromosome
    expected_significant_CpGs = total_CpGs * overall_ratio,  # Expected significant CpGs
    expected_non_significant_CpGs = total_CpGs - expected_significant_CpGs,  # Expected non-significant CpGs
    deviation = significant_CpGs - expected_significant_CpGs,  # Difference from expected
    direction = ifelse(deviation > 0, "Higher", "Lower")  # Higher or lower than expected
  ) %>%
  rowwise()



# Perform chi-squared test for each row and add p-value as a new column
chi_squared_results <- prepared_data %>%
  rowwise() %>%
  mutate(
    p_value = {
      # Observed values
      observed <- c(significant_CpGs, non_significant_CpGs)
      # Expected values
      expected <- c(expected_significant_CpGs, expected_non_significant_CpGs)
      # Chi-squared test
      chisq_result <- chisq.test(x = observed, p = expected / sum(expected))
      chisq_result$p.value
    }
  ) %>%
  ungroup()  # Remove rowwise structure





write.xlsx(chi_squared_results, "./plots/chi_squared-BW_dataset/chi_squared of significant CpGs over total CpGs across chromosomes.xlsx")


# plots

visual_data <- chi_squared_results %>%
  
  mutate(seqnames = factor(seqnames, levels = c(as.character(1:22), "X")))

# Bar Plot for Observed vs. Expected Significant CpGs
bar_plot <- ggplot(visual_data, aes(x = seqnames, fill = direction)) +
  geom_bar(aes(y = significant_CpGs), stat = "identity", alpha = 0.8, position = "dodge") +
  geom_point(aes(y = expected_significant_CpGs), shape = 21, size = 3, color = "black", fill = "yellow") +
  labs(
    title = "Observed vs. Expected Significant CpGs (FDR < 0.1)",
    x = "Chromosome",
    y = "Significant CpGs",
    fill = "Deviation Direction"
  ) +
  scale_fill_manual(values = c("Higher" = "red", "Lower" = "blue")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("./plots/chi_squared-BW_dataset/observed_vs_expected_barplot.png", bar_plot, width = 12, height = 6)
print(bar_plot)

# Scatter Plot for Observed vs. Expected Significant CpGs

scatter_plot <- ggplot(visual_data, aes(x = expected_significant_CpGs, y = significant_CpGs)) +
  geom_point(aes(color = direction, size = -log10(p_value)), alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  geom_text_repel(aes(label = seqnames), size = 3) +
  labs(
    title = "Scatter Plot: Observed vs. Expected Significant CpGs",
    x = "Expected Significant CpGs",
    y = "Observed Significant CpGs",
    color = "Deviation Direction",
    size = "Significance (-log10 p-value)"
  ) +
  scale_color_manual(values = c("Higher" = "red", "Lower" = "blue")) +
  theme_minimal()



scatter_plot <- ggplot(visual_data, aes(x = expected_significant_CpGs, y = significant_CpGs)) +
  geom_point(aes(color = direction, size = -log10(p_value), shape = p_value < 0.05), alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  geom_text_repel(
    aes(label = seqnames),  # Label all dots with chromosome names
    size = 3,
    box.padding = 0.8,      # Increase space around labels
    point.padding = 0.5,    # Add spacing between labels and points
    max.overlaps = Inf,     # Allow all labels to display
    segment.size = 0.3,     # Thinner line connecting label and point
    segment.color = "gray50"  # Lighter connecting lines
  ) +
  labs(
    title = "Scatter Plot: Observed vs. Expected Significant CpGs",
    x = "Expected Significant CpGs",
    y = "Observed Significant CpGs",
    color = "Deviation Direction",
    size = "-Log10 p-value",
    shape = "Significant (p < 0.05)"
  ) +
  scale_color_manual(values = c("Higher" = "red", "Lower" = "blue")) +
  scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 16)) +  # Triangle for significant, circle for non-significant
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis text for clarity
  )


ggsave("./plots/chi_squared-BW_dataset/scatter_plot_all_chromosomes_labeled.png", scatter_plot, width = 12, height = 6)
print(scatter_plot)





# Percentage Significant CpGs Plot
percentage_plot <- ggplot(visual_data, aes(x = seqnames, y = ratio * 100, fill = direction)) +
  geom_bar(stat = "identity", alpha = 0.8, position = "dodge") +
  labs(
    title = "Percentage of Significant CpGs (FDR < 0.1) per Chromosome",
    x = "Chromosome",
    y = "Percentage Significant CpGs (%)",
    fill = "Deviation Direction"
  ) +
  scale_fill_manual(values = c("Higher" = "red", "Lower" = "blue")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("./plots/chi_squared-BW_dataset/percentage_significant_CpGs.png", percentage_plot, width = 12, height = 6)
print(percentage_plot)

# -log10(P-value) Plot
log_p_plot <- ggplot(visual_data, aes(x = seqnames, y = -log10(p_value))) +
  geom_point(aes(color = direction, size = ratio), alpha = 0.8) +  # Map size to ratio
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_text_repel(aes(label = seqnames), size = 3) +
  labs(
    title = "-log10(P-value) of Chromosome Deviations",
    x = "Chromosome",
    y = "-log10(P-value)",
    color = "Deviation Direction",
    size = "Significant CpGs Ratio"
  ) +
  scale_color_manual(values = c("Higher" = "red", "Lower" = "blue")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("./plots/chi_squared-BW_dataset/log_pvalue_plot.png", log_p_plot, width = 12, height = 6)
print(log_p_plot)




