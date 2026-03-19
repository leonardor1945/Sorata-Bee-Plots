#### Alpha Diversity ####

library(phyloseq)
library(microbiome)
library(ggplot2)
library(reshape2)

# 1. Load and Clean
ps <- import_biom("/your/HiFi-16S-workflow/Results/results/feature-table-tax_vsearch.biom")
colnames(tax_table(ps)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax_table(ps) <- gsub("[a-z]__", "", tax_table(ps))
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

# 2. Plot Alpha Diversity with phyloseq
p_alpha <- plot_richness(ps, measures = c("Shannon", "Simpson")) + 
  geom_boxplot(fill="skyblue", alpha=0.5) +
  geom_point(size=3, color="darkblue") +
  theme_bw() +
  labs(title="Alpha Diversity - Sorata Bee Rectum", x="Samples")

# Show Plot
print(p_alpha)

# Save PDF
ggsave("Alpha_Diversity.pdf", plot=p_alpha, width=8, height=6)
