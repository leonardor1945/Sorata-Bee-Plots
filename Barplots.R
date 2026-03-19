#### Barplots ####

# 1. Load Packagess
library(phyloseq)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

# 2. Load and Prepare Data
# Make sure the path to your BIOM file is correct
ps <- import_biom("/your/HiFi-16S-workflow/Results/results/feature-table-tax_vsearch.biom")

# Map taxonomic ranks 
colnames(tax_table(ps)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Clean database prefixes (d__, s__, etc.)
tax_table(ps) <- gsub("[a-z]__", "", tax_table(ps))

# Remove ASVs with zero counts
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

# Transform to Relative Abundance (0 to 1 scale)
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

# --- PART 1: FULL COMMUNITY PLOTS (100% ABUNDANCE) ---

# Total Abundance: Family Level
p_fam_total <- plot_bar(ps_rel, fill="Family") + 
  geom_bar(aes(fill=Family), stat="identity", position="stack") +
  theme_minimal() + 
  labs(title="Full Microbiota Composition: Family Level", 
       y="Relative Abundance", x="Bee Samples (Sorata)") +
  theme(legend.text = element_text(size = 8))

ggsave("Full_Community_Family.pdf", plot=p_fam_total, width=11, height=8)

# Total Abundance: Genus Level
p_gen_total <- plot_bar(ps_rel, fill="Genus") + 
  geom_bar(aes(fill=Genus), stat="identity", position="stack") +
  theme_minimal() + 
  labs(title="Full Microbiota Composition: Genus Level", 
       y="Relative Abundance", x="Bee Samples (Sorata)") +
  theme(legend.text = element_text(size = 7))

ggsave("Full_Community_Genus.pdf", plot=p_gen_total, width=12, height=8)


# --- PART 2: TOP 10 GENERA ---

# Agglomerate by Genus to sum all ASVs of the same genus
ps_genus_glom <- tax_glom(ps_rel, taxrank = "Genus")
top10_gen_list <- names(sort(taxa_sums(ps_genus_glom), decreasing = TRUE)[1:10])
ps_top10_gen <- prune_taxa(top10_gen_list, ps_genus_glom)

p_top10_gen <- plot_bar(ps_top10_gen, fill="Genus") + 
  geom_bar(aes(fill=Genus), stat="identity", position="stack") +
  scale_fill_brewer(palette = "Paired") +
  theme_bw() + 
  labs(title="Top 10 Most Abundant Genera", 
       subtitle="Bee Rectum Samples - Sorata, Bolivia",
       y="Relative Abundance", x="Samples")

ggsave("Top10_Genera.pdf", plot=p_top10_gen, width=10, height=7)


# --- PART 3: TOP 10 SPECIES ---

# Agglomerate by Species
ps_spec_glom <- tax_glom(ps_rel, taxrank = "Species")
top10_spec_list <- names(sort(taxa_sums(ps_spec_glom), decreasing = TRUE)[1:10])
ps_top10_spec <- prune_taxa(top10_spec_list, ps_spec_glom)

p_top10_spec <- plot_bar(ps_top10_spec, fill="Species") + 
  geom_bar(aes(fill=Species), stat="identity", position="stack") +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() + 
  labs(title="Top 10 Most Abundant Species", 
       subtitle="Based on Full-Length 16S (V1-V9) HiFi sequencing",
       y="Relative Abundance", x="Samples") +
  theme(legend.text = element_text(face = "italic", size = 9))

ggsave("Top10_Species.pdf", plot=p_top10_spec, width=10, height=7)

print("Process finished 4 PDFs have been created.")
