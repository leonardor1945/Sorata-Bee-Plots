#### Rarefaction Curve ####

# 1. Upload Packages
library(phyloseq)
library(vegan)
library(RColorBrewer)

# 2. Upload and Prepare data
ps <- import_biom("/your/HiFi-16S-workflow/Results/results/feature-table-tax_vsearch.biom")

# 3. Obtain Abundance Matrix (OTU Table)
# Vegan needs the samples to be arranged in rows and species in columns
otu_tab <- as.data.frame(t(otu_table(ps)))

# 4. Define Colors for Samples
colores <- brewer.pal(n = min(nrow(otu_tab), 8), name = "Dark2")

# 5. Generate Rarefaction Plots
# 'step' defines the curve resolution. 100 is ideal for PacBio.
pdf("Rarefaction_Curves.pdf", width = 10, height = 7)

rare_res <- rarecurve(otu_tab, 
                      step = 100, 
                      col = colores, 
                      lwd = 2, 
                      label = TRUE,
                      main = "Rarefaction Curves: Bee Rectum Samples (Sorata)",
                      xlab = "Sequencing Depth (Number of Reads)", 
                      ylab = "Observed ASVs (Species Richness)")

# Add a vertical line to indicate the sample with the fewest readings (cutoff point)
abline(v = min(rowSums(otu_tab)), col = "red", lty = 2, lwd = 1.5)
legend("bottomright", legend = "Minimum depth cut-off", col = "red", lty = 2, bty = "n")

dev.off()

print("Rarefaction curve saved as 'Rarefaction_Curves_Sorata.pdf'")
