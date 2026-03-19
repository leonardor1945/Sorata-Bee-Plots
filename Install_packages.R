# 1. Install Packages 
install.packages(c("ggplot2", "vegan", "reshape2", "RColorBrewer", "dplyr", "readr", "ggpubr"), 
                 repos = "http://cran.us.r-project.org")

# 2. Install BiocManager (if needed)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# 3. Install Bioconductor Packages (for 16S)
BiocManager::install(c("phyloseq", "microbiome", "complexheatmap"))

# 4. Verify Packages
library(phyloseq)
library(microbiome)
library(vegan)
library(ggplot2)
library(dplyr)
library(reshape2)
library(readr)

print("R ready for analysis")

