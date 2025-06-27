#===============================================================================
# DEG Analysis using limma/voom
# Circadian Data

#-------------------------------------------------------------------------------

# To install packages use biocLite from bioconductor source
# source("http://bioconductor.org/biocLite.R")
# biocLite
# install.packages("BiocUpgrade")
# install.packages("Glimma")
# biocLite("Glimma")

## Load library
library(ggplot2)
library(edgeR)
library(limma)
library(dplyr)
library(factoextra)
library(FactoMineR)
library(gtools)
library(tidyverse)
library(DESeq2)

#-------------------------------------------------------------------------------

# Set working directory
setwd("/Users/judyabuel/Desktop/Xist/circadian_atlas")

##===========Import Data========================================================

# Read raw count file
counts <- read.delim("/Users/judyabuel/Desktop/Xist/circadian_atlas/GSE297702_circadian_atlas_rawcounts.txt", row.names = 1)

# Display the first few rows of the data frame
head(counts)


##===========Create Metadata============
# Read your metadata
metadata <- read.csv("/Users/judyabuel/Desktop/Xist/circadian_atlas/wt_circadian_traits.csv")

# Order metadata to match d exactly
metadata <- metadata[order(metadata$Sample_ID), ]

# Verify perfect alignment before proceeding
stopifnot(identical(colnames(d), metadata$Sample_ID))

# Create design matrix - now with guaranteed alignment
design <- model.matrix(~ Sex + Timepoint + Sex:Timepoint, data = metadata)
rownames(design) <- metadata$Sample_ID

# Final verification
stopifnot(
  nrow(design) == ncol(d),
  identical(rownames(design), colnames(d))
)


##===========Create a DGEList===================================================

d0 <- DGEList(counts)
d0 <- calcNormFactors(d0, method = "TMM")

# Filter lowly expressed genes (edgeR way)
keep <- rowSums(cpm(d0) >= 1) >= 3  # Keep genes with CPM ≥1 in at least 3 samples
d <- d0[keep, , keep.lib.sizes=FALSE]
dim(d) # Check how many number of genes left

design <- model.matrix(~ Sex + Timepoint + Sex:Timepoint, data = metadata)

v <- voom(d, design, plot = TRUE) # check the mean variance trend

# Fit linear model
fit <- lmFit(v, design) #standard limma pipeline
fit <- eBayes(fit) #Empirical Bayes moderation



#============Generate MDS plot Sex:Timepoint===================
mds <- plotMDS(d, 
               col = ifelse(metadata$Sex == "Male", "blue", "pink"),  # Colors by sex
               pch = ifelse(metadata$Timepoint == 0, 16,               # Shapes by timepoint
                            ifelse(metadata$Timepoint == 3, 17,
                                   ifelse(metadata$Timepoint == 6, 15, 18))),
               cex = 1.5,  # Point size
               main = "MDS Plot: Sex and Timepoint",
               xlab = "Dimension 1",
               ylab = "Dimension 2")

# Add legend for Sex
legend("topright", 
       legend = c("Male", "Female"), 
       col = c("blue", "pink"), 
       pch = 16, 
       title = "Sex")

# Add legend for Timepoint
legend("bottomright", 
       legend = metadata$Timepoint, 
       #pch = c(16, 17, 15, 18, 8, 9, 10, 12), 
       title = "Timepoint")

#-------MDS PLOT between Sexes-------------

# Convert MDS coordinates to data frame
mds_data <- data.frame(
  Dim1 = mds$x,
  Dim2 = mds$y,
  Sex = metadata$Sex
)

# Plot
ggplot(mds_data, aes(Dim1, Dim2, color = Sex)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Male" = "blue", "Female" = "pink")) +
  labs(x = "Dimension 1", y = "Dimension 2", title = "MDS Plot: Male vs. Female") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),  # Remove gray background
    axis.line = element_line(color = "black")  # Keep axis lines
  )


##========Make a contrast here==================================================

results_list <- list()

for (tp in timepoints) {
  # Subset data for this time point
  idx <- metadata$Timepoint == tp
  d_sub <- d[, idx]
  metadata_sub <- metadata[idx, ]
  
  # Redesign for simple Male vs Female comparison
  design_sub <- model.matrix(~ Sex, data = metadata_sub)
  d_sub <- estimateDisp(d_sub, design_sub)
  fit_sub <- glmQLFit(d_sub, design_sub)
  
  # Test Male vs Female (assumes "Male" is the numerator)
  qlf <- glmQLFTest(fit_sub, coef = 2)  # coef=2 is "sex_Male"
  res <- topTags(qlf, n = Inf)$table %>%
    rownames_to_column("gene") %>%
    mutate(timepoint = tp)
  
  results_list[[as.character(tp)]] <- res
}

all_results <- bind_rows(results_list)

write.csv(all_results, "Male_vs_Female_by_timepoint_edgeR.csv", row.names = FALSE)


# what are you comparing in your samples? Sex differences and timepoints
fit <- eBayes(lmFit(v, design))
results <- topTable(fit, coef = "SexMale", number = Inf)




# Classify genes
results <- results %>%
  mutate(Significance = case_when(
    adj.P.Val < sig_cutoff & logFC > fc_cutoff ~ "Male_up",
    adj.P.Val < sig_cutoff & logFC < -fc_cutoff ~ "Female_up",
    TRUE ~ "Not_sig"
  ))

# Custom color scale
color_scale <- c(
  "Male_up" = "#1f78b4",  # Blue
  "Female_up" = "#e31a1c", # Red
  "Not_sig" = "gray"
)

# Generate plot
volcano <- ggplot(results, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = Significance), alpha = 0.6, size = 2) +
  scale_color_manual(values = color_scale,
                     name = "Expression",
                     labels = c(paste("Female up (", sum(results$Significance == "Female_up"), ")"),
                                paste("Male up (", sum(results$Significance == "Male_up"), ")"),
                                "Not significant")) +
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(sig_cutoff), linetype = "dashed", color = "black") +
  labs(x = "log2 Fold Change (Male vs Female)",
       y = "-log10(p-value)",
       title = "DEGs Between Male and Female") +
  theme_minimal() +
  theme(legend.position = "bottom")

volcano


library(biomaRt)  # For gene ID conversion
library(org.Mm.eg.db)  # Mouse genome database (replace with org.Hs.eg.db for human)

# For mouse data (replace "mmusculus" with "hsapiens" for human)
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Get gene names
gene_names <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = rownames(results),
  mart = mart
)

# Merge with your results
results$gene_name <- gene_names$external_gene_name[match(rownames(results), gene_names$ensembl_gene_id)]

# Filter significant genes
sig_genes <- results %>%
  filter(Significance %in% c("Male_up", "Female_up"))

# Add labels to the existing plot
volcano_with_labels <- volcano +
  geom_text_repel(
    data = sig_genes,
    aes(label = gene_name),  # Use gene_name instead of rownames
    size = 3,
    box.padding = 0.5,
    max.overlaps = 50,
    segment.color = "grey50"  # Color of label lines
  )

# Display the labeled plot
print(volcano_with_labels)

volcano_clean <- volcano_with_labels +
  theme_classic() +  # Removes grid and sets white background
  theme(legend.position = "bottom")  # Optional legend adjustment

ggsave("volcano_with_labels_clean.png", plot = volcano_clean, width = 10, height = 8, dpi = 300)


##==============================================================================





