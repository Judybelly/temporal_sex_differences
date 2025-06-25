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

##===========Import Data==================

# Read raw count file
counts <- read.delim("/Users/judyabuel/Desktop/Xist/circadian_atlas/GSE297702_circadian_atlas_rawcounts.txt", row.names = 1)

# Display the first few rows of the data frame
head(counts)

##===========Normalize & Filter Data==============

d0 <- DGEList(counts)
d0 <- calcNormFactors(d0, method = "TMM")

# Filter lowly expressed genes (edgeR way)
keep <- rowSums(cpm(d0) >= 1) >= 3  # Keep genes with CPM ≥1 in at least 3 samples
d <- d0[keep, , keep.lib.sizes=FALSE]
dim(d) # Check how many number of genes left

##===========Create Proper Metadata============

# Create new metadata with EXACTLY the samples in d
metadata <- data.frame(
  Sample = colnames(d),
  Sex = factor(ifelse(grepl("^A", colnames(d)), "Male", "Female")),
  Timepoint = as.numeric(gsub("^A|^B", "", colnames(d))),
  stringsAsFactors = FALSE
)

# Order metadata to match d exactly
metadata <- metadata[order(metadata$Sample), ]

# Verify perfect alignment before proceeding
stopifnot(identical(colnames(d), metadata$Sample))

# Create design matrix - now with guaranteed alignment
design <- model.matrix(~ Sex + Timepoint + Sex:Timepoint, data = metadata)
rownames(design) <- metadata$Sample  # This MUST work now

# Final verification
stopifnot(
  nrow(design) == ncol(d),
  identical(rownames(design), colnames(d))
)

# Now run voom
v <- voom(d, design, plot = TRUE)


##==========Quality Control Plots=================
# MDS Plot
plotMDS(d0, 
        col = as.numeric(group),  # Colors: 1=Male, 2=Female
        pch = 16,                 # Solid circles
        cex = 1.5,                 # Point size
        main = "MDS Plot: Male vs Female",
        xlab = "Dimension 1",
        ylab = "Dimension 2")

# Add a legend to my plot
legend("topright", 
       legend = levels(group), 
       col = 1:2, 
       pch = 16, 
       title = "Sex")


##===========Create Exact Metadata============
# Manually specify timepoints for ALL samples in EXACT order they appear in your data
timepoints <- c(
  # Male samples (A) - first 46 samples
  rep(0, 7),  # A1-A3, A25-A27,A3
  rep(3, 6),   # A28-A30,A4-A6
  rep(6, 6),   # A31-A33,A7-A9
  rep(9, 6),   # A10-A12,A34-A36
  rep(12, 6),  # A13-A15,A37-A39
  rep(15, 6),  # A16-A18,A40-A42
  rep(18, 5),  # A19-A21,A43-A44
  rep(21, 4),  # A22-A24,A45-A46
  
  # Female samples (B) - next 46 samples
  rep(0, 7),   # B1-B3,B25-B27,B3
  rep(3, 6),    # B28-B30,B4-B6
  rep(6, 6),    # B31-B33,B7-B9
  rep(9, 6),    # B10-B12,B34-B36
  rep(12, 6),   # B13-B15,B37-B39
  rep(15, 6),   # B16-B18,B40-B42
  rep(18, 5),   # B19-B21,B43-B44
  rep(21, 4)    # B22-B24,B45-B46
)

# Create metadata - this CANNOT fail
metadata <- data.frame(
  Sample_ID = colnames(d),
  Sample_Name = ifelse(
    grepl("^A", colnames(d)),
    paste0("ZT", sprintf("%02d", timepoints[1:ncol(d)]), "_", colnames(d)),
    paste0("ZT", sprintf("%02d", timepoints[1:ncol(d)]), "_", colnames(d))
  ),
  Timepoint = factor(timepoints[1:ncol(d)], levels = c(0, 3, 6, 9, 12, 15, 18, 21)),
  Sex = ifelse(grepl("^A", colnames(d)), "Male", "Female"),
  stringsAsFactors = FALSE
)

# Verify
head(metadata, 20)
table(metadata$Sex, metadata$Timepoint)


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
       legend = levels(metadata$Timepoint), 
       pch = c(16, 17, 15, 18, 8, 9, 10, 12), 
       title = "Timepoint")




# Convert MDS coordinates to data frame
mds_data <- data.frame(
  Dimension1 = mds$x,
  Dimension2 = mds$y,
  Sex = metadata$Sex,
  Timepoint = metadata$Timepoint
)

# Create ggplot
ggplot(mds_data, aes(x = Dimension1, y = Dimension2, color = Sex, shape = Timepoint)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Male" = "blue", "Female" = "pink")) +
  scale_shape_manual(values = c(16, 17, 15, 18, 8, 9, 10, 12)) +
  labs(title = "MDS Plot: Sex and Timepoint",
       x = "Dimension 1",
       y = "Dimension 2") +
  theme_minimal() +  # Start with minimal theme
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),  # Remove gray background
    axis.line = element_line(color = "black")  # Keep axis lines
  )


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


##============================================
##______Visualize DEGs between Male and Female________________
design <- model.matrix(~ Sex, data = metadata)
v <- voom(d, design, plot = TRUE)
fit <- eBayes(lmFit(v, design))
results <- topTable(fit, coef = "SexMale", number = Inf)



# Ensure columns are numeric (fixes common errors)
results$logFC <- as.numeric(results$logFC)
results$P.Value <- as.numeric(results$P.Value)
results$adj.P.Val <- as.numeric(results$adj.P.Val)

# Define significance thresholds
sig_cutoff <- 0.01  # FDR cutoff
fc_cutoff <- 2      # |logFC| cutoff





# Define thresholds
sig_cutoff <- 0.05  # FDR cutoff
fc_cutoff <- 1      # Absolute logFC cutoff

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







