## Load library
library(ggplot2)
library(edgeR)
library(limma)
library(dplyr)
library(factoextra)
library(FactoMineR)
library(gtools)
library(tidyverse)
library(DESeq2)   # For DEG analysis

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicFeatures")  # For GTF/GFF import
BiocManager::install("rtracklayer")     # Alternative for GTF
BiocManager::install("AnnotationDbi")   # For gene database
BiocManager::install("txdbmaker")


## Load annotations GTF file
library(GenomicFeatures)

# Path to your downloaded GTF file (e.g., gencode.vM35.annotation.gtf.gz)
gtf_file <- "/Users/judyabuel/Desktop/Xist/circadian_atlas/gencode.vM35.annotation.gtf.gz"

# Create a TxDb object (database of genomic features)
txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")

# Extract genes, transcripts, or exons
genes <- genes(txdb)
transcripts <- transcripts(txdb)
exons <- exonsBy(txdb, by = "gene")  # Exons grouped by gene

# Read rawcount file
counts <- 
  read.delim("/Users/judyabuel/Desktop/Xist/circadian_atlas/GSE297702_circadian_atlas_rawcounts.txt", row.names = 1)
# Display the first few rows of the data frame
head(counts)


##===========Normalize & Filter Data==============
d0 <- DGEList(counts)
d0 <- calcNormFactors(d0, method = "TMM")


# Adding more information about your samples
# This prepares your data for your further analysis
traits <- 
  read.csv("/Users/judyabuel/Desktop/Xist/circadian_atlas/wt_circadian_traits.csv", header = TRUE)

# View the first few rows
head(traits)

##------SUBSET counts matrix----------------
# Get samples present in BOTH counts_data and sample_info
common_samples <- intersect(colnames(counts), traits$Sample_ID)
# Keep only columns (samples) that exist in metadata
counts <- counts[, common_samples]
counts_data <- counts[, common_samples]
View(counts_data)

# Ensure metadata rows are in the same order as counts columns
traits_info <- traits[match(common_samples, traits$Sample_ID), ]
View(traits_info)
# Now this should return TRUE
all(colnames(counts_data) == traits_info$Sample_ID)


# Convert categorical variables to factors (critical for statistical modeling)
traits_info$Timepoint <- factor(traits_info$Timepoint)
traits_info$Sex <- factor(traits_info$Sex)

# Create the DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = counts_data,       # Your aligned count matrix
  colData = traits_info,         # Your aligned metadata
  design = ~ Timepoint           # Design formula (adjust as needed)
)

keep <- rowSums(counts(dds) >= 10) >=3 # Keep genes with ≥10 counts in ≥3 samples
dds <- dds[keep, ]

# Normalize data and fit models
dds <- DESeq(dds)

# Extract results (e.g., Timepoint 3 vs 0)
results <- results(dds, contrast = c("Timepoint", "3", "0"))

# Get significant DEGs (padj < 0.05 and |log2FC| > 1)
degs <- results[!is.na(results$padj) & results$padj < 0.05 & abs(results$log2FoldChange) > 1, ]
degs <- degs[order(degs$padj), ]

# Filters low-expression gene on the dataset stored in DGEList object (d0)
# cutoff <- 1   # sets the threshold for filtering genes with low expression
#drop <- which(apply(cpm(d0),1,max) < cutoff)  # computes cpm for each gene in d0 and normalize the counts by library size
#d <- d0[-drop,]   # creates a new DGEList excluding low-expression genes
#dim(d) # Check how many number of genes left

# Install from Bioconductor (if not already installed)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("rtracklayer")

# Load the package
library(rtracklayer)

# Import GTF and filter for gene entries
gtf <- import("/Users/judyabuel/Desktop/Xist/circadian_atlas/gencode.vM35.annotation.gtf.gz")
genes_gtf <- gtf[gtf$type == "gene"]

# Create a gene ID-to-name mapping table
gene_annot <- data.frame(
  gene_id = sub("\\..*", "", genes_gtf$gene_id),  # Remove version numbers (e.g., ENSMUSG00000102693.5 → ENSMUSG00000102693)
  gene_name = genes_gtf$gene_name,
  gene_type = genes_gtf$gene_type
)

# Convert DEG results to data frame
degs_df <- as.data.frame(degs)
degs_df$gene_id <- rownames(degs_df)  # Add Ensembl IDs as a column

# Merge with annotations
degs_annotated <- merge(degs_df, gene_annot, by = "gene_id", all.x = TRUE)

# Sort by significance (padj)
degs_annotated <- degs_annotated[order(degs_annotated$padj), ]

# Save to file
write.csv(degs_annotated, "DEGs_annotated_with_GTF.csv", row.names = FALSE)


##==========CREATE PLOTS=================
# 1. Get variance-stabilized transformed (VST) counts
vst_counts <- vst(dds)  # 'dds' is your DESeqDataSet object

# 2. Calculate sample distances
sample_dist <- dist(t(assay(vst_counts)))

# 3. Perform MDS (PCA-like)
mds_data <- cmdscale(sample_dist, eig = TRUE, k = 2)  # k=2 dimensions
mds_plot <- data.frame(
  Sample = rownames(mds_data$points),
  Dim1 = mds_data$points[, 1],
  Dim2 = mds_data$points[, 2]
)

# 4. Merge with metadata (e.g., Timepoint, Sex)
mds_plot <- merge(mds_plot, traits_info, by.x = "Sample", by.y = "Sample_ID")

# 5. Plot
ggplot(mds_plot, aes(x = Dim1, y = Dim2, color = Timepoint, shape = Sex)) +
  geom_point(size = 4) +
  ggtitle("MDS Plot of VST-Normalized Counts") +
  xlab("MDS Dimension 1") +
  ylab("MDS Dimension 2") +
  theme_minimal()

## USING LIMMA========================================
# 1. Convert counts to logCPM (if not using DESeq2)
dge <- DGEList(counts = counts_data)  # Your raw counts
dge <- calcNormFactors(dge)
logcpm <- cpm(dge, log = TRUE)

keep <- filterByExpr(dge, group = traits_info$Timepoint)  # Or another condition
dge <- dge[keep, , keep.lib.sizes = FALSE]                # Apply filtering

# 2. Generate MDS plot
plotMDS(logcpm, 
        col = as.numeric(factor(traits_info$Timepoint)),
        pch = as.numeric(factor(traits_info$Sex)),
        main = "MDS Plot of logCPM-Normalized Counts")
legend("topright", 
       legend = levels(factor(traits_info$Timepoint)),
       col = 1:3, pch = 16, title = "Timepoint")

# Design matrix
design <- model.matrix(~ Timepoint + Sex, data = traits_info)

# Run voom
v <- voom(dge, design, plot = TRUE)  # Should show a nice Mean-Variance trend
fit <- lmFit(v, design)          # Fit linear model
fit <- eBayes(fit)               # Empirical Bayes moderation

##==========Compare each Timepoint against baseline (Timepoint0)================
# Find genes that change at each timepoint vs. the baseline
traits_info$Timepoint <- factor(traits_info$Timepoint, levels = c(0, 3, 6, 9, 12, 15, 18, 21))

design <- model.matrix(~ 0 + Timepoint, data = traits_info)
colnames(design) <- c("Timepoint0", "Timepoint3", "Timepoint6", "Timepoint9", 
                      "Timepoint12", "Timepoint15", "Timepoint18", "Timepoint21")

# Fit the model
fit <- lmFit(counts, design)

# Define contrasts (e.g., Timepoint3 vs Timepoint0)
contrasts <- makeContrasts(
  T3_vs_T0 = Timepoint3 - Timepoint0,
  T6_vs_T0 = Timepoint6 - Timepoint0,
  T9_vs_T0 = Timepoint9 - Timepoint0,
  T12_vs_T0 = Timepoint12 - Timepoint0,
  T15_vs_T0 = Timepoint15 - Timepoint0,
  T18_vs_T0 = Timepoint18 - Timepoint0,
  T21_vs_T0 = Timepoint21 - Timepoint0,
  # Add all other timepoints...
  levels = design
)

# Run differential expression
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

# Extract results results_T3  <- topTable(fit2, coef = "T3_vs_T0",  number = Inf, adjust.method = "BH")
results_T3  <- topTable(fit2, coef = "T3_vs_T0",  number = Inf, adjust.method = "BH")
results_T6  <- topTable(fit2, coef = "T6_vs_T0",  number = Inf, adjust.method = "BH")
results_T9  <- topTable(fit2, coef = "T9_vs_T0",  number = Inf, adjust.method = "BH")
results_T12  <- topTable(fit2, coef = "T12_vs_T0",  number = Inf, adjust.method = "BH")
results_T15  <- topTable(fit2, coef = "T15_vs_T0",  number = Inf, adjust.method = "BH")
results_T18  <- topTable(fit2, coef = "T18_vs_T0",  number = Inf, adjust.method = "BH")
results_T21  <- topTable(fit2, coef = "T21_vs_T0",  number = Inf, adjust.method = "BH")

# Save to CSV
write.csv(
  results_T3,
  file = "results_T3_vs_T0.csv",
  row.names = TRUE,    # Keep gene names
  quote = FALSE        # Avoid quotes
)

write.csv(
  results_T6,
  file = "results_T6_vs_T0.csv",
  row.names = TRUE,    # Keep gene names
  quote = FALSE        # Avoid quotes
)

write.csv(
  results_T9,
  file = "results_T9_vs_T0.csv",
  row.names = TRUE,    # Keep gene names
  quote = FALSE        # Avoid quotes
)

write.csv(
  results_T12,
  file = "results_T12_vs_T0.csv",
  row.names = TRUE,    # Keep gene names
  quote = FALSE        # Avoid quotes
)

write.csv(
  results_T15,
  file = "results_T15_vs_T0.csv",
  row.names = TRUE,    # Keep gene names
  quote = FALSE        # Avoid quotes
)

write.csv(
  results_T18,
  file = "results_T18_vs_T0.csv",
  row.names = TRUE,    # Keep gene names
  quote = FALSE        # Avoid quotes
)

write.csv(
  results_T21,
  file = "results_T21_vs_T0.csv",
  row.names = TRUE,    # Keep gene names
  quote = FALSE        # Avoid quotes
)


##=========GENERATE PLOTS===================
# Ensure data is clean
# Remove NA/Inf and clamp extreme values
volcano_data <- na.omit(volcano_data)
volcano_data$logFC[volcano_data$logFC > 10] <- 10    # Cap max logFC
volcano_data$logFC[volcano_data$logFC < -10] <- -10  # Cap min logFC
volcano_data$adj.P.Val[volcano_data$adj.P.Val < 1e-20] <- 1e-20  # Prevent -log10(0)


# Add significance thresholds (adjust as needed)
volcano_data$Significant <- ifelse(
  volcano_data$adj.P.Val < 0.05 & abs(volcano_data$logFC) > 1, 
  "DEG", 
  "Not DEG"
)

# Plot with refined aesthetics
ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = Significant), alpha = 0.7, size = 2) +
  scale_color_manual(values = c("DEG" = "red", "Not DEG" = "gray60")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", linewidth = 0.5) +
  labs(
    title = "Volcano Plot: T3 vs T0",
    x = "log2 Fold Change",
    y = "-log10(Adjusted P-value)",
    color = "Significance"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top"
  )

#plot(
  #volcano_data$logFC, 
  #-log10(volcano_data$adj.P.Val),
  #pch = 16, 
  #col = ifelse(volcano_data$adj.P.Val < 0.05 & abs(volcano_data$logFC) > 1, "red", "gray"),
  #xlab = "log2 Fold Change", 
  #ylab = "-log10(Adjusted P-value)",
  #main = "Volcano Plot (Base R)"
)
#abline(h = -log10(0.05), v = c(-1, 1), lty = 2)



