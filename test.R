#===============================================================================
# PURPOSE: DEG Analysis using limma/voom
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

#===============================================================================
#                        Import Data & Create Metadata
#===============================================================================

# Read raw count file
counts <- read.delim("/Users/judyabuel/Desktop/Xist/circadian_atlas/GSE297702_circadian_atlas_rawcounts.txt", row.names = 1)
# Display the first few rows of the data frame
head(counts)


# Read your metadata
metadata <- read.csv("/Users/judyabuel/Desktop/Xist/circadian_atlas/wt_circadian_traits.csv")


#===============================================================================
#                              Create DGEList
#===============================================================================

d0 <- DGEList(counts)
d0 <- calcNormFactors(d0, method = "TMM")

# Filter lowly expressed genes (edgeR way)
keep <- rowSums(cpm(d0) >= 1) >= 3  # Keep genes with CPM ≥1 in at least 3 samples
d <- d0[keep, , keep.lib.sizes=FALSE]
dim(d) # Check how many number of genes left

# design matrix
design <- model.matrix(~ Sex * Timepoint, data = metadata)
colnames(design) # Checks Check if your design matrix matches expectations

# Proceed with voom and lmFit as before
v <- voom(d, design, plot = TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit)


#===============================================================================
#                              Define Contrasts
#===============================================================================

coef_names <- colnames(fit$coefficients)
all_results <- list()

for (tp in c("00", "03", "06", "09", "12", "15", "18", "21")) {
  model_coefs <- colnames(fit$coefficients)
  contrast <- matrix(0, nrow = length(model_coefs), ncol = 1,
                     dimnames = list(model_coefs, paste0("ZT", tp)))
  
  if (tp == "00") {
    if ("SexMale" %in% model_coefs) {
      contrast["SexMale", 1] <- -1  # male vs female at ZT00
    } else {
      warning("SexMale coefficient not found for ZT00")
      next
    }
  } else {
    interaction_term <- paste0("SexMale:Timepoint", tp)
    if (all(c("SexMale", interaction_term) %in% model_coefs)) {
      contrast["SexMale", 1] <- -1
      contrast[interaction_term, 1] <- -1
    } else {
      warning(paste("Skipping ZT", tp, "- missing coefficient(s)"))
      next
    }
  }
  
  # Fit the model with the contrast
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  
  # Extract results with logFC included
  results <- topTable(fit2, coef = 1, number = Inf)
  results$Timepoint <- paste0("ZT", tp)
  
  all_results[[paste0("ZT", tp)]] <- results
  
  # Optional: print confirmation
  cat("Results stored for ZT", tp, "\n")
}


# Combine all the differential expression results into one data freame
combined_results <- do.call(rbind, all_results)


# Let's save it!
write.csv(combined_results, "DEG_Female_vs_Male_by_Timepoint.csv")



#===============================================================================
#                              Generate Plots
#===============================================================================

if (!requireNamespace("biomaRt", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("biomaRt")
}
library(biomaRt)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tibble)

#--------------------------
# Map Gene Names
#--------------------------

# Connect to Ensembl for mouse
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Get all unique gene IDs across your data
all_gene_ids <- unique(unlist(lapply(all_results, rownames)))
all_gene_ids <- sub("\\..*$", "", all_gene_ids)  # Remove version numbers like ".1"

# Query gene symbols from Ensembl
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = all_gene_ids,
  mart = ensembl
)


#-------------------------------
# Merge Gene Names Into Results
#-------------------------------
all_results <- lapply(all_results, function(df) {
  df <- df %>%
    tibble::rownames_to_column(var = "gene") %>%
    mutate(gene = sub("\\..*$", "", gene)) %>%
    left_join(gene_map, by = c("gene" = "ensembl_gene_id")) %>%
    rename(gene_symbol = external_gene_name)
})


#--------------------------
# Generate Volcane Plot
#--------------------------
make_volcano_plot <- function(df, timepoint, logFC_cutoff = 1, adjP_cutoff = 0.05, max_labels = 20) {
  
  # Define significance
  df <- df %>%
    mutate(
      Significance = case_when(
        adj.P.Val < adjP_cutoff & logFC > logFC_cutoff  ~ "Upregulated",
        adj.P.Val < adjP_cutoff & logFC < -logFC_cutoff ~ "Downregulated",
        TRUE                                             ~ "Not Significant"
      ),
      gene_label = ifelse(is.na(gene_symbol) | gene_symbol == "", gene, gene_symbol)
    )
  
  # Top genes to label
  label_df <- df %>%
    filter(Significance != "Not Significant") %>%
    arrange(adj.P.Val) %>%
    head(max_labels)
  
  # Create plot
  ggplot(df, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey70")) +
    geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed", color = "darkgrey") +
    geom_hline(yintercept = -log10(adjP_cutoff), linetype = "dashed", color = "darkgrey") +
    geom_text_repel(data = label_df, aes(label = gene_label), size = 3, max.overlaps = Inf) +
    theme_minimal() +
    labs(
      title = paste("Volcano Plot -", timepoint),
      x = "log2 Fold Change",
      y = "-log10 Adjusted p-value"
    )
}

#--------------------------
# SAVE as PNG file
#--------------------------
dir.create("volcano_plots", showWarnings = FALSE)

for (tp in names(all_results)) {
  df <- all_results[[tp]]
  p <- make_volcano_plot(df, tp)
  
  ggsave(
    filename = paste0("volcano_plots/volcano_", tp, ".png"),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
}



  
#===============================================================================
  ## HEATMAP: Top 15 DE Genes per Timepoint with Sex-specific Significance
#===============================================================================

library(dplyr)
library(tidyr)
library(pheatmap)

#--------------------------
# 1. Select Top Genes
#--------------------------

# Function to extract top 15 genes per timepoint
extract_top_genes <- function(df, top_n = 15) {
  df %>%
    arrange(adj.P.Val) %>%
    filter(!is.na(gene_symbol) & gene_symbol != "") %>%
    distinct(gene_symbol, .keep_all = TRUE) %>%
    head(top_n)
}

# Apply to all timepoints
top_genes_list <- lapply(all_results, extract_top_genes, top_n = 15)

# Combine all and keep unique gene symbols
top_genes <- bind_rows(top_genes_list) %>%
  pull(gene_symbol) %>%
  unique()

#--------------------------
# 2. Extract logFC for top genes across all timepoints
#--------------------------

# Build a matrix of logFC (rows = genes, cols = timepoints)
heatmap_matrix <- do.call(cbind, lapply(all_results, function(df) {
  df %>%
    filter(gene_symbol %in% top_genes) %>%
    select(gene_symbol, logFC) %>%
    tibble::column_to_rownames("gene_symbol") %>%
    arrange(match(rownames(.), top_genes))  # Preserve order
})) 

# Fix column names to match timepoints
colnames(heatmap_matrix) <- names(all_results)


#--------------------------
# 3. Plot the Heatmap
#--------------------------
# Desired order of timepoints
ordered_timepoints <- c("ZT00", "ZT03", "ZT06", "ZT09", "ZT12", "ZT15", "ZT18", "ZT21")

# Reorder columns
heatmap_matrix <- heatmap_matrix[, ordered_timepoints]

# Create and save the heatmap with adjusted cell size
png("Top_genes_heatmap.png", width = 8, height = 12, units = "in", res = 300)

pheatmap(
  heatmap_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "row",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  cellwidth = 20,     # width of each column (smaller = narrower boxes)
  cellheight = 10,    # height of each row (larger = taller boxes)
  fontsize_row = 7,
  fontsize_col = 10,
  main = "Top 15 DE Genes per Timepoint",
  filename = "Top_genes_heatmap.png",
  width = 8,          # in inches
  height = 12         # increase if labels get cut off
)

#--------------------------
# Save as PDF file
#--------------------------

pdf("Heatmap_top_genes_heatmap.pdf", width = 8, height = 12)

pheatmap(
  heatmap_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "row",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  cellwidth = 20,
  cellheight = 10,
  fontsize_row = 7,
  fontsize_col = 10,
  main = "Heatmap_Top 15 DE Genes per Timepoint"
)

dev.off()
  

  
  
  
#===============================================================================
## GO Analysis: Top 15 DE Genes per Timepoint with Sex-specific Significance
#===============================================================================  

if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler")
}
if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Mm.eg.db")
}
if (!requireNamespace("enrichplot", quietly = TRUE)) {
  BiocManager::install("enrichplot")
}

library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(dplyr)


# Collect top genes from each timepoint
top_genes <- lapply(all_results, function(df) {
  # Extract gene IDs from rownames before converting to tibble
  gene_ids <- rownames(df)
  gene_ids <- sub("\\..*$", "", gene_ids)  # remove version numbers
  
  # Add as a column (without converting to tibble yet)
  df$gene <- gene_ids
  
  df <- df %>%
    filter(adj.P.Val < 0.05) %>%
    arrange(adj.P.Val) %>%
    head(100)  # or top 15
  
  return(df$gene)
})

all_top_genes <- unique(unlist(top_genes))
head(all_top_genes)  # should show "ENSMUSG..." IDs now


library(clusterProfiler)
library(org.Mm.eg.db)

gene_df <- bitr(all_top_genes,
                fromType = "ENSEMBL",
                toType = "ENTREZID",
                OrgDb = org.Mm.eg.db)


# Convert Ensembl to Entrez
gene_df <- bitr(all_top_genes,
                fromType = "ENSEMBL",
                toType = "ENTREZID",
                OrgDb = org.Mm.eg.db)


ego <- enrichGO(
  gene = gene_df$ENTREZID,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",  # Can be "BP", "MF", or "CC"
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# Bar plot
barplot(ego, showCategory = 20, title = "GO Enrichment (BP)")

# Dot plot
dotplot(ego, showCategory = 20, title = "GO Enrichment (BP)")


# Save results to CSV
write.csv(as.data.frame(ego), "GO_enrichment_results.csv", row.names = FALSE)








  
  
#======================Generate MDS plot Sex:Timepoint==========================
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



  
  
  
  
  
  
  