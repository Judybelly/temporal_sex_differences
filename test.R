#===============================================================================
# PURPOSE: DEG Analysis
# Circadian Data
#===============================================================================

# Load libraries
library(edgeR)
library(limma)
library(dplyr)
library(tidyverse)
library(biomaRt)
library(ggplot2)
library(EnhancedVolcano)
library(ggrepel)
library(tibble)

# Set working directory
setwd("/Users/judyabuel/Desktop/Xist/circadian_atlas")

#===============================================================================
# Load Data
#===============================================================================

# Load count matrix
counts <- read.delim("GSE297702_circadian_atlas_rawcounts.txt", row.names = 1)

# Load metadata
metadata <- read.csv("wt_circadian_traits.csv")
metadata$Sex <- factor(metadata$Sex)
metadata$Timepoint <- factor(metadata$Timepoint)

# Relevel to set Male and ZT0 as references
metadata$Sex <- relevel(metadata$Sex, ref = "Male")
metadata$Timepoint <- relevel(metadata$Timepoint, ref = "0")
metadata$Timepoint <- factor(as.character(metadata$Timepoint))  # Ensure proper factor levels

#===============================================================================
# Preprocessing
#===============================================================================

# Create DGEList object and normalize
dge <- DGEList(counts)
dge <- calcNormFactors(dge, method = "TMM")

# Filter low-expressed genes
keep <- rowSums(cpm(dge) >= 1) >= 3
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Create design matrix
design <- model.matrix(~ Sex * Timepoint, data = metadata)
colnames(design) <- make.names(colnames(design))  # Clean column names

# Apply voom transformation and fit model
v <- voom(dge, design, plot = TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit)



#===============================================================================
# Define Contrasts
#===============================================================================

contrast.matrix <- makeContrasts(
  Female_vs_Male_T3  = SexFemale + SexFemale.Timepoint3,
  Female_vs_Male_T6  = SexFemale + SexFemale.Timepoint6,
  Female_vs_Male_T9  = SexFemale + SexFemale.Timepoint9,
  Female_vs_Male_T12 = SexFemale + SexFemale.Timepoint12,
  Female_vs_Male_T15 = SexFemale + SexFemale.Timepoint15,
  Female_vs_Male_T18 = SexFemale + SexFemale.Timepoint18,
  Female_vs_Male_T21 = SexFemale + SexFemale.Timepoint21,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

#===============================================================================
# Annotate Genes with MGI Symbols
#===============================================================================

# Get clean Ensembl IDs from DGE object
ensembl_ids <- rownames(dge)
ensembl_ids_clean <- sub("\\..*", "", ensembl_ids)

# Connect to Ensembl
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids_clean,
  mart = ensembl
)

#===============================================================================
# Generate Volcano Plots for Each Timepoint
#===============================================================================

# Create output directories
output_dir <- "New_volcano_plots"
results_dir <- "deg_results"
dir.create(output_dir, showWarnings = FALSE)
dir.create(results_dir, showWarnings = FALSE)

# Load required libraries
library(ggplot2)
library(dplyr)
library(ggrepel)


# Loop over all contrasts
for (contrast_name in colnames(contrast.matrix)) {
  message("Processing ", contrast_name, "...")
  
  # Get DEG results
  top_table <- topTable(fit2, coef = contrast_name, number = Inf, sort.by = "P")
  top_table$ensembl_gene_id <- sub("\\..*", "", rownames(top_table))
  top_table_annotated <- merge(top_table, gene_map, by = "ensembl_gene_id", all.x = TRUE)
  
  # Add fallback label
  top_table_annotated$label <- ifelse(
    is.na(top_table_annotated$mgi_symbol) | top_table_annotated$mgi_symbol == "",
    top_table_annotated$ensembl_gene_id,
    top_table_annotated$mgi_symbol
  )
  
  # Classify expression status
  top_table_annotated$Expression <- case_when(
    top_table_annotated$adj.P.Val < 0.05 & top_table_annotated$logFC > 1  ~ "Upregulated",
    top_table_annotated$adj.P.Val < 0.05 & top_table_annotated$logFC < -1 ~ "Downregulated",
    TRUE ~ "Not Significant"
  )
  
  # Save full DEG results
  write.csv(
    top_table_annotated,
    file = file.path(results_dir, paste0(contrast_name, "_DEG_results.csv")),
    row.names = FALSE
  )
  
  # Label ALL significant up/downregulated genes
  label_df <- top_table_annotated %>%
    filter(adj.P.Val < 0.05 & (logFC > 1 | logFC < -1)) %>%
    arrange(adj.P.Val) %>%
    head(30)
  
  
  # Plot
  p <- ggplot(top_table_annotated, aes(x = logFC, y = -log10(adj.P.Val), color = Expression)) +
    geom_point(alpha = 0.8, size = 1.5) +
    scale_color_manual(values = c(
      "Upregulated" = "red",
      "Downregulated" = "blue",
      "Not Significant" = "grey"
    )) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    labs(
      title = paste("Volcano Plot:", contrast_name),
      x = "Log2 Fold Change",
      y = expression(-log[10](adjusted~P~value)),
      color = "Expression"
    ) +
    theme_minimal(base_size = 14) +
    geom_text_repel(
      data = label_df,
      aes(label = label),
      size = 3,
      max.overlaps = 100,
      box.padding = 0.4,
      point.padding = 0.3,
      segment.color = "black"
    )
  
  # Save plot
  png(file.path(output_dir, paste0(contrast_name, "_volcano.png")), width = 2000, height = 1800, res = 300)
  print(p)
  dev.off()
}




#===============================================================================
# Generate Heatmap for Each Timepoint
#===============================================================================

library(pheatmap)
library(dplyr)
library(RColorBrewer)

install.packages("RColorBrewer")

# Load DEG CSV and gene_map
deg_all <- read.csv("DEG_Female_vs_Male_by_Timepoint.csv")
deg_all <- deg_all %>%
  mutate(ensembl_gene_id = sub(".*\\.", "", X))  # Extract Ensembl ID

# Join with gene_map to add gene symbols
deg_annot <- left_join(deg_all, gene_map, by = "ensembl_gene_id")

# Use gene symbol for labels; fall back to Ensembl if missing
deg_annot <- deg_annot %>%
  mutate(label = ifelse(is.na(mgi_symbol) | mgi_symbol == "", ensembl_gene_id, mgi_symbol))

# Get top 15 DEGs per timepoint (adj.P.Val < 0.05)
deg_top <- deg_annot %>%
  filter(adj.P.Val < 0.05) %>%
  group_by(Timepoint) %>%
  arrange(adj.P.Val, .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

# Get unique top genes
top_genes <- unique(deg_top$ensembl_gene_id)

# Subset voom-normalized expression matrix
expr_top <- v$E[rownames(v$E) %in% top_genes, ]

# Match rownames to gene symbols
gene_labels <- deg_top %>%
  dplyr::select(ensembl_gene_id, label) %>%
  distinct()


rownames(expr_top) <- gene_labels$label[match(rownames(expr_top), gene_labels$ensembl_gene_id)]

# Optional: scale per gene (z-score)
expr_scaled <- t(scale(t(expr_top)))

# Annotation for columns (samples)
ann_colors <- list(Sex = c(Female = "firebrick", Male = "steelblue"))
annot_col <- metadata %>%
  dplyr::select(Sample_ID, Sex, Timepoint) %>%
  filter(Sample_ID %in% colnames(expr_scaled)) %>%
  tibble::column_to_rownames("Sample_ID")


# Save heatmap
output_dir <- "heatmaps"
dir.create(output_dir, showWarnings = FALSE)

heatmap_file <- file.path(output_dir, "Heatmap_Top10_DEGs_AllTimepoints.png")
print(paste("Saving heatmap to:", heatmap_file))

# Plot heatmap
png("heatmaps/Heatmap_DEGs_AllTimepoints.png", width = 2200, height = 1600, res = 300)
pheatmap(
  expr_scaled,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = annot_col,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("navy", "white", "firebrick"))(100),
  fontsize_row = 6,
  fontsize_col = 6,
  angle_col = 45,
  main = "Top DEGs at Each Timepoint: Female vs Male",
  fontsize = 10
)

dev.off()



#===============================================================================
## HEATMAP: Top 15 DE Genes per Timepoint with Sex-specific Significance
#===============================================================================

library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(tibble)

# Load DEG CSV and gene_map
deg_all <- read.csv("DEG_Female_vs_Male_by_Timepoint.csv")
deg_all <- deg_all %>%
  mutate(ensembl_gene_id = sub(".*\\.", "", X))  # Extract Ensembl ID

# Join with gene_map to add gene symbols
deg_annot <- left_join(deg_all, gene_map, by = "ensembl_gene_id")

# Use gene symbol for labels; fall back to Ensembl ID if missing
deg_annot <- deg_annot %>%
  mutate(label = ifelse(is.na(mgi_symbol) | mgi_symbol == "", ensembl_gene_id, mgi_symbol))

# Get top 10 DEGs per timepoint (adj.P.Val < 0.05)
deg_top <- deg_annot %>%
  filter(adj.P.Val < 0.05) %>%
  group_by(Timepoint) %>%
  arrange(adj.P.Val, .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

# Get unique top genes
top_genes <- unique(deg_top$ensembl_gene_id)

# Subset voom-normalized expression matrix
expr_top <- v$E[rownames(v$E) %in% top_genes, ]

# Match rownames to gene symbols
gene_labels <- deg_top %>%
  dplyr::select(ensembl_gene_id, label) %>%
  distinct()
rownames(expr_top) <- gene_labels$label[match(rownames(expr_top), gene_labels$ensembl_gene_id)]

# Optional: scale per gene (z-score)
expr_scaled <- t(scale(t(expr_top)))

# Column annotation using metadata
annot_col <- metadata %>%
  dplyr::select(Sample_ID, Sex, Timepoint) %>%
  filter(Sample_ID %in% colnames(expr_scaled)) %>%
  tibble::column_to_rownames("Sample_ID")

# Update colnames to "Timepoint_SampleID" for cleaner display
colnames(expr_scaled) <- paste0(annot_col[colnames(expr_scaled), "Timepoint"], "_", colnames(expr_scaled))

# Define colors for annotation
ann_colors <- list(
  Sex = c(Female = "firebrick", Male = "steelblue"),
  Timepoint = setNames(RColorBrewer::brewer.pal(length(unique(annot_col$Timepoint)), "Set3"),
                       unique(annot_col$Timepoint))
)

# Match updated colnames back to annot_col
annot_col <- annot_col[colnames(expr_scaled), , drop = FALSE]
rownames(annot_col) <- colnames(expr_scaled)  # match to expr_scaled

# Create output directory
output_dir <- "Heatmap"
dir.create(output_dir, showWarnings = FALSE)

# Save heatmap
png(file.path(output_dir, "Heatmap_Top10_DEGs_AllTimepoints.png"), width = 1600, height = 1200, res = 300)

pheatmap(
  expr_scaled,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annot_col,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("navy", "white", "firebrick"))(100),
  fontsize_row = 9,
  fontsize_col = 6,
  angle_col = 45,
  main = "Top DEGs per Timepoint: Female vs Male"
)

dev.off()


#===============================================================================
#              Line graph of the Top 10 genes across timepoints
#===============================================================================

# Load required packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggforce)  # For pagination

# Read data and select genes
data <- read.csv("Heatmap_Top_DEGs_logFC_matrix.csv", row.names = 1)
specific_genes <- c("Uty", "Ddx3y", "Kdm5d", "Eif2s3y", "Xist", "Eif2s3x", "Kdm6a", "Glp2r")  # Your genes

# Reshape data
data_long <- data %>%
  filter(rownames(.) %in% specific_genes) %>%
  mutate(Gene = factor(rownames(.))) %>%  # Ensure proper ordering
  pivot_longer(cols = -Gene, names_to = "Timepoint", values_to = "logFC")

# Order timepoints
data_long$Timepoint <- factor(data_long$Timepoint, levels = paste0("ZT", sprintf("%02d", seq(0, 21, 3))))

# Create a function to plot in batches
plot_genes_paginated <- function(data, ncol = 3, nrow = 2, per_page = ncol * nrow) {
  n_pages <- ceiling(length(unique(data$Gene)) / per_page)
  
  for (i in 1:n_pages) {
    p <- ggplot(data, aes(x = Timepoint, y = logFC, group = Gene, color = Gene)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      ggforce::facet_wrap_paginate(~ Gene, ncol = ncol, nrow = nrow, page = i) +
      scale_y_continuous(limits = c(-10, 15), breaks = seq(-10, 15, by = 5)) +
      labs(title = paste("Changes in Gene Expression per Timepoint", i), x = "Timepoint", y = "log2 Fold Change") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            strip.background = element_rect(fill = "lightgrey"))
    
    ggsave(
      filename = paste0("Selected_Genes_Page_", i, ".pdf"),
      plot = p,
      device = "pdf",
      width = 10,
      height = 7,
      dpi = 300
    )
  }
}

# Run the function (adjust nrow/ncol as needed)
plot_genes_paginated(data_long, ncol = 3, nrow = 2)  # 6 genes per page (3 columns x 2 rows)


#===============================================================================
#                         Generate UpSet Plot
#===============================================================================

# Load libraries
install.packages("ComplexUpset")
install.packages("ggplot2")  # needed for plotting

library(ComplexUpset)
library(ggplot2)


# Define timepoints
all_tp <- c("ZT00", "ZT03", "ZT06", "ZT09", "ZT12", "ZT15", "ZT18", "ZT21")

# Basic UpSet plot
upset_plot <- upset(
  deg_df,
  intersect = all_tp,
  name = "DEG Timepoints",
  base_annotations = list(
    'Intersection size' = intersection_size(
      text = list(size = 3)
    )
  ),
  set_sizes = upset_set_size(),
  sort_intersections_by = "cardinality",
  n_intersections = 15,
  min_size = 1,
  width_ratio = 0.3,
  themes = upset_modify_themes(list(
    'intersect_size' = theme(text = element_text(size = 10)),
    'overall_sizes' = theme(axis.text.y = element_text(size = 10)),
    'intersections_matrix' = theme(axis.text.x = element_text(size = 8, angle = 90))
  ))
) +
  labs(
    title = "Top DEG Intersections Across Circadian Timepoints",
    subtitle = NULL
  )

# Save to PNG
ggsave("UpSet_DEGs_Top15_Basic.png", plot = upset_plot, width = 10, height = 6, dpi = 300)






#===============================================================================
#                           Generate GO Analysis
#===============================================================================

# Load required packages
library(clusterProfiler)
library(org.Mm.eg.db) # For mouse genes (change if different organism)
library(enrichplot)
library(ggplot2)
library(dplyr)

# 1. Prepare significant genes from your DEG data
sig_genes <- deg_annot %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1) %>% # Adjust thresholds as needed
  pull(ensembl_gene_id) %>%
  unique()

# 2. Convert ENSEMBL IDs to ENTREZID
gene_list <- bitr(sig_genes, 
                  fromType = "ENSEMBL",
                  toType = "ENTREZID",
                  OrgDb = org.Mm.eg.db)

# 3. Perform KEGG enrichment
kegg_enrich <- enrichKEGG(
  gene = gene_list$ENTREZID,
  organism = "mmu", # "mmu" for mouse, "hsa" for human
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2
)

# Simplify redundant terms
kegg_enrich <- clusterProfiler::simplify(kegg_enrich)

# 4. Visualize results (multiple options - choose your preferred one)

# Option A: Dot plot
dot_plot <- dotplot(kegg_enrich, showCategory = 20) +
  ggtitle("KEGG Pathway Enrichment") +
  theme(axis.text.y = element_text(size = 8))

# Option B: Bar plot
bar_plot <- barplot(kegg_enrich, showCategory = 15) +
  ggtitle("KEGG Pathway Enrichment") +
  xlab("Gene Count") +
  theme(axis.text.y = element_text(size = 8))

# Option C: Enrichment Map
enrich_map <- emapplot(pairwise_termsim(kegg_enrich), showCategory = 20)

# Option D: Pathway-specific visualization (example)
# First check available pathways:
View(kegg_enrich@result)
# Then visualize specific pathway (example with mmu05200 - Cancer pathways)
# browseKEGG(kegg_enrich, "mmu05200") # Uncomment to view in browser

# 5. Save plots
ggsave("KEGG_dotplot.png", plot = dot_plot, width = 10, height = 8, dpi = 300)
ggsave("KEGG_barplot.png", plot = bar_plot, width = 10, height = 8, dpi = 300)
ggsave("KEGG_enrichmap.png", plot = enrich_map, width = 12, height = 10, dpi = 300)

# Show one of the plots
print(dot_plot)








#===============================================================================
#                 Generate MDS plot Sex:Timepoint
#===============================================================================

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

