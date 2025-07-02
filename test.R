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

# 1. Connect to the Ensembl BioMart for mouse
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# 2. Extract all unique Ensembl gene IDs from your results
all_gene_ids <- unique(unlist(lapply(all_results, rownames)))

# 3. Query biomart to get gene symbols
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = all_gene_ids,  # ⬅️ use this instead
  mart = ensembl
)

# Ensure character types
gene_map$ensembl_gene_id <- as.character(gene_map$ensembl_gene_id)
gene_map$external_gene_name <- as.character(gene_map$external_gene_name)

# Join gene names to each df in all_results
all_results <- lapply(all_results, function(df) {

  library(dplyr)
  library(tibble)
  
  # Pick one dataset from all_results
  df <- all_results[[1]]
  
  # Add gene column
  df <- tibble::rownames_to_column(df, var = "gene")
  
  # Strip version numbers (if any)
  df$gene <- sub("\\..*$", "", df$gene)
  
  # Ensure types match
  df$gene <- as.character(df$gene)
  gene_map$ensembl_gene_id <- as.character(gene_map$ensembl_gene_id)
  
  # Do the join and inspect it
  df_joined <- left_join(df, gene_map, by = c("gene" = "ensembl_gene_id"))
  
  # Check result
  print(head(df_joined[, c("gene", "external_gene_name")], 10))
  




# Enhanced Volcano Plot with Gene Labels
  make_volcano_plot <- function(df, timepoint, logFC_cutoff = 1, adjP_cutoff = 0.05, max_labels = 20) {
    
    df <- df %>%
      mutate(
        Significance = case_when(
          adj.P.Val < adjP_cutoff & logFC > logFC_cutoff  ~ "Upregulated",
          adj.P.Val < adjP_cutoff & logFC < -logFC_cutoff ~ "Downregulated",
          TRUE                                             ~ "Not Significant"
        ),
        # Create label: use gene symbol if available, else fall back to Ensembl ID
        gene_label = ifelse(is.na(gene_symbol) | gene_symbol == "", gene, gene_symbol)
      )
    
    # Select top DEGs to label
    label_df <- df %>%
      filter(Significance != "Not Significant") %>%
      arrange(adj.P.Val) %>%
      head(max_labels)
    
    # Make plot
    p <- ggplot(df, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
      geom_point(alpha = 0.6) +
      scale_color_manual(values = c(
        "Upregulated" = "firebrick",
        "Downregulated" = "royalblue",
        "Not Significant" = "grey70"
      )) +
      geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed", color = "darkgrey") +
      geom_hline(yintercept = -log10(adjP_cutoff), linetype = "dashed", color = "darkgrey") +
      geom_text_repel(data = label_df, aes(label = gene_label), size = 3, max.overlaps = Inf) +
      theme_minimal() +
      labs(
        title = paste("Volcano Plot -", timepoint),
        x = "log2 Fold Change",
        y = "-log10 Adjusted p-value"
      )
    
    return(p)
  }
  
  for (tp in names(all_results)) {
    df <- all_results[[tp]]
    dir.create("volcano_plots", showWarnings = FALSE)
    p <- make_volcano_plot(df, tp)
    
    ggsave(
      filename = paste0("volcano_plots/volcano_", tp, ".png"),
      plot = p,
      width = 8,
      height = 6
    )
  }
  
  

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



  
  
  
  
  
  
  