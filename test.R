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


##===========Create Metadata====================================================

# Read your metadata
metadata <- read.csv("/Users/judyabuel/Desktop/Xist/circadian_atlas/wt_circadian_traits.csv")

# Ensure factors are properly set
metadata$Timepoint <- factor(metadata$Timepoint)
metadata$Sex <- factor(metadata$Sex, levels = c("Female", "Male")) # Female as reference


##===========Create a DGEList===================================================

d0 <- DGEList(counts)
d0 <- calcNormFactors(d0, method = "TMM")

# Filter lowly expressed genes (edgeR way)
keep <- rowSums(cpm(d0) >= 1) >= 3  # Keep genes with CPM ≥1 in at least 3 samples
d <- d0[keep, , keep.lib.sizes=FALSE]
dim(d) # Check how many number of genes left

# design matrix
design <- model.matrix(~ Sex + Timepoint + Sex:Timepoint, data = metadata)

# Proceed with voom and lmFit as before
v <- voom(d, design, plot = TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit)


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

# Set up timepoints
timepoints <- c("0", "3", "6", "9", "12", "15", "18", "21") # Adjust based on your actual timepoints

# Initialize list to store results
all_results <- list()

# Loop through each timepoint
for (tp in timepoints) {
  # Create contrast name
  contrast_name <- paste0("Female_vs_Male_ZT", tp)
  
  # Create the contrast dynamically
  if (tp == "0") {
    # Baseline comparison (just SexMale)
    contrast <- makeContrasts(SexMale, levels = design)
    colnames(contrast) <- contrast_name
  } else {
    # Timepoint-specific comparison (SexMale + interaction)
    contrast <- makeContrasts(
      paste0("SexMale + SexMale:Timepoint", tp),
      levels = design
    )
    colnames(contrast) <- contrast_name
  }
  
  # Fit contrasts
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  
  # Get results
  results <- topTable(fit2, coef = contrast_name, number = Inf, sort.by = "p")
  results$Timepoint <- paste0("ZT", tp)
  results$Contrast <- contrast_name
  
  # Add gene symbols
  results$gene_symbol <- mapIds(org.Mm.eg.db,
                                keys = rownames(results),
                                column = "SYMBOL",
                                keytype = "ENSEMBL",
                                multiVals = "first")
  
  # Store results
  all_results[[contrast_name]] <- results
  
  # Create and save volcano plot for this timepoint
  volcano_plot <- create_volcano_plot(results, paste0("ZT", tp))
  print(volcano_plot)
  
  # Optional: Save plot to file
  ggsave(paste0("volcano_plot_ZT", tp, ".png"), plot = volcano_plot, 
         width = 8, height = 6, dpi = 300)
}

# Combine all results into one dataframe
combined_results <- do.call(rbind, all_results)

# Save combined results
write.csv(combined_results, "female_vs_male_all_timepoints_results.csv", row.names = TRUE)


##====================Create Volcano Plot=======================================

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Mm.eg.db")  # Mouse genome
library(org.Mm.eg.db)

create_volcano_plot <- function(results, timepoint, title_suffix = "")
  # Define significance thresholds
  logFC_threshold <- 1
  pval_threshold <- 0.05
  
  # Prepare data for plotting
  plot_data <- results %>%
    mutate(
      significance = case_when(
        abs(logFC) > logFC_threshold & adj.P.Val < pval_threshold ~ "Significant (FDR < 0.05)",
        abs(logFC) > logFC_threshold & P.Value < pval_threshold ~ "Significant (P < 0.05)",
        TRUE ~ "Not significant"
      )
     
  
  # Create volcano plot
  create_volcano_plot <- function(results, timepoint, title_suffix = "")
    # Define significance thresholds
    logFC_threshold <- 1
    pval_threshold <- 0.05
    
    # Prepare data for plotting
    plot_data <- results %>%
      mutate(
        significance = case_when(
          abs(logFC) > logFC_threshold & adj.P.Val < pval_threshold ~ "Significant (FDR < 0.05)",
          abs(logFC) > logFC_threshold & P.Value < pval_threshold ~ "Significant (P < 0.05)",
          TRUE ~ "Not significant"
        ),
        gene_label = ifelse(
          (abs(logFC) > logFC_threshold & P.Value < pval_threshold) |
            gene_symbol %in% c("Xist", "Tsix", "Eif2s3y", "Ddx3y", "Uty", "Kdm5d"), # Highlight sex chromosome genes
          gene_symbol, ""
        )
      )
    
    # Create volcano plot
    volcano_plot <- ggplot(plot_data, aes(x = logFC, y = -log10(P.Value), color = significance)) +
      geom_point(alpha = 0.6, size = 2) +
      scale_color_manual(values = c(
        "Significant (FDR < 0.05)" = "red",
        "Significant (P < 0.05)" = "blue",
        "Not significant" = "gray60"
      )) +
      geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "grey") +
      geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "grey") +
      labs(
        x = "Log2 fold change (Female/Male)",
        y = "-Log10 p-value",
        title = paste("Female vs Male at", timepoint, title_suffix),
        color = "Significance"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        legend.position = "bottom",
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blan)
      






