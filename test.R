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
design <- model.matrix(~ Sex + Timepoint + Sex:Timepoint, data = metadata)
colnames(design) # Checks Check if your design matrix matches expectations

# Proceed with voom and lmFit as before
v <- voom(d, design, plot = TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit)


#===============================================================================
#                              Define Contrasts
#===============================================================================
# Check fitted model coefficients
model_coefs <- colnames(coef(fit))
print(model_coefs)

# Check design matrix
print(colnames(design))

# Check metadata structure
print(levels(metadata$Sex))
print(levels(metadata$Timepoint))

# Clean setup - ensure proper factors
metadata$Timepoint <- factor(metadata$Timepoint)

# For mixed format data (some with ZT prefix, some without)
metadata$Timepoint <- factor(
  formatC(as.numeric(gsub("ZT", "", metadata$Timepoint)), width = 2, flag = "0"),
  levels = c("00", "03", "06", "09", "12", "15", "18", "21")
)

# Create design matrix with explicit contrasts
design <- model.matrix(
  ~ Sex + Timepoint + Sex:Timepoint,
  data = metadata,
  contrasts.arg = list(
    Sex = contr.treatment(2, base = 1, contrasts = FALSE),  # Keep both levels
    Timepoint = contr.treatment(8, base = 1, contrasts = FALSE)  # Keep all timepoints
  )
)

# Verify all timepoints are represented
design_cols <- colnames(design)
print(design_cols)
metadata$Sex <- factor(metadata$Sex, levels = c("Female", "Male")) # Female as reference

all_results <- list()

model_coefs <- rownames(fit$coefficients)



###This is the part that is not working for me.. I am trying to loop my timepoints using the contrast that I created.
for (tp in c("00", "03", "06", "09", "12", "15", "18", "21")) {
  contrast <- matrix(0, nrow = length(model_coefs), ncol = 1,
                     dimnames = list(model_coefs, paste0("ZT", tp)))
  
  if (tp == "00") {
    if (!(female_coef %in% model_coefs)) {
      stop(paste("female_coef", female_coef, "not found"))
    }
    contrast[female_coef, 1] <- 1
  } else {
    interaction_term <- paste0(female_coef, ":Timepoint", tp)
    if (interaction_term %in% model_coefs) {
      contrast[female_coef, 1] <- 1
      contrast[interaction_term, 1] <- 1
    } else {
      warning(paste("Skipping ZT", tp, "- missing interaction term"))
      next
    }
  }
  
  cat("\nContrast for ZT", tp, ":\n")
  print(contrast[contrast[,1] != 0, , drop = FALSE])
  
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  
  results <- topTable(fit2, coef = 1, number = Inf)
  results$Timepoint <- paste0("ZT", tp)
  all_results[[paste0("ZT", tp)]] <- results
}















#===============================================================================
#                              Generate Plots
#===============================================================================

for (tp_name in names(results)) {
  df <- results[[tp_name]]
  
  # Create significance categories
  df$Expression <- ifelse(
    df$adj.P.Val < 0.05 & df$logFC > 0.5, "Upregulated",
    ifelse(df$adj.P.Val < 0.05 & df$logFC < -0.5, "Downregulated", "Not significant")
  )
  
  # Get top 5 up/down genes for labeling
  top_up <- head(df[df$Expression == "Upregulated" & !is.na(df$gene_symbol), ], 5)
  top_down <- head(df[df$Expression == "Downregulated" & !is.na(df$gene_symbol), ], 5)
  
  p <- ggplot(df, aes(x = logFC, y = -log10(adj.P.Val), color = Expression)) +
    geom_point(aes(color = Expression), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Downregulated" = "blue", 
                                  "Upregulated" = "red", 
                                  "Not significant" = "grey")) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    ggtitle(paste("Female vs Male at", gsub("Female_vs_Male_", "", tp_name))) +
    labs(x = "log2 Fold Change", y = "-log10(Adjusted p-value)") +
    theme_minimal() +
    theme(legend.position = "right") +
    
    # Add labels for top genes
    geom_text_repel(
      data = top_up,
      aes(label = gene_symbol),
      color = "red",
      size = 3,
      box.padding = 0.5,
      max.overlaps = Inf
    ) +
    geom_text_repel(
      data = top_down,
      aes(label = gene_symbol),
      color = "blue",
      size = 3,
      box.padding = 0.5,
      max.overlaps = Inf
    )
  
  # Save plot
  ggsave(
    paste0("volcano_", tp_name, ".png"),
    plot = p,
    width = 10,  # Slightly wider to accommodate labels
    height = 8,
    dpi = 300
  )
  
  # Print plot to screen
  print(p)
}




library(ggplot2)
library(ggrepel)

generate_volcano_plot <- function(df, timepoint) {
  # Classify genes
  df$Expression <- ifelse(
    df$adj.P.Val < 0.05 & df$logFC > 0.5, "Upregulated",
    ifelse(df$adj.P.Val < 0.05 & df$logFC < -0.5, "Downregulated", "Not significant")
  )
  
  # Get top 5 up/downregulated genes for labeling
  top_up <- df %>%
    filter(Expression == "Upregulated", !is.na(gene_symbol)) %>%
    arrange(adj.P.Val) %>%
    head(5)
  
  top_down <- df %>%
    filter(Expression == "Downregulated", !is.na(gene_symbol)) %>%
    arrange(adj.P.Val) %>%
    head(5)
  
  # Generate plot
  p <- ggplot(df, aes(x = logFC, y = -log10(adj.P.Val), color = Expression)) +
    geom_point(aes(color = Expression), alpha = 0.6, size = 2) +
    scale_color_manual(
      values = c("Upregulated" = "red", "Downregulated" = "blue", "Not significant" = "grey50"),
      name = "Expression"
    ) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black", linewidth = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.5) +
    labs(
      title = paste("Female vs Male at", timepoint),
      x = "log2 Fold Change (Female / Male)",
      y = "-log10(Adjusted p-value)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    # Label top genes
    geom_text_repel(
      data = top_up,
      aes(label = gene_symbol),
      color = "red",
      size = 3.5,
      box.padding = 0.5,
      max.overlaps = 20
    ) +
    geom_text_repel(
      data = top_down,
      aes(label = gene_symbol),
      color = "blue",
      size = 3.5,
      box.padding = 0.5,
      max.overlaps = 20
    )
  
  return(p)
}


# Create plots for each timepoint
for (tp_name in names(all_results)) {
  # Extract timepoint (e.g., "ZT0" from "Female_vs_Male_ZT0")
  timepoint <- gsub("Female_vs_Male_", "", tp_name)
  
  # Generate plot
  volcano_plot <- generate_volcano_plot(all_results[[tp_name]], timepoint)
  
  # Save plot (high resolution, 300 dpi)
  ggsave(
    filename = paste0("volcano_plot_", timepoint, ".png"),
    plot = volcano_plot,
    width = 10,
    height = 8,
    dpi = 300,
    bg = "white"
  )
  
  # Print plot to RStudio
  print(volcano_plot)
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



  
  
  
  
  
  
  