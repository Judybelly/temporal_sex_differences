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

##===========Import Data==================
# Read raw count file
counts <- 
  read.delim("/Users/judyabuel/Desktop/Xist/circadian_atlas/GSE297702_circadian_atlas_rawcounts.txt", row.names = 1)
# Display the first few rows of the data frame
head(counts)

##===========Normalize & Filter Data==============
# edgeR approach
d0 <- DGEList(counts)
d0 <- calcNormFactors(d0, method = "TMM")

# Filter lowly expressed genes (edgeR way)
keep <- rowSums(cpm(d0) >= 1) >= 3  # Keep genes with CPM ≥1 in at least 3 samples
d <- d0[keep, , keep.lib.sizes=FALSE]
dim(d) # Check how many number of genes left

##==========Create a group============
# Extract group labels from column names (e.g., "A1", "A10", etc.)
group <- ifelse(grepl("^A", colnames(d0$counts)), "Male", "Female")
group <- factor(group)  # Convert to factor


##==========FIT THE MODEL + Plot=================
plotMDS(d0, 
        col = as.numeric(group),  # Colors: 1=Male, 2=Female
        pch = 16,                 # Solid circles
        cex = 1.5,                 # Point size
        main = "MDS Plot: Male vs Female Samples",
        xlab = "Dimension 1",
        ylab = "Dimension 2")

# Add a legend to my plot
legend("topright", 
       legend = levels(group), 
       col = 1:2, 
       pch = 16, 
       title = "Sex")

# Convert counts to VOOM-transformed values
design <- model.matrix(~ group)  # Replace 'group' with your experimental factor
v <- voom(d0, design, plot = TRUE)  # The 'plot=TRUE' generates the voom plot

# Fit the model
fit <- glmQLFit(d0, design, robust = TRUE)

#==========================================
# 1. Create metadata 
metadata <- data.frame(
  sample = colnames(d0),
  sex = factor(ifelse(grepl("^A", colnames(d0)), "Male", "Female")),
  timepoint = factor(gsub("[A-Z]", "", colnames(d0))),
  stringsAsFactors = FALSE
)

# Add group column separately
metadata$group <- factor(paste0("group", metadata$sex, "_", metadata$timepoint))

# 2. Filter and normalize
keep <- filterByExpr(d0, group = metadata$group)
d0 <- d0[keep, , keep.lib.sizes = FALSE]
d0 <- calcNormFactors(d0, method = "TMM")

# 3. Create design matrix
design <- model.matrix(~0 + group, data = metadata)
colnames(design) <- gsub("group", "", colnames(design))  # Clean column names


# If you can combine some timepoints biologically:
metadata$period <- ifelse(metadata$timepoint %in% 1:10, "early", "late")
design <- model.matrix(~sex * period, data = metadata)

# MA plot for a timepoint
plot(rowMeans(logCPM), timepoint_results[["1"]]$logFC, 
     main = "MA Plot (Timepoint 1)",
     xlab = "Average expression", ylab = "logFC (Male/Female)")

##=========GENERATE PLOTS===================




