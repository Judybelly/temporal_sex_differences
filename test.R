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
# Full edgeR + voom pipeline
d0 <- calcNormFactors(d0, method = "TMM")
design <- model.matrix(~ group)
v <- voom(d0, design, plot = TRUE)  # Generates voom plot

fit <- lmFit(v, design)
fit <- eBayes(fit)
topTable(fit, coef = 2)  # Extract DEGs


#========Compare Male vs Female at each timpoints========
# First, let's check your actual sample count
total_samples <- length(colnames(counts))  # Should be 92

# Create metadata with proper lengths
metadata <- data.frame(
  Sample = colnames(counts),
  Sex = rep(c("Male", "Female"), each = 46),
  Timepoint = c(
    rep(c(0, 3, 6, 9, 12, 15, 18, 21), length.out = 46),  # Male timepoints
    rep(c(0, 3, 6, 9, 12, 15, 18, 21), length.out = 46)   # Female timepoints
  ),
  row.names = colnames(counts)
)

# Verify
dim(metadata)  # Should show 92 rows, 3 columns
table(metadata$Sex, metadata$Timepoint)


# Create DGEList
dge <- DGEList(counts = counts, 
               group = interaction(metadata$Sex, metadata$Timepoint))

# Filter low-expressed genes (keep genes with ≥1 CPM in ≥3 samples)
keep <- filterByExpr(dge, group = dge$samples$group)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)  # TMM normalization

# Design matrix (account for Sex + Timepoint + Sex:Timepoint interaction)
# Create a more robust design matrix
design <- model.matrix(~ 0 + Sex + Timepoint + Sex:Timepoint, data = metadata)
colnames(design) <- gsub(":", ".", colnames(design))  # Replace : with .

# Fit the model
fit <- lmFit(v, design)
fit <- eBayes(fit)



##=========GENERATE PLOTS===================




