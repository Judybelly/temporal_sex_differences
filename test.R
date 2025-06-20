## Load library
library(ggplot2)
library(edgeR)
library(limma)
library(dplyr)
library(factoextra)
library(FactoMineR)
library(gtools)


# Read a file in table format and assign a name to the variable: counts
counts <- 
  read.delim("/Users/judyabuel/Desktop/Xist/circadian_atlas/GSE297702_circadian_atlas_rawcounts.txt")

# Display the first few rows of the data frame
head(counts)


##===========Normalize Data==============
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Mm.eg.db")      # For gene annotations (symbols, Entrez IDs, etc.)
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene") # For genomic coordinates (optional)
BiocManager::install("AnnotationDbi")     # For mapping gene IDs

library(org.Mm.eg.db)
library(AnnotationDbi)
library(biomaRt)
library(edgeR)

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Create a DGEList, a special data structure in edgeR designed for storing RNA-seq count data.
# Calculate noramlization factors to adjust for differences in library sizes and RNA composition between samples.
d0 <- DGEList(counts)
d0 <- calcNormFactors(d0, method = "TMM") #TMM Trimmed Mean of M-values

annotations <- getBM(
  attributes = c(
    "ensembl_transcript_id",   # Ensembl Transcript ID
    "external_gene_name",      # Gene Symbol
    "entrezgene_id",           # Entrez ID
    "description"              # Gene description
  ),
  filters = "ensembl_transcript_id", # Filter by Ensembl Transcript ID
  values = gene_ids,           # Your transcript IDs
  mart = mart
)

# Merge with DGEList
d0$genes <- annotations[match(gene_ids, annotations$ensembl_transcript_id), ]

head(d0$genes)


# Filters low-expression gene on the dataset stored in DGEList object (d0)
cutoff <- 1   # sets the threshold for filtering genes with low expression
drop <- which(apply(cpm(d0),1,max) < cutoff)  # computes cpm for each gene in d0 and normalize the counts by library size
d <- d0[-drop,]   # creates a new DGEList excluding low-expression genes
dim(d) # Check how many number of genes left

# Get the genes that passed filtering
remaining_genes <- rownames(d$counts)

# View the first few genes
head(remaining_genes)

# Count how many genes remain
length(remaining_genes)

# Access the gene annotation data frame
gene_info <- d$genes

# View the first few rows
head(gene_info)

# To save the remaining genes to a file:
write.table(remaining_genes, file = "remaining_genes.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# Adding more information about your samples
# This prepares your data for your further analysis
traits <- 
  read.csv("/Users/judyabuel/Desktop/Xist/circadian_atlas/wt_circadian_traits.csv")

# View the first few rows
head(traits)


# Read in Annotation to link diff expression results to gene names (Not sure where I'd get this)
#anno <- read.delim("esembl_mm_.tsv", as.is=T)
#dim(anno)

##==========CREATE PLOTS=================
# Get sample names from counts data (replace 'counts' with your dataframe)
sample_names <- colnames(counts)  # e.g., "A1", "A2", "B1", "B2"

# Create a 'group' vector where A* = Male, B* = Female
group <- ifelse(grepl("^A", sample_names), "Male", "Female")

# Convert to factor (useful for plotting)
group <- factor(group, levels = c("Male", "Female"))

# Generate an MDS plot to visualize sample similarities/differences (e.g., gene expression)
# Open PNG device FIRST
png(
  filename = "limma_voom/newplotMDS_plot.png",
  width    = 2000,          
  height   = 1500,          
  res      = 300,           
  units    = "px"           
)

# Move gene IDs from column 1 to row names
rownames(counts) <- counts[, 1]    # Set row names = gene IDs
counts <- counts[, -1]  

# Convert counts to DGEList (edgeR/limma format)
dge <- DGEList(counts = counts)

# Normalize data (voom for RNA-seq)
v <- voom(dge, plot = FALSE)

# MDS plot colored by group
plotMDS(v, col = as.numeric(group), pch = 16, cex = 1.5)
legend("topright", legend = levels(group), col = 1:2, pch = 16)

# Close the device (saves the file)
dev.off()


##=========Combine Data=========
# I want to create a new data set that contains the first 3 columns from my traits data set.
# I want to 
library(dplyr)
new_dataset <- traits %>% select(1:3)



















