#### Diffretential Transcript Expression Analysis : Seizure Type
## Edited by RASHI VERMA
## Date: 07/23/2024
## REFS https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html

# Load Libraries
library(statmod)
library(edgeR)
library(limma)
library(pheatmap)
library(ggplot2)

# set Seed 
set.seed(99)

# Upload Count and Pheno Data
countsTable <- read.csv(file='transcript_count_matrix.csv', row.names="transcript_id")
head(countsTable)
targets <- read.csv(file= 'pheno_data.csv', row.names="sample")
head(targets)

all(colnames(countsTable) %in% rownames(targets))
all(colnames(countsTable) == rownames(targets))

sex <- as.factor(targets$sex)
age <- as.numeric(targets$age) # Assuming age is numeric
order <- as.numeric(targets$order)
race <- as.factor(targets$race)
sampleID <- as.factor(targets$subject)
type <- as.factor(targets$seizure_type)

# Create DGEList object
d0 <- DGEList(countsTable)
class(d0)
dim(d0)
d0

# Set Factors
group <- factor(paste(targets$type,targets$order,sep="."))
cbind(targets,Group=group)
d0$samples$group <- group

# Calculate Normalization Factors
d1 <- calcNormFactors(d0, method = "TMM")

# Filter Low-expressed Genes
cutoff <- 1
drop <- which(apply(cpm(d1), 1, max) < cutoff)
d2 <- d1[-drop,] 
dim(d2)

# Calculate Library Size
size <- d2$samples$norm.factors
max <- max(size)
min <- min(size)
ratio <- max/min
ratio

# Make Design
design <- model.matrix(~ 0+ group + sex + age + race)
design

# use this chunk when ratio of largest and smallest library size more than 3 fold
v <- voom(d2, design, plot=TRUE, normalize="quantile")
cor <- duplicateCorrelation(v, design, block=sampleID)
cor$consensus.correlation
vfit <- lmFit(object = v, design=design, block=sampleID, correlation=cor$consensus.correlation)

# Contrast Matrix
con.YESvsNO <- makeContrasts(
  YvsN = (group2.1 - group1.1),
  levels=design)

con.YESvsNO

# Estimate Contrast
vfit <- contrasts.fit(vfit, contrasts=con.YESvsNO)
efit <- eBayes(vfit)
plotSA(efit)
summary(decideTests(efit))

# Differentially expressed Transcripts
top.table<- topTable(efit, sort.by = "P", n = Inf, adjust.method = "BH")
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05))
length(which(top.table$logFC > 1))

# Filter and Write top.table to a file with adj.P.Val < 0.05 & Log2FC
filtered_genes <- top.table[top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(2), ]
output_file <- "Base_PNES_vs_Seizure.csv"
write.csv(filtered_genes, file = output_file, row.names = TRUE)

############
# pheatmap
new_column_name <- "ID"
filtered_genes[, new_column_name] <- rownames(filtered_genes)
str(filtered_genes)
top_genes <- head(filtered_genes, 203)
gene_names <- top_genes$ID
top_gene_expression <- cpm(d2)[gene_names, ]
baseline_columns <- grep("_Base$", colnames(cpm(d2)))
top_gene_expression_baseline <- top_gene_expression[, baseline_columns]
library(pheatmap)
pheatmap(top_gene_expression_baseline,
         scale = "row",        # Scale rows (genes)
         col = colorRampPalette(c("blue", "white", "red"))(100),  # Color scheme
         breaks = seq(-2, 2, length.out = 101),
         margins = c(5, 10),   # Adjust margins for labels
         main = "Heatmap of Top Genes",
         xlab = "Sample ID",
         ylab = "Genes",
         cex.axis = 1,
         #cluster_rows = FALSE,
         fontsize_row = 8)

############
# volcano plot
library(ggplot2)

# Your existing code ...

# Create a column in top.table to specify point color
top.table$Color <- ifelse(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > 1, "orange",
                          ifelse(top.table$logFC < -1, "red",
                                 ifelse(top.table$logFC > 1, "blue", "black")))

# Create a volcano plot
ggplot(top.table, aes(x = logFC, y = -log10(adj.P.Val), color = Color)) +
  geom_point(size = 1) + coord_cartesian(xlim = c(-10,10), ylim = c(0,6)) +
  scale_color_manual(values = c("red" = "red", "blue" = "blue", "orange" = "orange", "black" = "black")) +
  theme_minimal() +
  labs(x = "log2(FC)", y = "-log10(Adj P-Value)") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "green") +
  geom_vline(xintercept = c(-1, 1), linetype = "solid", color = "green") +  # Add green vertical lines
  ggtitle("Baseline") +
  theme(legend.position = "none", text=element_text(size=32,  family="TT Arial"))  # Set the background to white

############
# PCA Plot
filtered_counts <- countsTable[rownames(filtered_genes), ]

# Perform PCA
pca <- prcomp(t(filtered_counts), scale. = TRUE)

# Create a PCA plot
library(ggplot2)
library(scales)  # Load the scales package

# Create a dataframe with the PCA results
pca_data <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2])

# Add the "change" variable as a factor (assuming it's in 'targets' dataframe)
pca_data$Change <- factor(targets$type)

# Create the PCA plot with colors based on "change"
baseline_rows <- grep("_Base$", rownames(pca_data))
pca_data_baseline <- pca_data[baseline_rows,]
ggplot(pca_data_baseline, aes(x = PC1, y = PC2, color = Change)) +
  geom_point(size = 4, stroke = 1) +  # Increase the size of points
  labs(title = "PCA Plot for Filtered Genes") +
  scale_color_manual(values = c("2" = "aquamarine4", "3" = "deeppink4"),
                     labels = c("Focal", "Generalized")) +
  theme_minimal() +
  xlab(paste("PC1 (", percent(pca$sdev[1]^2 / sum(pca$sdev^2)), " variance)")) +
  ylab(paste("PC2 (", percent(pca$sdev[2]^2 / sum(pca$sdev^2)), " variance)"))
  
############################################




