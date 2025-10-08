####Set Working Directory####
setwd("./GSE37448")

####Install and Load Packages####
######General Bioconductor packages#####
library(Biobase)
library(oligoClasses)
library(affy)

#####Annotation and data import packages#####
library(ArrayExpress)
library(pd.mogene.1.0.st.v1)
library(mogene10sttranscriptcluster.db)
library(GEOquery)
library(biomaRt)

#####Quality control and pre-processing packages#####
library(oligo)
library(arrayQualityMetrics)
library(affyPLM)

#####Analysis and statistics packages#####
library(limma)
library(topGO)
library(ggkegg) # required to run ReactomePA
library(ReactomePA)
library(clusterProfiler)

#####Plotting and color options packages#####
library(gplots)
library(ggplot2)
library(geneplotter)
library(RColorBrewer)
library(pheatmap)
library(enrichplot)
library(ggrepel)

#####Formatting/documentation packages#####
library(BiocStyle)
library(dplyr)
library(tidyr)
install.packages("knitr")

#####Helpers#####
library(stringr)
library(matrixStats)
library(genefilter)
library(openxlsx)
library(devtools)

#### Series Matrix and Analysis ####
# Following : https://gtk-teaching.github.io/Microarrays-R/04-MetaData/index.html
##### Retrieving as expression set - will download series matrix file ####
gsm <- getGEO('GSE37448')
length(gsm) #  1
class(gsm[[1]]) # "ExpressionSet" attr(,"package") "Biobase"
gse37448 <- gsm[[1]]

varLabels(gse37448)

gse37448$supplementary_file

pd <- pData(gse37448)

##### Parse through series matrix but does not parse through the soft file ####
gsm2 <- getGEO('GSE37448',GSEMatrix=FALSE)
names(GSMList(gsm2))

gsm2@gsms$GSM920648@dataTable@table

##### Retrieve raw data ####
filePaths <- getGEOSuppFiles('GSE37448')

##### Reading CEL data ####
  # Download the RAW file on command line -  curl -O ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE37nnn/GSE37448/suppl/GSE37448_RAW.tar
  # Unzip the RAW file - tar -xvf GSE37448_RAW.tar

cel_dir <- "./GSE37448/GSE37448_CEL"

cel_files <- list.celfiles(cel_dir, full.names = TRUE)
length(cel_files)
head(cel_files) 

gse37448_celdata <- read.celfiles(cel_files)

##### Prepare phenotype data so that it correctly points to the actual CEL files #####
# Extract file names
pd$cel_file <- basename(pd$supplementary_file)
head(pd$cel_file)

# Check files exist
all(file.exists(file.path(cel_dir, pd$cel_file)))

# Read CEL file with phenotype data
gse37448_celdata<- read.celfiles(
  file.path(cel_dir, pd$cel_file),
  phenoData = phenoData(gse37448)
)

# Will get warning because some extra information about detection channels is added

# Access information on the data sets
pData(gse37448_celdata)[,c("geo_accession", "cell type:ch1", "tissue:ch1")]


##### RMA on raw data ####
# Summarisation - reduces the size of the data for each sample to the number of measured transcripts
oligo_summarised <- rma(gse37448_celdata,background=FALSE,normalize=FALSE)
hist(oligo_summarised)


# Quantile normalization puts the data on a common empirical distribution
gse37448_eset <- rma(gse37448_celdata)

# Plot PCA
expr_matrix <- exprs(gse37448_eset)

PCA_eset <- prcomp(t(expr_matrix), scale = FALSE)

percentVar <- round(100*PCA_eset$sdev^2/sum(PCA_eset$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA_eset$x[,1], PC2 = PCA_eset$x[,2],
                     Source = 
                       Biobase::pData(gse37448_eset)$source_name_ch1)

my_colors <- c(
  "red", "blue", "green", "orange", "purple", "pink", "brown",
  "cyan", "magenta", "yellow", "grey", "darkgreen", "gold",
  "tan", "plum4", "gray19", "burlywood4"
)

ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(colour = Source), size = 3) +
  ggtitle("PCA plot of the calibrated, summarized data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio) +
  scale_color_manual(values = my_colors)

##### Statistical Tests####
# https://dputhier.github.io/ASG/practicals/microarrays_student_test/DenBoer_Student_test.html#student_test_and_p-value

###### Mean Expression ####

# calculate mean expression level
mean.per.gene <- apply(expr_matrix, 1, mean, na.rm=TRUE)
head(mean.per.gene)

g1.columns <- which(cell_type == "abT Cell")
g2.columns <- which(cell_type == "Macrophage, mixed background")

g1.means <- apply(expr_matrix[, g1.columns], 1, mean, na.rm = TRUE)
g2.means <- apply(expr_matrix[, g2.columns], 1, mean, na.rm = TRUE)

logFC_means <- log2(g1.means + 1) - log2(g2.means + 1) 

gene_means <- data.frame(
  Mean_per_gene = mean.per.gene,
  abT_Cell_Mean = g1.means,
  Macrophage_mixed_Mean = g2.means,
  log2FC = logFC_means
)

  # inspect top fold change
head(gene_means[order(abs(gene_means$log2FC), decreasing = TRUE), ])
kable(head(gene_means),
      caption="First rows of the result table, used to collect different statistics per gene.",
      digits=2)

  # plot mean comparsion plot
plot(gene_means$abT_Cell_Mean, gene_means$Macrophage_mixed_Mean)

g1 <- "abT_Cell_Mean"
g2 <- "Macrophage_mixed_Mean"
axis.range <- c(0, ceiling(max(gene_means[,c("abT_Cell_Mean","Macrophage_mixed_Mean")])))
plot(gene_means$abT_Cell_Mean, 
     gene_means$Macrophage_mixed_Mean,
     main="Comparison between abT and macrophage from mixed background means",
     xlab=paste(g1, "mean (RMA-norm log2 values)"),
     ylab=paste(g2, "mean (RMA-norm log2 values)"),
     xlim=axis.range, ylim=axis.range,
     pch=20,
     col = densCols(x=gene_means$abT_Cell_Mean,
                    y=gene_means$Macrophage_mixed_Mean)
)
grid(lty="solid", col="darkgray")

abline(h=0, col="black")
abline(v=0, col="black")
abline(a=0, b=1)
abline(a=1, b=1, col="red")
abline(a=-1, b=1, col="red")


###### Variance and Standard Deviation ####
# Variance
n <- ncol(expr_matrix)
gene_means$variance<- (n-1)/n *  apply(expr_matrix, 1, var, na.rm=TRUE)

n1 <- length(g1.columns) 
gene_means$variance_abT  <- (n1-1)/n1*apply(expr_matrix[,g1.columns], 1, var, na.rm=TRUE)

n2 <- length(g2.columns)
gene_means$variance_mf  <- (n2-1)/n2*apply(expr_matrix[,g2.columns], 1, var, na.rm=TRUE)

# Standard Deviation
gene_means$sd <- sqrt(gene_means$variance)

gene_means$sd_abT <- sqrt(gene_means$variance_abT)
gene_means$sd_mf <- sqrt(gene_means$variance_mf)


  # plot
# VARIANCE 
max.var <- max(gene_means[, c("variance_abT","variance_mf")])
var.range <- c(0, ceiling(max.var))
plot(gene_means$variance_abT,
     gene_means$variance_mf,
     main="Group-wise variances",
     xlab=paste("Sample variance for", g1),
     ylab=paste("Sample variance for", g2),
     xlim=var.range, ylim=var.range, 
     col=densCols(gene_means$variance_abT, gene_means$variance_mf)
)
grid(lty="solid", col="grey")
abline(col="black", a=0, b=1) 

# DEVIATION
max.sd <- max(gene_means[, c("sd_abT","sd_mf")])
sd.range <- c(0, ceiling(max.sd))
plot(gene_means$sd_abT,
     gene_means$sd_mf,
     main="Group-wise standard dev.",
     xlab=paste(g1, "sample standard deviation"),
     ylab=paste(g2, "sample standard deviation"),
     xlim=sd.range, ylim=sd.range, 
     col=densCols(gene_means$sd_abT, gene_means$sd_mf)
)
grid(lty="solid", col="grey")
abline(col="black", a=0, b=1) 

###### Mean vs Variance ####
plot(gene_means$Mean_per_gene, gene_means$sd, 
     main=paste("All samples: mean vs var"),
     xlab=paste("Sample mean"),
     ylab=paste("Sample standard deviation"),
     col=densCols(gene_means$Mean_per_gene, gene_means$sd)
)
grid(lty="solid", col="gray")

###### Effect Size - MA Plot ####
 # Effect Size = difference of the two mean grps = D (log2fc)
gene_means$D <- gene_means$abT_Cell_Mean - gene_means$Macrophage_mixed_Mean

  # Average expression across two grps = A
gene_means$A <- (gene_means$abT_Cell_Mean + gene_means$Macrophage_mixed_Mean)/2

plot(gene_means$A,
     gene_means$D,
     main="MA plot",
     xlab= "Average expression: A = (m1 + m2)/2",
     ylab = "Effect size: d = m2 - m1",
     col=densCols(gene_means$D, gene_means$A)
)
grid(lty="solid", col="grey")
abline(col="black", h=0)
abline(col="red", h=c(1, -1))

###### Student T-Test ####
gene_means$var.pooled <- (n1* gene_means$variance_abT + n2* gene_means$variance_mf) / (n1+n2-2)
gene_means$sd.pooled <- sqrt(gene_means$var.pooled)

gene_means$t.Student <- gene_means$D / (gene_means$sd.pooled * sqrt(1/n1+1/n2))

  # Plot t-test vs effect size
plot(gene_means$D, xlab="Effect size (d = m2 - m1)",
     gene_means$t.Student, ylab="Student t statistics",
     main = "Student t statistics vs. effect size",
     col=densCols(gene_means$D, gene_means$t.Student)
)
grid(col="grey", lty="solid")
abline(v=0, col="black") 
abline(h=0, col="black") 
abline(v=c(-1,1), col="red") 

  # up and down regulated genes
gene_means$sign <- "null"
gene_means$sign[gene_means$D > 0] <- "up"
gene_means$sign[gene_means$D < 0] <- "down"

print(table(gene_means$sign))

  # p value
gene_means$p.value <- 2*pt(abs(gene_means$t.Student), 
                             df=n1 + n2 - 2, lower.tail = FALSE)

  # plot pooled deviation vs p value 
plot(gene_means$D, xlab="Effect size (d = m2 - m1)",
     gene_means$p.value, ylab="Student p-value",
     main = "Student p-value vs. effect size",
     col=densCols(gene_means$D, gene_means$p.value)
)
grid(col="grey", lty="solid")
abline(v=0, col="black") 
abline(v=c(-1,1), col="red") 
abline(h=0.01, col="red") # 5% p value threshold 

  # volcano plot - log10(p-value) as a function of the effect size
plot(gene_means$D, xlab="Effect size (d = m2 - m1)",
     -log10(gene_means$p.value), ylab="-log10(p-value)",
     main = "Volcano plot",
     col=densCols(gene_means$D, -log10(gene_means$p.value))
)
grid(col="grey", lty="solid")
abline(v=0, col="black") 
abline(v=c(-1,1), col="red") 
abline(h=-log10(0.05), col="red") # 5% p value threshold = Genes above this line are declared significant

hist(gene_means$p.value, breaks=20, main="P-value histogram", xlab="Pvalue", ylab="Number of probesets", col="grey")

  # expected number of false positives
alpha <- 0.05 
n.positive <- sum(gene_means$p.value < alpha)
# number of probes
n.tests <- nrow(gene_means) 
# proportions of genes positive
f.positive <-  n.positive / n.tests
exp.n.positive <- alpha * n.tests

  # FDR
gene_means$fdr <- stats::p.adjust(gene_means$p.value, method="fdr")
plot(gene_means$p.value, gene_means$fdr, log="xy"); abline(a=0, b=1); grid()

plot(gene_means$D, xlab="Effect size (d = m2 - m1)",
     -log10(gene_means$fdr), ylab="sig = -log10(FDR)",
     main = "FDR volcano plot",
     col=densCols(gene_means$D,-log10(gene_means$fdr)),
     ylim=c(min(-log10(gene_means$fdr)), max(-log10(gene_means$fdr))*1.2)
)
grid(col="grey", lty="solid")
abline(v=0, col="black") 
abline(v=c(-1,1), col="red") 
abline(h=-log10(alpha), col="red") 


##### Differential Genes Expression Analysis ####
###### 1. Generate model matrix - want to compare cell types ####
# check the cell types 
varLabels(gse37448)
head(gse37448$`cell type:ch1`)
pData(gse37448_celdata)[,c("cell type:ch1")]

# cell type information stored in description
pData(gse37448_celdata)[,c("description")]
  # check the number of cell types 
levels(factor(pData(gse37448_celdata)[, "description"]))

  # count how many samples fall into each category
cell_type <- pData(gse37448_eset)$description
print(data.frame(sort(table(cell_type), decreasing = TRUE)))


  # remove "" - step not performed - if done will have to :
    # remove samples which correspond to "" in eset obj 
    # re run rma 
pData(gse37448_celdata)$description[pData(gse37448_celdata)$description == ""] <- NA

###### Build Design without contrast ####
design <- model.matrix( ~ gse37448_eset[['description']])
design

  # Assign new column names
colnames(design)
new_names <- c("abT_Cell",
  "abT_C57BL6J", "B_Cell", "B_Cell_C57BL6J", "DC", "DC_C57BL6J", "DC_NA",
  "Eosinophil_C57BL6", "Eosinophil_C57BL6J", "Macrophage_C57BL6J", "Macrophage_mixed",
  "Neutrophil_C57BL6J", "Stromal_C57BL6J", "Stromal_C57BL6J_129SJL_WT",
  "Tgd_vg2_Th1", "Tgd_vg2_Th2", "Tgd_vg2_Th3", 
  "Tgd_vg4_Th1", "Tgd_vg4_Th2", "Tgd_vg4_Th3",
  "Tgd_vg4_lo1", "Tgd_vg4_lo2"
)
colnames(design) <- c("(Intercept)", new_names)

###### 2. Fit the model to the design ####
fit <- lmFit(gse37448_eset,design)

  # Empirical Bayes correction
fitted.ebayes <- eBayes(fit)

# extract DEGs
topTable(fitted.ebayes)

# lfc = 1, p value = 0.05 in abT cells - example
summary(decideTests(fitted.ebayes[,"abT_Cell"],lfc=1, p.value = 0.05))


###### Build Design with contrast ####
  # use to compare specific conditions, for this project : abT_Cell vs Macrophage_mixed
contrast_matrix <- makeContrasts(
  abT_vs_mf_mixed = abT_Cell - Macrophage_mixed,
  levels = design
)
contrast_matrix

# Fit Model - abT vs Macrophage mixed background
fit2 <- contrasts.fit(fit,contrasts=contrast_matrix)
fit2 <- eBayes(fit2)
summary(decideTests(fit2,lfc=1))

topTable(fit2)


##### Annotation ####
# AnnotationDbi interface for mouse
columns(mogene10sttranscriptcluster.db)
keytypes(mogene10sttranscriptcluster.db) # columns can be used as keys
head(keys(mogene10sttranscriptcluster.db,keytype="PROBEID"))

# extract probe ids from fitted model without contrast
top_res <- topTable(fitted.ebayes, number=Inf, p.value = 0.05, lfc=1)
top_res$PROBEID <- rownames(top_res)
pids <- rownames(topTable(fitted.ebayes,number=Inf,p.value = 0.05,lfc=1))

# add annotations
annot_df <- AnnotationDbi::select(mogene10sttranscriptcluster.db,pids,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
head(annot_df)
annotated_res <- merge(top_res, annot_df, by="PROBEID")
sum(is.na(annotated_res$SYMBOL))

##### Visualization ####

###### 1. Without Contrast ####
names(fitted.ebayes)

volcanoplot(fitted.ebayes, coef=1)
limma::plotMA(fitted.ebayes, coef=1, main="MA Plot: Fitted eBayes")
hist(fitted.ebayes$p.value[,1], breaks=50, main="P-value Distribution")
plotSA(fitted.ebayes, main="Mean-variance trend after eBayes")
qqt(fitted.ebayes$t[,1], df=fitted.ebayes$df.prior + fitted.ebayes$df.residual)

top_res_2 <- topTable(fitted.ebayes, number=Inf, p.value = 0.05, lfc=2)
head(top_res_2)

# add annotations
gene_symbols <- mapIds(mogene10sttranscriptcluster.db,
                       keys = rownames(top_res_2),
                       column = "SYMBOL",
                       keytype = "PROBEID",   # may need to change depending on probe type
                       multiVals = "first")
head(gene_symbols)
top_res_2$Gene <- gene_symbols
head(top_res_2[, c("Gene", "B_Cell")])

volcano_df <- top_res_2
volcano_df$logFC <- volcano_df$abT_Cell
volcano_df$negLogP <- -log10(volcano_df$adj.P.Val)

# Thresholds
lfc_cutoff <- 2
padj_cutoff <- 0.05

# Significance categories
volcano_df$significance <- "Not Sig"
volcano_df$significance[volcano_df$logFC >= lfc_cutoff & volcano_df$adj.P.Val < padj_cutoff] <- "Up"
volcano_df$significance[volcano_df$logFC <= -lfc_cutoff & volcano_df$adj.P.Val < padj_cutoff] <- "Down"

head(volcano_df[, c("Gene", "significance")])
table(volcano_df$significance)

top_genes <- volcano_df[volcano_df$significance != "Not Sig", ]
top_genes <- top_genes[order(top_genes$adj.P.Val), ]
top_genes <- head(top_genes, 10)

# Plot
ggplot(volcano_df, aes(x = logFC, y = negLogP, color = significance)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not Sig" = "grey")) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "black") +
  geom_text_repel(data = top_genes, aes(label = Gene), size = 3, max.overlaps = 10) +
  labs(title = "abT_Cell vs Baseline",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-Value") +
  theme_minimal()


# MA Plot with Significance
ggplot(volcano_df, aes(x = AveExpr, y = logFC, color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not Sig" = "grey")) +
  geom_text_repel(data = top_genes, aes(label = Gene), size = 3, max.overlaps = 10) +
  labs(title = "MA Plot with Significance",
       x = "Average Expression",
       y = "Log2 Fold Change") +
  theme_minimal()

###### 2. With Contrast ####

names(fit2)

volcanoplot(fit2, coef=1)
limma::plotMA(fit2, coef=1, main="MA Plot: Fitted eBayes with contrast - abT vs Macrophage mixed background")
plotSA(fit2, main="Mean-variance trend after eBayes with contrast - abT vs Macrophage mixed background")
qqt(fit2$t[,1], df=fit2$df.prior + fit2$df.residual)


top_res <- topTable(fit2, number=Inf, p.value = 0.05, lfc=2)
head(top_res)

# add annotations
gene_symbols_contrast <- mapIds(mogene10sttranscriptcluster.db,
                       keys = rownames(top_res),
                       column = "SYMBOL",
                       keytype = "PROBEID",   # may need to change depending on probe type
                       multiVals = "first")
head(gene_symbols_contrast)
top_res$Gene <- gene_symbols_contrast
head(top_res[, c("Gene")])

volcano_df_C <- top_res
volcano_df_C$negLogP <- -log10(volcano_df_C$adj.P.Val)

# Thresholds
lfc_cutoff <- 2
padj_cutoff <- 0.05

# Significance categories
volcano_df_C$significance <- "Not Sig"
volcano_df_C$significance[volcano_df_C$logFC >= lfc_cutoff & volcano_df_C$adj.P.Val < padj_cutoff] <- "Up"
volcano_df_C$significance[volcano_df_C$logFC <= -lfc_cutoff & volcano_df_C$adj.P.Val < padj_cutoff] <- "Down"

head(volcano_df_C[, c("Gene", "significance")])
table(volcano_df_C$significance)

top_genes_C <- volcano_df_C[volcano_df_C$significance != "Not Sig", ]
top_genes_C <- top_genes_C[order(top_genes_C$adj.P.Val), ]
top_genes_C <- head(top_genes_C, 10)

# Plot
ggplot(volcano_df_C, aes(x = logFC, y = negLogP, color = significance)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not Sig" = "grey")) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "black") +
  geom_text_repel(data = top_genes_C, aes(label = Gene), size = 3, max.overlaps = 10) +
  labs(title = "abT_Cell vs Macrophages Mixed Background",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-Value") +
  theme_minimal()

# MA Plot with Significance
ggplot(volcano_df_C, aes(x = AveExpr, y = logFC, color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not Sig" = "grey")) +
  geom_text_repel(data = top_genes_C, aes(label = Gene), size = 3, max.overlaps = 10) +
  labs(title = "MA Plot with Significance",
       x = "Average Expression",
       y = "Log2 Fold Change") +
  theme_minimal()

# find H2-D1, Actb, Gapdh, or Lrg1
volcano_df_C[volcano_df_C$Gene %in% c("H2-D1", "Actb", "Gapdh", "Lrg1"), 
             c("Gene", "logFC", "adj.P.Val", "significance")]



##### Multi Density Plot ####
keep <- cell_type %in% c("abT Cell", "Macrophage, mixed background")
expr_subset <- expr_matrix[, keep]
cell_subset <- cell_type[keep]


plotDensities(expr_subset, group = cell_subset,
              main = "Multi-density plot: abT Cell vs Macrophage_mixed",
              legend = "topright")

expr_long <- as.data.frame(t(expr_subset))
expr_long$CellType <- cell_subset

expr_long <- pivot_longer(expr_long, 
                          cols = -CellType, 
                          names_to = "Gene", 
                          values_to = "Expression")

ggplot(expr_long, aes(x = Expression, fill = CellType)) +
  geom_density(alpha = 0.4) +
  theme_bw() +
  labs(title = "Expression distribution: abT Cell vs Macrophage_mixed",
       x = "Expression (log2 RMA normalized)",
       y = "Density")

##### Running TopGO ####
head(expr)
head(top_res)

# pull out genes
  # DEGs from top_res (adjusted p < 0.05 and |logFC| >= 2)
deg_genes <- top_res %>%
  filter(adj.P.Val < 0.05 & abs(logFC) >= 1) %>%
  pull(Gene)
length(deg_genes)

 # genes from expr_matrix
expr_GOmat <- expr_matrix
# annotate - using different method for practice = BiomaRt
mart <- useMart("ensembl")
ensembl <- useDataset("mmusculus_gene_ensembl", mart)
values <- rownames(expr_GOmat)
attributes <- listAttributes(ensembl)

biomart <- getBM(attributes = c("ensembl_gene_id", "affy_mogene_1_0_st_v1"),
                 filters = "affy_mogene_1_0_st_v1", values = values, mart = ensembl)
head(biomart)

mart_1 <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
genes <- biomart$ensembl_gene_id
G_list <- getBM(filters= "ensembl_gene_id", 
                attributes= c("ensembl_gene_id","external_gene_name"),
                values=genes,mart= mart_1)
head(G_list)

Probe_data <- full_join(biomart, G_list , by = c("ensembl_gene_id" = "ensembl_gene_id"), relationship = "many-to-many")
View(Probe_data)
head(Probe_data)

Probe_data <- Probe_data %>%
  mutate(affy_mogene_1_0_st_v1 = as.character(affy_mogene_1_0_st_v1))

expr_GOmat <- expr_matrix %>% 
  as.data.frame() %>% 
  mutate(PROBEID = rownames(expr_matrix)) %>% 
  left_join(Probe_data %>% select(affy_mogene_1_0_st_v1, external_gene_name), 
            by = c("PROBEID" = "affy_mogene_1_0_st_v1")) %>%
  filter(!is.na(external_gene_name)) %>%
  rename(Gene = external_gene_name)
head(expr_GOmat)

# extract genes
all_genes <- expr_GOmat$Gene

length(deg_genes)
length(all_genes)

geneList <- factor(ifelse(all_genes %in% deg_genes, 1, 0))
names(geneList) <- all_genes
levels(geneList)
table(geneList)

run_topGO <- function(ontology){
  GOdata <- new(
    "topGOdata",
    ontology = ontology,
    allGenes = geneList,
    geneSelectionFun = function(x) x == 1,
    annot = annFUN.org,
    mapping = "org.Mm.eg.db",
    ID = "symbol"
  )
  
  resultClassic <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  resultElim    <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
  
  topGOres <- GenTable(
    GOdata,
    classicFisher = resultClassic,
    elimFisher = resultElim,
    orderBy = "classicFisher",
    ranksOf = "classicFisher",
    topNodes = 20
  )
  
  # Barplot
  ggplot(topGOres, aes(x = reorder(Term, -as.numeric(classicFisher)),
                       y = -log10(as.numeric(classicFisher)))) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    ylab("-log10(p-value)") +
    xlab(paste("GO Term (", ontology, ")", sep="")) +
    ggtitle(paste("TopGO Enrichment -", ontology)) +
    theme_minimal()
  
  return(topGOres)
}

topGO_BP <- run_topGO("BP")
head(topGO_BP)

pvals_GO <- topGO_BP$classicFisher
pvals_GO <- ifelse(pvals_GO == "<1e-30", 1e-30, pvals_GO)
pvals_GO <- as.numeric(pvals_GO)
names(pvals_GO) <- topGO_BP$GO.ID
head(pvals_GO)

go_df <- data.frame(
  GO_ID = names(pvals_GO),
  pvalue = pvals_GO
)

go_df$Term <- topGO_BP$Term
go_df$Significant <- topGO_BP$Significant
go_df$Annotated <- topGO_BP$Annotated

# plot dot plot
ggplot(go_df, aes(x = -log10(pvalue), y = Term, size = Significant)) +
  geom_point(color = "darkred") +
  xlab("-log10(p-value)") +
  ylab("GO Term") +
  ggtitle("GO Enrichment - TopGO") +
  theme_minimal()

# plot bubble plot
go_df$ratio <- go_df$Significant / go_df$Annotated

ggplot(go_df, aes(x = ratio, y = Term, size = Significant, color = -log10(pvalue))) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  xlab("DE genes / total annotated genes") +
  ylab("GO Term") +
  ggtitle("GO Enrichment Bubble Plot") +
  theme_minimal()

# optional
topGO_MF <- run_topGO("MF")
topGO_CC <- run_topGO("CC")
  
  
#### Immgen Dataset Analysis ####
# Following : https://bioconductor.org/packages/release/workflows/vignettes/maEndToEnd/inst/doc/MA-Workflow.html
###### Load data #####
GSE37448_Normalized_Data <- read.csv("~GSE37448/GSE37448_Normalized_Data.csv")
View(GSE37448_Normalized_Data)

head(GSE37448_Normalized_Data)

###### Remove multiple mappings #####
# Count unique ProbeSetID mappings for each GeneSymbol
mapping_summary <- GSE37448_Normalized_Data %>%
  group_by(GeneSymbol) %>%
  summarise(gene_count = n_distinct(ProbeSetID))

# Count the number of GeneSymbol with more than one ProbeSetIDs
multi_gene_count <- mapping_summary %>%
  filter(gene_count > 1) %>%
  summarise(total = n())
print(multi_gene_count) # 1854

# Remove multiple mappings
immgen_data <- GSE37448_Normalized_Data %>%
  distinct(GeneSymbol, .keep_all = TRUE)

# Set the rownames as the gene symbols
rownames(immgen_data) <- immgen_data$GeneSymbol
immgen_data <- subset(immgen_data, select=-GeneSymbol)
immgen_data <- subset(immgen_data, select=-ProbeSetID)
dim(immgen_data)


###### PCA #####
pca_results <- prcomp(t(immgen_data), scale. = FALSE) 
pca_df <- as.data.frame(pca_results$x)

ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(size = 1) +
  xlab(paste0("PC1 (", round(100 * pca_results$sdev[1]^2 / sum(pca_results$sdev^2), 2), "%)")) +
  ylab(paste0("PC2 (", round(100 * pca_results$sdev[2]^2 / sum(pca_results$sdev^2), 2), "%)")) +
  theme_minimal() +
  ggtitle("PCA for RMA log-normalized signal intensity of GSE37448 Immgen Dataset") 

###### UMAP #####
umap <- umap(t(immgen_data), n_neighbors = 5, random_state = 123)
umap_df <- as.data.frame(umap$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")

ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(size = 1) +                                  # Plot points
  theme_minimal() +
  labs(title = "UMAP for RMA log-normalized signal intensity of GSE37448 Immgen Dataset")


sample_names <- as.data.frame(colnames(immgen_data))
colnames(sample_names) <- 'Samples'
head(sample_names, 20)
