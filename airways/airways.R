# Install required packages ----

# Use renv and rig for controlling versions of packages and r.
# Use Bioconductor packages for various bioinformatics functions.
# Use SummarizedExperiment to access an example dataset.

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("airway", "biomaRt", "edgeR", "limma", "SummarizedExperiment"))

# Load libraries ----

install.packages("dyplyr")
library(airway)
library(biomaRt)
library(edgeR)
library(ggplot2)
library(limma)
library(ggrepel)
library(SummarizedExperiment)

install.packages(c("lintr", "styler"))

# Analysis variables ----

# Create an output directory
save_dir <- "results/"
dir.create(save_dir)

# Gene filtering
min.cpm <- 1
min.cpm.fraction <- 1 / 4

# Data preparation ----
## Load and inspect experiment data ----

# The 'airway' package containes an example RNAseq experiment.
# Himes et al. 'RNA-Seq Transcriptome Profiling Identifies CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates Cytokine Function in Airway Smooth Muscle Cells.'
# PLoS One. 2014 Jun 13;9(6):e99625. PMID: 24926665. GEO: GSE52778.

data(airway, package = "airway")
str(airway)
str(airway@assays@data@listData[["counts"]])

# Inspect the samples and experiment groupings (dex treatment); 4 treated vs 4 untreated
table(airway$Sample, airway$dex)
treatment.groups <- factor(airway$dex)

# Could edit here to rename the treatment groups "Control" and "DEX" to match the paper, then see below that later code needs to be updated?
# treatment.groups <- vector()
# treatment.groups[airway$dex=="untrt"] <- "Control"
# treatment.groups[airway$dex=="trt"] <- "DEX"
# treatment.groups <- factor(treatment.groups)

# Obtain the raw counts as a data matrix: genes(rows) x sample(columns)
raw.counts.matrix <- assay(airway, "counts")
str(raw.counts.matrix)
# Rows are genes, columns are samples, values are integer expression counts


## Filter out low-expression counts ----

# Counts per million (cpm)
cpm <- cpm(raw.counts.matrix)
dim(cpm)

# Require genes have at least 'min.cpm' in at least 'min.cpm.fraction' of the samples
keepGenes <- rowSums(cpm >= min.cpm) >= ncol(cpm) * min.cpm.fraction
sum(keepGenes)
filteredCountsMat <- data.matrix(raw.counts.matrix[keepGenes, ])
dim(filteredCountsMat)
print(paste0("Keeping ", sum(keepGenes), " of ", nrow(cpm), " genes."))


# Convert Ensembl ID numbers to Gene Symbols ----
# This section covers a typical issue, where genes are referred to by an accession number (in this case Ensembl ID).
# It is easier to interpret them if they are named with their gene symbol. The biomaRt package can be used.

# get an object that references the relevant biomart database (human gene names)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# # Could specify which version of the ensembl gene names for better reproducibility?
# # Some gene names changed over time (with new research,
# # or fixing problematic gene names eg. "MARCH7" text name autocorrecting to a
# # date value if you ever view it in Excel... https://doi.org/10.1371/journal.pcbi.1008984)
# ensembl_110 <- useEnsembl(biomart = 'genes',
#                           dataset = 'hsapiens_gene_ensembl',
#                           version = 110)

ensemble.ids <- row.names(filteredCountsMat)
listAttributes(mart = ensembl)
gene.table <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name"),
  filters = "ensembl_gene_id", values = ensemble.ids,
  mart = ensembl, uniqueRows = TRUE
)

# an easy mistake would be to mismatch the rows, or not check for 1:1 relationship. (manual checking steps)
dim(gene.table)
table(is.na(gene.table$hgnc_symbol))
table(duplicated(gene.table$hgnc_symbol))

# inspect some of the problem rows
which(duplicated(gene.table$hgnc_symbol))
gene.table[2809, ]
gene.table$hgnc_symbol[which(duplicated(gene.table$hgnc_symbol))]

# There is a blank value for gene symbol in a few hundred cases.
# They can fall back to Ensemble Id.
gene.table$hgnc_symbol[gene.table$hgnc_symbol == ""] <- gene.table$ensembl_gene_id[gene.table$hgnc_symbol == ""]
table(is.na(gene.table$hgnc_symbol)) # ok
table(duplicated(gene.table$hgnc_symbol)) # ok

gene.table$hgnc_symbol[grep(x = gene.table$hgnc_symbol, "MARCH")]
# just showing at some point they inserted an "F" in the MARCH* genes

# Using the gene.table to convert gene names
row.names(filteredCountsMat) <- gene.table$hgnc_symbol[match(
  row.names(filteredCountsMat),
  table = gene.table$ensembl_gene_id
)]

## Normalise counts ----

# Upper-quartile normalise function
# This is a method to correct or differences in library depth between samples.
# The library-normalised counts can then be used to compare expression of a gene between the samples.
quartileNormalise <- function(data, q = 4) {
  # Get the upper quartile (if q=4) value of each sample (samples as the rows of 'data')
  upperquartiles <- vector(length = nrow(data))
  for (i in 1:nrow(data)) {
    upperquartiles[i] <- quantile(data[i, ])[q]
  }
  # now scale the data.. divide by uq and multiply by the average-uq
  normalised.data <- data
  for (i in 1:nrow(data)) {
    normalised.data[i, ] <- normalised.data[i, ] / upperquartiles[i] * median(upperquartiles)
  }
  # give back the normalised data result
  return(normalised.data)
}
filteredCountsMatUQ <- t(quartileNormalise(t(filteredCountsMat)))

# View distribution before normalising
boxplot(log2(1 + filteredCountsMat),
  las = 2, cex.axis = 0.5, cex = 0.5,
  ylab = "log2(1+counts)"
)

# View distribution after normalising
boxplot(log2(1 + filteredCountsMatUQ),
  las = 2, cex.axis = 0.5, cex = 0.5,
  ylab = "log2(1+counts)"
)

# Output: table of prepared gene expression data (counts x samples)
write.csv(filteredCountsMatUQ, file = paste0(save_dir, format(Sys.time(), "%Y-%m-%dT%H-%M-%S"), "_", "filteredCountsMatUQ.csv"))

# Exploratory plots to view the amount of data collected in each sample:
# Plot the depth of sequencing per sample
plot(colSums(filteredCountsMat) / 10^6,
  ylab = "Counts (millions)",
  ylim = c(0, max(colSums(filteredCountsMat) / 10^6)),
  xlab = "Sample #"
)
# Add experiment groupings
plot(
  y = colSums(filteredCountsMat) / 10^6, ylab = "Counts (millions)",
  ylim = c(0, max(colSums(filteredCountsMat) / 10^6)),
  x = factor(airway$dex)
)

# DEG analysis ----
# (Differential Gene Expression tests)

## limma ----

# 'dgelist' object with groups
X.dgelist <- DGEList(filteredCountsMat, group = treatment.groups)
dim(X.dgelist)

# Normalise libraries
X.dgelist <- calcNormFactors(X.dgelist, method = "TMM")
X.dgelist$samples$norm.factors

# limma test design and contrast objects
design <- model.matrix(~ 0 + treatment.groups)
colnames(design) <- c("trt", "untrt")

contr <- makeContrasts(
  dex    = trt - untrt,
  levels = colnames(design)
)

v <- voom(counts = X.dgelist, design = design, plot = TRUE)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts = contr)
efit <- eBayes(vfit)
plotSA(efit)

topTable.dex <- topTable(efit, coef = "dex", number = Inf, sort.by = "p")
head(topTable.dex)
hist(topTable.dex$P.Value)

# Check the published DEG: CRISPLD2 (ENSG00000103196)
gene.of.interest <- "ENSG00000103196" # (i.e. CRISPLD2, if not convering from the provided Ensembl ID)
gene.of.interest <- "CRISPLD2" # after converting to gene symbols

topTable.dex[which(topTable.dex$ID == gene.of.interest), ]
plot(log2(1 + filteredCountsMatUQ[gene.of.interest, ]),
  x = treatment.groups,
  ylab = "UQ-normalised mRNA counts", xlab = "", main = gene.of.interest
)

# Could continue with typical visualisations (heatmap, volcano plot), PCA analysis, etc.

## Volcano plot ----

# basic volcano plot (x-axis log2 fold change, y-axis Significance)
volcanoplot(efit,
  coef = 1, style = "p-value", highlight = 10,
  names = row.names(efit), hl.col = "blue",
  xlab = "Log2 Fold Change", ylab = NULL, pch = 16, cex = 0.35
)
# issues seen (overlapping labels, unequal x-axis range)

# custom volcano plot
volcano.plot.data.frame <- data.frame(topTable.dex, EnsembleId = row.names(topTable.dex))

n.highlighted <- 30 # number of genes to be labelled, in order of p-value
volcano.plot.data.frame$label <- NA
volcano.plot.data.frame$label[1:n.highlighted] <-
  volcano.plot.data.frame$ID[1:n.highlighted]

# basic test version
ggplot(volcano.plot.data.frame, aes(x = logFC, y = (-log10(adj.P.Val)))) +
  geom_point() +
  geom_text_repel(aes(label = label), max.overlaps = Inf, size = 3, colour = "blue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_light()

# (apply customisations from manual testing)
ggplot(volcano.plot.data.frame, aes(x = logFC, y = (-log10(adj.P.Val)))) +
  xlim(c(-8, 8)) +
  ylim(c(0, 4.5)) +
  geom_point(size = 1, alpha = 0.4) +
  geom_text_repel(aes(label = label),
    max.overlaps = Inf, size = 3, colour = "blue",
    force = 5, min.segment.length = 0, segment.alpha = 0.2
  ) +
  geom_hline(yintercept = -log10(0.05), alpha = 0.5, colour = "darkgreen", linetype = "dashed") +
  annotate("label", y = -log10(0.05), x = -7, label = "adj.p=0.05", colour = "darkgreen") +
  theme_light()
