BiocManager::install("monocle")
library(monocle)
library(tidyverse)
library(MASS)
library(reshape)
library(ggplot2)


expr_mtx <- read_csv("./MDA_adata_R/MDA_adata.csv", col_names=FALSE)#expression matrix
expr_mtx <- t(as.matrix(expr_mtx))
sample_sheet <- read_csv("./MDA_adata_R/obs.csv" ) #pheno data
row.names(sample_sheet) <- colnames(expr_mtx)
genes_ann <- as.data.frame(read_csv("./MDA_adata_R/var.csv")) #feature data
row.names(genes_ann) <- row.names(expr_mtx)

pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = genes_ann)
HCC_MK <- newCellDataSet(as.matrix(expr_mtx),
                         phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())

HCC_MK <- estimateSizeFactors(HCC_MK)
HCC_MK <- estimateDispersions(HCC_MK)

HCC_MK <- detectGenes(HCC_MK, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(HCC_MK),
                                    num_cells_expressed >= 10))

# Log-transform each value in the expression matrix.
L <- log(exprs(HCC_MK[expressed_genes,]))

# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))

# Plot the distribution of the standardized gene expression values.
qplot(value, geom = "density", data = melted_dens_df) +
  stat_function(fun = dnorm, size = 0.5, color = 'red') +
  xlab("Standardized log(FPKM)") +
  ylab("Density")

#Proceed with Trajectory Analysis: technically parental vs LM2 are beginning vs. end of experimental selection
#can use as timepoints for trajectory analysis #might need to adapt for MDA dataset
fData(HCC_MK)$use_for_ordering <-
  fData(HCC_MK)$num_cells_expressed > 0.05 * ncol(HCC_MK)
plot_pc_variance_explained(HCC_MK, return_all = F)

HCC_MK <- reduceDimension(HCC_MK,
                          max_components = 3,
                          norm_method = 'log',
                          num_dim = 3,
                          reduction_method = 'tSNE', 
                          residualModelFormulaStr = "~num_genes_expressed",
                          verbose = T)
HCC_MK <- clusterCells(HCC_MK, verbose=F)
#sanity check for clustering
plot_cell_clusters(HCC_MK, color_by = 'Cluster')
plot_cell_clusters(HCC_MK, color_by = 'SampleType')

#DEG
clustering_DEG_genes <-
  differentialGeneTest(HCC_MK[expressed_genes,],
                       fullModelFormulaStr = '~Cluster',
                       cores = 2)
#Selecting top 1000 DE genes for ordering
HCC_MK_ordering_genes <-
  row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
HCC_MK_ordering_genes <-setOrderingFilter(HCC_MK,ordering_genes = HCC_MK_ordering_genes)
HCC_MK <-reduceDimension(HCC_MK, method = 'DDRTree')
#order along pseudotime
HCC_MK <-orderCells(HCC_MK)
#write function to define initial state; otherwise monocle chooses arbitrarily
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$SampleType)[,"TGL"]
    return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
HCC_AA <-orderCells(HCC_MK, root_state = GM_state(HCC_MK))
plot_cell_trajectory(HCC_MK, color_by = "SampleType")
plot_cell_trajectory(HCC_MK, color_by = "Pseudotime")
plot_cell_trajectory(HCC_MK, color_by = "SampleType") + facet_wrap("~State")

#Finding genes that change as a function of branching
#1 first critical branch along pseudotime between parental and LM2
BEAM_res <- BEAM(HCC_MK, branch_point = 1, cores = 2)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

plot_genes_branched_heatmap(HCC_AA[row.names(subset(BEAM_res,
                                                    qval < 1e-4)),],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)
#gotta figure out cluster optimizations: what should it correspond to?

#repeat BEAM for each critical branch point
