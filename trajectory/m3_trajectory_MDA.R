library(monocle3)
library(data.table)

R_data_directory_to_cds <- function(R_data_dir_path){
  sample_sheet <- as.data.frame(fread(paste0(R_data_dir_path,'/obs.csv')))
  colnames(sample_sheet)[1] <- 'sampleNames'
  rownames(sample_sheet) <- sample_sheet$sampleNames
  
  gene_annotation <- as.data.frame(fread(paste0(R_data_dir_path,'/var.csv')))
  colnames(gene_annotation)[1:2] <- c('gene_short_name','featureNames')
  rownames(gene_annotation) <- gene_annotation$gene_short_name
  
  expr <- t(as.matrix(fread(paste0(R_data_dir_path,str_sub(basename(R_data_dir_path),1,-3),'.csv'),header = F)))
  colnames(expr) <- sample_sheet$sampleNames
  rownames(expr) <- gene_annotation$gene_short_name
  
  cds <- new_cell_data_set(expr,sample_sheet,gene_annotation)
  return(cds)
  
}

mda_cds <- R_data_directory_to_cds('/wynton/scratch/bp205/processed/MDA_adata_R//')

mda_cds <- preprocess_cds(mda_cds,method = 'PCA')
plot_pc_variance_explained(mda_cds)

mda_cds <- reduce_dimension(mda_cds,max_components = 3)

mda_cds <- cluster_cells(mda_cds,cluster_method = 'leiden')

plot_cells(mda_cds,color_cells_by = 'SampleType')

mda_cds <- learn_graph(mda_cds,verbose = T)

mda_cds <- order_cells(mda_cds)

plot_cells(mda_cds,color_cells_by = "SampleType")

