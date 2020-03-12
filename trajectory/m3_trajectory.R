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

hcc_mk_cds <- R_data_directory_to_cds('/wynton/scratch/bp205/processed/HCC_MK_adata_R/')

hcc_mk_cds <- preprocess_cds(hcc_mk_cds)
plot_pc_variance_explained(hcc_mk_cds)

hcc_mk_cds <- reduce_dimension(hcc_mk_cds,max_components = 3)
hcc_mk_cds <- cluster_cells(hcc_mk_cds,cluster_method = 'leiden')

hcc_mk_cds <- learn_graph(hcc_mk_cds,verbose = T)

plot_cells(hcc_mk_cds,color_cells_by = 'partition')
plot_cells(hcc_mk_cds,color_cells_by = 'SampleType')

hcc_mk_cds <- order_cells(hcc_mk_cds)

plot_cells_3d(hcc_mk_cds,
           color_cells_by = "pseudotime")
plot_cells(hcc_mk_cds,color_cells_by = 'SampleType')


## Graph testing
closest_vertex <-hcc_mk_cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex <- as.matrix(closest_vertex[colnames(hcc_mk_cds), ])

hcc_mk_cds_pr_test_res <- graph_test(hcc_mk_cds,neighbor_graph = "principal_graph",cores=20)
pr_deg_ids <- hcc_mk_cds_pr_test_res %>% filter(q_value < 0.05) %>% arrange(q_value)


plot_cells(hcc_mk_cds,genes=pr_deg_ids$gene_short_name[1000],
           show_trajectory_graph = T,label_roots = T)

gene_module_df <- find_gene_modules(hcc_mk_cds[pr_deg_ids$gene_short_name,],cores = 10)

cell_group_df <- tibble(cell=rownames(colData(hcc_mk_cds)),cell_group=colData(hcc_mk_cds)$SampleType)

agg_mat <- aggregate_gene_expression(hcc_mk_cds,gene_module_df,cell_group_df)

plot_cells(hcc_mk_cds,genes = gene_module_df)

plot_genes_in_pseudotime(hcc_mk_cds[rowData(hcc_mk_cds)$gene_short_name %in% pr_deg_ids$gene_short_name[1:4],],color_cells_by = 'SampleType')




