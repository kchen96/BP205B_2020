library(VISION)
library(tidyverse)
library(data.table)
options(mc.cores=25)

create_signature_set <- function(dataset,DE_results,ref_gene_set,abs=T){
  filenames <- stringr::str_sub(basename(DE_results),1,-5)
  
  allSigSets <- unlist(lapply(1:length(DE_results),function(f){
    gmt <- data.table::fread(DE_results[f],sep=',')
    sigData <- sign(gmt$logfoldchange)
    names(sigData) <- gmt$genes
    matched <- names(sigData) %in% ref_gene_set
    
    print(paste0('Unmatched ',f,' ',filenames[f],' ',sum(!matched)/nrow(gmt)))
    
    sigSets <- unlist(lapply(c(0,1,2),function(thresh){
      fold_data <- sigData[abs(gmt$logfoldchange)>=thresh & matched]
      if(abs){
        fold_data <- abs(fold_data)
      }
      sig_set_name <- paste0(dataset,'_',filenames[f],'_FC_',thresh)
      
      a <- VISION::createGeneSignature(name = sig_set_name,sigData = fold_data)

      retlist <- list(a)
      names(retlist) <- sig_set_name
  
      
      return(retlist)
    }))
    return(sigSets)
  }))
  return(allSigSets)
}

h5read <- function(h5path){
  hfile <- hdf5r::h5file(filename = h5path, mode = 'r')
  
  m <- as.matrix(Matrix::sparseMatrix(i=hfile[['X']][['indices']][]+1,
                             p=hfile[['X']][['indptr']][],
                             x=hfile[['X']][['data']][]
  ))
  hdf5r::h5close(hfile)
  return(m)
}

h5read_umap <- function(h5path){
  hfile <- hdf5r::h5file(h5path,mode = 'r')

  umap <- as.data.frame(t(hfile[['obsm']][['X_umap']][,]))
  
  hdf5r::h5close(hfile)
  return(umap)

}

## read in the data
pdx_exprmat <- h5read('/wynton/scratch/bp205/processed//PDX_normalized_adata.h5ad')
# norm.factor = median(colSums(pdx_exprmat))
# pdx_exprmat <- t( t(pdx_exprmat) / colSums(pdx_exprmat)) * norm.factor

pdx_meta <- as.data.frame(data.table::fread('/wynton/scratch/bp205/processed/PDX_adata_R/obs.csv',sep = ','))
pdx_meta$SampleType <- factor(pdx_meta$SampleType,levels = paste0('Sample',1:26))
pdx_genes <- as.data.frame(data.table::fread('/wynton/scratch/bp205/processed/PDX_adata_R/var.csv',sep = ','))
rownames(pdx_genes) <- pdx_genes$V1
## Rename the exprmat
rownames(pdx_exprmat) <- rownames(pdx_genes)
colnames(pdx_exprmat) <- 1:ncol(pdx_exprmat)
###
pdx_umap <- h5read_umap('/wynton/scratch/bp205/processed/PDX_adata.h5ad')
rownames(pdx_umap) <- 1:ncol(pdx_exprmat)

## Load and generate the signature files
hcc_mk_DE_results <- list.files('/wynton/scratch/bp205/signatures/BP205B_2020/DE/HCC_MK/FINAL/',
                                recursive = T,full.names = T,pattern = '.csv')
mda_DE_results <- list.files('/wynton/scratch/bp205/signatures/BP205B_2020/DE_results/Final/',
                             recursive = T,full.names = T,pattern = '.csv')
mda_DE_signatures <- create_signature_set('MDA',mda_DE_results,pdx_genes$V1,abs = F)
hcc_mk_DE_signatures <- create_signature_set('HCC_MK',hcc_mk_DE_results,pdx_genes$V1,abs = F)

all_DE_signatures <- c(mda_DE_signatures,hcc_mk_DE_signatures)
## Load signature files by adding a sign to them 

pdx_vis <- VISION::Vision(data=pdx_exprmat,
                          signatures=all_DE_signatures,
                          sig_gene_threshold=0.005,
                          meta=pdx_meta,projection_methods=c('UMAP','tSNE30'),
                          pool=FALSE,name='PDX_agg_sizefactor_norm')

pdx_vis_result <- VISION::analyze(pdx_vis)
pdx_vis_result<- addProjection(pdx_vis_result,"UMAP",pdx_umap)
save(pdx_vis_result,file = '/wynton/home/students/snanda/rds/bp205/analysis/pdx/PDX_agg_sizefactor_norm.rda',version = 2)


# load('/wynton/home/students/snanda/rds/bp205/analysis/pdx/PDX_agg_sizefactor_norm.rda')
VISION::viewResults(pdx_vis_result)


