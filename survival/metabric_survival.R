library(data.table)
library(tidyverse)
library(MASS)
library(ggplot2)
library(gridExtra)
library(readxl)
library(survival)
library(survminer)
library(ranger)
library(HGNChelper)

load_signature_set <- function(dataset,DE_results,ref_gene_set,abs=T){
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
      
      a <- fold_data
      
      retlist <- list(a)
      names(retlist) <- sig_set_name
      
      
      return(retlist)
    }),recursive = F)
    return(sigSets)
  }),recursive = F)
  return(allSigSets)
}

read_metabric_data <- function(fix_names=F){
  data <- fread('/wynton/scratch/bp205/METABRIC/METABRIC_MedCen_Collapsed_GEOannot(n=1992).merge.txt')
  
  if(fix_names==TRUE){
    checked <- suppressWarnings(HGNChelper::checkGeneSymbols(data$V1)) %>% 
      mutate(to_keep = !(str_detect(Suggested.Symbol,'//') | is.na(Suggested.Symbol)))
    
    data2  <- data[checked$to_keep,]
    data2$V1 <- checked$Suggested.Symbol[checked$to_keep]
    
    duplicated_names <- names(which(table(data2$V1) > 1))
    duplicated_to_drop <- unlist(lapply(duplicated_names,function(name){
      which(data2$V1 %in% name)[-1]
    }))
    data2 <- data2[-duplicated_to_drop,,drop=FALSE]
  }else{
    data2 <- data
  }
  
  data3 <- dcast(melt(data2,id.vars = 'V1'),formula = as.character(variable)~V1)
  setDF(data3)
  
  colnames(data3)[1] <- 'patient'
  
  return(data3)
}

read_metabric_metadata <- function(censor = 120){
  metadata <- readxl::read_excel('/wynton/scratch/bp205/METABRIC/brca_metabric_clinical_data.xlsx') %>%
    rename(patient=`Patient ID`,sampleID=`Sample ID`,
           vital_status = `Patient's Vital Status`,
           overall_survival_months=`Overall Survival (Months)`) %>%
    mutate(patient= str_replace(patient,'-',''),
           sampleID = str_replace(sampleID,'-',''),
           vital_status = ifelse(vital_status=='Died of Disease',1,0),
           vital_status = ifelse(overall_survival_months>censor,0,vital_status),
           overall_survival_months = ifelse(overall_survival_months>censor,censor,overall_survival_months),
           overall_survival_years = overall_survival_months/12)
  return(metadata)
}

select_and_score_metagene <- function(df,metagene,facets=c()){
  default_facets<- c('patient','overall_survival_months','overall_survival_years','vital_status')
  added_facets <- facets
  facets <- c(default_facets,added_facets)

  if(is.character(metagene)){
    genes <- metagene[metagene %in% colnames(df)]
    sign <- rep(1,length(genes))
    names(sign) <- genes
  }else{
    genes <- names(metagene)[names(metagene) %in% colnames(df)]
    sign <- metagene[names(metagene) %in% colnames(df)]
  }
  scored <- df %>%
    select_at(vars(c(genes,facets))) %>%
    pivot_longer(all_of(genes)) %>%
    group_by_at(vars(facets)) %>%
    summarize(activation=sum(value)) %>%
    
    #summarize(activation=sum(value*sign[name])) %>%
    ungroup %>% mutate_at(vars(added_facets),as.factor)
  return(scored)
}

stratify_activation <- function(df,ntiles=2,threshold=T){
  
  if(ntiles == 2 && threshold==TRUE){
    ## Do the iteration over the tresholds to optimize Up/Down partiioning
    thresholds <- sort(df$activation)
    thresholds <- thresholds[10:(length(thresholds)-10)]
    
    thresh_surv <- Surv(df$overall_survival_years,df$vital_status)
    message('Running threshold optimization')
    pvals <- unlist(lapply(thresholds,function(t){
      updn <- as.numeric(df$activation>t)
      diff <- survdiff(thresh_surv~updn)
      pv <- (pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE))
      return(pv)
    }))
    final_threshold <-thresholds[which.max(-log10(pvals))]
    strata <- factor(ifelse(df$activation>final_threshold,'High','Low'),levels=c('High','Low'))
    
  }else if(ntiles==3){
    strata <- factor(c('Hi','Med','Lo')[dplyr::ntile(df$activation,n_tiles)],levels=c('Lo','Med','Hi'))
  }else{
    ## For ntiles==2, this is median partitioning
    strata <- as.factor(dplyr::ntile(df$activation,ntiles))
  }
  df$strata <- strata
  return(df)
}

#########################################################################
### Implement wrapper function for running all signatures and scoring them 
score_stratify_fit_metagene <- function(df,metagene_name,metagene,ntiles=2,threshold=T,other_facets,outdir){
  print(metagene_name)
  # Create the directories if not present
  lapply(c('km_plots','factor_plots','hazard'),function(path){
    dir.create(paste0(outdir,'/',path,'/'),showWarnings = F)
    return(NULL)
  })

  km_plot_path <- paste0(outdir,'/km_plots/',metagene_name,'_km_plot.pdf')
  factor_plot_path <- paste0(outdir,'/factor_plots/',metagene_name,'_factor_plot.pdf')
  hazard_plot_path <- paste0(outdir,'/hazard_plots/',metagene_name,'_hazard_plot.pdf')
  
  ## Generate the dataset by scoring the metagene
  test <- df %>%
    select_and_score_metagene(metagene,facets=other_facets) %>%
    stratify_activation(ntiles = ntiles,threshold = threshold)
  rm(df)
  ## Run Kaplan-Meir analysis
  
  message('Running surviavl analysis')
  km_fit <- survminer::surv_fit(survival::Surv(overall_survival_years, vital_status) ~ strata, data=test)
  km_pval <- c(survminer::surv_pvalue(km_fit)$pval,
               summary(coxph(Surv(overall_survival_years, vital_status) ~ strata,data=test))$coef[2],
               km_fit$n)
  names(km_pval) <- c('pval_KM','HZ_KM','n_low','n_high')
  
  km_plt <- arrange_ggsurvplots(list(ggsurvplot(km_fit,
                                                pval=T,conf.int = T,risk.table = T,palette = 'Dark2',risk.table.col='strata')),
                                print=FALSE,ncol=1,nrow=1
  )
  ggsave(plot = km_plt,filename = km_plot_path,width = 8,height = 6)
  
  # message('Running cox fit')
  cox_fit <- coxph(Surv(overall_survival_years, vital_status) ~ strata +`Chemotherapy` + `Pam50 + Claudin-low subtype`, data=test,model = T)
  cox_pval <- c(summary(cox_fit)$coef[,c(2,5)],recursive=T)
  names(cox_pval) <- c(paste0('HZ_',rownames(summary(cox_fit)$coef)),paste0('pval_',rownames(summary(cox_fit)$coef)))

  cox_plt <- suppressWarnings(ggforest(cox_fit,data=test))
  ggsave(plot = cox_plt,filename = hazard_plot_path, width=10,height=7)

  result <- c(km_pval,cox_pval)
  rm(test)
  ## Run Factor analysis
  return(result)
}

##################
## Run the analysis
##################
data <- read_metabric_data(fix_names = T)
metadata <- read_metabric_metadata()
df <- inner_join(data,metadata,"patient")


hcc_mk_DE_results <- list.files('/wynton/scratch/bp205/signatures/BP205B_2020/DE/HCC_MK/FINAL/',
                                recursive = T,full.names = T,pattern = '.csv')
mda_DE_results <- list.files('/wynton/scratch/bp205/signatures/BP205B_2020/DE_results/Final/',
                             recursive = T,full.names = T,pattern = '.csv')

hcc_sig_set <- load_signature_set('HCC',hcc_mk_DE_results,colnames(df),abs = F)
mda_sig_set <- load_signature_set('MDA',mda_DE_results,colnames(df),abs = T)
all_sig_set <- c(mda_sig_set,hcc_sig_set)
all_sig_set <- all_sig_set[sapply(all_sig_set,length) != 0]

other_facets <- c('Cancer Type' , 'Cellularity' , 'Chemotherapy' ,'ER Status' , 'HER2 Status' , 'PR Status','Tumor Stage' , 'Age at Diagnosis' , 'Pam50 + Claudin-low subtype')

outdir <- '/wynton/home/students/snanda/rds/bp205/analysis/survival/results/'

rm(data,hcc_sig_set,mda_sig_set,metadata)

models <-parallel::mcmapply(FUN=score_stratify_fit_metagene,
       metagene_name=names(all_sig_set),
       metagene = all_sig_set,
       ntiles=2,threshold=T,outdir = outdir,SIMPLIFY = F,
       MoreArgs = list(df=df,other_facets = other_facets),mc.cores = 25)


models_df <- as_tibble(rownames_to_column(as.data.frame(do.call(rbind,models)))) %>% arrange(desc(HZ_KM))
colnames(models_df)[1] <- 'signature'

data.table::fwrite(models_df,paste0(outdir,'survival_results.csv'),sep = ',')



# facet_boxplot <- function(df,facet){
#   f <- ensym(facet)
#   colors <- c("#DA6162","#0C9D8B","#98BE37","#005C9A","#C291EE","#CD9849","#2FAD29","#002B95","#E579C6","#00899D")
#   active_colors <- colors[1:length(unique(df[[facet]]))]
#   
#   test <- kruskal.test(as.formula(paste0('activation ~ `',facet,'`')),data = df)
#   
#   plt <- df %>% filter(!is.na(!! f)) %>%
#     ggplot(aes(x=!! f,y=activation,fill=!! f))+geom_boxplot()+
#     scale_fill_manual(values = active_colors)+
#     annotate('text',-Inf, Inf,label=paste0('p = ',round(test$p.value,5)),hjust=-0.2,vjust=2)+
#     theme_bw(30) + 
#     labs(x=facet,y='Activation')+ggtitle(paste0(facet))+
#     theme(text = element_text(size=11), panel.background = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank())
#   
#   return(plt)
# }
# 
# 
# active_facets <- names(which(sapply(test,is.factor)))
# plts <- lapply(active_facets,facet_boxplot,df=test)
# 
# grobs_arranged <- do.call(grid.arrange, c(plts, ncol=4))
# 


# process_metagene <- function(metagene,refset){
#   if(is.character(metagene)){
#     genes <- metagene[metagene %in% refset]
#     sign <- rep(1,length(genes))
#     names(sign) <- genes
#   }else{
#     genes <- names(metagene)[names(metagene) %in% refset]
#     sign <- metagene[names(metagene) %in% refset]
#   }
#   
#   return(list(sign=sign,genes=genes))
# }
# 
# get_facets <- function(facets){
#   default_facets<- c('patient','overall_survival_months','overall_survival_years','vital_status')
#   return(list(facets=c(default_facets,facets),added_facets=facets))
# }
# 
# select_metagene <- function(df,metagene,facets=c()){
#   f <- get_facets(facets)
#   mg <- process_metagene(metagene,colnames(df))
#   
#   selected <- df %>% 
#     select_at(vars(c(mg$genes,f$facets)))
#   
#   return(selected)
# }
# 
# score_metagene <- function(df,metagene,facets=c()){
#   f <- get_facets(facets)
#   mg <- process_metagene(metagene,colnames(df))
#   
#   scored <- df %>%
#     pivot_longer(all_of(mg$genes)) %>%
#     group_by_at(vars(f$facets)) %>%
#     summarize(activation=sum(value*mg$sign[name])) %>%
#     ungroup %>% mutate_at(vars(f$added_facets),as.factor)
#   return(scored)
# }
