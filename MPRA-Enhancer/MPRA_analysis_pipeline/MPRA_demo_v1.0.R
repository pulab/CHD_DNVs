library(jmuOutlier) 
library(distr6)
library(ggplot2)
library(DESeq2)
library(Rfast)
library(dplyr)
library(reshape2)
library(fitdistrplus)
library(compositions)
library(stats)
library(scrime)
library(corrplot)
library(pheatmap)
library(fdrtool)
library(coin)
library(patchwork)
library(ggpubr)
library(cowplot)
library(grid)
library(stringr)
library(patchwork)
Enrich_score <- function(datasets, groups_col , value_col , order_col){
  ALL <- datasets[,value_col]
  Max <- length(ALL)
  datasets$Enrich_score <-0
  datasets$group_order <-0
  pvalue <- c()
  groups_names <-c()
  rows_number <-length(unique(datasets[,groups_col]))
  datasets2 <- datasets[FALSE,]
  #datasets2 <- datasets
  datasets3 <- datasets
  for (cols_name in unique(datasets[,groups_col])) {
    datasets3 <- datasets
    PC1 <- datasets3[datasets3[,groups_col]==cols_name,]
    p_PC1 <- perm.test(PC1[,value_col],ALL)$p.value
    #PC1$group_order <-0
    PC1$group_order <- rank(-as.numeric(as.character(PC1[,value_col])),ties.method = "first")
    Max_order <- max(PC1[,"group_order"])
    datasets3[rownames(PC1),"group_order"]<- PC1$group_order
    Max_pc1 <- length(PC1[,order_col])
    Max_pc1_order <- max(PC1[,order_col])
    for (order_name in c(1:Max_order)){
      left_value <- 0
      order_last <- order_name-1
      right_value<- datasets3[datasets3$group_order==order_name,order_col]
      if (order_name ==1){
        left_value <- 0
      }
      if (order_name >=2) { 
        left_value <- datasets3[datasets3$group_order==order_last,order_col]
      }
      datasets3$group_order[datasets3[,order_col]<right_value & datasets3[,order_col]> left_value]<- order_last
    }
    if (order_name==Max_order){
      right_value<- datasets3[datasets3$group_order==Max_order,order_col]
      datasets3$group_order[datasets3[,order_col]> right_value]<- Max_order
    }
    Max <- length(datasets3[,order_col])
    datasets3$Enrich_score <- datasets3[,"group_order"]/Max_pc1 - datasets3[,order_col]/Max

    pvalue <- c(pvalue, p_PC1)
    groups_names <- c(groups_names,cols_name)
    datasets3$score_group <- cols_name
    datasets3$final_ES <- datasets3$Enrich_score

    datasets2 <- rbind(datasets2,datasets3)
    
  }  
  dat_text <- data.frame(
    label = pvalue,
    score_group = groups_names
  )
  
  newlist <-list("dat_text"=dat_text, "datasets"=datasets2)
  return(newlist)
}
MPRA_SizeFactors <- function(datasets,groups,cols){
  filter_norm_data_nc <- datasets[groups,cols] %>% as.matrix()
  filter_norm_data_full_nc <- datasets[,cols] %>% as.matrix()
  filter_norm_data_nc <- filter_norm_data_nc[rowMins(filter_norm_data_nc,value = T)>0,]
  filter_norm_data_nc_ratio_pseudo_exp <- geometricmeanRow(filter_norm_data_nc)
  filter_norm_data_nc_ratio_scaled <- t(t(filter_norm_data_nc)/filter_norm_data_nc_ratio_pseudo_exp)
  filter_norm_data_nc_normalize_factor <- colMedians(filter_norm_data_nc_ratio_scaled)
  filter_norm_data_nc_normalize_exp <- t(t(filter_norm_data_full_nc)/filter_norm_data_nc_normalize_factor)
  return (filter_norm_data_nc_normalize_exp)
}
##################PCC plot##############################
panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits=2)
  txt <- paste0("R = ", r)
  cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor*1 )
}

# Customize upper panel
upper.panel<-function(x, y){
  points(x,y, pch = 19, col = my_cols[iris$Species],cex=0.1)
}
# Create the plots
pairs(filter_norm_data_nc_normalize_exp, 
      lower.panel = panel.cor,
      upper.panel = upper.panel)
#################DEseq2 active enhancers##########################################
library(DESeq2)
Active_enhancer<-function(filter_for_degs,coldata,condition="condition",RNA="RNA",DNA="DNA"){
  genes=as.matrix(filter_for_degs)
  dds <- DESeqDataSetFromMatrix(countData = genes, colData=coldata, design=~condition)
  dds <- DESeq(dds)
  res <- results( dds, contrast = c(condition,RNA,DNA) )
  #res2 <- results( dds, contrast = c("condition","d17_DNA_2","d17_RNA_2") )
  res$qvalue_0_05 <- 0
  res$qvalue_0_05[res$padj < 0.05 & res$log2FoldChange >0] <- 1
  return (res)
}
####################FDR 0.95 for negative control##################################
FDR_negative <- function(Negative_C,colname,distribution="norm"){
  fitg<-fitdist(Negative_C[,colname] ,distribution)
  
  FDR_0.95 <- qnorm(0.95, mean = fitg$estimate["mean"], sd = fitg$estimate["sd"])
  return(FDR_0.95)
}

##################Different enhancers #########################
Diff_enhancers<- function(data,conditon_1,conditon_2){
  pvalue <- c()
  for (i in 1:dim(data[,conditon_1])[1]) {
    ref_value <- as.numeric(data[,conditon_1][i,])
    alt_value <- as.numeric(data[,conditon_2][i,])
    pvalue_i <- t.test(ref_value, alt_value)
    pvalue <- c(pvalue, pvalue_i$p.value)
  }
  qvalue <- pvalue * length(pvalue) /rank(pvalue)
  data$C1 <- rowMeans(data[,conditon_1])
  data$C2 <- rowMeans(data[,conditon_1])
  data$qvalue <- as.numeric(as.character(qvalue))
  data$pvalue <- as.numeric(as.character(pvalue))
  return(data)
}
###############odds ratio for motifs##########################
get_fisher_p <- function(df){
  mat <- matrix(as.numeric(df[c(1:4)]), ncol=2)
  f <- fisher.test(as.table(mat), alt="two.sided")
  return(f$p.value)
}
get_fisher_odds <- function(df){
  mat <- matrix(as.numeric(df[c(1:4)]), ncol=2)
  f <- fisher.test(as.table(mat), alt="two.sided")
  return(f$estimate)
}



Get_odds_ratio <-function(data, guess_cols, background_cols) {
  Odd_result<-data[,c(guess_cols,background_cols)]
  SSA_matrix_merge=data.frame(row.names = rownames(data))
  for (cols_name in guess_cols) {
    total_number <-sum(data[,cols_name])
    for (bgc_name in background_cols){
      total_bg_number <- sum(data[,bgc_name])
      Odd_result$guss_res <-total_number- Odd_result[,cols_name]
      Odd_result$bdg_res <-total_bg_number- Odd_result[,bgc_name]
      Odd_input<-as.data.frame(Odd_result[,c(cols_name, "guss_res", bgc_name, "bdg_res")])
      res_list_p <- as.data.frame(apply(Odd_input, 1,  get_fisher_p))
      colnames(res_list_p)<-paste0(cols_name,"_",bgc_name,"_","Pvalue")
      res_list_odds <- as.data.frame(apply(Odd_input, 1,  get_fisher_odds))
      res_list_odds[rownames(Odd_input[Odd_input[,cols_name]==0 & Odd_input[,bgc_name]==0,]),1]<-NA
      colnames(res_list_odds)<-paste0(cols_name,"_",bgc_name,"_","odds")
      #odds_col_name <-paste0(cols_name,"_",bgc_name,"_","odds")
      #f <- fisher.test(as.table(mat), alt="two.sided")
      SSA_matrix_merge<-cbind(SSA_matrix_merge,res_list_odds)
      SSA_matrix_merge<-cbind(SSA_matrix_merge,res_list_p)
      ######rm both 0 value####################################
      
      
    }
  }
  return(SSA_matrix_merge)
}
#####################################################


#filter_norm_data_scaled_grouped_test_A <-Enrich_score(filter_norm_data_scaled_grouped, "group", "log_A_ratio", "order_log_A_mean_ratio")
#filter_norm_data_scaled_grouped_test <-Enrich_score(filter_norm_data_scaled_grouped, "group2", "log_d24_ratio", "order_d24_ratio")
#dat_text_A <-filter_norm_data_scaled_grouped_test$dat_text
#filter_norm_data_for_plot_A <- filter_norm_data_scaled_grouped_test$datasets