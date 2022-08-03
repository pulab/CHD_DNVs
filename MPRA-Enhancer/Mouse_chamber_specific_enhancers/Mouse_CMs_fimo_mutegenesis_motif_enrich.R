library(slingshot)
library(BUSpaRse)
library(tidyverse)
library(tidymodels)
library(Seurat)
library(scales)
library(viridis)
library(Matrix)
library(rsample)
library(dplyr)
library(Rfast)
library(pheatmap)
library(ranger)
library(reshape2)
setwd("~/Desktop/work/Bill/yangpo/MPRA_4")
#data1 <- read.table(file = "mutagenesis_WT-ALT_fimo_new.txt",header=F)
#data1 <- read.table(file = "common_mutagenesis_library_MPRA_result_motif_region_result.txt",header=T)
data1 <- read.table(file = "mutagenesis_WT-ALT_fimo_new_2_with_expressed_new.txt",header=T)
data2 <- read.table(file = "AV_mutagenesis_MPRA_4_normalized_log2_ratio_fold_pvalue_add1_with_AV_diff_new_region.txt",header=T)
colnames(data1) <-c("motif_name","motif_region_name","region_name","nlog_qvalue_in_WT","nlog_pvalue_in_WT","nlog_qvalue_in_ALT","nlog_pvalue_in_ALT")
filter_norm_data1 <- data1[as.numeric(data1[,"nlog_pvalue_in_WT"]) >=4 |
                                   as.numeric(data1[,"nlog_pvalue_in_ALT"]) >=4 ,]
rownames(filter_norm_data1) <- filter_norm_data1$motif_region_name
filter_norm_data1_by <- filter_norm_data1 %>% group_by(motif_region_name)
filter_norm_data1_by_sum <- filter_norm_data1_by %>% summarise(
  nlog_qvalue_in_WT = sum(nlog_qvalue_in_WT),
  nlog_pvalue_in_WT = sum(nlog_pvalue_in_WT),
  nlog_qvalue_in_ALT = sum(nlog_qvalue_in_ALT),
  nlog_pvalue_in_ALT = sum(nlog_pvalue_in_ALT)
  )
filter_norm_data1_by_sum$motif_name <- filter_norm_data1_by[filter_norm_data1_by_sum$motif_region_name %in% filter_norm_data1_by$motif_region_name,]$motif_name
filter_norm_data1_by_sum$region_name <- filter_norm_data1_by[filter_norm_data1_by_sum$motif_region_name %in% filter_norm_data1_by$motif_region_name,]$region_name
data2$region_name <- rownames(data2)
data_u <-merge(filter_norm_data1_by_sum,data2, by="region_name")

#list_file <- read.table(file = "mutagenesis_fold_final_sig.txt",header = T,row.names = 1)
#list_file <- read.table(file = "mutagenesis_fold_final_sig.txt",header = T,row.names = 1)
#rownames(data1) = data1[,1]
#data_u <- data1
########################TBX and ERR #########################
data_TBX <- data_u[data_u$motif_name=="TBX3_HUMAN.H11MO.0.C",]
data_TBX <- data_u[data_u$motif_name=="ERR2_MOUSE.H11MO.0.A",]
#data_TBX_for_plot <- data_TBX[,c("REF_A","ALT_A","REF_V","ALT_V","nlog_pvalue_in_WT","nlog_pvalue_in_ALT")]
filter_norm_data_tbx <- data_TBX[as.numeric(data_TBX[,"nlog_pvalue_in_WT"]) >=4 |
                                       as.numeric(data_TBX[,"nlog_pvalue_in_ALT"]) >=4 ,]
data_TBX_for_plot <- filter_norm_data_tbx[,c("region_name","REF_A","ALT_A","REF_V","ALT_V","final_sig_A","final_sig_V","nlog_pvalue_in_WT","nlog_pvalue_in_ALT")]
data_TBX_for_plot$final_sig_WT <- data2[data_TBX_for_plot$region,]$final_sig_WT_AV
data_TBX_for_plot$final_sig_ALT <- data2[data_TBX_for_plot$region,]$final_sig_ALT_AV

data_TBX_for_plot <- data_TBX_for_plot[as.numeric(data_TBX_for_plot[,"final_sig_A"]) >0 |
                                   as.numeric(data_TBX_for_plot[,"final_sig_V"]) >0 |
                                ! is.na(data_TBX_for_plot[,"final_sig_WT"])|
                                  ! is.na(data_TBX_for_plot[,"final_sig_ALT"]) ,]
rownames(data_TBX_for_plot) = data_TBX_for_plot[,1]
data_TBX_for_plot = data_TBX_for_plot[,2:length(data_TBX_for_plot[1,])]
pheatmap(data_TBX_for_plot, cluster_rows=T,cluster_cols=F,
         
         breaks=seq(-5,5,0.1),color = colorRampPalette(c("navy", "white","gold"))(100),
         show_colnames = T, fontsize=7)
#####################################################################################
filter_norm_data <- data_u
dim(data_u)
colnames(data_u)
filter_norm_data <- filter_norm_data[as.numeric(filter_norm_data[,"nlog_pvalue_in_WT"]) >=4 |
                                       as.numeric(filter_norm_data[,"nlog_pvalue_in_ALT"]) >=4 ,]
filter_norm_data <- filter_norm_data[(as.numeric(filter_norm_data[,"nlog_pvalue_in_WT"]) -
                                       as.numeric(filter_norm_data[,"nlog_pvalue_in_ALT"]) >=2) | as.numeric(filter_norm_data[,"nlog_pvalue_in_WT"]) -
                                       as.numeric(filter_norm_data[,"nlog_pvalue_in_ALT"]) <=-2,]
data_motif_A <- read.table(file="seedmotif_gene_expressed_in_A_tpm5_cluster.txt",header=F,row.names = 1)
data_motif_V <- read.table(file="seedmotif_gene_expressed_in_V_tpm5_cluster.txt",header=F,row.names = 1)
colnames(data_motif_A)<-c("gene","Cluster")
colnames(data_motif_V)<-c("gene","Cluster")
data_motif_A$motif_name<-rownames(data_motif_A)
data_motif_V$motif_name<-rownames(data_motif_V)
data_motif_A<-data_motif_A[,c("motif_name","gene")]
data_motif_V<-data_motif_V[,c("motif_name","gene")]
#data_motif_V <- read.table(file="data_motif_V.txt",header=T,row.names = 1)
filter_norm_data_genes <- merge(filter_norm_data, data_motif_A, by="motif_name", all.x=T)
filter_norm_data_genes_2 <- merge(filter_norm_data_genes, data_motif_V, by="motif_name", all.x=T)
filter_norm_data$motif_name_1<-str_split_fixed(filter_norm_data$motif_name, pattern="_HUMAN",2)[,1]
filter_norm_data$motif_name_2<-str_split_fixed(filter_norm_data$motif_name_1, pattern="_MOUSE",2)[,1]
region_motif <- aggregate(gene.x~region_name, data = filter_norm_data_genes_2, paste0, collapse=",")
region_motif <- aggregate(gene.y~region_name, data = filter_norm_data_genes_2, paste0, collapse=",")
#region_motif <- aggregate(c(gene.x, gene.y)~region_name, data = filter_norm_data_genes_2, paste0, collapse=",")
region_motif <- aggregate(motif_name_2~region_name, data = filter_norm_data, paste0, collapse=",")
write.table(region_motif,file="region_changed_motifs_gene_A_new.list",sep = "\t",quote = FALSE)
write.table(region_motif,file="region_changed_motifs_gene_V_new.list",sep = "\t",quote = FALSE)
write.table(region_motif,file="region_changed_motifs_new.list",sep = "\t",quote = FALSE)
#filter_norm_data <- filter_norm_data[as.numeric(filter_norm_data[,4]) >=1 |
#                                       as.numeric(filter_norm_data[,6]) >=1 ,]
filter_norm_data$foldpvalue <- filter_norm_data$nlog_pvalue_in_ALT - filter_norm_data$nlog_pvalue_in_WT
filter_norm_data <-filter_norm_data[filter_norm_data$foldpvalue>=4 | filter_norm_data$foldpvalue<=-4,]
dim(filter_norm_data[filter_norm_data$foldpvalue <= -4,])
head(filter_norm_data)
#filter_norm_data$fold <- list_file[as.character(filter_norm_data$V2),"fold"]
#dim(filter_norm_data[!is.na(filter_norm_data$fold),])
#filter_norm_data <-filter_norm_data[!is.na(filter_norm_data$fold),]
###########################################################################3
#Namelist <- read.table(file="~/Desktop/work/Bill/Fengxiao/motifs/db.list",header=T,row.names = 1)
#final_sig_list <- read.table(file = "../new_deepbind/CHD_variant_fold_pvalue_final_sig.txt",header=T,row.names = 1)
for_heat <- as.data.frame(filter_norm_data)
#final_sig_list[final_sig_list$final_sig ==2,]
for_heat_G_A <- filter_norm_data[filter_norm_data$final_sig_A ==2,]
for_heat_L_A <- filter_norm_data[filter_norm_data$final_sig_A ==1,]
for_heat_no_A <- filter_norm_data[filter_norm_data$final_sig_A ==0,]
for_heat_G_V <- filter_norm_data[filter_norm_data$final_sig_V ==2,]
for_heat_L_V <- filter_norm_data[filter_norm_data$final_sig_V ==1,]
for_heat_no_V <- filter_norm_data[filter_norm_data$final_sig_V ==0,]
for_heat_all <- filter_norm_data
##############################################################################################################

for_heat_G_M_A <- as.matrix(for_heat_G_A[,c("motif_region_name","motif_name","region_name","final_sig_A",
                                            "nlog_pvalue_in_WT","nlog_pvalue_in_WT","foldpvalue")])
for_heat_L_M_A <- as.matrix(for_heat_L_A[,c("motif_region_name","motif_name","region_name","final_sig_A",
                                            "nlog_pvalue_in_WT","nlog_pvalue_in_WT","foldpvalue")])
for_heat_no_M_A <- as.matrix(for_heat_no_A[,c("motif_region_name","motif_name","region_name","final_sig_A",
                                              "nlog_pvalue_in_WT","nlog_pvalue_in_WT","foldpvalue")])
for_heat_all_M_A <- as.matrix(for_heat_all[,c("motif_region_name","motif_name","region_name","final_sig_A",
                                              "nlog_pvalue_in_WT","nlog_pvalue_in_WT","foldpvalue")])

rownames(for_heat_G_M_A) = for_heat_G_A$motif_name
rownames(for_heat_L_M_A) = for_heat_L_A$motif_name
rownames(for_heat_no_M_A) = for_heat_no_A$motif_name
rownames(for_heat_all_M_A) = for_heat_all$motif_name
melt_for_heat_G_A <- unique(melt(for_heat_G_A[,c("region_name","motif_name","foldpvalue")]))
melt_for_heat_G_M_A <- acast(melt_for_heat_G_A[,c(1,2,4)],formula =motif_name~region_name, fill = 0 )
melt_for_heat_L_A <- unique(melt(for_heat_L_A[,c("region_name","motif_name","foldpvalue")]))
melt_for_heat_L_M_A <- acast(melt_for_heat_L_A[,c(1,2,4)],formula =motif_name~region_name, fill = 0 )
melt_for_heat_no_A <- unique(melt(for_heat_no_A[,c("region_name","motif_name","foldpvalue")]))
melt_for_heat_no_M_A <- acast(melt_for_heat_no_A[,c(1,2,4)],formula =motif_name~region_name, fill = 0 )
melt_for_heat_all_A <- unique(melt(for_heat_all[,c("region_name","motif_name","foldpvalue")]))
melt_for_heat_all_M_A <- acast(melt_for_heat_all_A[,c(1,2,4)],formula =motif_name~region_name, fill = 0 )
##################################################################################################
for_heat_G_M_V <- as.matrix(for_heat_G_V[,c("motif_region_name","motif_name","region_name","final_sig_V",
                                            "nlog_pvalue_in_WT","nlog_pvalue_in_WT","foldpvalue")])
for_heat_L_M_V <- as.matrix(for_heat_L_V[,c("motif_region_name","motif_name","region_name","final_sig_V",
                                            "nlog_pvalue_in_WT","nlog_pvalue_in_WT","foldpvalue")])
for_heat_no_M_V <- as.matrix(for_heat_no_V[,c("motif_region_name","motif_name","region_name","final_sig_V",
                                              "nlog_pvalue_in_WT","nlog_pvalue_in_WT","foldpvalue")])
for_heat_all_M_V <- as.matrix(for_heat_all[,c("motif_region_name","motif_name","region_name","final_sig_V",
                                              "nlog_pvalue_in_WT","nlog_pvalue_in_WT","foldpvalue")])

rownames(for_heat_G_M_V) = for_heat_G_V$motif_name
rownames(for_heat_L_M_V) = for_heat_L_V$motif_name
rownames(for_heat_no_M_V) = for_heat_no_V$motif_name
rownames(for_heat_all_M_V) = for_heat_all$motif_name
melt_for_heat_G_V <- unique(melt(for_heat_G_V[,c("region_name","motif_name","foldpvalue")]))
melt_for_heat_G_M_V <- acast(melt_for_heat_G_V[,c(1,2,4)],formula =motif_name~region_name, fill = 0 )
melt_for_heat_L_V <- unique(melt(for_heat_L_V[,c("region_name","motif_name","foldpvalue")]))
melt_for_heat_L_M_V <- acast(melt_for_heat_L_V[,c(1,2,4)],formula =motif_name~region_name, fill = 0 )
melt_for_heat_no_V <- unique(melt(for_heat_no_V[,c("region_name","motif_name","foldpvalue")]))
melt_for_heat_no_M_V <- acast(melt_for_heat_no_V[,c(1,2,4)],formula =motif_name~region_name, fill = 0 )
melt_for_heat_all_V <- unique(melt(for_heat_all[,c("region_name","motif_name","foldpvalue")]))
melt_for_heat_all_M_V <- acast(melt_for_heat_all_V[,c(1,2,4)],formula =motif_name~region_name, fill = 0)

###################################################################



##################################################################
dim(for_heat_G_M_A)
dim(for_heat_L_M_A)
dim(for_heat_no_M_A)
dim(for_heat_all_M_V)

dim(for_heat_G_M_V)
dim(for_heat_L_M_V)
dim(for_heat_no_M_V)
dim(for_heat_all_M_V)

dim(melt_for_heat_G_M_A)
dim(melt_for_heat_L_M_A)
dim(melt_for_heat_no_M_A)
dim(melt_for_heat_all_M_V)
#annotation <- as.data.frame(for_heat_L_M[,6])
#rownames(annotation)<-rownames(for_heat_L_M)

#pheatmap(for_heat_G_M[,c(2,4)], cluster_rows=T,cluster_cols=F, annotation_row = annotation,
         
#         breaks=seq(-10,10,0.2),color = colorRampPalette(c("navy", "white","gold"))(100), show_colnames = F)

pheatmap(melt_for_heat_G_M_A, cluster_rows=T,cluster_cols=T,
         
         breaks=seq(-5,5,0.1),color = colorRampPalette(c("navy", "white","gold"))(100), show_colnames = F, fontsize=7)

pheatmap(melt_for_heat_L_M_A, cluster_rows=T,cluster_cols=T,
         
         breaks=seq(-5,5,0.1),color = colorRampPalette(c("navy", "white","gold"))(100), show_colnames = F, fontsize=7)
pheatmap(melt_for_heat_no_M_A, cluster_rows=T,cluster_cols=T,
         
         breaks=seq(-5,5,0.1),color = colorRampPalette(c("navy", "white","gold"))(100), show_colnames = F, fontsize=5)

pheatmap(melt_for_heat_all_M_V, cluster_rows=T,cluster_cols=T,
         
         breaks=seq(-5,5,0.1),color = colorRampPalette(c("navy", "white","gold"))(100), show_colnames = F, fontsize=5)
################################################################################################3
plot_gt_ls <- data.frame(gain_motif = rowSums(melt_for_heat_G_M_A >=4), lost_motif = rowSums(melt_for_heat_G_M_A <=-4))
plot_gt_ls <- as.matrix(plot_gt_ls)
rownames(plot_gt_ls) <- rownames(melt_for_heat_G_M_A)
pheatmap(plot_gt_ls, cluster_rows=T,cluster_cols=F, 
         color = colorRampPalette(c("white","gold"))(100), show_colnames = T)

plot_gt_ls2 <- data.frame(gain_motif = rowSums(melt_for_heat_L_M_A >=4), lost_motif = rowSums(melt_for_heat_L_M_A <=-4))
plot_gt_ls2 <- as.matrix(plot_gt_ls2)
rownames(plot_gt_ls2) <- rownames(melt_for_heat_L_M_A)
pheatmap(plot_gt_ls2, cluster_rows=T,cluster_cols=F, 
         color = colorRampPalette(c( "white","gold"))(100), show_colnames = T)
plot_gt_ls <-as.data.frame(plot_gt_ls)
plot_gt_ls2 <-as.data.frame(plot_gt_ls2)
plot_gt_ls3 <- merge(plot_gt_ls,plot_gt_ls2,by=0,all=T)
plot_gt_ls3[is.na(plot_gt_ls3)] <- 0
rownames(plot_gt_ls3) = plot_gt_ls3[,1]
plot_gt_ls3 = plot_gt_ls3[,2:length(plot_gt_ls3[1,])]
pheatmap(as.matrix(plot_gt_ls3), cluster_rows=T,cluster_cols=F, 
         color = colorRampPalette(c( "white","gold"))(100), show_colnames = T,fontsize = 8)
#######################################################################
data_motif_A <- read.table(file="data_motif_A.txt",header=T,row.names = 1)

plot_gt_ls4 <- merge(plot_gt_ls3, data_motif_A, by.x="row.names",by.y="motif_name") 
plot_gt_ls4_2 <-as.matrix(plot_gt_ls4[,c(2,3,4,5)])
rownames(plot_gt_ls4_2) <- plot_gt_ls4[,"gene"]
pdf("motif_number_lef_GoF_rigth_LoF_for_A.pdf",width =8,height = 6)
pheatmap(as.matrix(plot_gt_ls4_2), cluster_rows=T,cluster_cols=F,color = colorRampPalette(c( "blue","white","gold"))(100),breaks = seq(-5,5,0.1), 
         show_colnames = T,fontsize = 6)
dev.off()
#############################################################

plot_gt_ls5 <- data.frame(gain_motif = rowSums(melt_for_heat_no_M_A >=4), lost_motif = rowSums(melt_for_heat_no_M_A <=-4))
plot_gt_ls5 <- as.matrix(plot_gt_ls5)
rownames(plot_gt_ls5) <- rownames(melt_for_heat_no_M_A)
pheatmap(plot_gt_ls5, cluster_rows=T,cluster_cols=F, 
         color = colorRampPalette(c("white","gold"))(100), show_colnames = T)

plot_gt_ls6 <- merge(plot_gt_ls3, plot_gt_ls5,by=0,all=T)
plot_gt_ls6[is.na(plot_gt_ls6)] <- 0
plot_gt_ls6<-plot_gt_ls6[rowMaxs(as.matrix(plot_gt_ls6[,c(2:7)]),value = T)>1,]
rownames(plot_gt_ls6) = plot_gt_ls6[,1]
#plot_gt_ls6<-plot_gt_ls6[rowMaxs(plot_gt_ls6)>1,]
plot_gt_ls6 = plot_gt_ls6[,2:length(plot_gt_ls6[1,])]
################################################################################

plot_gt_ls46<-plot_gt_ls6

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
plot_gt_ls46_lost<-Get_odds_ratio(plot_gt_ls46,c("lost_motif.x","lost_motif.y"),"lost_motif")
plot_gt_ls46_gain<-Get_odds_ratio(plot_gt_ls46,c("gain_motif.x","gain_motif.y"),"gain_motif")
plot_gt_ls47<- cbind(plot_gt_ls46,plot_gt_ls46_gain)
plot_gt_ls47<- cbind(plot_gt_ls47,plot_gt_ls46_lost)
head(plot_gt_ls47[is.na(plot_gt_ls47$gain_motif.x_gain_motif_odds),])
Melt_plot_gt_ls47_G_G <-plot_gt_ls47[,c("gain_motif.x_gain_motif_odds","gain_motif.x_gain_motif_Pvalue")]
Melt_plot_gt_ls47_G_G$type <-"GoF_enhancers_Gain_motif"
colnames(Melt_plot_gt_ls47_G_G)<- c("Odds","Pvalue","type")
Melt_plot_gt_ls47_G_G$motif<-rownames(Melt_plot_gt_ls47_G_G)
Melt_plot_gt_ls47_G_L <-plot_gt_ls47[,c("lost_motif.x_lost_motif_odds","lost_motif.x_lost_motif_Pvalue")]
Melt_plot_gt_ls47_G_L$type <-"GoF_enhancers_lost_motif"
colnames(Melt_plot_gt_ls47_G_L)<- c("Odds","Pvalue","type")
Melt_plot_gt_ls47_G_L$motif<-rownames(Melt_plot_gt_ls47_G_L)
Melt_plot_gt_ls47_L_G <-plot_gt_ls47[,c("gain_motif.y_gain_motif_odds","gain_motif.y_gain_motif_Pvalue")]
Melt_plot_gt_ls47_L_G$type <-"LoF_enhancers_Gain_motif"
colnames(Melt_plot_gt_ls47_L_G)<- c("Odds","Pvalue","type")
Melt_plot_gt_ls47_L_G$motif<-rownames(Melt_plot_gt_ls47_L_G)
Melt_plot_gt_ls47_L_L <-plot_gt_ls47[,c("lost_motif.y_lost_motif_odds","lost_motif.y_lost_motif_Pvalue")]
Melt_plot_gt_ls47_L_L$type <-"LoF_enhancers_lost_motif"
colnames(Melt_plot_gt_ls47_L_L)<- c("Odds","Pvalue","type")
Melt_plot_gt_ls47_L_L$motif<-rownames(Melt_plot_gt_ls47_L_L)

for_plot_gt_ls47 <-rbind(Melt_plot_gt_ls47_G_G,Melt_plot_gt_ls47_G_L)
for_plot_gt_ls47<-rbind(for_plot_gt_ls47,Melt_plot_gt_ls47_L_G)
for_plot_gt_ls47<-rbind(for_plot_gt_ls47,Melt_plot_gt_ls47_L_L)
#for_plot_gt_ls47$motif<-rownames(for_plot_gt_ls47)
for_plot_gt_ls48<-for_plot_gt_ls47
for_plot_gt_ls47<-for_plot_gt_ls48
for_plot_gt_ls47[is.infinite(for_plot_gt_ls47$Odds),]$Odds<-10
for_plot_gt_ls47 <- for_plot_gt_ls47[!is.na(for_plot_gt_ls47$Odds),]
for_plot_gt_ls47[for_plot_gt_ls47$Odds>10,]$Odds<-10
#for_plot_gt_ls47<-for_plot_gt_ls47[which(for_plot_gt_ls47$Odds>1),]
#for_plot_gt_ls47<-for_plot_gt_ls47[which(for_plot_gt_ls47$Pvalue<=0.1),]
ggplot(for_plot_gt_ls47,aes(x=type,y = motif,size=Odds, color=-log10(Pvalue)))+ geom_point()+
  scale_color_gradient2(midpoint=1, low="white", mid="blue",
                        high="red", space ="Lab" )+
  theme_bw()
for_plot_gt_ls47_A <-for_plot_gt_ls47
################################################################################
dim(for_heat_G_M_A)
dim(for_heat_L_M_A)
dim(for_heat_no_M_A)
dim(for_heat_all_M_V)
#length(unique(for_heat_G_M_A[,"motif_region_name"]))

plot_gt_ls7 <-plot_gt_ls6
plot_gt_ls7[,c(1:2)]<- plot_gt_ls7[,c(1:2)]/24
plot_gt_ls7[,c(3:4)]<- plot_gt_ls7[,c(3:4)]/74
plot_gt_ls7[,c(5:6)]<- plot_gt_ls7[,c(5:6)]/1302

pheatmap(as.matrix(plot_gt_ls6), cluster_rows=T,cluster_cols=F, 
         color = colorRampPalette(c( "white","gold"))(100), show_colnames = T,fontsize = 8)
pheatmap(as.matrix(plot_gt_ls7), cluster_rows=T,cluster_cols=F, 
         color = colorRampPalette(c( "white","gold"))(100), show_colnames = T,fontsize = 8)

plot_gt_ls8 <- as.data.frame(plot_gt_ls7)
plot_gt_ls8$GoF_diff <- plot_gt_ls8$gain_motif.x-plot_gt_ls8$lost_motif.x
plot_gt_ls8$LoF_diff <- plot_gt_ls8$lost_motif.y-plot_gt_ls8$gain_motif.y
plot_gt_ls8$Nochange <- abs(plot_gt_ls8$gain_motif-plot_gt_ls8$lost_motif)
##################
plot_gt_ls8$GoF_diff <- plot_gt_ls8$gain_motif.x+plot_gt_ls8$lost_motif.x
plot_gt_ls8$LoF_diff <- plot_gt_ls8$lost_motif.y+plot_gt_ls8$gain_motif.y
plot_gt_ls8$Nochange <- abs(plot_gt_ls8$gain_motif+plot_gt_ls8$lost_motif)
##########
#plot_gt_ls8$Nochange <- abs(plot_gt_ls8$gain_motif-plot_gt_ls8$lost_motif)
plot_gt_ls8$G_fold_to_backgroud <- log2(plot_gt_ls8$GoF_diff/(plot_gt_ls8$Nochange+0.0000001)+1)
plot_gt_ls8$L_fold_to_backgroud <- log2(plot_gt_ls8$LoF_diff/(plot_gt_ls8$Nochange+0.0000001)+1)
plot_gt_ls9 <- plot_gt_ls8[,10:11]
plot_gt_ls9_A <-plot_gt_ls9
plot_gt_ls9[is.na(plot_gt_ls9)] <- 0
pheatmap(as.matrix(plot_gt_ls9), cluster_rows=T,cluster_cols=F,color = colorRampPalette(c( "blue","white","gold"))(100),breaks = seq(-10,10,0.2), show_colnames = T,fontsize = 8)

plot_gt_ls10 <- plot_gt_ls9[plot_gt_ls9$G_fold_to_backgroud !=0 | plot_gt_ls9$L_fold_to_backgroud !=0 , ]
pheatmap(as.matrix(plot_gt_ls10), cluster_rows=T,cluster_cols=F,color = colorRampPalette(c( "blue","white","gold"))(100),breaks = seq(-10,10,0.2), 
         show_colnames = T,fontsize = 6)

data_motif_A_cluster <- read.table(file="seedmotif_name_genenames_mouse_with_cluster.txt",header=F,row.names = 1)
plot_gt_ls11 <- merge(plot_gt_ls10, data_motif_A, by.x="row.names",by.y="motif_name") 
plot_gt_ls12 <-as.matrix(plot_gt_ls11[,c(2,3)])
rownames(plot_gt_ls12) <- plot_gt_ls11[,"gene"]
plot_gt_ls_full <- merge(plot_gt_ls11, data_motif_A_cluster, by.x="Row.names",by.y="row.names")
plot_gt_ls_full$full_name <- paste(plot_gt_ls_full$Row.names,"~",plot_gt_ls_full$V3,"~",plot_gt_ls_full$gene)
plot_gt_ls12_full <-as.matrix(plot_gt_ls_full[,c(2,3)])
rownames(plot_gt_ls12_full) <- plot_gt_ls_full[,"full_name"]
pdf("motif_function_score_for_A_add.pdf",width =9,height = 6)
#pheatmap(as.matrix(plot_gt_ls12), cluster_rows=T,cluster_cols=F,color = colorRampPalette(c( "blue","white","gold"))(100),breaks = seq(-10,10,0.2), 
#         show_colnames = T,fontsize = 6)
row_names<-pheatmap(as.matrix(plot_gt_ls12_full), cluster_rows=T,cluster_cols=F,color = colorRampPalette(c( "blue","white","gold"))(100),breaks = seq(-10,10,0.2), 
         show_colnames = T,fontsize = 6)

## Extract the right grob
grob_classes <- purrr::map(row_names$gtable$grobs, class)
idx_grob <- which(purrr::map_lgl(grob_classes, function(cl) 'gTree' %in% cl))[1]
grob_names <- names(row_names$gtable$grobs[[idx_grob]]$children)
idx_rect <- grob_names[grep('rect', grob_names)][1]

## Remove borders around cells
row_names$gtable$grobs[[idx_grob]]$children[[idx_rect]]$gp$col <- "black"
row_names$gtable$grobs[[idx_grob]]$children[[idx_rect]]$gp$lwd <- 0.1

## Plot result
graphics::plot.new()
print(row_names)
dev.off()
dev.off()

plot_gt_ls15 <-plot_gt_ls6
plot_gt_ls15[,c(1)]<- (plot_gt_ls15[,c(1)])/(plot_gt_ls15[,c(5)]+plot_gt_ls15[,c(1)])
plot_gt_ls15[,c(2)]<- (plot_gt_ls15[,c(2)])/(plot_gt_ls15[,c(6)]+plot_gt_ls15[,c(2)])
plot_gt_ls15[,c(3)]<- (plot_gt_ls15[,c(3)])/(plot_gt_ls15[,c(5)]+plot_gt_ls15[,c(3)])
plot_gt_ls15[,c(4)]<- (plot_gt_ls15[,c(4)])/(plot_gt_ls15[,c(6)]+plot_gt_ls15[,c(4)])
plot_gt_ls15[is.na(plot_gt_ls15)] <- 0
plot_gt_ls16 <-plot_gt_ls15[plot_gt_ls15$gain_motif.x >=0.15 | plot_gt_ls15$lost_motif.x >=0.15 |plot_gt_ls15$gain_motif.y >=0.15 | plot_gt_ls15$lost_motif.y >=0.15,
                            c(1:4)]
plot_gt_ls17 <- merge(plot_gt_ls16, data_motif_A, by.x="row.names",by.y="motif_name") 
plot_gt_ls18 <-as.matrix(plot_gt_ls17[,c(2,3,4,5)])
rownames(plot_gt_ls18) <- plot_gt_ls17[,"gene"]
pdf("motif_persentage_left_GoF_right_LoF_for_A.pdf",width =5,height = 6)
pheatmap(as.matrix(plot_gt_ls18), cluster_rows=T,cluster_cols=F,color = colorRampPalette(c( "blue","white","gold"))(100),breaks = seq(-1,1,0.02), 
         show_colnames = T,fontsize = 6)
dev.off()

#########################################################################################
################################################################################################3
plot_gt_ls <- data.frame(gain_motif = rowSums(melt_for_heat_G_M_V >=4), lost_motif = rowSums(melt_for_heat_G_M_V <=-4))
plot_gt_ls <- as.matrix(plot_gt_ls)
rownames(plot_gt_ls) <- rownames(melt_for_heat_G_M_V)
pheatmap(plot_gt_ls, cluster_rows=T,cluster_cols=F, 
         color = colorRampPalette(c("white","gold"))(100), show_colnames = T)

plot_gt_ls2 <- data.frame(gain_motif = rowSums(melt_for_heat_L_M_V >=4), lost_motif = rowSums(melt_for_heat_L_M_V <=-4))
plot_gt_ls2 <- as.matrix(plot_gt_ls2)
rownames(plot_gt_ls2) <- rownames(melt_for_heat_L_M_V)
pheatmap(plot_gt_ls2, cluster_rows=T,cluster_cols=F, 
         color = colorRampPalette(c( "white","gold"))(100), show_colnames = T)
plot_gt_ls <-as.data.frame(plot_gt_ls)
plot_gt_ls2 <-as.data.frame(plot_gt_ls2)
plot_gt_ls3 <- merge(plot_gt_ls,plot_gt_ls2,by=0,all=T)
plot_gt_ls3[is.na(plot_gt_ls3)] <- 0
rownames(plot_gt_ls3) = plot_gt_ls3[,1]
plot_gt_ls3 = plot_gt_ls3[,2:length(plot_gt_ls3[1,])]
pheatmap(as.matrix(plot_gt_ls3), cluster_rows=T,cluster_cols=F, 
         color = colorRampPalette(c( "white","gold"))(100), show_colnames = T,fontsize = 8)
#######################################################################
data_motif_V <- read.table(file="data_motif_V.txt",header=T,row.names = 1)
plot_gt_ls4 <- merge(plot_gt_ls3, data_motif_V, by.x="row.names",by.y="motif_name") 
plot_gt_ls4_2 <-as.matrix(plot_gt_ls4[,c(2,3,4,5)])
rownames(plot_gt_ls4_2) <- plot_gt_ls4[,"gene"]
pdf("motif_number_lef_GoF_rigth_LoF_for_V.pdf",width =8,height = 6)
pheatmap(as.matrix(plot_gt_ls4_2), cluster_rows=T,cluster_cols=F,color = colorRampPalette(c( "blue","white","gold"))(100),breaks = seq(-5,5,0.1), 
         show_colnames = T,fontsize = 6)
dev.off()
#############################################################
plot_gt_ls5 <- data.frame(gain_motif = rowSums(melt_for_heat_no_M_V >=4), lost_motif = rowSums(melt_for_heat_no_M_V <=-4))
plot_gt_ls5 <- as.matrix(plot_gt_ls5)
rownames(plot_gt_ls5) <- rownames(melt_for_heat_no_M_V)
pheatmap(plot_gt_ls5, cluster_rows=T,cluster_cols=F, 
         color = colorRampPalette(c("white","gold"))(100), show_colnames = T)

plot_gt_ls6 <- merge(plot_gt_ls3, plot_gt_ls5,by=0,all=T)
plot_gt_ls6[is.na(plot_gt_ls6)] <- 0
plot_gt_ls6<-plot_gt_ls6[rowMaxs(as.matrix(plot_gt_ls6[,c(2:7)]),value = T)>1,]
rownames(plot_gt_ls6) = plot_gt_ls6[,1]
plot_gt_ls6 = plot_gt_ls6[,2:length(plot_gt_ls6[1,])]
################################################################################
total_motif_change_matrix<-colSums(plot_gt_ls6)
plot_gt_ls46<-plot_gt_ls6

plot_gt_ls46_lost<-Get_odds_ratio(plot_gt_ls46,c("lost_motif.x","lost_motif.y"),"lost_motif")
plot_gt_ls46_gain<-Get_odds_ratio(plot_gt_ls46,c("gain_motif.x","gain_motif.y"),"gain_motif")
plot_gt_ls47<- cbind(plot_gt_ls46,plot_gt_ls46_gain)
plot_gt_ls47<- cbind(plot_gt_ls47,plot_gt_ls46_lost)
head(plot_gt_ls47[is.na(plot_gt_ls47$gain_motif.x_gain_motif_odds),])
Melt_plot_gt_ls47_G_G <-plot_gt_ls47[,c("gain_motif.x_gain_motif_odds","gain_motif.x_gain_motif_Pvalue")]
Melt_plot_gt_ls47_G_G$type <-"GoF_enhancers_Gain_motif"
colnames(Melt_plot_gt_ls47_G_G)<- c("Odds","Pvalue","type")
Melt_plot_gt_ls47_G_G$motif<-rownames(Melt_plot_gt_ls47_G_G)
Melt_plot_gt_ls47_G_L <-plot_gt_ls47[,c("lost_motif.x_lost_motif_odds","lost_motif.x_lost_motif_Pvalue")]
Melt_plot_gt_ls47_G_L$type <-"GoF_enhancers_lost_motif"
colnames(Melt_plot_gt_ls47_G_L)<- c("Odds","Pvalue","type")
Melt_plot_gt_ls47_G_L$motif<-rownames(Melt_plot_gt_ls47_G_L)
Melt_plot_gt_ls47_L_G <-plot_gt_ls47[,c("gain_motif.y_gain_motif_odds","gain_motif.y_gain_motif_Pvalue")]
Melt_plot_gt_ls47_L_G$type <-"LoF_enhancers_Gain_motif"
colnames(Melt_plot_gt_ls47_L_G)<- c("Odds","Pvalue","type")
Melt_plot_gt_ls47_L_G$motif<-rownames(Melt_plot_gt_ls47_L_G)
Melt_plot_gt_ls47_L_L <-plot_gt_ls47[,c("lost_motif.y_lost_motif_odds","lost_motif.y_lost_motif_Pvalue")]
Melt_plot_gt_ls47_L_L$type <-"LoF_enhancers_lost_motif"
colnames(Melt_plot_gt_ls47_L_L)<- c("Odds","Pvalue","type")
Melt_plot_gt_ls47_L_L$motif<-rownames(Melt_plot_gt_ls47_L_L)


for_plot_gt_ls47 <-rbind(Melt_plot_gt_ls47_G_G,Melt_plot_gt_ls47_G_L)
for_plot_gt_ls47<-rbind(for_plot_gt_ls47,Melt_plot_gt_ls47_L_G)
for_plot_gt_ls47<-rbind(for_plot_gt_ls47,Melt_plot_gt_ls47_L_L)
#for_plot_gt_ls47$motif<-rownames(for_plot_gt_ls47)
for_plot_gt_ls48<-for_plot_gt_ls47
for_plot_gt_ls47<-for_plot_gt_ls48
for_plot_gt_ls47[is.infinite(for_plot_gt_ls47$Odds),]$Odds<-10
for_plot_gt_ls47 <- for_plot_gt_ls47[!is.na(for_plot_gt_ls47$Odds),]
for_plot_gt_ls47[for_plot_gt_ls47$Odds>10,]$Odds<-10
#for_plot_gt_ls47<-for_plot_gt_ls47[which(for_plot_gt_ls47$Odds>1),]
#for_plot_gt_ls47<-for_plot_gt_ls47[which(for_plot_gt_ls47$Pvalue<=0.1),]
ggplot(for_plot_gt_ls47,aes(x=type,y = motif,size=Odds, color=-log10(Pvalue)))+ geom_point()+
  scale_color_gradient2(midpoint=1, low="white", mid="blue",
                        high="red", space ="Lab" )+
  theme_bw()
for_plot_gt_ls47_V <-for_plot_gt_ls47
############A and V############################
for_plot_gt_ls47_V$celltype <- "V"
for_plot_gt_ls47_A$celltype <- "A"
for_plot_gt_ls47_AV <-rbind(for_plot_gt_ls47_A,for_plot_gt_ls47_V)
#for_plot_gt_ls47_AV<-for_plot_gt_ls47_AV[which(for_plot_gt_ls47_AV$Pvalue<0.05),]
rows_chose<-unique(for_plot_gt_ls47_AV[which(for_plot_gt_ls47_AV$Pvalue<0.05),"motif"])
for_plot_gt_ls47_AV_chose<-for_plot_gt_ls47_AV[for_plot_gt_ls47_AV$motif %in% rows_chose,]
ggplot(for_plot_gt_ls47_AV,aes(x=type,y = motif,size=Odds, color=-log10(Pvalue)))+ geom_point()+
  scale_color_gradient2(midpoint=1, low="white", mid="blue",
                        high="red", space ="Lab" )+
  theme_bw()+facet_wrap(~celltype)
ggplot(for_plot_gt_ls47_AV_chose,aes(x=type,y = motif,size=Odds, color=-log10(Pvalue)))+ geom_point()+
  scale_color_gradient2(midpoint=1, low="white", mid="blue",
                        high="red", space ="Lab" )+
  theme_bw()+facet_wrap(~celltype)
head(for_plot_gt_ls47_AV_chose)

for_plot_gt_ls48_AV <- merge(for_plot_gt_ls47_AV_chose, data_motif_V, by.x="motif",by.y="motif_name")
#for_plot_gt_ls48_AV <- merge(for_plot_gt_ls48_AV, data_motif_A, by.x="motif",by.y="motif_name",all.x=T)
#for_plot_gt_ls48_AV <-for_plot_gt_ls48_AV[!is.na(for_plot_gt_ls48_AV$gene.x)|!is.na(for_plot_gt_ls48_AV$gene.y),]
asces_for_plot_gt_ls48_AV_chose_A <-acast(for_plot_gt_ls48_AV[for_plot_gt_ls48_AV$celltype=="A",c("gene","type","Odds")],formula =gene~type, fill = 0 )
asces_for_plot_gt_ls48_AV_chose_V <-acast(for_plot_gt_ls48_AV[for_plot_gt_ls48_AV$celltype=="V",c("gene","type","Odds")],formula =gene~type, fill = 0 )
asces_for_plot_gt_ls48_AV_chose<- cbind(asces_for_plot_gt_ls48_AV_chose_A,asces_for_plot_gt_ls48_AV_chose_V)
rownames_sorted<-pheatmap(as.matrix(asces_for_plot_gt_ls48_AV_chose), cluster_rows=T,cluster_cols=F,color = colorRampPalette(c( "white","yellow","orange","red"))(100),breaks = seq(0,20,0.2), 
         show_colnames = T,fontsize = 6)
#for_plot_gt_ls48_AV <-as.data.frame(for_plot_gt_ls48_AV[,c(2:6)])
#rownames(for_plot_gt_ls48_AV) <- for_plot_gt_ls48_AV[,"gene"]
row_sorted<-rownames(asces_for_plot_gt_ls48_AV_chose[rownames_sorted$tree_row[["order"]],])
row_sorted <-c("ESRRA,ESRRG","GATA4,GATA6","ZFX","PRDM5","SREBF1,SREBF2","TCF3","MEF2A,MEF2C,MEF2D","ZBTB6",
               "RXRA","NFYA,NFYB,NFYC,PBX3","MGA,TBX20,TBX5","PLAG1","RARA,RARB,RXRA,RXRB,THRA","STAT3,STAT5B")
#for_plot_gt_ls48_AV <-for_plot_gt_ls48_AV[for_plot_gt_ls48_AV$Pvalue<=0.05,]
for_plot_gt_ls48_AV <- for_plot_gt_ls48_AV[for_plot_gt_ls48_AV$Odds>1,]
for_plot_gt_ls48_AV$sig<-0
for_plot_gt_ls48_AV[for_plot_gt_ls48_AV$Pvalue<=0.05,]$sig<-1
for_plot_gt_ls48_AV[for_plot_gt_ls48_AV$type=="GoF_enhancers_lost_motif" | for_plot_gt_ls48_AV$type=="LoF_enhancers_Gain_motif",]$Odds <- -1*for_plot_gt_ls48_AV[for_plot_gt_ls48_AV$type=="GoF_enhancers_lost_motif" | for_plot_gt_ls48_AV$type=="loF_enhancers_Gain_motif",]$Odds
#for_plot_gt_ls48_AV[for_plot_gt_ls48_AV$Pvalue<=0.001,]$Pvalue<-0.001
pdf("AV_motif_odds_ratio_pvalue_odd_with_negative.pdf")
ggplot(for_plot_gt_ls48_AV,aes(x=type,y = factor(gene,levels =rev(row_sorted)),color=Odds, size=-1*log(Pvalue)))+ geom_point()+
  #scale_color_viridis_c(option = "geyser")+
  scale_color_gradientn(colours = (colorspace::divergingx_hcl(7,palette = "Fall")))+
  #scale_color_gradient2(midpoint=2.995732, low="blue", mid="white",
  #                      high="darkred", space ="Lab" )+
  facet_grid(~celltype)+scale_size(range = c(0, 10))+
  theme_bw()
ggplot(for_plot_gt_ls48_AV,aes(x=type,y = factor(gene,levels =rev(row_sorted)),color=factor(sig), size=-1*log(Pvalue)))+ geom_point()+
 # scale_color_gradient2(midpoint=2.995732, low="blue", mid="white",
#                        high="darkred", space ="Lab" )+
  facet_grid(~celltype)+scale_color_manual(breaks = c("1", "0"), values=c("red","lightgray"))+
  scale_size(range = c(0, 10))+
  theme_bw()
ggplot(for_plot_gt_ls48_AV,aes(x=type,y = factor(gene,levels = rev(row_sorted)),color=-log10(Pvalue), size=Odds))+ geom_point()+
  #scale_color_gradient2(midpoint=1.3, low="white", mid="blue",
  #                      high="red", space ="Lab" )+
  scale_color_viridis_c(option="C")+
  facet_grid(~celltype)+scale_size(range = c(0, 10))+
  theme_bw()
dev.off()
asces_for_plot_gt_ls48_AV_chose<-asces_for_plot_gt_ls48_AV_chose[row_sorted,]
pdf("AV_motif_odds_ratio_heatmap.pdf")

pheatmap(as.matrix(asces_for_plot_gt_ls48_AV_chose), cluster_rows=F,cluster_cols=F,color = colorRampPalette(c( "white","yellow","orange","red"))(100),breaks = seq(0,20,0.2), 
         show_colnames = T,fontsize = 6)
dev.off()
for_plot_gt_ls47_AV[for_plot_gt_ls47_AV$type=="GoF_enhancers_lost_motif" | for_plot_gt_ls47_AV$type=="LoF_enhancers_Gain_motif",]$Odds <- 
  -1*for_plot_gt_ls47_AV[for_plot_gt_ls47_AV$type=="GoF_enhancers_lost_motif" | for_plot_gt_ls47_AV$type=="LoF_enhancers_Gain_motif",]$Odds
for_plot_gt_ls49_AV <- merge(for_plot_gt_ls47_AV, data_motif_V, by.x="motif",by.y="motif_name",all.x=T)
colnames(for_plot_gt_ls49_AV)<-c("motif","Odds","Pvalue","type","celltype","gene_V")
for_plot_gt_ls49_AV <- merge(for_plot_gt_ls49_AV, data_motif_A, by.x="motif",by.y="motif_name",all.x=T)
colnames(for_plot_gt_ls49_AV)<-c("motif","Odds","Pvalue","type","celltype","gene_A")
write.table(for_plot_gt_ls49_AV,"AV_motif_odds_ratio_all.txt",sep="\t",quote = F)
for_plot_gt_ls47_AV_chose<-for_plot_gt_ls47_AV[for_plot_gt_ls47_AV$motif %in% rows_chose,]
head(for_plot_gt_ls47_AV_chose)
for_plot_gt_ls49_AV_chose <- merge(for_plot_gt_ls47_AV_chose, data_motif_V, by.x="motif",by.y="motif_name",all.x=T)
colnames(for_plot_gt_ls49_AV_chose)<-c("motif","Odds","Pvalue","type","celltype","gene_V")
for_plot_gt_ls49_AV_chose <- merge(for_plot_gt_ls49_AV_chose, data_motif_A, by.x="motif",by.y="motif_name",all.x=T)
colnames(for_plot_gt_ls49_AV_chose)<-c("motif","Odds","Pvalue","type","celltype","gene_V","gene_A")
write.table(for_plot_gt_ls49_AV_chose,"AV_motif_odds_ratio_pvalue_0.05.txt",sep="\t",quote = F)





#write.table(for_plot_gt_ls48_AV,"AV_motif_odds_ratio_all.txt",sep="\t",quote = F)
################################################################################
plot_gt_ls7 <-plot_gt_ls6
dim(for_heat_G_M_V)
dim(for_heat_L_M_V)
dim(for_heat_no_M_V)
dim(for_heat_all_M_V)
plot_gt_ls7[,c(1:2)]<- plot_gt_ls7[,c(1:2)]/15
plot_gt_ls7[,c(3:4)]<- plot_gt_ls7[,c(3:4)]/49
plot_gt_ls7[,c(5:6)]<- plot_gt_ls7[,c(5:6)]/1336
pheatmap(as.matrix(plot_gt_ls6), cluster_rows=T,cluster_cols=F, 
         color = colorRampPalette(c( "white","gold"))(100), show_colnames = T,fontsize = 8)
pheatmap(as.matrix(plot_gt_ls7), cluster_rows=T,cluster_cols=F, 
         color = colorRampPalette(c( "white","gold"))(100), show_colnames = T,fontsize = 8)


plot_gt_ls8 <- as.data.frame(plot_gt_ls7)
plot_gt_ls8$GoF_diff <- plot_gt_ls8$gain_motif.x-plot_gt_ls8$lost_motif.x
plot_gt_ls8$LoF_diff <- plot_gt_ls8$lost_motif.y-plot_gt_ls8$gain_motif.y
plot_gt_ls8$Nochange <- abs(plot_gt_ls8$gain_motif-plot_gt_ls8$lost_motif)
###############3
plot_gt_ls8$GoF_diff <- plot_gt_ls8$gain_motif.x+plot_gt_ls8$lost_motif.x
plot_gt_ls8$LoF_diff <- plot_gt_ls8$lost_motif.y+plot_gt_ls8$gain_motif.y
plot_gt_ls8$Nochange <- abs(plot_gt_ls8$gain_motif+plot_gt_ls8$lost_motif)
############################
plot_gt_ls8$G_fold_to_backgroud <- log2(plot_gt_ls8$GoF_diff/(plot_gt_ls8$Nochange+0.0000001)+1)
plot_gt_ls8$L_fold_to_backgroud <- log2(plot_gt_ls8$LoF_diff/(plot_gt_ls8$Nochange+0.0000001)+1)
plot_gt_ls9 <- plot_gt_ls8[,10:11]
plot_gt_ls9_V <-plot_gt_ls9
plot_gt_ls9[is.na(plot_gt_ls9)] <- 0
pheatmap(as.matrix(plot_gt_ls9), cluster_rows=T,cluster_cols=F,color = colorRampPalette(c( "blue","white","gold"))(100),breaks = seq(-10,10,0.2), show_colnames = T,fontsize = 8)

plot_gt_ls10 <- plot_gt_ls9[plot_gt_ls9$G_fold_to_backgroud !=0 | plot_gt_ls9$L_fold_to_backgroud !=0 , ]
pheatmap(as.matrix(plot_gt_ls10), cluster_rows=T,cluster_cols=F,color = colorRampPalette(c( "blue","white","gold"))(100),breaks = seq(-10,10,0.2), 
         show_colnames = T,fontsize = 6)

data_motif_V <- read.table(file="data_motif_V.txt",header=T,row.names = 1)
plot_gt_ls11 <- merge(plot_gt_ls10, data_motif_V, by.x="row.names",by.y="motif_name") 
plot_gt_ls12 <-as.matrix(plot_gt_ls11[,c(2,3)])
rownames(plot_gt_ls12) <- plot_gt_ls11[,"gene"]
###############################################
data_motif_V_cluster <- read.table(file="seedmotif_name_genenames_mouse_with_cluster.txt",header=F,row.names = 1)

plot_gt_ls_full <- merge(plot_gt_ls11, data_motif_V_cluster, by.x="Row.names",by.y="row.names")
plot_gt_ls_full$full_name <- paste(plot_gt_ls_full$Row.names,"~",plot_gt_ls_full$V3,"~",plot_gt_ls_full$gene)
plot_gt_ls12_full <-as.matrix(plot_gt_ls_full[,c(2,3)])
rownames(plot_gt_ls12_full) <- plot_gt_ls_full[,"full_name"]
pdf("motif_function_score_for_V.pdf",width =5,height = 6)
#pheatmap(as.matrix(plot_gt_ls12), cluster_rows=T,cluster_cols=F,color = colorRampPalette(c( "blue","white","gold"))(100),breaks = seq(-10,10,0.2), 
#         show_colnames = T,fontsize = 6)
row_names<-pheatmap(as.matrix(plot_gt_ls12_full), cluster_rows=T,cluster_cols=F,color = colorRampPalette(c( "blue","white","gold"))(100),breaks = seq(-5,5,0.1), 
         show_colnames = T,fontsize = 6)
grob_classes <- purrr::map(row_names$gtable$grobs, class)
idx_grob <- which(purrr::map_lgl(grob_classes, function(cl) 'gTree' %in% cl))[1]
grob_names <- names(row_names$gtable$grobs[[idx_grob]]$children)
idx_rect <- grob_names[grep('rect', grob_names)][1]

## Remove borders around cells
row_names$gtable$grobs[[idx_grob]]$children[[idx_rect]]$gp$col <- "black"
row_names$gtable$grobs[[idx_grob]]$children[[idx_rect]]$gp$lwd <- 0.1

## Plot result
graphics::plot.new()
print(row_names)
dev.off()

###################################################
pdf("motif_function_score_for_V.pdf",width =5,height = 6)
pheatmap(as.matrix(plot_gt_ls12), cluster_rows=T,cluster_cols=F,color = colorRampPalette(c( "blue","white","gold"))(100),breaks = seq(-10,10,0.2), 
         show_colnames = T,fontsize = 6)
dev.off()

plot_gt_ls15 <-plot_gt_ls6
plot_gt_ls15[,c(1)]<- (plot_gt_ls15[,c(1)])/(plot_gt_ls15[,c(5)]+plot_gt_ls15[,c(1)])
plot_gt_ls15[,c(2)]<- (plot_gt_ls15[,c(2)])/(plot_gt_ls15[,c(6)]+plot_gt_ls15[,c(2)])
plot_gt_ls15[,c(3)]<- (plot_gt_ls15[,c(3)])/(plot_gt_ls15[,c(5)]+plot_gt_ls15[,c(3)])
plot_gt_ls15[,c(4)]<- (plot_gt_ls15[,c(4)])/(plot_gt_ls15[,c(6)]+plot_gt_ls15[,c(4)])
plot_gt_ls15[is.na(plot_gt_ls15)] <- 0
plot_gt_ls16 <-plot_gt_ls15[plot_gt_ls15$gain_motif.x >=0.15 | plot_gt_ls15$lost_motif.x >=0.15 |plot_gt_ls15$gain_motif.y >=0.15 | plot_gt_ls15$lost_motif.y >=0.15,
                            c(1:4)]
plot_gt_ls17 <- merge(plot_gt_ls16, data_motif_V, by.x="row.names",by.y="motif_name") 
plot_gt_ls18 <-as.matrix(plot_gt_ls17[,c(2,3,4,5)])
rownames(plot_gt_ls18) <- plot_gt_ls17[,"gene"]
pdf("motif_persentage_left_GoF_right_LoF_for_V.pdf",width =5,height = 6)
pheatmap(as.matrix(plot_gt_ls18), cluster_rows=T,cluster_cols=F,color = colorRampPalette(c( "blue","white","gold"))(100),breaks = seq(-1,1,0.02), 
         show_colnames = T,fontsize = 6)
dev.off()

#############################A V plot###############################
plot_gt_ls9_merge <- merge(plot_gt_ls9_A, plot_gt_ls9_V, by="row.names")
plot_gt_ls9_merge[is.na(plot_gt_ls9_merge)] <- 0
plot_gt_ls10_merge <- plot_gt_ls9_merge[plot_gt_ls9_merge$G_fold_to_backgroud.x !=0 | plot_gt_ls9_merge$L_fold_to_backgroud.x !=0 |
                                    plot_gt_ls9_merge$G_fold_to_backgroud.y !=0 | plot_gt_ls9_merge$L_fold_to_backgroud.y !=0 , ]
rownames(plot_gt_ls10_merge) <- plot_gt_ls10_merge$Row.names
plot_gt_ls11_merge <- merge(plot_gt_ls10_merge, data_motif_V, by.x="Row.names",by.y="motif_name") 
plot_gt_ls12_merge <-as.matrix(plot_gt_ls11_merge[,c(2,3,4,5)])
rownames(plot_gt_ls12_merge) <- plot_gt_ls11_merge[,"gene"]
pdf("motif_fuction_score_left_A_right_V_for_A_andV.pdf",width =5,height = 6)
pheatmap(as.matrix(plot_gt_ls12_merge), cluster_rows=T,cluster_cols=F,color = colorRampPalette(c( "blue","white","gold"))(100),breaks = seq(-5,5,0.1), 
         show_colnames = T,fontsize = 6)
dev.off()
write.table(plot_gt_ls12_merge, "mutagenesis_A_V_function_score.txt",sep="\t",quote=F)
#########################################################################################3
pheatmap(for_heat_L_M[for_heat_L$Specics=="Homo" | for_heat_L$Specics=="Mus" , ], cluster_rows=T,cluster_cols=T, 
         breaks=seq(-10,10,0.2),color = colorRampPalette(c("navy", "white","gold"))(100), show_colnames = F)
pheatmap(for_heat_G_M[for_heat_G$Specics=="Homo", ], cluster_rows=T,cluster_cols=T, color = colorRampPalette(c("navy", "white","gold"))(100))
pheatmap(for_heat_G_M[for_heat_G$Specics=="Homo" | for_heat_G$Specics=="Mus" , ], cluster_rows=T,cluster_cols=T, color = colorRampPalette(c("navy", "white","gold"))(100))
for_heat_G_M$factor_name <- row.names(for_heat_G_M)
for_heat_G$change ="Gain"
for_heat_L$change ="Loss"


ggplot(for_heat_G,aes(y=name))+ geom_bin(fill = after_stat(count),binwidth = c(0.1, 0.1))
for_heat_G$gt_4 <- rowSums(for_heat_G[,1:149] >= 4)
for_heat_G %>% dplyr::group_by(name) %>% dplyr::summarise(sum(gt_4)) %>% as.data.frame()

for_heat_G$ls_m4 <- rowSums(for_heat_G[,1:149] <= -4)
for_heat_G %>% dplyr::group_by(name) %>% dplyr::summarise(sum(ls_m4)) %>% as.data.frame()


###################################################################################
#install.packages("neuralnet")
library(neuralnet)
library("ISLR")
#df <- data.frame(t(melt_for_heat_all_M_V))

dat_use_df <- data.frame(t(melt_for_heat_all_M_A))
dat_use_df[1:5,1:5]

fold_info <- unique(filter_norm_data[,c("region_name","fold_A")]) %>% as.data.frame()
rownames(fold_info) <- fold_info$region
dat_use_df$fold <- fold_info[rownames(dat_use_df),]$fold_A
#####################################################################
mean_data <- apply(dat_use_df[1:110], 2, mean)
sd_data <- apply(dat_use_df[1:110], 2, sd)
data_scaled <- as.data.frame(scale(dat_use_df[,1:110],center = mean_data, scale = sd_data))
head(data_scaled, n=20)
#index = sample(1:nrow(data),round(0.70*nrow(data)))
#train_data <- as.data.frame(data_scaled[index,])
#test_data <- as.data.frame(data_scaled[-index,])
#####################################################################
dat_split <- initial_split(data_scaled)
dat_split <- initial_split(dat_use_df)
dat_train <- training(dat_split)
dat_val <- testing(dat_split)
n = colnames(dat_train)
f = as.formula(paste("fold ~", paste(n[!n %in% "fold"], collapse = " + ")))
#nn=neuralnet(f,data=dat_train, hidden=3,act.fct = "logistic",
#            linear.output = FALSE)
nn=neuralnet(f,data=dat_train, hidden=4,act.fct="tanh",
             linear.output = FALSE)
predict_net_test <- compute(nn,dat_val[,1:(length(dat_val[1,])-1)])
predict_net_test <- compute(nn,dat_val[,1:155])
MSE.net <- sum((dat_val$fold - predict_net_test$net.result)^2)/nrow(dat_val)
Lm_Mod <- lm(fold~., data=dat_train)
AA<-summary(Lm_Mod)
for_plot2 <- as.data.frame(AA$coefficients)
for_plot <-  for_plot2[for_plot2$`Pr(>|t|)`<0.05,4]
names(for_plot) <- rownames(for_plot2[for_plot2$`Pr(>|t|)`<0.05,])
for_plot_new<- merge(for_plot,data_motif_A,by.x="row.names",by.y="motif_name")
for_plot_new_2 <- as.matrix(for_plot_new[,2])
rownames(for_plot_new_2)<-for_plot_new$gene

pheatmap(for_plot_new_2,cluster_rows=F,cluster_cols=F, color = colorRampPalette(c("gold", "white","navy"))(100))
predict_lm <- predict(Lm_Mod,dat_val)
MSE.lm <- sum((predict_lm - dat_val$fold)^2)/nrow(dat_val)
Lm_Mod <- lm(fold~., data=dat_val)
summary(Lm_Mod)
#predict_lm <- predict(Lm_Mod,dat_val)
#MSE.lm <- sum((predict_lm - dat_val$fold)^2)/nrow(dat_val)
par(mfrow=c(1,2))
plot(dat_val$fold,predict_net_test$net.result,col='blue',main='Real vs predicted for neural network',pch=18,
     cex=2,xlim=c(-3, 3),ylim=c(-3,3))
abline(0,1,lwd=5)

plot(dat_val$fold,predict_lm,col='blue',main='Real vs predicted for linear regression',pch=18,
     cex=2,xlim=c(-3, 3),ylim=c(-3,3))
abline(0,1,lwd=5)
#################################
mplot_density(dat_val$fold,predict_net_test$net.result )
library(ROCR)
pred <- prediction(predict_net_test$net.result, dat_val$fold)


#################################
#dat_train <- as.matrix(dat_train)
dat_use_df <- data.frame(t(melt_for_heat_all_M_A))
dat_use_df[1:5,1:5]

fold_info <- unique(filter_norm_data[,c("region_name","fold_A","REF_A","ALT_A")]) %>% as.data.frame()
rownames(fold_info) <- fold_info$region_name
library(kmer)
woodmouse.kdist <- dist(fold_info[,2:4])
set.seed(1)
eg4.kmeans <- kmeans(fold_info[ , 2:4], 8)
fold_info$cluster <- eg4.kmeans$cluster
ggplot(fold_info,aes(x=REF_A, y=ALT_A, color=factor(cluster)))+geom_point()

dat_use_df$fold <- fold_info[rownames(dat_use_df),]$cluster
dat_split <- initial_split(dat_use_df)
dat_train <- training(dat_split)
dat_val <- testing(dat_split)


dim(dat_train)
model <- rand_forest(mtry = 100, trees = 2000, min_n = 15, mode = "regression") %>%
  set_engine("ranger", importance = "impurity", num.threads = 3) %>%
  fit(fold ~ ., data = dat_train)
plot(model)
var_imp <- sort(model$fit$variable.importance, decreasing = TRUE)
top_genes <- names(var_imp)
#Namelist <- read.table(file="~/Desktop/work/Bill/Fengxiao/motifs/db.list",header=T,row.names = 1)
var_imp_t <- as.data.frame(var_imp) 


var_imp_t$name <- rownames(var_imp_t)
#var_imp_t$Specics <- Namelist[rownames(var_imp_t),2]
head(var_imp_t)
data_for_plot <- as.matrix(var_imp_t[1:50,1])
rownames(data_for_plot) <- var_imp_t$name[1:50]
data_for_plot_new<- merge(data_for_plot,data_motif_A,by.x="row.names",by.y="motif_name")
data_for_plot_new_22 <- data_for_plot_new[order(-data_for_plot_new$V1),]
data_for_plot_new_2 <- as.matrix(data_for_plot_new_22[,2])
rownames(data_for_plot_new_2)<-data_for_plot_new_22$gene

pheatmap(data_for_plot_new_2,cluster_rows=F,cluster_cols=F, color = colorRampPalette(c("navy", "white","gold"))(100))

pheatmap(data_for_plot,cluster_rows=F,cluster_cols=F, color = colorRampPalette(c("navy", "white","gold"))(100))


prediction <-predict(model, dat_val[,1:length(dat_val[1,]-1)])
dat_val$pre <- prediction$.pred
#confusionMatrix(prediction, dat_val$fold)

ggplot(dat_val,aes(x= pre,y=fold))+geom_point(color="darkblue")+xlim(-3,3)+ylim(-3,3)+
  geom_abline(slope=1)
ggplot(dat_val,aes(x= pre,y=fold))+geom_point(color="darkblue")+xlim(1,8)+ylim(1,8)+
  geom_abline(slope=1)
##########################################################################3

