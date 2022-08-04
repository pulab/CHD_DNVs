setwd("~/Desktop/work/Bill/Fengxiao/muatiaon_MPRA/raw_counts")
library(ggplot2)
library(DESeq2)
library(Rfast)
library(dplyr)
library(reshape2)
library(fitdistrplus)
#library(sensitivity)
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
#data<-read.table(file="format_table.txt",header=T,row.names = 1)
data<-read.table(file="Mutagenesis_raw_counts_DNA1-4_RNA1-4.txt",header=T,row.names = 1)
###################################DESeq2########################
head(data)
dim(data)
filter_for_degs<- data[,c(1:8)]
#filter_for_degs<- data_sub[rowMins(as.matrix(norm_data[,13:18]), value = T) >=20,13:24]
coldata <-as.data.frame(colnames(filter_for_degs))

coldata$type <- coldata$`colnames(filter_for_degs)`
coldata$condition <-gsub("\\d","",coldata$type)
genes <-as.matrix((filter_for_degs))
dds <- DESeqDataSetFromMatrix(countData = genes, colData=coldata, design=~condition)
dds <- DESeq(dds)
res <- results( dds, contrast = c("condition","RNA","DNA") )

res$qvalue_0_05 <- 0
res$qvalue_0_05[res$padj < 0.05 & res$log2FoldChange >0] <- 1
sum(res$qvalue_0_05)
#res2$qvalue_0_05 <- 0
#res2$qvalue_0_05[res2$padj < 0.05 & res2$log2FoldChange >0] <- 1
#sum(res2$qvalue_0_05)
######################to FPM#####################################
data_sub <- as.matrix(data[,1:8])
norm_data <- t(t(data_sub)/colSums(data_sub)*1000000)
#filter_norm_data <- cbind(norm_data, data[,6:7]) %>% dplyr::filter(DNA >= 5 &rowSums(.[,1:4]>0))
test <- dplyr::select(as.data.frame(norm_data), contains("DNA"))
filter_norm_data <- norm_data[rowMins(as.matrix(test), value = T) >=20,]
head(test)
new <- data.frame(matrix(ncol = 0, nrow = length(rownames(filter_norm_data))))
rownames(new) <- rownames(filter_norm_data)

for_plot <- as.data.frame(test) 
ggplot(for_plot, aes(x = DNA3))+ geom_histogram(binwidth = 5,color="skyblue")+theme_bw()+theme_classic()+
  geom_vline(xintercept = 20,color="skyblue",size=0.4,linetype="dashed")+xlim(0,500)

for (i in 1:4){
  j <- i +4
  ratio <- (filter_norm_data[,j]+1)/(filter_norm_data[,i]+1)
  ratio <- log2(ratio)
  new <- cbind(new,ratio)
}

new_colname <- paste0(colnames(filter_norm_data)[grep("RNA", colnames(filter_norm_data))], "_vs_DNA_log")
colnames(new) <- new_colname
head(new)
######################size factor#####################################
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
Negative_C <- t(dplyr::select(as.data.frame(t(new)), !contains("::")))
Positive_C <- t(dplyr::select(as.data.frame(t(new)), contains("::F")))
#both_C <- rownames(Positive_C)
both_C <- rownames(rbind(Negative_C,Positive_C))
filter_norm_data_nc_normalize_exp <-MPRA_SizeFactors(new,both_C,c(1:4))
filter_norm_data_nc_normalize_exp <- as.data.frame(filter_norm_data_nc_normalize_exp)

#filter_norm_data_nc <- Positive_C[,c(1:4)] %>% as.matrix()
#filter_norm_data_full_nc <- new[,c(1:4)] %>% as.matrix()
#filter_norm_data_nc <- filter_norm_data_nc[rowMins(filter_norm_data_nc,value = T)>0,]
#filter_norm_data_nc_ratio_pseudo_exp <- geometricmeanRow(filter_norm_data_nc)
#filter_norm_data_nc_ratio_scaled <- t(t(filter_norm_data_nc)/filter_norm_data_nc_ratio_pseudo_exp)

#filter_norm_data_nc_normalize_factor <- colMedians(filter_norm_data_nc_ratio_scaled)
#filter_norm_data_nc_normalize_exp <- t(t(new[,1:4])/filter_norm_data_nc_normalize_factor)
#filter_norm_data_scaled <- cbind(filter_norm_data[,c(5:7)], filter_norm_data_nc_normalize_exp)
#filter_norm_data_scaled <- cbind(new[,c(1:4)], filter_norm_data_nc_normalize_exp)
#filter_norm_data_scaled$A_ratio <- (filter_norm_data_scaled$A1 + filter_norm_data_scaled$A2)/2/filter_norm_data_scaled$DNA
#filter_norm_data_scaled$V_ratio <- (filter_norm_data_scaled$V1 + filter_norm_data_scaled$V2)/2/filter_norm_data_scaled$DNA
#filter_norm_data_nc_normalize_exp <- as.data.frame(filter_norm_data_nc_normalize_exp)
A_V_cor <- cor(filter_norm_data_nc_normalize_exp)
A_V_cor2 <- cor(new)
corrplot(A_V_cor, method  = "number",tl.cex=0.5,hclust.method="median")
corrplot(A_V_cor2, method  = "number",tl.cex=0.5,hclust.method="median")
ggplot(filter_norm_data_nc_normalize_exp,aes(x=RNA1_vs_DNA_log,y=RNA2_vs_DNA_log))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(0,3)+ylim(0,3)
ggplot(filter_norm_data_nc_normalize_exp,aes(x=RNA1_vs_DNA_log,y=RNA3_vs_DNA_log))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(0,3)+ylim(0,3)
ggplot(filter_norm_data_nc_normalize_exp,aes(x=RNA1_vs_DNA_log,y=RNA4_vs_DNA_log))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(0,3)+ylim(0,3)
ggplot(filter_norm_data_nc_normalize_exp,aes(x=RNA2_vs_DNA_log,y=RNA3_vs_DNA_log))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(0,3)+ylim(0,3)

write.table(new,file="10bp_tiled_log2_ratio_add1.txt",sep="\t")
write.table(filter_norm_data,file="10bp_tiled_log2_FPM_add1.txt",sep="\t")
write.table(filter_norm_data_nc_normalize_exp,file="10bp_tiled_log2_normalized_log2_ratio_add1.txt",sep="\t")

#################################rank figure#####################################
filter_norm_data_nc_normalize_exp <- read.table(file="10bp_tiled_log2_normalized_log2_ratio_add1.txt",header=T,row.names = 1)

mean_log_ratio <- rowMeans(filter_norm_data_nc_normalize_exp)
filter_norm_data_scaled <- cbind(filter_norm_data_nc_normalize_exp,mean_log_ratio)
filter_norm_data_scaled<-as.data.frame(filter_norm_data_scaled)
filter_norm_data_scaled$rank_value <- rank(-as.numeric(as.character(filter_norm_data_scaled$mean_log_ratio)))

select_ALT <- t(filter_norm_data_scaled) %>% as.data.frame() %>% dplyr::select(contains("::") & !contains("_17")) %>% t()
select_ALT <- as.data.frame(select_ALT)
select_ALT$group <- "ALT"
select_REF <- t(filter_norm_data_scaled) %>% as.data.frame() %>% dplyr::select(contains("_17")) %>% t()
select_REF <- as.data.frame(select_REF)
select_REF$group <- "WT"

select_NC <- t(filter_norm_data_scaled) %>% as.data.frame() %>% dplyr::select(!contains("::")) %>% t()
select_NC <- as.data.frame(select_NC)
select_NC$group <- "NegativeControl"
#select_PC <- t(filter_norm_data_scaled) %>% as.data.frame() %>% dplyr::select(contains("::F")) %>% t()
#select_PC <- as.data.frame(select_PC)
#select_PC$group <- "PostiveControl"
filter_norm_data_scaled_grouped <-rbind(select_REF,select_ALT,select_NC)
filter_norm_data_scaled_grouped <- as.data.frame(filter_norm_data_scaled_grouped)
##########################FDR 0.05################################################
fitg<-fitdist(select_NC$mean_log_ratio ,"norm")

FDR_0.95 <- qnorm(0.95, mean = fitg$estimate["mean"], sd = fitg$estimate["sd"])
FDR_0.95
ggplot(filter_norm_data_scaled_grouped[filter_norm_data_scaled_grouped$group=="NegativeControl",], aes(x = mean_log_ratio, color=group))+geom_density()+theme_bw()+theme_classic()+
  geom_vline(xintercept = FDR_0.95,color="skyblue",size=0.4,linetype="dashed")


############################## heatmap #############################

filter_norm_data_scaled_grouped$order <- rank(-as.numeric(as.character(filter_norm_data_scaled_grouped$mean_log_ratio)))
filter_norm_data_scaled_grouped$active <- res[rownames(filter_norm_data_scaled_grouped),]$qvalue_0_05
#ggplot(filter_norm_data_scaled_grouped,aes(x=order,y=mean_log_ratio))+geom_point(size=2,color="darkblue")+
#  geom_hline(yintercept=FDR_0.95,color="red",linetype="dashed")+theme_bw()+theme_classic()+geom_text(aes(0,round(FDR_0.95,3),label = round(FDR_0.95,3), vjust = -1))

ggplot(filter_norm_data_scaled_grouped,aes(x=order,y=mean_log_ratio,color=factor(active)))+geom_point(size=2)+
  scale_color_manual(breaks = c("1", "0"), values=c("red", "darkblue"))+
  geom_hline(yintercept=FDR_0.95,color="red",linetype="dashed")+theme_bw()+theme_classic()+geom_text(aes(0,round(FDR_0.95,3),label = round(FDR_0.95,3), vjust = -1))

ggplot(filter_norm_data_scaled_grouped[filter_norm_data_scaled_grouped$group == "WT",],aes(x=order,y=mean_log_ratio,color=factor(active)))+geom_point(size=2)+
  scale_color_manual(breaks = c("1", "0"), values=c("red", "darkblue"))+
  geom_hline(yintercept=FDR_0.95,color="red",linetype="dashed")+theme_bw()+theme_classic()+geom_text(aes(0,round(FDR_0.95,3),label = round(FDR_0.95,3), vjust = -1))


ggplot(filter_norm_data_scaled_grouped, aes(x = order, y = group))+stat_bin2d(aes(fill = after_stat(count)), binwidth = c(3,1))
write.table(filter_norm_data_scaled_grouped,file="10bp_tiled_log2_normalized_log2_ratio_with_active_add1.txt",sep="\t",quote = F)
#########################################################

#####################################residual ######################3
filter_norm_data <- data
filter_norm_data$order_A_residuals <- rank(-as.numeric(as.character(filter_norm_data$A_residuals)))
filter_norm_data$order_V_residuals <- rank(-as.numeric(as.character(filter_norm_data$V_residuals)))
filter_norm_data$AVsV <- log2(filter_norm_data$A_ratio+1) /log2(filter_norm_data$V_ratio+1)
head(filter_norm_data)
ggplot(filter_norm_data,aes(x=order_A_residuals ,y=A_residuals))+geom_point(size=2,color="darkblue")+theme_bw()+theme_classic()
ggplot(filter_norm_data,aes(x=order_V_residuals ,y=V_residuals))+geom_point(size=2,color="darkblue")+theme_bw()+theme_classic()
filter_norm_data$Group6 <- paste0(filter_norm_data$aMTF,filter_norm_data$vMTF,filter_norm_data$aP300,filter_norm_data$vP300,filter_norm_data$aGene,filter_norm_data$vGene)
ggplot(filter_norm_data, aes(x = order_A_residuals, y = as.factor(Group6)))+stat_bin2d(aes(fill = after_stat(count)), binwidth = c(3,1))
filter_norm_data$order_log_V_ratio <- rank(-as.numeric(as.character(filter_norm_data$log_V_ratio)))
ggplot(filter_norm_data,aes(x=order_V_residuals ,y=V_residuals))+geom_point(size=2,color="darkblue")+theme_bw()+theme_classic()
ggplot(filter_norm_data, aes(x = order_V_residuals, y = as.factor(Group6)))+stat_bin2d(aes(fill = after_stat(count)), binwidth = c(3,1))

####################################################################3


################################diff enhancers####################3
table <- read.table("10bp_tiled_log2_normalized_log2_ratio_with_active_add1.txt", header = T, row.names = 1)
table <- filter_norm_data_scaled_grouped
head(table)
#table <- filter_norm_data_scaled_grouped
table_RNA <- table %>% dplyr::select(contains("DNA_log"))
table_active <- table %>% dplyr::select(contains("active"))
select_REF <- t(table_RNA) %>% as.data.frame() %>% dplyr::select(contains("_17")) %>% t()
select_REF_active <- t(table_active) %>% as.data.frame() %>% dplyr::select(contains("_17")) %>% t()
sum(as.data.frame(select_REF_active)$active)
dim(select_REF_active)
i = 0
WT_vs_KO <- NULL
for (i in 0:16){
  rownames_KO <- gsub("_17", paste0("_", i), rownames(select_REF))
  for (j in 1:length(rownames_KO)){
    if (rownames_KO[j] %in% rownames(table_RNA)){
      KO_value <- as.numeric(as.character(table_RNA[rownames_KO[j],]))
      WT_value <- as.numeric(as.character(table_RNA[rownames(select_REF)[j],]))
      mean_KO_j <- mean(KO_value)
      mean_WT_j <- mean(WT_value)
      pvalue_KO_j <- t.test(KO_value, WT_value)
      WT_vs_KO <- rbind(WT_vs_KO, c(rownames_KO[j], rownames(select_REF)[j], 
                                    KO_value, WT_value,mean_KO_j,mean_WT_j , pvalue_KO_j$p.value))      
    }
  }
}
colnames(WT_vs_KO) <- c("KO","WT","KO_1","KO_2","KO_3","KO_4","WT_1","WT_2","WT_3","WT_4","mean_KO","mean_WT","pvalue")
WT_vs_KO <- as.matrix(WT_vs_KO)
pvalue <- as.numeric(as.character(WT_vs_KO[,"pvalue"]))
qvalue <- pvalue * length(pvalue) /rank(pvalue)
WT_vs_KO <- as.data.frame(WT_vs_KO)
WT_vs_KO$qvalue <- as.numeric(as.character(qvalue))

######################################################################3
inter_join <-WT_vs_KO

#inter_join$active_KO <- table[inter_join$KO,]$active
#inter_join$active_WT <- table[inter_join$WT,]$active
inter_join$active_WT <- res[inter_join$WT,]$qvalue_0_05
inter_join$active_KO <- res[inter_join$KO,]$qvalue_0_05
inter_join$fold <- as.numeric(as.character(inter_join$mean_KO))-as.numeric(as.character(inter_join$mean_WT))
inter_join$qvalue_0_05 <- 0
inter_join$qvalue_0_05[inter_join$qvalue < 0.05] <- 1
inter_join$fold_1_5 <- 0
inter_join$fold_1_5 [inter_join$fold <= -0.58] <- 1
inter_join$fold_1_5 [inter_join$fold >= 0.58] <- 2
inter_join$active <- 0
#inter_join$active[as.numeric(as.character(inter_join$mean_WT))>FDR_0.95] <-1
#inter_join$active[as.numeric(as.character(inter_join$mean_KO))>FDR_0.95] <-1
#inter_join$active <- 0
inter_join$active[as.numeric(as.character(inter_join$active_WT))>0] <-1
inter_join$active[as.numeric(as.character(inter_join$active_KO))>0] <-1
inter_join$final_sig <-0
inter_join$final_sig[inter_join$qvalue_0_05 >0 & inter_join$fold <= -0.58 & inter_join$active >0] <-1
inter_join$final_sig[inter_join$qvalue_0_05 >0 & inter_join$fold >= 0.58 & inter_join$active >0] <-2
inter_join$final_sig_2 <-0
inter_join$final_sig_2[inter_join$qvalue_0_05 >0 & inter_join$fold >= 0.58] <-2
inter_join$final_sig_2[inter_join$qvalue_0_05 >0 & inter_join$fold <= -0.58] <-1
write.table(inter_join,file="10bp_tiled_log2_normalized_log2_ratio_fold_pvalue_add_1.txt",sep="\t")

####################################################################
inter_join<-read.table(file="10bp_tiled_log2_normalized_log2_ratio_fold_pvalue_add_1.txt",header = T, row.names = 1)

inter_join$mean_KO <- as.numeric(as.character(inter_join$mean_KO))
inter_join$mean_WT <- as.numeric(as.character(inter_join$mean_WT))

ggplot(inter_join,aes(log2(fold), -log10(qvalue),color=factor(fold_1_5,levels=c('0','1','2')),alpha =factor(active,levels=c("0","1"))))+
  geom_point()+theme_bw()+theme_classic()+geom_vline(xintercept = 0.58,color="black",size=0.2,linetype="dashed")+
  geom_vline(xintercept = -0.58,color="black",size=0.2,linetype="dashed")+geom_hline(yintercept =-log10(0.05),color="black",size=0.2,linetype="dashed")+scale_color_manual(breaks = c("0", "1", "2"), values=c("gray", "red", "blue"))

ggplot(inter_join,aes(mean_WT, mean_KO,color=factor(final_sig_2,levels=c('0','1','2')),alpha =factor(active,levels=c("0","1"))))+
  geom_point()+theme_bw()+theme_classic()+geom_abline(slope = 1,color="black",size=0.2,linetype="dashed")+
  geom_abline(slope = 1.5,color="black",size=0.2,linetype="dashed")+
  geom_abline(slope=0.667,color="black",size=0.2,linetype="dashed")+
  scale_color_manual(breaks = c("0", "1", "2"), values=c("gray", "red", "blue"))+
  scale_alpha_manual(breaks = c("0", "1"), values=c(0.2, 1))+
  geom_vline(xintercept = FDR_0.95,color="skyblue",size=0.4,linetype="dashed")+
  geom_hline(yintercept = FDR_0.95,color="skyblue",size=0.4,linetype="dashed")+
  xlim(0,4)+ylim(0,4)
#########################DEE plot ########################
ggplot(inter_join,aes(mean_WT, mean_KO,color=factor(final_sig,levels=c('0','1','2')),alpha =factor(qvalue_0_05,levels=c("0","1"))))+
  geom_point()+theme_bw()+theme_classic()+geom_abline(slope = 1,color="black",size=0.2,linetype="dashed")+
  geom_abline(slope = 1, intercept =-0.58, color="black",size=0.2,linetype="dashed")+
  geom_abline(slope= 1,intercept =0.58, color="black",size=0.2,linetype="dashed")+
  scale_color_manual(breaks = c("0", "1", "2"), values=c("gray", "red", "blue"))+
  scale_alpha_manual(breaks = c("0", "1"), values=c(0.2, 1))+
  geom_vline(xintercept = FDR_0.95,color="skyblue",size=0.4,linetype="dashed")+
  geom_hline(yintercept = FDR_0.95,color="skyblue",size=0.4,linetype="dashed")+
  xlim(-5,5)+ylim(-5,5)

ggplot(inter_join,aes(mean_WT, mean_KO,color=factor(final_sig,levels=c('0','1','2')), alpha =factor(qvalue_0_05,levels=c("0","1"))))+
  geom_point()+theme_bw()+theme_classic()+geom_abline(slope = 1,color="black",size=0.2,linetype="dashed")+
  geom_abline(slope = 1.5,color="black",size=0.2,linetype="dashed")+
  geom_abline(slope=0.667,color="black",size=0.2,linetype="dashed")+
  scale_color_manual(breaks = c("0", "1", "2"), values=c("gray", "red", "blue"))+
  scale_alpha_manual(breaks = c("0", "1"), values=c(0.2, 1))
sub <-inter_join[inter_join$final_sig==1,1:8]      
pheatmap(sub,cluster_rows=T,cluster_cols=F,scale="row",color = colorRampPalette(c("navy", "white","gold"))(50))
sum(inter_join$final_sig==1)
sum(inter_join$final_sig==2)
###########################################################################################

###########################################################################################
groups_split<-str_split_fixed(inter_join$KO, pattern="::|_",3)
inter_join_location <- inter_join[,c(1:2,11:20)]
inter_join_location$group_name <-groups_split[,1]
inter_join_location$group_number <-groups_split[,3]
inter_join_location$group_F <-groups_split[,2]
head(inter_join_location)
#inter_join_location$active_WT <- 0
#inter_join_location$active_WT[as.numeric(as.character(inter_join_location$mean_WT))>FDR_0.95] <-1
###################################################################################3
#################################################################################3
##################################################################################3
#groups_split<-str_split_fixed(rownames(filter_norm_data_scaled_grouped), pattern="::|_",3)
#filter_norm_data_scaled_grouped_2 <- filter_norm_data_scaled_grouped[,c(5:9)]
#filter_norm_data_scaled_grouped_2 <- filter_norm_data_scaled_grouped[,c(5:6)]
#filter_norm_data_scaled_grouped_2$group_name <-groups_split[,1]
#filter_norm_data_scaled_grouped_2$group_number <-groups_split[,3]
#filter_norm_data_scaled_grouped_2$group_F <-groups_split[,2]
#filter_norm_data_scaled_grouped_2$active_WT <- 0
#filter_norm_data_scaled_grouped_2$active_WT[as.numeric(as.character(filter_norm_data_scaled_grouped_2$active))>0 &filter_norm_data_scaled_grouped_2$group=="WT"] <-1

#head(filter_norm_data_scaled_grouped_2)

data2<-read.table("library3_171_full_names.txt",header= F, row.names = 1)
#filter_norm_data_scaled_grouped_2$full_names <- data2[as.character(filter_norm_data_scaled_grouped_2$group_name),]
inter_join_location$full_names <- data2[as.character(inter_join_location$group_name),]
#head(filter_norm_data_scaled_grouped_2)
head(inter_join_location)
####################################################

length(unique(inter_join_location[inter_join_location$active_WT >0,"full_names"]))
length(unique(inter_join_location[,"full_names"]))
B<-unique(inter_join_location[,"full_names"])
B

A<-unique(inter_join_location[inter_join_location$active_WT >0,"full_names"])
C<- intersect(A,B)


active_paired <- unique(inter_join_location[inter_join_location$active_WT >0,c("full_names","group_F")])
active_F1 <- active_paired[active_paired$group_F == "F1",]
active_F2 <- active_paired[active_paired$group_F == "F2",]
active_F3 <- active_paired[active_paired$group_F == "F3",]
ggplot(active_paired, aes(x=group_F,y=full_names))+geom_bin_2d()+
  theme_bw()+theme_classic()
active_paired$value=1
active_paire_matrix <-acast(active_paired, full_names~group_F,value=value, fill =0)
pheatmap(active_paire_matrix,cluster_cols = F,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

dim(active_F1)
dim(active_F2)
dim(active_F3)
active_F1_F2 <- intersect(active_F1$full_names,active_F2$full_names)
active_F1_F3 <- intersect(active_F1$full_names,active_F3$full_names)
active_F2_F3 <- intersect(active_F2$full_names,active_F3$full_names)
length(active_F1_F2)
length(active_F2_F3)
length(active_F1_F3)
active_F1_F2_F3 <- intersect(active_F1_F2,active_F3$full_names)
length(active_F1_F2_F3)
##################################################3
inter_join_location$final_sig <- inter_join[rownames(inter_join_location),]$final_sig
length(unique(inter_join_location[inter_join_location$final_sig==1,]$full_names))
length(unique(inter_join_location[inter_join_location$final_sig==2,]$full_names))
inter_join_location_1<-inter_join_location[inter_join_location$group_name=="chr20:49039671-49039842",]
ggplot(inter_join_location_1,aes(as.numeric (group_number), log2(fold),color=factor(fold_1_5,levels=c('0','1','2')),alpha =factor(qvalue_0_05,levels=c("0","1"))))+
  geom_point()+theme_bw()+theme_classic()+
  geom_hline(yintercept = 0,color="black",size=0.4,linetype="dashed")+
  scale_color_manual(breaks = c("0", "1", "2"), values=c("gray", "red", "blue"))+
  scale_alpha_manual(breaks = c("0", "1"), values=c(0.2, 1))+
  geom_hline(yintercept = log2(1.5),color="skyblue",size=0.4,linetype="dashed")+
  geom_hline(yintercept = log2(0.677),color="skyblue",size=0.4,linetype="dashed")+
  facet_grid(rows = vars(group_name),cols = vars(group_F),scales="free")
####################################################################################
head(inter_join_location)
KO_median<- as.data.frame(base::tapply( inter_join_location$mean_KO,inter_join_location$WT, FUN =median))
colnames(KO_median)<-"KO_median"
head(KO_median)
inter_join_location$KO_median <- KO_median[inter_join_location$WT,]
###################################################################################

inter_join_location_1<-inter_join_location[inter_join_location$full_names=="chr20:49530503-49530903",]
inter_join_location_1<-inter_join_location[inter_join_location$full_names=="chr10:35047312-35047712",]
inter_join_location_1<-inter_join_location[inter_join_location$full_names=="chr6:39745032-39745432",]
inter_join_location_1<-inter_join_location[inter_join_location$full_names=="chr3:34053580-34053980",]
full_names_1 <- inter_join_location[inter_join_location$final_sig>0,]

table(unique(full_names_1)$full_names)[as.character(full_names_1$full_names)]

inter_join_location_1<-inter_join_location[inter_join_location$full_names %in% full_names_1,]
inter_join_location_names <- head(inter_join_location[inter_join_location$final_sig>0,]$full_names)
inter_join_location_names <- c("chr7:93870775-93871175","chr14:59536249-59536649","chr4:148231410-148231810","chr1:245748592-245748992")
inter_join_location_names <- c("chr15:96808175-96808575","chr2:47925801-47926201","chr18:22819712-22820112","chr13:39332360-39332760")
inter_join_location_names <- c("chr3:82149507-82149907",
                               "chr20:49157407-49157807",
                               "chr9:137494462-137494862",
                               "chr9:101837589-101837989",
                               "chr4:77586893-77587293")
inter_join_location_1<-inter_join_location[inter_join_location$full_names %in% inter_join_location_names,]
#inter_join_location_1$mean_KO=as.numeric(levels(inter_join_location_1$mean_KO))[inter_join_location_1$mean_KO]

ggplot(inter_join_location_1,aes(as.numeric (group_number), log2(fold),color=factor(fold_1_5,levels=c('0','1','2')),alpha =factor(qvalue_0_05,levels=c("0","1"))))+
  geom_point()+theme_bw()+theme_classic()+
  geom_hline(yintercept = 0,color="black",size=0.4,linetype="dashed")+
  scale_color_manual(breaks = c("0", "1", "2"), values=c("gray", "red", "blue"))+
  scale_alpha_manual(breaks = c("0", "1"), values=c(0.2, 1))+
  geom_hline(yintercept = log2(1.5),color="skyblue",size=0.4,linetype="dashed")+
  geom_hline(yintercept = log2(0.677),color="skyblue",size=0.4,linetype="dashed")+
  facet_grid(cols = vars(group_F),scales="free")+ggtitle(inter_join_location_1$full_names)

ggplot(inter_join_location_1,aes(as.numeric (group_number), mean_KO,color=factor(final_sig,levels=c('0','1','2')),alpha =factor(qvalue_0_05,levels=c("0","1"))))+
  geom_point()+theme_bw()+theme_classic()+
  geom_hline(yintercept = 0,color="black",size=0.4,linetype="dashed")+
  scale_color_manual(breaks = c("0", "1", "2"), values=c("gray", "red", "blue"))+
  scale_alpha_manual(breaks = c("0", "1"), values=c(0.2, 1))+
  geom_line(aes(as.numeric (group_number), mean_WT),color="blue",linetype="dashed")+
  facet_grid(cols = vars(group_F),scales="free")+ggtitle(inter_join_location_1$full_names)

ggplot(inter_join_location_1,aes(as.numeric (group_number), as.numeric(as.character(mean_KO)),color=factor(final_sig,levels=c('0','1','2')), alpha =factor(qvalue_0_05,levels=c("0","1"))))+
  
  geom_point(aes(size=factor(active,levels=c('0','1'))))+
  theme_bw()+theme_classic()+
  geom_hline(yintercept = 0,color="black",size=0.4,linetype="dashed")+
  scale_color_manual(breaks = c("0", "1", "2"), values=c("gray", "red", "blue"))+
  scale_alpha_manual(breaks = c("0", "1"), values=c(0.2, 1))+
  scale_size_manual(breaks = c("0", "1"), values=c(2, 4))+
  geom_line(aes(as.numeric (group_number), as.numeric(as.character(mean_WT)),alpha =factor(active_WT,levels=c("0","1"))),color="blue",linetype="dashed")+
  facet_grid(rows = vars(full_names),cols = vars(group_F),scales="free")
 
ggplot(inter_join_location_1,aes(as.numeric (group_number), as.numeric(as.character(mean_KO))))+
  
  geom_point(aes(color=factor(final_sig,levels=c('0','1','2')), alpha =factor(qvalue_0_05,levels=c("0","1")),size=factor(active,levels=c('0','1'))))+
  theme_bw()+theme_classic()+
  geom_hline(yintercept = 0,color="black",size=0.4,linetype="dashed")+
  scale_color_manual(breaks = c("0", "1", "2"), values=c("gray", "red", "blue"))+
  scale_alpha_manual(breaks = c("0", "1"), values=c(0.2, 1))+
  scale_size_manual(breaks = c("0", "1"), values=c(2, 4))+
  geom_line(aes(as.numeric (group_number), as.numeric(as.character(mean_WT)),alpha =factor(active_WT,levels=c("0","1"))),color="blue",linetype="dashed")+
  facet_grid(rows = vars(full_names),cols = vars(group_F),scales="free")+
  geom_line(aes(x=as.numeric (group_number), y=as.numeric(as.character(mean_KO))), alpha=0.2)

ggplot(inter_join_location_1,aes(as.numeric (group_number), as.numeric(as.character(mean_KO))))+
  
  geom_point(aes(color=factor(final_sig,levels=c('0','1','2')), alpha =factor(qvalue_0_05,levels=c("0","1")),size=factor(active,levels=c('0','1'))))+
  theme_bw()+theme_classic()+
  geom_hline(yintercept = 0,color="black",size=0.4,linetype="dashed")+
  scale_color_manual(breaks = c("0", "1", "2"), values=c("gray", "red", "blue"))+
  scale_alpha_manual(breaks = c("0", "1"), values=c(0.2, 1))+
  scale_size_manual(breaks = c("0", "1"), values=c(2, 4))+
  geom_line(aes(as.numeric (group_number), as.numeric(as.character(mean_WT)),alpha =factor(active_WT,levels=c("0","1"))),color="blue",linetype="dashed")+
  geom_line(aes(as.numeric (group_number), as.numeric(as.character(KO_median)),alpha =factor(active_KO,levels=c("0","1"))),color="coral",linetype="dotdash")+
  facet_grid(rows = vars(full_names),cols = vars(group_F),scales="free")+
  geom_line(aes(x=as.numeric (group_number), y=as.numeric(as.character(mean_KO))), alpha=0.2)










ggplot(filter_norm_data_scaled_grouped_2[filter_norm_data_scaled_grouped_2$active_WT>0,], aes(x = group_F))+ geom_bar(color="skyblue")+theme_bw()+theme_classic()
geom_bar(aes(log2(fold),color=factor(fold_1_5,levels=c('0','1','2')),alpha =factor(qvalue_0_05,levels=c("0","1"))))



write.table(inter_join_location,file="10bp_tiled_log2_normalized_log2_ratio_fold_pvalue_final.txt",sep="\t")  








inter_join_location[inter_join_location$mean_WT == 1.82784864351185,]

head(inter_join)
WT_selct<- unique(inter_join[,c(2,7:10,12,15)])
row.names(WT_selct)<-WT_selct$WT
WT_selct <- merge(filter_norm_data,WT_selct,by=0)
write.table(WT_selct,file="10bp_Mutagenesis_WT_final.txt",sep="\t",quote = T)

