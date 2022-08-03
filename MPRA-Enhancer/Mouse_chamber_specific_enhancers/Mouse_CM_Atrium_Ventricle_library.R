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
library(fuzzySim)
library(pheatmap)
library(fdrtool)
library(data.table)
library(coin)
library(stringr)
setwd("~/Desktop/work/Bill/yangpo/MPRA_5/")
#data1<-read.table(file="MPRA_5_barcode_number.txt",header=T,row.names = 1)
#data1<-read.table(file="MPRA_5_UMI_number_header.txt",header=T,row.names = 1)
#data2<-read.table(file="MPRA_2_A_3-6_V3-6_g3-6.txt",header=T,row.names = 1)
#data2<-read.table(file="AV_MPRA_batch2_batch3.txt",header=T,row.names = 1)
data1<-read.table(file="MPRA_7/A-V_DNA5-A5-V5_UMI_counts_merged_header.txt",header=T,row.names = 1)
rownames(data1) <- gsub("__", "_", rownames(data1))

head(data1)
dim(data1)

data_sub <- data1
######################################################################
annotate_files <- read.table(file="MPRA_7/MPRA_AV_300bp_annotation_final_with_header.txt",header=T, row.names = 1)
annotate_files$Group6 <- paste0(annotate_files$ATF,annotate_files$VTF,annotate_files$AP300,annotate_files$VP300)
Neg_sub <- annotate_files[annotate_files$Group6=="0000",]
Negative_C_annotation <- as.data.frame(t(dplyr::select(as.data.frame(t(Neg_sub)), contains("ESC"))))
dim(Negative_C_annotation)
Negative_C_annotation$group_new <-"Neg"
annotate_files$group<-"new"
annotate_files[rownames(Negative_C_annotation),]$group <-"Neg"



annotate_files$MTF <- paste0(annotate_files$AMTF,annotate_files$VMTF)
annotate_files$P300 <- paste0(annotate_files$AP300,annotate_files$VP300)
annotate_files[(annotate_files$ATF==1 |annotate_files$ATF==2 |annotate_files$VTF==1| annotate_files$VTF==2) &
                 annotate_files$P300=="00",]$group <-"TF1-2"
annotate_files[(annotate_files$ATF==3 |annotate_files$ATF==4 |annotate_files$VTF==3| annotate_files$VTF==4) &
                 annotate_files$P300=="00",]$group <-"TF3-4"
annotate_files[(annotate_files$MTF != "00" | annotate_files$P300 != "00"),]$group<-"MTF or P300"


annotate_files$spMTF <- paste0(annotate_files$AspMTF,annotate_files$VspMTF)
annotate_files$spP300 <- paste0(annotate_files$AspP300,annotate_files$VspP300)
#annotate_files$Gene <- paste0(annotate_files$aGene,annotate_files$vGene)
annotate_files$spMTF_spP300 <- paste0(annotate_files$spMTF,annotate_files$spP300)
annotate_files[annotate_files$spMTF_spP300 != "0000" ,]$group <-"spMTF or spP300"

table(annotate_files$group)
annotate_files_sub <- annotate_files[annotate_files$group != "new",]

write.table(annotate_files_sub,"AVlibrary_400bp_final_annotation_3897.txt",sep="\t",quote=F)
data1$group <- annotate_files_sub[rownames(data1),"group"]
data_sub <- data1[!is.na(data1$group),]
##########################################################################
#data_sub <- merge(data1, annotate_files_sub, by=rownames)
data_sub <- data_sub[,c(1:15)]
#######################################################################3
norm_data <- t(t(data_sub)/colSums(data_sub)*1000000)
total_sum <-as.data.frame(colSums(data_sub))
ggplot(total_sum, aes(x=row.names(total_sum),y=colSums(data_sub)))+ geom_bar(color="skyblue",stat="identity")+theme_bw()+theme_classic()
for_plot <- as.data.frame(norm_data) 
ggplot(for_plot, aes(x = DNA1))+ geom_histogram(binwidth = 5,color="skyblue")+theme_bw()+theme_classic()+
  geom_vline(xintercept = 5,color="skyblue",size=0.4,linetype="dashed")+xlim(0,500)
ggplot(for_plot, aes(x = DNA2))+ geom_histogram(binwidth = 5,color="skyblue")+theme_bw()+theme_classic()+
  geom_vline(xintercept = 10,color="skyblue",size=0.4,linetype="dashed")+xlim(0,500)
ggplot(for_plot, aes(x = DNA3))+ geom_histogram(binwidth = 5,color="skyblue")+theme_bw()+theme_classic()+
  geom_vline(xintercept = 5,color="skyblue",size=0.4,linetype="dashed")+xlim(0,500)
ggplot(for_plot, aes(x = DNA4))+ geom_histogram(binwidth = 5,color="skyblue")+theme_bw()+theme_classic()+
  geom_vline(xintercept = 10,color="skyblue",size=0.4,linetype="dashed")+xlim(0,500)
ggplot(for_plot, aes(x = DNA5))+ geom_histogram(binwidth = 5,color="skyblue")+theme_bw()+theme_classic()+
  geom_vline(xintercept = 10,color="skyblue",size=0.4,linetype="dashed")+xlim(0,500)
#filter_norm_data <- cbind(norm_data, data_annotation) %>% dplyr::filter(rowMins(as.matrix(norm_data[,c(5,14:21)]), value = T) >=10)

#DNA_total <-as.data.frame(rowMin(norm_data[,c(11:15)]))
DNA_total <-as.data.frame(rowMax(norm_data[,c(11:15)]))
colnames(DNA_total)<-c("DNA_min_value")
pdf("AV_400bp_DNA_counts_distribution.pdf")
ggplot(DNA_total, aes(x = DNA_min_value))+ geom_histogram(binwidth = 5,color="skyblue")+theme_bw()+theme_classic()+
  geom_vline(xintercept = 5,color="skyblue",size=0.4,linetype="dashed")+xlim(0,500)
dev.off()

grep("chr9:109730912-109731312",rownames(norm_data))
#filter_norm_data <- norm_data[rowMins(as.matrix(norm_data[,c(11:15)]), value = T) >=5,]
filter_norm_data <- norm_data[rowMaxs(as.matrix(norm_data[,c(11:15)]), value = T) >=5,]
dim(filter_norm_data)
filter_norm_data <- as.data.frame(filter_norm_data)
#filter_norm_data$mean_DNA <-(filter_norm_data$g1+filter_norm_data$g3+filter_norm_data$g4+filter_norm_data$g5+filter_norm_data$g6)/5
#filter_norm_data$mean_DNA <-(filter_norm_data$DNA1+filter_norm_data$DNA2+filter_norm_data$DNA3+filter_norm_data$DNA4+filter_norm_data$DNA5)/5
####################################################################3
total_number <- length(DNA_total[,1])
kept_number <- length(filter_norm_data[,1])
keptratio <- kept_number/total_number
DNA_counts <- data.frame(
  label=c("total_number","kept_number","keptratio"),
  value=c(as.character(total_number),as.character(kept_number),as.character(keptratio))
)
filter_norm_data$group <- annotate_files_sub[rownames(filter_norm_data),"group"]
write.table(DNA_counts,"AVlibrary_DNA_distribution_filter_log.txt",sep="\t",quote=F)
numbers<-as.data.frame(table(filter_norm_data$group))
write.table(numbers,"AVlibrary_DNA_distribution_filter_log_diff_group.txt",sep="\t",quote=F)
####################################################################
filter_for_degs<-data_sub[,1:15]
#filter_for_degs<- data_sub[rowMins(as.matrix(norm_data[,13:18]), value = T) >=20,13:24]
coldata <-as.data.frame(colnames(filter_for_degs))

coldata$type <- coldata$`colnames(filter_for_degs)`
coldata$condition <-gsub("\\d","",coldata$type)
genes <- as.matrix(data_sub)
#genes <-as.matrix((filter_for_degs))
#mode(genes) <-"integer"
dds <- DESeqDataSetFromMatrix(countData = genes, colData=coldata, design=~condition)
dds <- DESeq(dds)
res <- results( dds, contrast = c("condition","A","DNA") )
res2 <- results( dds, contrast = c("condition","V","DNA") )
res2 <- res2[!is.na(res$padj),]
res <- res[!is.na(res$padj),]
res$qvalue_0_05 <- 0
res$qvalue_0_05[res$padj < 0.05 & res$log2FoldChange >0] <- 1
sum(res$qvalue_0_05)
res2$qvalue_0_05 <- 0
res2$qvalue_0_05[res2$padj < 0.05 & res2$log2FoldChange >0] <- 1
sum(res2$qvalue_0_05)

############################################################################
new <- data.frame(matrix(ncol = 0, nrow = length(rownames(filter_norm_data))))
rownames(new) <- rownames(filter_norm_data)
filter_norm_data <- as.data.frame(filter_norm_data)
filter_norm_data$DNA <-rowMeans(filter_norm_data[,11:15])
#for (i in 1:10){
#  ratio <- filter_norm_data[,i]+1/filter_norm_data$DNA
#  ratio <- log2(ratio+1)
#  new <- cbind(new,ratio)
#}
for (i in 1:10){
  ratio <- (filter_norm_data[,i]+1)/(filter_norm_data$DNA+1)
  ratio <- log2(ratio)
  new <- cbind(new,ratio)
}
new_colname <- paste0(colnames(filter_norm_data[,c(1:10)]), "_vs_DNA_log")
colnames(new) <- new_colname
head(new)
#filter_norm_data["aMTF+p300-101-0::chr11:86916688-86917088",]
new$log_A_ratio <- rowMeans(new[,1:5])
new$log_V_ratio <- rowMeans(new[,6:10])
as.data.frame(colnames(filter_norm_data))
as.data.frame(colnames(new))

######################size factor#####################################
ggplot(new,aes(x=A1_vs_DNA_log,A3_vs_DNA_log))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(-8,8)+ylim(-8,8)
ggplot(new,aes(x=A1_vs_DNA_log,A2_vs_DNA_log))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(-8,8)+ylim(-8,8)
ggplot(new,aes(x=A4_vs_DNA_log,A3_vs_DNA_log))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(-8,8)+ylim(-8,8)
ggplot(new,aes(x=V1_vs_DNA_log,V3_vs_DNA_log))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(-8,8)+ylim(-8,8)
ggplot(new,aes(x=V3_vs_DNA_log,V2_vs_DNA_log))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(-8,8)+ylim(-8,8)
ggplot(new,aes(x=V4_vs_DNA_log,V3_vs_DNA_log))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(-8,8)+ylim(-8,8)
##################pcc figures##############################
my_cols2 <-  colorRampPalette(c("blue","white", "orange"))(n = 299)
my_cols <- c("black")  
panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits=2)
  txt <- paste0("R = ", r)
  cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor*2 )
}

# Customize upper panel
upper.panel<-function(x, y){
  points(x,y, pch = 19, col = my_cols[iris$Species],cex=0.1)
}
# Create the plots
pairs(new[,1:10], 
      lower.panel = panel.cor,
      upper.panel = upper.panel)


data_test<-new[,c(1:10)]
pdf("MPRA_5_pcc.pdf",width = 11, height = 11)
A_V_cor <- cor(data_test)
p<-corrplot(A_V_cor, method  = "number",order = "alphabet")
pairs(new[,1:10], 
      lower.panel = panel.cor,
      upper.panel = upper.panel)
dev.off()

################################################################



#########################Barcode and UMI###########################################3
#barcode<-data_test
#compare <- intersect(rownames(barcode),rownames(data_test))
#A_V_cor <- cbind(data_test[compare,],barcode[compare,])
#colnames(A_V_cor) <- c("UMI_A1","UMI_A2","UMI_A3","UMI_A4","UMI_A5","UMI_V1","UMI_V2","UMI_V3","UMI_V4","UMI_V5","raw_A1","raw_A2","raw_A3","raw_A4","raw_A5","raw_V1","raw_V2","raw_V3","raw_V4","raw_V5")
#A_V_cor2 <- cor(A_V_cor)
#library(RColorBrewer)
#corrplot(A_V_cor2, method  = "number",order = "alphabet",number.font = 2,number.cex = 0.7)
###################################################################
write.table(filter_norm_data,file="MPRA_7/A_V_mouse_MPRA_7_rawcount_version_add_1_filtered_new.txt",sep="\t")
write.table(new,file="MPRA_7/A_V_mouse_MPRA_7_normlized_count_version_add_1_filtered_new.txt",sep="\t")
write.table(filter_norm_data,file="MPRA_7/A_V_mouse_MPRA_7_UMI_rawcount_version_add_1_filtered_new.txt",sep="\t")
write.table(new,file="MPRA_7/A_V_mouse_MPRA_7_normlized_UMIcount_version_add_1_filtered.txt_new",sep="\t")
new <-read.table(file="MPRA_7/A_V_mouse_MPRA_7_normlized_UMIcount_version_add_1_filtered.txt_new",header = T,row.names = 1)
#new <-cbind(new,filter_norm_data[,c(22:37)])

head(new)
filter_norm_data <- new
filter_norm_data$log_A_ratio <- rowMeans(filter_norm_data[,1:5])
filter_norm_data$log_V_ratio <- rowMeans(filter_norm_data[,6:10])
filter_norm_data$order_log_A_ratio <- rank(-as.numeric(as.character(filter_norm_data$log_A_ratio)))
filter_norm_data$order_log_V_ratio <- rank(-as.numeric(as.character(filter_norm_data$log_A_ratio)))
head(filter_norm_data)
################################Size factor###################################
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

new <- filter_norm_data
Negative_C <- t(dplyr::select(as.data.frame(t(new)), contains("ESC")))
#Positive_C <- t(dplyr::select(as.data.frame(t(new)), order_log_A_ratio<=1000 &order_log_V_ratio<=1000))
Positive_C <-new[new$order_log_A_ratio<=1000 &new$order_log_V_ratio<=1000,]
#Positive_C <-t(dplyr::select(as.data.frame(t(new)), contains("BNA")))
both_C <- unique(c(rownames(Negative_C),rownames(Positive_C)))
filter_norm_data_nc_normalize_exp <-MPRA_SizeFactors(new,both_C,c(1:10))
filter_norm_data_nc_normalize_exp <- as.data.frame(filter_norm_data_nc_normalize_exp)

#################################################################################################
head(filter_norm_data_nc_normalize_exp)
#filter_norm_data_scaled <- cbind(new[,11:14],filter_norm_data_nc_normalize_exp[,1:10])
filter_norm_data_scaled <- as.data.frame(as.matrix(filter_norm_data_nc_normalize_exp))
head(filter_norm_data_scaled)
filter_norm_data_scaled$log_A_ratio <- rowMeans(filter_norm_data_scaled[,1:5])
filter_norm_data_scaled$log_V_ratio <- rowMeans(filter_norm_data_scaled[,6:10])

filter_norm_data_scaled$order_log_A_mean_ratio <- rank(-as.numeric(as.character(filter_norm_data_scaled$log_A_ratio)),ties.method = "first")
filter_norm_data_scaled$order_log_V_mean_ratio <- rank(-as.numeric(as.character(filter_norm_data_scaled$log_V_ratio)),ties.method = "first")
#data_test<-filter_norm_data_scaled[,c(1:12)]
#A_V_cor <- cor(data_test)
#corrplot(A_V_cor, method  = "number")
#ggplot(filter_norm_data_scaled,aes(x=log_A1_ratio,y=log_A2_ratio))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(0,8)+ylim(0,8)
#ggplot(filter_norm_data_scaled,aes(x=log_V1_ratio,y=log_V2_ratio))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(0,8)+ylim(0,8)
write.table(filter_norm_data_scaled,file="MPRA_7/A_V_mouse_MPRA_7_rawcounts_scaled_version_add_1_filtered_new.txt",sep="\t")
#write.table(filter_norm_data_scaled,file="MPRA_7/A_V_mouse_MPRA_7_rawcounts_scaled_version_add_1_BNA_normalize.txt",sep="\t")
#write.table(filter_norm_data_scaled,file="MPRA_7/A_V_mouse_MPRA_5_UMIscaled_version_add_1_filtered_new.txt",sep="\t")
filter_norm_data_scaled <- read.table(file = "MPRA_7/A_V_mouse_MPRA_7_rawcounts_scaled_version_add_1_filtered_new.txt",header = T, row.names = 1)
head(filter_norm_data_scaled)
#########################################FDR 0.05#######################
###############################annotation ########################################
#annotate_files <- read.table(file="MPRA_7/MPRA_new_annotation.txt",header=T, row.names = 1)
#annotate_files <- read.table(file="MPRA_7/MRAP_annotation_add_MTF_12_P300_2cols_final_new.txt",header=T, row.names = 1)
#annotate_files <- read.table(file="MPRA_annotation_add_P300_2_final_newest_with_header.txt",header=T, row.names = 1)
filter_norm_data_scaled$group <- annotate_files[rownames(filter_norm_data_scaled),"group"]
filter_norm_data_scaled$ATF <- annotate_files[rownames(filter_norm_data_scaled),"ATF"]
filter_norm_data_scaled$VTF <- annotate_files[rownames(filter_norm_data_scaled),"VTF"]
filter_norm_data_scaled$ATF1_2 <-0
#filter_norm_data_scaled[(filter_norm_data_scaled$ATF==1 |filter_norm_data_scaled$ATF==2) &filter_norm_data_scaled$group=="TF1-2",]$ATF1_2 <-1
filter_norm_data_scaled[(filter_norm_data_scaled$ATF==1 |filter_norm_data_scaled$ATF==2) ,]$ATF1_2 <-1
filter_norm_data_scaled$VTF1_2 <-0
filter_norm_data_scaled[(filter_norm_data_scaled$VTF==1 |filter_norm_data_scaled$VTF==2) ,]$VTF1_2 <-1
filter_norm_data_scaled$ATF3_4 <-0
filter_norm_data_scaled[(filter_norm_data_scaled$ATF==3 |filter_norm_data_scaled$ATF==4) ,]$ATF3_4 <-1
filter_norm_data_scaled$VTF3_4 <-0
filter_norm_data_scaled[(filter_norm_data_scaled$VTF==3 |filter_norm_data_scaled$VTF==4),]$VTF3_4 <-1
filter_norm_data_scaled$MTF <- annotate_files[rownames(filter_norm_data_scaled),"MTF"]
filter_norm_data_scaled$P300 <- annotate_files[rownames(filter_norm_data_scaled),"P300"]
filter_norm_data_scaled$spMTF <- annotate_files[rownames(filter_norm_data_scaled),"spMTF"]
filter_norm_data_scaled$spP300 <- annotate_files[rownames(filter_norm_data_scaled),"spP300"]
filter_norm_data_scaled$MTF_P300 <- annotate_files[rownames(filter_norm_data_scaled),"MTF_P300"]
filter_norm_data_scaled$spMTF_spP300 <- annotate_files[rownames(filter_norm_data_scaled),"spMTF_spP300"]
filter_norm_data_scaled$AMTF <- annotate_files[rownames(filter_norm_data_scaled),"AMTF"]
filter_norm_data_scaled$VMTF <- annotate_files[rownames(filter_norm_data_scaled),"VMTF"]
filter_norm_data_scaled$AspMTF <- annotate_files[rownames(filter_norm_data_scaled),"AspMTF"]
filter_norm_data_scaled$VspMTF <- annotate_files[rownames(filter_norm_data_scaled),"VspMTF"]

filter_norm_data_scaled$AP300 <- annotate_files[rownames(filter_norm_data_scaled),"AP300"]
filter_norm_data_scaled$VP300 <- annotate_files[rownames(filter_norm_data_scaled),"VP300"]
filter_norm_data_scaled$AspP300 <- annotate_files[rownames(filter_norm_data_scaled),"AspP300"]
filter_norm_data_scaled$VspP300 <- annotate_files[rownames(filter_norm_data_scaled),"VspP300"]

filter_norm_data_scaled$Group6 <- annotate_files[rownames(filter_norm_data_scaled),"Group6"]

filter_norm_data_scaled$order_log_A_mean_ratio <- rank(-as.numeric(as.character(filter_norm_data_scaled$log_A_ratio)),ties.method = "first")
filter_norm_data_scaled$order_log_V_mean_ratio <- rank(-as.numeric(as.character(filter_norm_data_scaled$log_V_ratio)),ties.method = "first")
#Negative_C <- t(dplyr::select(as.data.frame(t(filter_norm_data_scaled)), contains("ESC")))
Negative_C <- as.data.frame(filter_norm_data_scaled[filter_norm_data_scaled$Group6=="0000",])
Negative_C <-as.data.frame(Negative_C)
fitg<-fitdist(as.numeric(Negative_C$log_A_ratio) ,"norm")
FDR_A_0.95 <- qnorm(0.95, mean = fitg$estimate["mean"], sd = fitg$estimate["sd"])

fitg2<-fitdist(as.numeric(Negative_C$log_V_ratio) ,"norm")
FDR_V_0.95 <- qnorm(0.95, mean = fitg2$estimate["mean"], sd = fitg$estimate["sd"])

FDR_A_0.95
FDR_V_0.95
FDR_A_0.95<-2.318062
FDR_V_0.95<-2.480459
#ggplot(filter_norm_data_scaled, aes(x = log_V_mean, color=group))+geom_density()+theme_bw()+theme_classic()
ggplot(Negative_C, aes(x = as.numeric(log_A_ratio)))+geom_density()+theme_bw()+theme_classic()
ggplot(Negative_C, aes(x = as.numeric(log_V_ratio)))+geom_density()+theme_bw()+theme_classic()

################################################################################################
filter_norm_data_scaled$A_active <- res[rownames(filter_norm_data_scaled),]$qvalue_0_05
filter_norm_data_scaled$V_active <- res2[rownames(filter_norm_data_scaled),]$qvalue_0_05
table(filter_norm_data_scaled$A_active)
table(filter_norm_data_scaled$V_active)
#new <- rename(new, A_active=`res$qvalue_0_05`, V_active=`res2$qvalue_0_05`)
W1<-as.data.frame(table(filter_norm_data_scaled$A_active,filter_norm_data_scaled$group))
A1 <-as.data.frame(table(filter_norm_data_scaled$V_active,filter_norm_data_scaled$group))
dat_text_3 <- cbind(W1, A1) 
colnames(dat_text_3) <-c("A_active","group","number","V_active","group","number")
write.table(dat_text_3, "MPRA_A_V_active_regions_number.txt",quote=F, sep = "\t")
#####################################################group marker#################
write.table(filter_norm_data_scaled,file="MPRA_7/A_V_mouse_MPRA_7_UMIscaled_version_add_1_filtered_new_2.txt",sep="\t")
filter_norm_data_scaled <- read.table("MPRA_7/A_V_mouse_MPRA_7_UMIscaled_version_add_1_filtered_new_2.txt",header=T,row.names = 1)
library(caret)
SSA<-function(confmatrix){
  Acu<-as.data.frame(confusionMatrix(confmatrix)$overall["Accuracy"])[,1]
  Sen<-sensitivity(confmatrix)
  Spe<-specificity(confmatrix)
  SSA_matrix<-data.frame(label=c("Sensitivity","Specificiy","Accuracy"),
                         value=c(as.character(Sen),as.character(Spe),as.character(Acu)))
  return(SSA_matrix)
}
SSA_merge<-function(data,cols_P,cols_R){
  SSA_matrix_merge=as.data.frame(matrix(nrow = 0,ncol = 4))
  for (i in cols_P){
    for (j in cols_R) {
      predicted_values<-data[,i]
      actual_values<-data[,j]
      confmatrix<-table(predicted_values,actual_values)
      Acu<-as.data.frame(confusionMatrix(confmatrix)$overall["Accuracy"])[,1]
      Sen<-sensitivity(confmatrix)
      Spe<-specificity(confmatrix)
      SSA_matrix<-data.frame(label=c("Sensitivity","Specificity","Accuracy"),
                             value=c(as.character(Spe),as.character(Sen),as.character(Acu)))
      SSA_matrix$Factor <-as.character(i)
      SSA_matrix$cell <-as.character(j)
      SSA_matrix_merge<-rbind(SSA_matrix_merge,SSA_matrix)
    }
    }
  return (SSA_matrix_merge)
}
SSA_merge_matrix<-function(data,cols_P,cols_R){
  SSA_matrix_merge=as.data.frame(matrix(nrow = 0,ncol = 4))
  for (i in cols_P){
    for (j in cols_R) {
      predicted_values<-data[,i]
      actual_values<-data[,j]
      confmatrix<-table(predicted_values,actual_values)
      
      SSA_matrix<-as.matrix.data.frame(confmatrix)
      rownames(SSA_matrix) <-paste0(c(0,1),"_",as.character(i))
      colnames(SSA_matrix) <-paste0(c(0,1),"_",as.character(j))
      SSA_matrix_merge<-rbind(SSA_matrix_merge,SSA_matrix)
    }
  }
  return (SSA_matrix_merge)
}

model_filter_norm_data_scaled <- filter_norm_data_scaled[,c("AMTF","VMTF","A_active","V_active","AP300",
                                                            "VP300","ATF1_2","VTF1_2","ATF3_4","VTF3_4","Group6")]

ATAC_anntation<-read.table(file = "MPRA_7/new_annotation_sensitive/MPRA_with_A-V_ATAC_As_Vs_ATAC_full_list.txt")
loops_anntation<-read.table(file = "MPRA_7/new_annotation_sensitive/MPRA_with_A-V_loops_As_Vs_loops_full_list.txt")
ATAC_anntation$V1 <- gsub("__", "_", ATAC_anntation$V1)
loops_anntation$V1 <- gsub("__", "_", loops_anntation$V1)
TFs_anntation<-read.table(file = "MPRA_7/new_annotation_sensitive/MPRA_with_A-V_GACNSTD_full_list.txt")
TFs_anntation$V1 <- gsub("__", "_", TFs_anntation$V1)
dim(ATAC_anntation)

dim(filter_norm_data_scaled)
model_filter_norm_data_scaled$ATAC_A <-ATAC_anntation[ATAC_anntation$V1 %in% rownames(model_filter_norm_data_scaled),"V2"]
model_filter_norm_data_scaled$ATAC_V <-ATAC_anntation[ATAC_anntation$V1 %in% rownames(model_filter_norm_data_scaled),"V3"]
model_filter_norm_data_scaled$loops_A <-loops_anntation[loops_anntation$V1 %in% rownames(model_filter_norm_data_scaled),"V2"]
model_filter_norm_data_scaled$loops_V <-loops_anntation[loops_anntation$V1 %in% rownames(model_filter_norm_data_scaled),"V3"]
#########################################################################
model_filter_norm_data_scaled$Gata4_A <-TFs_anntation[TFs_anntation$V1 %in% rownames(model_filter_norm_data_scaled),"V2"]
model_filter_norm_data_scaled$Mef2a_A <-TFs_anntation[TFs_anntation$V1 %in% rownames(model_filter_norm_data_scaled),"V3"]
model_filter_norm_data_scaled$Mef2c_A <-TFs_anntation[TFs_anntation$V1 %in% rownames(model_filter_norm_data_scaled),"V4"]
model_filter_norm_data_scaled$Nkx2.5_A <-TFs_anntation[TFs_anntation$V1 %in% rownames(model_filter_norm_data_scaled),"V5"]
model_filter_norm_data_scaled$Srf_A <-TFs_anntation[TFs_anntation$V1 %in% rownames(model_filter_norm_data_scaled),"V6"]
model_filter_norm_data_scaled$Tbx5_A <-TFs_anntation[TFs_anntation$V1 %in% rownames(model_filter_norm_data_scaled),"V7"]
model_filter_norm_data_scaled$Tead1_A <-TFs_anntation[TFs_anntation$V1 %in% rownames(model_filter_norm_data_scaled),"V8"]
model_filter_norm_data_scaled$Gata4_V <-TFs_anntation[TFs_anntation$V1 %in% rownames(model_filter_norm_data_scaled),"V9"]
model_filter_norm_data_scaled$Mef2a_V <-TFs_anntation[TFs_anntation$V1 %in% rownames(model_filter_norm_data_scaled),"V10"]
model_filter_norm_data_scaled$Mef2c_V <-TFs_anntation[TFs_anntation$V1 %in% rownames(model_filter_norm_data_scaled),"V11"]
model_filter_norm_data_scaled$Nkx2.5_V <-TFs_anntation[TFs_anntation$V1 %in% rownames(model_filter_norm_data_scaled),"V12"]
model_filter_norm_data_scaled$Srf_V <-TFs_anntation[TFs_anntation$V1 %in% rownames(model_filter_norm_data_scaled),"V13"]
model_filter_norm_data_scaled$Tbx5_V <-TFs_anntation[TFs_anntation$V1 %in% rownames(model_filter_norm_data_scaled),"V14"]
model_filter_norm_data_scaled$Tead1_V <-TFs_anntation[TFs_anntation$V1 %in% rownames(model_filter_norm_data_scaled),"V15"]
SSA_matrix_merge_A_Tfs <- SSA_merge(model_filter_norm_data_scaled,c("Gata4_A","Mef2a_A","Mef2c_A","Nkx2.5_A",
                                                                "Srf_A","Tbx5_A","Tead1_A"),c("A_active"))
SSA_matrix_merge_V_Tfs <- SSA_merge(model_filter_norm_data_scaled,c("Gata4_V","Mef2a_V","Mef2c_V","Nkx2.5_V",
                                                                    "Srf_V","Tbx5_V","Tead1_V"),c("V_active"))
library(cowplot)
library(grid)
library(patchwork)
pdf("single_MTF_400bp_SSA.pdf",width = 7,height = 7)
A<-ggplot(SSA_matrix_merge_A_Tfs,aes(y=factor(Factor,levels=rev(c("Gata4_A","Mef2a_A","Mef2c_A","Nkx2.5_A",
                                                                  "Srf_A","Tbx5_A","Tead1_A"))),
                                 x=as.numeric(value),group=label,fill=factor(label)))+
  geom_bar(stat="identity")+ facet_wrap(vars(factor(label,levels=c("Specificity","Sensitivity","Accuracy"))),ncol = 3)+
  theme_classic()
V<-ggplot(SSA_matrix_merge_V_Tfs,aes(y=factor(Factor,levels =rev(c("Gata4_V","Mef2a_V","Mef2c_V","Nkx2.5_V",
                                                               "Srf_V","Tbx5_V","Tead1_V"))),
                                 x=as.numeric(value),group=label,fill=factor(label)))+
  geom_bar(stat="identity")+ facet_wrap(vars(factor(label,levels=c("Specificity","Sensitivity","Accuracy"))),ncol = 3)+
  theme_classic()

A/V
dev.off()
##########################################################################
model_filter_norm_data_scaled$A_merged <-0
model_filter_norm_data_scaled[rowSums(model_filter_norm_data_scaled[,c("AMTF","AP300",
                                                                       "ATAC_A","loops_A")])>0,]$A_merged<- 1
model_filter_norm_data_scaled$V_merged <-0
model_filter_norm_data_scaled[rowSums(model_filter_norm_data_scaled[,c("VMTF","VP300",
                                                                       "ATAC_V","loops_V")])>0,]$V_merged<- 1

SSA_matrix_merge_A <- SSA_merge(model_filter_norm_data_scaled,c("AMTF","AP300",
                                                          "ATAC_A","loops_A","A_merged"),c("A_active"))
SSA_matrix_merge_V <- SSA_merge(model_filter_norm_data_scaled,c("VMTF","VP300","ATAC_V","loops_V","V_merged"),
                                c("V_active"))
SSA_matrix_merge_A_no_ng <- SSA_merge(model_filter_norm_data_scaled[model_filter_norm_data_scaled$Group6 !="0",],
                                      c("AMTF","AP300", "ATAC_A","loops_A","A_merged"),c("A_active"))
SSA_matrix_merge_V_no_bg <- SSA_merge(model_filter_norm_data_scaled[model_filter_norm_data_scaled$Group6 !="0",],c("VMTF","VP300",
                                           "ATAC_V","loops_V","V_merged"),
                                c("V_active"))
SSA_matrix_merged <- rbind(SSA_matrix_merge_A,SSA_matrix_merge_V)
SSA_matrix_merged <- rbind(SSA_matrix_merge_A_no_ng,SSA_matrix_merge_V_no_bg)
SSA_matrix_merge_A_matrix <- SSA_merge_matrix(model_filter_norm_data_scaled,c("AMTF","AP300",
                                                                "ATAC_A","loops_A","ATF1_2","ATF3_4","A_merged"),c("A_active"))
SSA_matrix_merge_V_matrix <- SSA_merge_matrix(model_filter_norm_data_scaled,c("VMTF","VP300","ATAC_V","loops_V",
                                                                              "V_merged"),
                                c("V_active"))
write.table(SSA_matrix_merge_A_matrix, file="A_confmatrix.txt",sep = "\t",quote = F)
write.table(SSA_matrix_merge_V_matrix, file="V_confmatrix.txt",sep = "\t",quote = F)
library(cowplot)
library(grid)
library(patchwork)
pdf("SSA_A_MTF_P300_gene_ATAC_loops_merged.pdf",width = 7,height = 6)
A<-ggplot(SSA_matrix_merge_A,aes(y=factor(Factor,levels=rev(c("AMTF","AP300",
                                                          "Agene","ATAC_A","loops_A","ATF1_2","ATF3_4","A_merged"))),
                                 x=as.numeric(value),group=label,fill=factor(label)))+
  geom_bar(stat="identity")+ facet_wrap(vars(factor(label,levels=c("Specificity","Sensitivity","Accuracy"))),ncol = 3)+
  theme_classic()
V<-ggplot(SSA_matrix_merge_V,aes(y=factor(Factor,levels =rev(c("VMTF","VP300","Vgene","ATAC_V","loops_V","VTF1_2","VTF3_4","V_merged"))),
                                 x=as.numeric(value),group=label,fill=factor(label)))+
  geom_bar(stat="identity")+ facet_wrap(vars(factor(label,levels=c("Specificity","Sensitivity","Accuracy"))),ncol = 3)+
  theme_classic()

A/V
A_no<-ggplot(SSA_matrix_merge_A_no_ng,aes(y=factor(Factor,levels=rev(c("AMTF","AP300",
                                                              "Agene","ATAC_A","loops_A","ATF1_2","ATF3_4","A_merged"))),
                                 x=as.numeric(value),group=label,fill=factor(label)))+
  geom_bar(stat="identity")+ facet_wrap(vars(factor(label,levels=c("Specificity","Sensitivity","Accuracy"))),ncol = 3)+
  theme_classic()
V_no<-ggplot(SSA_matrix_merge_V_no_bg,aes(y=factor(Factor,levels =rev(c("VMTF","VP300","Vgene","ATAC_V","loops_V","VTF1_2","VTF3_4","V_merged"))),
                                 x=as.numeric(value),group=label,fill=factor(label)))+
  geom_bar(stat="identity")+ facet_wrap(vars(factor(label,levels=c("Specificity","Sensitivity","Accuracy"))),ncol = 3)+
  theme_classic()

A_no/V_no
dev.off()

confmatrix<-table(model_filter_norm_data_scaled$Agene,model_filter_norm_data_scaled$A_active)
write.table(SSA_matrix_merged,"AVlibrary_Sensitivity_specitivity_Accuracy_merged.txt",sep="\t",quote=F)
write.table(SSA_matrix_merged,"AVlibrary_Sensitivity_specitivity_Accuracy_merged.txt",sep="\t",quote=F)
#################################################################################

filter_norm_data_scaled$rawname<-str_split_fixed(rownames(filter_norm_data_scaled), pattern="::",2)[,1]
filter_norm_data_scaled$region_locs<-str_split_fixed(rownames(filter_norm_data_scaled), pattern="::",2)[,2]
################################################################################
#filter_norm_data_scaled_sub <-filter_norm_data_scaled[(filter_norm_data_scaled$Group6=="0" | filter_norm_data_scaled$AMTF==1 |
#                                                         filter_norm_data_scaled$AP300==1 |
#                                                         filter_norm_data_scaled$ATF1_2==1 |
#                                                         filter_norm_data_scaled$ATF3_4==1),
#                                                          c("Group6","AMTF","AP300","ATF1_2","ATF3_4","order_log_A_mean_ratio")]
filter_norm_data_scaled_sub1<-filter_norm_data_scaled[filter_norm_data_scaled$Group6=="0",] 
filter_norm_data_scaled_sub1$marker<-"NC"
filter_norm_data_scaled_sub2<-filter_norm_data_scaled[filter_norm_data_scaled$AMTF==1,] 
filter_norm_data_scaled_sub2$marker<-"AMTF"
filter_norm_data_scaled_sub3<-filter_norm_data_scaled[filter_norm_data_scaled$AP300==1,] 
filter_norm_data_scaled_sub3$marker<-"AP300"
filter_norm_data_scaled_sub4<-filter_norm_data_scaled[filter_norm_data_scaled$ATF1_2==1,] 
filter_norm_data_scaled_sub4$marker<-"ATF1_2"
filter_norm_data_scaled_sub5<-filter_norm_data_scaled[filter_norm_data_scaled$ATF3_4==1,] 
filter_norm_data_scaled_sub5$marker<-"ATF3_4"
#filter_norm_data_scaled_sub4<-filter_norm_data_scaled[filter_norm_data_scaled$Agene==1,] 
#filter_norm_data_scaled_sub4$marker<-"AGene"
filter_norm_data_scaled_sub_merged_A <- rbind(filter_norm_data_scaled_sub1[,c("order_log_A_mean_ratio","marker","log_A_ratio")],
                                            filter_norm_data_scaled_sub2[,c("order_log_A_mean_ratio","marker","log_A_ratio")],
                                            filter_norm_data_scaled_sub3[,c("order_log_A_mean_ratio","marker","log_A_ratio")],
                                            filter_norm_data_scaled_sub4[,c("order_log_A_mean_ratio","marker","log_A_ratio")],
                                            filter_norm_data_scaled_sub5[,c("order_log_A_mean_ratio","marker","log_A_ratio")])
#filter_norm_data_scaled_sub <-filter_norm_data_scaled[(filter_norm_data_scaled$Group6=="000000" | filter_norm_data_scaled$vMTF_2==2 |filter_norm_data_scaled$vGene_3==1 |
#                                                         filter_norm_data_scaled$vP300==1 ),
#                                                      c("Group6","vMTF_2","vGene_3","vP300","order_log_V_mean_ratio")]
filter_norm_data_scaled_sub1<-filter_norm_data_scaled[filter_norm_data_scaled$Group6=="0",] 
filter_norm_data_scaled_sub1$marker<-"NC"
filter_norm_data_scaled_sub2<-filter_norm_data_scaled[filter_norm_data_scaled$VMTF==1,] 
filter_norm_data_scaled_sub2$marker<-"VMTF"
filter_norm_data_scaled_sub3<-filter_norm_data_scaled[filter_norm_data_scaled$VP300==1,] 
filter_norm_data_scaled_sub3$marker<-"VP300"
filter_norm_data_scaled_sub4<-filter_norm_data_scaled[filter_norm_data_scaled$VTF1_2==1,] 
filter_norm_data_scaled_sub4$marker<-"VTF1_2"
filter_norm_data_scaled_sub5<-filter_norm_data_scaled[filter_norm_data_scaled$VTF3_4==1,] 
filter_norm_data_scaled_sub5$marker<-"VTF3_4"
#filter_norm_data_scaled_sub4$marker<-"VGene"
filter_norm_data_scaled_sub_merged_V <- rbind(filter_norm_data_scaled_sub1[,c("order_log_V_mean_ratio","marker","log_V_ratio")],
                                            filter_norm_data_scaled_sub2[,c("order_log_V_mean_ratio","marker","log_V_ratio")],
                                            filter_norm_data_scaled_sub3[,c("order_log_V_mean_ratio","marker","log_V_ratio")],
                                            filter_norm_data_scaled_sub4[,c("order_log_V_mean_ratio","marker","log_V_ratio")],
                                            filter_norm_data_scaled_sub5[,c("order_log_V_mean_ratio","marker","log_V_ratio")])
library(patchwork)

pdf("new_MTF_single_AV_final.pdf", width =12,height = 8)
A<-ggplot(filter_norm_data_scaled,aes(x=order_log_A_mean_ratio,y=log_A_ratio,color=factor(A_active)))+
  geom_point(size=2)+
  #geom_hline(yintercept=FDR_A_0.95,color="red",linetype="dashed")+
  theme_bw()+
  scale_color_manual(breaks = c("1", "0"), values=c("red", "gray"))+

  theme_classic()
  #+geom_text(aes(0,round(FDR_A_0.95,3),label = round(FDR_A_0.95,3), vjust = -1))

B<-ggplot(filter_norm_data_scaled_sub_merged_A, aes(x = order_log_A_mean_ratio, y = factor(marker,levels=c("NC","AGene",
                                                                                                         "AP300","AMTF","ATF1_2","ATF3_4"))))+
         stat_bin2d(aes(fill = after_stat(count)), binwidth = c(1,1),color="Black")+theme_bw()
A/B
A<-ggplot(filter_norm_data_scaled,aes(x=order_log_A_mean_ratio,y=log_A_ratio,color=factor(A_active)))+
  geom_point(size=2)+
  #geom_hline(yintercept=FDR_A_0.95,color="red",linetype="dashed")+
  theme_bw()+
  scale_color_manual(breaks = c("1", "0"), values=c("red", "gray"))+
  theme_classic()
  #geom_text(aes(0,round(FDR_A_0.95,3),label = round(FDR_A_0.95,3), vjust = -1))

B<-ggplot(filter_norm_data_scaled_sub_merged_A, aes(x = order_log_A_mean_ratio, y = factor(marker,levels=c("NC",
                                                                                                           "AP300","AMTF","ATF1_2","ATF3_4"))))+
  stat_bin2d(aes(fill = after_stat(count)), binwidth = c(3,1))+theme_bw()
A/B
C<-ggplot(filter_norm_data_scaled,aes(x=order_log_V_mean_ratio,y=log_V_ratio,color=factor(V_active)))+
  geom_point(size=2)+
  geom_hline(yintercept=FDR_V_0.95,color="red",linetype="dashed")+theme_bw()+
  scale_color_manual(breaks = c("1", "0"), values=c("red", "gray"))+
  theme_classic()+geom_text(aes(0,round(FDR_V_0.95,3),label = round(FDR_V_0.95,3), vjust = -1))
D<-ggplot(filter_norm_data_scaled_sub_merged_V, aes(x = order_log_V_mean_ratio, y = factor(marker,levels=c("NC","VGene",
                                                                                                         "VP300","VMTF","VTF1_2","VTF3_4"))))+
  stat_bin2d(aes(fill = after_stat(count)), binwidth = c(1,1),color="Black")+theme_bw()
C/D
C<-ggplot(filter_norm_data_scaled,aes(x=order_log_V_mean_ratio,y=log_V_ratio,color=factor(V_active)))+
  geom_point(size=2)+
  #geom_hline(yintercept=FDR_V_0.95,color="red",linetype="dashed")+
  theme_bw()+
  scale_color_manual(breaks = c("1", "0"), values=c("red", "gray"))+
  theme_classic()
  #geom_text(aes(0,round(FDR_V_0.95,3),label = round(FDR_V_0.95,3), vjust = -1))
D<-ggplot(filter_norm_data_scaled_sub_merged_V, aes(x = order_log_V_mean_ratio, y = factor(marker,levels=c("NC","VGene",
                                                                                                           "VP300","VMTF","VTF1_2","VTF3_4"))))+
  stat_bin2d(aes(fill = after_stat(count)), binwidth = c(3,1))+theme_bw()
C/D
dev.off()
############################################
library(jmuOutlier) 
library(distr6)
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
    #mean_value_p <- mean(datasets3$Enrich_score[datasets3$Enrich_score >=0])
    #mean_value_n <- mean(datasets3$Enrich_score[datasets3$Enrich_score <0])
    #datasets3$NES[datasets3$Enrich_score >=0] <- 
    #  datasets3$Enrich_score[datasets3$Enrich_score >=0]/mean_value_p
    #datasets3$NES[datasets3$Enrich_score <0] <- 
    #  datasets3$Enrich_score[datasets3$Enrich_score <0]/mean_value_n
    #datasets3$negative_Enrich_score <- (Max_pc1-datasets3[,"group_order"]+1)/Max_pc1 - (Max-datasets3[,order_col]+1)/Max
    #mean_value_negative_p <- mean(datasets3$negative_Enrich_score[datasets3$negative_Enrich_score >=0])
    #mean_value_negative_n <- mean(datasets3$negative_Enrich_score[datasets3$negative_Enrich_score <0])
    ##datasets3$negative_NES[datasets3$negative_Enrich_score >=0] <-
    ##  -1*datasets3$negative_Enrich_score[datasets3$negative_Enrich_score >=0]/mean_value_negative_p
    ##datasets3$negative_NES[datasets3$negative_Enrich_score <0] <-
    ##  datasets3$negative_Enrich_score[datasets3$negative_Enrich_score <0]/mean_value_negative_n
    pvalue <- c(pvalue, p_PC1)
    groups_names <- c(groups_names,cols_name)
    datasets3$score_group <- cols_name
    #datasets3$final_ES <- datasets3$Enrich_score - datasets3$negative_Enrich_score
    datasets3$final_ES <- datasets3$Enrich_score
    #if (max(abs(datasets3$negative_Enrich_score)) > max(abs(datasets3$Enrich_score))){
    #  datasets3$final_ES <-
    #  datasets3$negative_NES
    #}
    #else {
    #  datasets3$final_ES <-datasets3$NES
    #}
    datasets2 <- rbind(datasets2,datasets3)
    
  }  
  dat_text <- data.frame(
    label = pvalue,
    score_group = groups_names
  )
  
  newlist <-list("dat_text"=dat_text, "datasets"=datasets2)
  return(newlist)
}
#############################################3
filter_norm_data_scaled_grouped_test <-Enrich_score(filter_norm_data_scaled, "AMTF", "log_A_ratio", "order_log_A_mean_ratio")
dat_text_A <-filter_norm_data_scaled_grouped_test$dat_text
filter_norm_data_for_plot_A <- filter_norm_data_scaled_grouped_test$datasets
filter_norm_data_for_plot_A_1 <-filter_norm_data_for_plot_A[,c("order_log_A_mean_ratio","final_ES","score_group")]
filter_norm_data_for_plot_A_1$score_group_2<-0
filter_norm_data_for_plot_A_1[filter_norm_data_for_plot_A_1$score_group==1,]$score_group_2<-"AMTF"
filter_norm_data_scaled_grouped_test_2 <-Enrich_score(filter_norm_data_scaled, "AP300", "log_A_ratio", "order_log_A_mean_ratio")
dat_text_A_2 <-filter_norm_data_scaled_grouped_test_2$dat_text
filter_norm_data_for_plot_A_2 <- filter_norm_data_scaled_grouped_test_2$datasets
filter_norm_data_for_plot_A_2 <-filter_norm_data_for_plot_A_2[,c("order_log_A_mean_ratio","final_ES","score_group")]
filter_norm_data_for_plot_A_2$score_group_2<-0
filter_norm_data_for_plot_A_2[filter_norm_data_for_plot_A_2$score_group==1,]$score_group_2<-"AP300"

filter_norm_data_scaled$Neg <-0
filter_norm_data_scaled[filter_norm_data_scaled$Group6=="0000",]$Neg <-1

data_for_plot_e_A <-as.data.frame(rbind(filter_norm_data_for_plot_A_1,filter_norm_data_for_plot_A_2))
##################################################
#################################################
##################################################
################################################
pdf("new_MTF_single_AV_final.pdf", width =12,height = 8)
p_A<-ggplot(data_for_plot_e_A, aes(order_log_A_mean_ratio, final_ES,color=score_group_2))+geom_line()+theme_bw()+theme_classic()+
  #facet_wrap(~score_group,scales = "free",nrow =8)+ 
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  scale_color_manual(breaks = c("AMTF","AP300","0"), values=c("darkred", "purple","black"))+
  geom_hline(yintercept = 0,color="skyblue",size=0.4,linetype="dashed")+
  geom_text(
    data    = dat_text_A,
    mapping = aes(x = -Inf, y = -Inf, label = label),
    hjust   = -2,
    vjust   = -0.25) 
A<-ggplot(filter_norm_data_scaled,aes(x=order_log_A_mean_ratio,y=log_A_ratio,color=factor(A_active)))+
  geom_point(size=2)+
  geom_hline(yintercept=FDR_A_0.95,color="red",linetype="dashed")+theme_bw()+
  scale_color_manual(breaks = c("1", "0"), values=c("red", "gray"))+
  theme_classic()+geom_text(aes(0,round(FDR_A_0.95,3),label = round(FDR_A_0.95,3), vjust = -1))

B<-ggplot(filter_norm_data_scaled_sub_merged_A, aes(x = order_log_A_mean_ratio, y = factor(marker,levels=c("NC","ATF1_2","ATF3_4",
                                                                                                           "AP300","AMTF"))))+
  
  geom_bin2d(aes(alpha = ..density..), fill = "slategray",binwidth = c(5,1)) +
  scale_fill_continuous(guide = "none")+
  theme_bw()
A/B
C<-ggplot(filter_norm_data_scaled,aes(x=order_log_V_mean_ratio,y=log_V_ratio,color=factor(V_active)))+
  geom_point(size=2)+
  geom_hline(yintercept=FDR_V_0.95,color="red",linetype="dashed")+theme_bw()+
  scale_color_manual(breaks = c("1", "0"), values=c("red", "gray"))+
  theme_classic()+geom_text(aes(0,round(FDR_V_0.95,3),label = round(FDR_V_0.95,3), vjust = -1))

D<-ggplot(filter_norm_data_scaled_sub_merged_V, aes(x = order_log_V_mean_ratio, y = factor(marker,levels=c("NC",
                                                                                                           "VTF1_2","VTF3_4","VP300","VMTF"))))+
  
  geom_bin2d(aes(alpha = ..density..), fill = "slategray",binwidth = c(10,1)) +
  scale_fill_continuous(guide = "none")+
  theme_bw()
C/D
dev.off()
A/p_A/B
##############################################
pdf("new_MTF_single_AV_final_boxplot.pdf")
ggplot(filter_norm_data_scaled_sub_merged_A, aes(y = log_A_ratio, x = factor(marker, levels = c("NC","ATF1_2","ATF3_4",
                                                                                                "AP300","AMTF")),group=marker,
                                  color=marker))+
  geom_boxplot(notch = TRUE)+
  #geom_hline(yintercept = 0,color="skyblue",size=0.4,linetype="dashed")+
  scale_color_manual(breaks = c("NC","ATF1_2","ATF3_4", "AP300","AMTF"), values=c("gray","red", "red","red","red"))+
  theme_classic()
ggplot(filter_norm_data_scaled_sub_merged_A, aes(y = log_A_ratio, x = factor(marker, levels = c("NC","ATF1_2","ATF3_4",
                                                                                                "AP300","AMTF")),group=marker,
                                                 fill=marker))+
  geom_boxplot(notch = TRUE)+
  #geom_hline(yintercept = 0,color="skyblue",size=0.4,linetype="dashed")+
  scale_fill_manual(breaks = c("NC","ATF1_2","ATF3_4", "AP300","AMTF"), values=c("gray","red", "red","red","red"))+
  theme_classic()

ggplot(filter_norm_data_scaled_sub_merged_V, aes(y = log_V_ratio, x = factor(marker, levels = c("NC","VTF1_2","VTF3_4",
                                                                                                "VP300","VMTF")),group=marker,
                                                 color=marker))+
  geom_boxplot(notch = TRUE)+
  #geom_hline(yintercept = 0,color="skyblue",size=0.4,linetype="dashed")+
  scale_color_manual(breaks = c("NC","VTF1_2","VTF3_4", "VP300","VMTF"), values=c("gray","blue", "blue", "blue","blue"))+
  theme_classic()
ggplot(filter_norm_data_scaled_sub_merged_V, aes(y = log_V_ratio, x = factor(marker, levels = c("NC","VTF1_2","VTF3_4",
                                                                                                "VP300","VMTF")),group=marker,
                                                 fill=marker))+
  geom_boxplot(notch = TRUE)+
  #geom_hline(yintercept = 0,color="skyblue",size=0.4,linetype="dashed")+
  scale_fill_manual(breaks = c("NC","VTF1_2","VTF3_4", "VP300","VMTF"), values=c("gray", "blue","blue", "blue","blue"))+
  theme_classic()

dev.off()

########################################################################################

write.table(filter_norm_data_scaled,file="MPRA_7/A_V_mouse_MPRA_7_UMI_scaled__active_version1.txt",sep="\t")
write.table(filter_norm_data_scaled,file="MPRA_7/A_V_mouse_MPRA_7_UMI_scaled__active_version1_annotated.txt",sep="\t")
#write.table(filter_norm_data_scaled,file="MPRA_7/A_V_mouse_MPRA_7_UMI_scaled__active_version1_normalized_by_BNA.txt",sep="\t")
######################################################################################3

##############################A plot ###############################################33

filter_norm_data_scaled_grouped <-filter_norm_data_scaled

Enrich_score <- function(datasets, groups_col , value_col , order_col){
  ALL <- datasets[,value_col]
  Max <- length(ALL)
  datasets$Enrich_score <-0
  datasets$group_order <-0
  pvalue <- c()
  groups_names <-c()
  rows_number <-length(unique(datasets[,groups_col]))
  datasets2 <- datasets[FALSE,]
  datasets3 <- datasets
  for (cols_name in unique(datasets[,groups_col])) {
    datasets3 <- datasets
    PC1 <- datasets3[datasets3[,groups_col]==cols_name,]
    p_PC1 <- perm.test(PC1[,value_col],ALL)$p.value
    PC1$group_order <-0
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
        left_value <- -1
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
    mean_value_p <- mean(datasets3$Enrich_score[datasets3$Enrich_score >=0])
    mean_value_n <- mean(datasets3$Enrich_score[datasets3$Enrich_score <0])
    datasets3$NES[datasets3$Enrich_score >=0] <- 
      datasets3$Enrich_score[datasets3$Enrich_score >=0]/mean_value_p
    datasets3$NES[datasets3$Enrich_score <0] <- 
      datasets3$Enrich_score[datasets3$Enrich_score <0]/mean_value_n
    datasets3$negative_Enrich_score <- (Max_pc1-datasets3[,"group_order"]+1)/Max_pc1 - (Max-datasets3[,order_col]+1)/Max
    mean_value_negative_p <- mean(datasets3$negative_Enrich_score[datasets3$negative_Enrich_score >=0])
    mean_value_negative_n <- mean(datasets3$negative_Enrich_score[datasets3$negative_Enrich_score <0])
    #datasets3$negative_NES[datasets3$negative_Enrich_score >=0] <-
    #  -1*datasets3$negative_Enrich_score[datasets3$negative_Enrich_score >=0]/mean_value_negative_p
    #datasets3$negative_NES[datasets3$negative_Enrich_score <0] <-
    #  datasets3$negative_Enrich_score[datasets3$negative_Enrich_score <0]/mean_value_negative_n
    pvalue <- c(pvalue, p_PC1)
    groups_names <- c(groups_names,cols_name)
    datasets3$score_group <- cols_name
    datasets3$final_ES <- datasets3$Enrich_score - datasets3$negative_Enrich_score
    #if (max(abs(datasets3$negative_Enrich_score)) > max(abs(datasets3$Enrich_score))){
    #  datasets3$final_ES <-
    #  datasets3$negative_NES
    #}
    #else {
    #  datasets3$final_ES <-datasets3$NES
    #}
    datasets2 <- rbind(datasets2,datasets3)
    
  }  
  dat_text <- data.frame(
    label = pvalue,
    score_group = groups_names
  )
  
  newlist <-list("dat_text"=dat_text, "datasets"=datasets2)
  return(newlist)
}
################################3A###########################
#filter_norm_data_scaled_grouped_test_A <-Enrich_score(filter_norm_data_scaled_grouped, "group", "log_A_ratio", "order_log_A_mean_ratio")
filter_norm_data_scaled_grouped_test_A <-Enrich_score(filter_norm_data_scaled_grouped, "AMTF", "log_A_ratio", "order_log_A_mean_ratio")
dat_text_A <-filter_norm_data_scaled_grouped_test_A$dat_text
filter_norm_data_for_plot_A <- filter_norm_data_scaled_grouped_test_A$datasets
pdf("A_MTF_2_ES.pdf")
p_A<-ggplot(filter_norm_data_for_plot_A, aes(order_log_A_mean_ratio, final_ES,color=factor(score_group)))+geom_line()+theme_bw()+theme_classic()+
  facet_wrap(~score_group,scales = "free",nrow =8)+ theme(panel.grid.major.y = element_blank(),
                                    panel.grid.minor.y = element_blank(),
                                    panel.grid.major.x = element_blank(),
                                    panel.grid.minor.x = element_blank())+
  geom_hline(yintercept = 0,color="skyblue",size=0.4,linetype="dashed")
p_A+geom_text(
  data    = dat_text_A,
  mapping = aes(x = -Inf, y = -Inf, label = label),
  hjust   = -2,
  vjust   = -0.25) 
#########################################################3
filter_norm_data_scaled_grouped_test_V <-Enrich_score(filter_norm_data_scaled_grouped, "MTF", "log_V_ratio", "order_log_V_mean_ratio")
dat_text_V <-filter_norm_data_scaled_grouped_test_V$dat_text
filter_norm_data_for_plot_V <- filter_norm_data_scaled_grouped_test_V$datasets

p_V<-ggplot(filter_norm_data_for_plot_V, aes(order_log_V_mean_ratio, final_ES,color=factor(score_group)))+geom_line()+theme_bw()+theme_classic()+
  facet_wrap(~score_group,scales = "free",nrow =8)+ theme(panel.grid.major.y = element_blank(),
                                                          panel.grid.minor.y = element_blank(),
                                                          panel.grid.major.x = element_blank(),
                                                          panel.grid.minor.x = element_blank())+
  geom_hline(yintercept = 0,color="skyblue",size=0.4,linetype="dashed")
p_V+geom_text(
  data    = dat_text_V,
  mapping = aes(x = -Inf, y = -Inf, label = label),
  hjust   = -2,
  vjust   = -0.25) 
dev.off()
#######################################normlized to DNA counts################
fit <- lm(log2(A_ratio+1) ~ log2(DNA+1), data = filter_norm_data_scaled)
filter_norm_data_scaled$A_residuals <- residuals(fit)

fit <- lm(log2(V_ratio+1) ~ log(DNA), data = filter_norm_data_scaled)
filter_norm_data_scaled$V_residuals <- residuals(fit)

melt_filter_norm_data <- melt(filter_norm_data_scaled[,c(3,8,9)])
ggplot(melt_filter_norm_data, aes(x = variable, y=log2(value+1), fill = Group))+geom_boxplot()+theme_bw()+theme_classic()
ggplot(melt_filter_norm_data, aes(x = Group, y=log2(value+1), fill = variable))+geom_boxplot()+theme_bw()+theme_classic()+theme(axis.text=element_text(size=7),
                                                                                                                                axis.title=element_text(size=7,face="bold"))

melt_filter_norm_data2 <- melt(filter_norm_data_scaled[,c(3,10,11)])
ggplot(melt_filter_norm_data2, aes(x = variable, y=value, fill = Group))+geom_boxplot()+theme_bw()+theme_classic()
ggplot(melt_filter_norm_data2, aes(x = Group, y=value, fill = variable))+geom_boxplot()+theme_bw()+theme_classic()+theme(axis.text=element_text(size=7),
                                                                                                                         axis.title=element_text(size=7,face="bold"))


melt_filter_norm_data2 <- melt(filter_norm_data[,c(7,16,17)])
ggplot(melt_filter_norm_data2, aes(x = variable, y=log2(value+1), fill = Group))+geom_boxplot()+theme_bw()+theme_classic()
ggplot(melt_filter_norm_data2, aes(x = Group, y=log2(value+1), fill = variable))+geom_boxplot()+theme_bw()+theme_classic()+theme(axis.text=element_text(size=7),
                                                                                                                                 axis.title=element_text(size=7,face="bold"))
write.table(filter_norm_data_scaled,file="A_V_mouse_MPRA_normalized_count.txt",sep="\t")

head(melt_filter_norm_data2)
########################################################################
################################diff enhancers####################
#table <- read.table("A_V_mouse_MPRA_normalized_count.txt", header = T, row.names = 1)
#head(table)
#noTSS <- read.table(file = "yangpo_av_mpra_library_information_no_TSS_1kb.bed",header = F)
#head(filter_norm_data_scaled)
#head(noTSS)
#rownames(noTSS) <- noTSS[,4]
#noTSS_regions <- intersect(rownames(noTSS),rownames(filter_norm_data_scaled))
#table_RNA <- filter_norm_data_scaled[noTSS_regions,]
table_RNA <- filter_norm_data_scaled
head(table_RNA)
pvalue <- c()
for (i in 1:dim(table_RNA)[1]) {
  ref_value <- as.numeric(table_RNA[i,c(1:5)])
  alt_value <- as.numeric(table_RNA[i,c(6:10)])
  pvalue_i <- t.test(ref_value, alt_value)
  pvalue <- c(pvalue, pvalue_i$p.value)
}
qvalue <- pvalue * length(pvalue) /rank(pvalue)
inter_join<-table_RNA
inter_join$A <- inter_join$log_A_ratio 
inter_join$V <- inter_join$log_V_ratio
inter_join$qvalue <- as.numeric(as.character(qvalue))
inter_join$pvalue <- as.numeric(as.character(pvalue))
inter_join$fold <- inter_join$A-inter_join$V

inter_join$qvalue_0_05 <- 0
inter_join$qvalue_0_05[inter_join$qvalue < 0.05] <- 1
inter_join$pvalue_0_05 <- 0
inter_join$pvalue_0_05[inter_join$pvalue < 0.05] <- 1
inter_join$fold_1_5 <- 0
inter_join$fold_1_5 [inter_join$fold <= -0.58] <- 2
inter_join$fold_1_5 [inter_join$fold >= 0.58] <- 1
inter_join$active <- 0
#inter_join$active[inter_join$A_active==1] <-1
inter_join$active[inter_join$A_active==1] <-1
inter_join$active[inter_join$V_active==1] <-1
#inter_join$active[inter_join$A_active==1] <-1
#inter_join$active[inter_join$V_active==1] <-1
inter_join$final_sig <-0
inter_join$final_sig[inter_join$qvalue_0_05 >0 & inter_join$fold >= 0.58 & inter_join$active >0] <-2
inter_join$final_sig[inter_join$qvalue_0_05 >0 &  inter_join$fold <= -0.58 & inter_join$active >0] <-1
inter_join$order_fold <- rank(-as.numeric(as.character(inter_join$fold)))
inter_join$Asp<-0
inter_join[inter_join$final_sig==2,]$Asp<-1
inter_join$Vsp<-0
inter_join[inter_join$final_sig==1,]$Vsp<-1
write.table(inter_join,file="MPRA_7/A_V_mouse_MPRA_7_UMI_normalized_count_fd_version_add_1_new_annotation.txt",sep="\t",quote = F)
#write.table(inter_join,file="MPRA_7/A_V_mouse_MPRA_7_UMI_normalized_count_fd_version_add_1_norm_by_BNA.txt",sep="\t")
inter_join<-read.table(file="MPRA_7/A_V_mouse_MPRA_7_UMI_normalized_count_fd_version_add_1_new_annotation.txt",header = T, row.names = 1,sep = "\t")
##################################3
mutagenesis_list<-read.table(file="mutagenesis_42_regions.txt",header=F)
over_with_mutagenesis_list <-inter_join[inter_join$region_locs %in% mutagenesis_list[,1],]
over_with_mutagenesis_list
write.table(over_with_mutagenesis_list, "mutagenesis_over_max_5.txt",quote=F,sep="\t")
#######################################################################################
###########################################real figure
pdf("A_V_400bp_DEE.pdf")
ggplot(inter_join,aes(A, V,color=factor(final_sig,levels=c('0','1','2')),alpha =factor(qvalue_0_05,levels=c("0","1"))))+
  geom_point()+theme_bw()+theme_classic()+geom_abline(slope = 1,color="black",size=0.2,linetype="dashed")+
  geom_abline(slope = 1,intercept=0.58,color="black",size=0.2,linetype="dashed")+
  geom_abline(slope=1,intercept=-0.58,color="black",size=0.2,linetype="dashed")+
  scale_color_manual(breaks = c("0", "1", "2"), values=c("gray", "blue", "red"))+
  scale_alpha_manual(breaks = c("0", "1"), values=c(0.2, 1))+
  geom_vline(xintercept = FDR_A_0.95,color="skyblue",size=0.4,linetype="dashed")+
  geom_hline(yintercept = FDR_V_0.95,color="skyblue",size=0.4,linetype="dashed")+
  xlim(-6,9)+ylim(-6,9)
dev.off()
#############################################################################
number <- as.data.frame(table(inter_join$final_sig))
dat_text_2 <- data.frame(
  label = c("nochange","V_specific_enhancer","A_specific_enhancers"),
  score_group = number$Freq
)
write.table(dat_text_2,file = "Diff_enhancers_number_mutagenesis.txt",quote=F, sep = "\t")
###########################################################################333

head(inter_join)
table(inter_join$active)
sum(inter_join$sig==1)
sum(inter_join$sig==2)
sum(inter_join$final_sig==1)
sum(inter_join$final_sig==2)
#inter_join <- as.data.frame(inter_join)
############################ fold figure
#####################################
#####################################
inter_join_sub8<-inter_join[((inter_join$MTF=="11" & inter_join$P300=="00")|(inter_join$MTF=="11"& inter_join$P300=="11") |
                              (inter_join$MTF=="00" & inter_join$P300=="11")),] 
inter_join_sub8<-inter_join[((inter_join$MTF=="11" & inter_join$P300=="0")|(inter_join$MTF=="11"& inter_join$P300=="11") |
                               (inter_join$MTF=="0" & inter_join$P300=="11")),] 
inter_join_sub8$marker<-"PC"
inter_join_sub8$marker2<-"PC"
inter_join_sub1<-inter_join[(inter_join$ATF1_2==1 & inter_join$VTF<1),]
inter_join_sub1$marker<-"AspTF1_2"
inter_join_sub1$marker2<-"A"
inter_join_sub5<-inter_join[(inter_join$VTF1_2==1 & inter_join$ATF<1),]
inter_join_sub5$marker<-"VspTF1_2"
inter_join_sub5$marker2<-"V"
inter_join_sub6<-inter_join[(inter_join$ATF3_4==1 & inter_join$VTF<3),]
inter_join_sub6$marker<-"AspTF3_4"
inter_join_sub6$marker2<-"A"
inter_join_sub7<-inter_join[(inter_join$VTF3_4==1 & inter_join$ATF<3),]
inter_join_sub7$marker<-"VspTF3_4"
inter_join_sub7$marker2<-"V"

inter_join_sub2<-inter_join[inter_join$VspMTF==1,] 
inter_join_sub2$marker<-"V-sMTF"
inter_join_sub4<-inter_join[inter_join$AspMTF==1,]
inter_join_sub4$marker <- "A-sMTF"
inter_join_sub2$marker2<-"V"
inter_join_sub4$marker2 <- "A"
inter_join_sub3<-inter_join[(inter_join$AspP300==1 | inter_join$VspP300==1),] 
inter_join_sub3$marker<-"V-sP300"
inter_join_sub3[inter_join_sub3$AspP300==1,]$marker <- "A-sP300"
inter_join_sub3$marker2<-"V"
inter_join_sub3[inter_join_sub3$AspP300==1,]$marker2 <- "A"
#inter_join_sub4<-inter_join[inter_join$Gene_3=="10" |inter_join$Gene_3=="01",] 
#inter_join_sub4$marker<-"V-sGene"
#inter_join_sub4[inter_join_sub4$Gene_3=="10",]$marker<-"A-sGene"
#inter_join_sub4$marker2<-"V"
#inter_join_sub4[inter_join_sub4$Gene_3=="10",]$marker2<-"A"
inter_join_sub_merged <- rbind(#inter_join_sub1[,c("order_fold","fold","marker","marker2")],
                               inter_join_sub4[,c("order_fold","fold","marker","marker2")],
                               inter_join_sub2[,c("order_fold","fold","marker","marker2")],
                               inter_join_sub3[,c("order_fold","fold","marker","marker2")],
                               inter_join_sub4[,c("order_fold","fold","marker","marker2")],
                               #inter_join_sub5[,c("order_fold","fold","marker","marker2")],
                               #inter_join_sub6[,c("order_fold","fold","marker","marker2")],
                               #inter_join_sub7[,c("order_fold","fold","marker","marker2")],
                               inter_join_sub8[,c("order_fold","fold","marker","marker2")])
                               
                               #inter_join_sub4[,c("order_fold","fold","marker","marker2")])

library(patchwork)
pdf("fold_AV_rank_figures_new.pdf", width =12,height = 8)
A<-ggplot(inter_join,aes(x=order_fold,y=fold,color=factor(final_sig),size=factor(active)))+
  geom_point()+
  geom_hline(yintercept=0,color="black",linetype="dashed")+theme_bw()+
  scale_color_manual(breaks = c("2","1","0"), values=c("red", "blue","gray"))+
  scale_size_manual(breaks = c("1","0"), values=c(2,0.01))+
  theme_classic()

B<-ggplot(inter_join_sub_merged, aes(x = order_fold, y = factor(marker,levels=c("PC","VspTF1_2","AspTF1_2","VspTF3_4","AspTF3_4",
                                                                                "V-sP300","A-sP300",
                                                                                "V-sMTF","A-sMTF"))))+
  geom_bin2d(aes(alpha = ..density..), fill = "slategray",binwidth = c(3,1)) +
  scale_fill_continuous(guide = "none")+
  theme_bw()
  #stat_bin2d(aes(fill = after_stat(count)), binwidth = c(1,1),color="Black")+theme_bw()
A/B
A<-ggplot(inter_join,aes(x=order_fold,y=fold,color=factor(final_sig),size=factor(active)))+
  geom_point()+
  geom_hline(yintercept=0,color="black",linetype="dashed")+theme_bw()+
  scale_color_manual(breaks = c("2","1","0"), values=c("red", "blue","gray"))+
  scale_size_manual(breaks = c("1","0"), values=c(2,0.01))+
  theme_classic()

B<-ggplot(inter_join_sub_merged, aes(x = order_fold, y = factor(marker,levels=c("PC",
                                                                                "V-sP300","A-sP300",
                                                                                "V-sMTF","A-sMTF"))))+
  geom_bin2d(aes(alpha = ..density..), fill = "slategray",binwidth = c(3,1)) +
  scale_fill_continuous(guide = "none")+
  theme_bw()
  #stat_bin2d(aes(fill = after_stat(count)), binwidth = c(3,1))+theme_bw()
A/B
dev.off()
pdf("fold_AV_boxplot_new.pdf")
ggplot(inter_join_sub_merged, aes(y = fold, x = factor(marker, levels = c("PC","VspTF1_2","AspTF1_2","VspTF3_4","AspTF3_4",
                                                                          "V-sP300","A-sP300",
                                                                          "V-sMTF","A-sMTF")),group=marker,
                                  color=factor(marker2, levels = c("NC","A","V"))))+
  geom_boxplot(notch = TRUE)+
  geom_hline(yintercept = 0,color="skyblue",size=0.4,linetype="dashed")+
  scale_color_manual(breaks = c("NC","A","V"), values=c("gray","red", "blue"))+
  theme_classic()
ggplot(inter_join_sub_merged, aes(y = fold, x = factor(marker, levels = c("PC","VspTF1_2","AspTF1_2","VspTF3_4","AspTF3_4",
                                                                          "V-sP300","A-sP300",
                                                                          "V-sMTF","A-sMTF")),group=marker,
                                  fill=factor(marker2, levels = c("PC","A","V"))))+
  geom_boxplot(notch = TRUE)+
  geom_hline(yintercept = 0,color="black",size=0.4,linetype="dashed")+
  scale_fill_manual(breaks = c("NC","A","V"), values=c("gray","red", "blue"))+
  theme_classic()
dev.off()
inter_join_sub<-inter_join[,c("A1_vs_DNA_log","A2_vs_DNA_log","A3_vs_DNA_log","A4_vs_DNA_log",
                                                        "A5_vs_DNA_log", "V1_vs_DNA_log"  ,   "V2_vs_DNA_log" ,   "V3_vs_DNA_log",
                                                        "V4_vs_DNA_log" ,  "V5_vs_DNA_log",   "log_A_ratio" ,  "log_V_ratio",
                                                        "order_log_A_mean_ratio" ,"order_log_V_mean_ratio", "Group6", "MTF","spMTF","ATF","VTF","AMTF","VMTF","AspMTF","VspMTF",
                                                        "P300","spP300","AP300","VP300","AspP300","VspP300","ATF1-2","VTF1-2","ATF3-4","VTF3-4",
                              "A_active",   "V_active",   "rawname",  "group",
                                                        "A", "V",  "qvalue", "pvalue","fold","active","final_sig","order_fold")]
write.table(inter_join_sub, file="A-V_400bp_MPRA_final_group_annotation.txt",sep = "\t",quote = F)
###################################################################################3

#inter_join <-as.data.frame(as.character(inter_join))
#inter_join_with_bed <- cbind(data[rownames(inter_join), 7:9],inter_join)
#write.table(inter_join_with_bed,file="MPRA_7/A_V_mouse_MPRA_normalized_count_fd_with_bed.txt",sep="\t",quote = F)
#write.table(inter_join_with_bed,file="MPRA_7/A_V_mouse_MPRA_normalized_count_fd_with_bed_norm_BNA.txt",sep="\t",quote = F)


################################3A###########################
filter_norm_data_scaled_grouped_test_f <-Enrich_score(inter_join, "P300", "fold", "order_fold")
dat_text_f <-filter_norm_data_scaled_grouped_test_f$dat_text
filter_norm_data_for_plot_f<- filter_norm_data_scaled_grouped_test_f$datasets
pdf("fold_AV_MTF_2_ES.pdf")
pdf("fold_AV_P300_ES.pdf")
p_f<-ggplot(filter_norm_data_for_plot_f, aes(order_fold, final_ES,color=factor(score_group)))+geom_line()+theme_bw()+theme_classic()+
  facet_wrap(~score_group,scales = "free",nrow =8)+ theme(panel.grid.major.y = element_blank(),
                                                          panel.grid.minor.y = element_blank(),
                                                          panel.grid.major.x = element_blank(),
                                                          panel.grid.minor.x = element_blank())+
  geom_hline(yintercept = 0,color="skyblue",size=0.4,linetype="dashed")
p_f+geom_text(
  data    = dat_text_f,
  mapping = aes(x = -Inf, y = -Inf, label = label),
  hjust   = -2,
  vjust   = -0.25) 
dev.off()


######################################not need
ggplot(inter_join,aes(x=order_fold,y=log2(fold),color = factor(final_sig,levels=c("0","1"))))+geom_point(size=2)+
  geom_hline(yintercept=log2(1.5),color="red",linetype="dashed")+
  geom_hline(yintercept=log2(0.667),color="red",linetype="dashed")+
  theme_bw()+theme_classic()+scale_color_manual(breaks = c("0", "1"), values=c("gray", "blue"))
inter_join <-as.data.frame(as.character(inter_join))
ggplot(inter_join, aes(x = order_fold, y = formatC(group, width = 6, format = "d", flag = "0")))+stat_bin2d(aes(fill = after_stat(count)), binwidth = c(3,1))
       
sub <-inter_join[inter_join$pvalue_0_05==1,1:5]      
pheatmap(sub[2:5],cluster_rows=T,cluster_cols=F,scale="row",color = colorRampPalette(c("navy", "white","gold"))(50))
ggplot(inter_join[inter_join$pvalue_0_05==1,], aes(x = order_fold, y = formatC(Group6, width = 6, format = "d", flag = "0"),color = factor(final_sig_2,levels=c("1","2"))))+
stat_bin2d(aes(fill = after_stat(count)), binwidth = c(3,1))+scale_color_manual(breaks = c("1", "2"), values=c("blue", "red")) 
       
ggplot(inter_join[inter_join$final_sig_2 >0,],aes(x=order_fold,y=log2(fold),color = factor(final_sig_2,levels=c("1","2"))))+geom_point(size=2)+
  geom_hline(yintercept=log2(1.5),color="red",linetype="dashed")+
  geom_hline(yintercept=log2(0.667),color="red",linetype="dashed")+
  theme_bw()+theme_classic()+scale_color_manual(breaks = c("1", "2"), values=c("blue", "red"))       
#inter_join_with_bed <- cbind(data[rownames(inter_join), 7:9],inter_join)
#write.table(inter_join_with_bed,file="A_V_mouse_MPRA_normalized_count_fd_with_bed.txt",sep="\t",quote = F)
write.table(inter_join,file="A_V_mouse_MPRA_normalized_count_fd_with_bed_new.txt",sep="\t",quote = F)
###################################################################
################################################################################
inter_join<- read.table("A_V_mouse_MPRA_normalized_count_fd_with_bed_new.txt",header=T,row.names = 1)
inter_join_2<-inter_join
dim(inter_join)
#inter_join_2$Aspf<-0
#inter_join_2[inter_join_2$fold>=0.58,]$Asp<-1
#inter_join_2$Vspf<-0
#inter_join_2[inter_join_2$fold<=-0.58,]$Vsp<-1
#inter_join_2$AspMTF <-0
#inter_join_2[inter_join_2$ATF>=5 & inter_join_2$VTF<=3,]$AspMTF <-1
#inter_join_2$VspMTF <-0
#inter_join_2[inter_join_2$VTF>=5 & inter_join_2$ATF<=3,]$VspMTF <-1
model_filter_norm_data_scaled_2 <-inter_join_2[,c("AspMTF","VspMTF","Asp","Vsp","AspP300",
                                                "VspP300","AP300","VP300","ATF","VTF")]

ATAC_anntation<-read.table(file = "MPRA_7/new_annotation_sensitive/MPRA_with_A-V_ATAC_As_Vs_ATAC_full_list.txt")
loops_anntation<-read.table(file = "MPRA_7/new_annotation_sensitive/MPRA_with_A-V_loops_As_Vs_loops_full_list.txt")
annotate_files_2 <- read.table(file="MPRA_7/new_annotation_sensitive/MPRA_list_TFs_and_ATAC_new.txt",header=T, row.names = 1)
annotate_files_3 <- read.table(file="MPRA_7/new_annotation_sensitive/MPRA_hichp_A-V_counts_distance.txt",header=F, row.names = 1)

dim(ATAC_anntation)
dim(model_filter_norm_data_scaled_2)
model_filter_norm_data_scaled_2$ATAC_Asp <-ATAC_anntation[ATAC_anntation$V1 %in% rownames(model_filter_norm_data_scaled_2),"V4"]
model_filter_norm_data_scaled_2$ATAC_Vsp <-ATAC_anntation[ATAC_anntation$V1 %in% rownames(model_filter_norm_data_scaled_2),"V5"]
model_filter_norm_data_scaled_2$loops_Asp <-loops_anntation[loops_anntation$V1 %in% rownames(model_filter_norm_data_scaled_2),"V4"]
model_filter_norm_data_scaled_2$loops_Vsp <-loops_anntation[loops_anntation$V1 %in% rownames(model_filter_norm_data_scaled_2),"V5"]

model_filter_norm_data_scaled_2$count_A <-annotate_files_3[rownames(model_filter_norm_data_scaled_2),"V2"]
model_filter_norm_data_scaled_2$count_V <-annotate_files_3[rownames(model_filter_norm_data_scaled_2),"V3"]
model_filter_norm_data_scaled_2$dis_A <-annotate_files_3[rownames(model_filter_norm_data_scaled_2),"V4"]
model_filter_norm_data_scaled_2$dis_V <-annotate_files_3[rownames(model_filter_norm_data_scaled_2),"V5"]
model_filter_norm_data_scaled_2$count_A_2<-0
model_filter_norm_data_scaled_2[(model_filter_norm_data_scaled_2$count_A>=5 & (model_filter_norm_data_scaled_2$count_A+1)/(model_filter_norm_data_scaled_2$count_V+1)>=1.5),]$count_A_2 <- 1
model_filter_norm_data_scaled_2$count_V_2<-0
model_filter_norm_data_scaled_2[(model_filter_norm_data_scaled_2$count_V>=5 & (model_filter_norm_data_scaled_2$count_V+1)/(model_filter_norm_data_scaled_2$count_A+1)>=1.5),]$count_V_2 <- 1
model_filter_norm_data_scaled_2$dis_V_2<-0
model_filter_norm_data_scaled_2[(model_filter_norm_data_scaled_2$count_V>=5 & (model_filter_norm_data_scaled_2$dis_V+1)-(model_filter_norm_data_scaled_2$dis_A+1)>=10000),]$dis_V_2 <- 1
model_filter_norm_data_scaled_2$dis_A_2<-0
model_filter_norm_data_scaled_2[(model_filter_norm_data_scaled_2$count_A>=5 & (model_filter_norm_data_scaled_2$dis_A+1)-(model_filter_norm_data_scaled_2$dis_V+1)>=10000),]$dis_A_2 <- 1

#model_filter_norm_data_scaled_2$dis_V_2<-0
#model_filter_norm_data_scaled_2[(model_filter_norm_data_scaled_2$count_V>=5 & model_filter_norm_data_scaled_2$dis_V>=100000),]$dis_V_2 <- 1
#model_filter_norm_data_scaled_2$dis_A_2<-0
#model_filter_norm_data_scaled_2[(model_filter_norm_data_scaled_2$count_A>=5 & model_filter_norm_data_scaled_2$dis_A>=100000),]$dis_A_2 <- 1

model_filter_norm_data_scaled_2$ATAC_A <-annotate_files_2[rownames(model_filter_norm_data_scaled_2),"ATAC_A"]
model_filter_norm_data_scaled_2$ATAC_V <-annotate_files_2[rownames(model_filter_norm_data_scaled_2),"ATAC_V"]
model_filter_norm_data_scaled_2$P300_A <-annotate_files_2[rownames(model_filter_norm_data_scaled_2),"P300_A"]
model_filter_norm_data_scaled_2$P300_V <-annotate_files_2[rownames(model_filter_norm_data_scaled_2),"P300_V"]
model_filter_norm_data_scaled_2$MTF_A <-annotate_files_2[rownames(model_filter_norm_data_scaled_2),"MTF_A"]
model_filter_norm_data_scaled_2$MTF_V <-annotate_files_2[rownames(model_filter_norm_data_scaled_2),"MTF_V"]
model_filter_norm_data_scaled_2$P300_2_A<-0
model_filter_norm_data_scaled_2[model_filter_norm_data_scaled_2$AP300==1 &( model_filter_norm_data_scaled_2$P300_A/model_filter_norm_data_scaled_2$P300_V>=1.5),]$P300_2_A<-1
model_filter_norm_data_scaled_2$P300_2_V<-0
model_filter_norm_data_scaled_2[model_filter_norm_data_scaled_2$VP300==1 &( model_filter_norm_data_scaled_2$P300_V/model_filter_norm_data_scaled_2$P300_A>=1.5),]$P300_2_V<-1
model_filter_norm_data_scaled_2$MTF_2_A<-0
model_filter_norm_data_scaled_2[model_filter_norm_data_scaled_2$ATF>=4 & model_filter_norm_data_scaled_2$VTF<=3,]$MTF_2_A<-1
model_filter_norm_data_scaled_2$MTF_2_V<-0
model_filter_norm_data_scaled_2[model_filter_norm_data_scaled_2$VTF>=4 & model_filter_norm_data_scaled_2$ATF<=3,]$MTF_2_V<-1
model_filter_norm_data_scaled_2$MTF_3_A<-0
model_filter_norm_data_scaled_2[(model_filter_norm_data_scaled_2$MTF_A+1)/(model_filter_norm_data_scaled_2$MTF_V+1)>=1.5,]$MTF_3_A<-1
model_filter_norm_data_scaled_2$MTF_3_V<-0
model_filter_norm_data_scaled_2[(model_filter_norm_data_scaled_2$MTF_V+1)/(model_filter_norm_data_scaled_2$MTF_A+1)>=1.5,]$MTF_3_V<-1
model_filter_norm_data_scaled_2$ATAC_2_A<-0
model_filter_norm_data_scaled_2[(model_filter_norm_data_scaled_2$ATAC_A+1)/(model_filter_norm_data_scaled_2$ATAC_V+1)>=1.5,]$ATAC_2_A<-1
model_filter_norm_data_scaled_2$ATAC_2_V<-0
model_filter_norm_data_scaled_2[(model_filter_norm_data_scaled_2$ATAC_V+1)/(model_filter_norm_data_scaled_2$ATAC_A+1)>=1.5,]$ATAC_2_V<-1
model_filter_norm_data_scaled_2$merged_V<-0
model_filter_norm_data_scaled_2$merged_A<-0
model_filter_norm_data_scaled_2[(model_filter_norm_data_scaled_2$AspMTF ==1 |model_filter_norm_data_scaled_2$AspP300==1|
                                   model_filter_norm_data_scaled_2$ATAC_Asp ==1 |model_filter_norm_data_scaled_2$loops_Asp==1) ,]$merged_A<-1
model_filter_norm_data_scaled_2[(model_filter_norm_data_scaled_2$VspMTF ==1 |model_filter_norm_data_scaled_2$VspP300==1|
                                   model_filter_norm_data_scaled_2$ATAC_Vsp ==1 |model_filter_norm_data_scaled_2$loops_Vsp==1) ,]$merged_V<-1
SSA_matrix_merge_A_2 <- SSA_merge(model_filter_norm_data_scaled_2,c("AspMTF","AspP300",
                                                                "ATAC_Asp","loops_Asp","merged_A"),c("Asp"))
SSA_matrix_merge_V_2 <- SSA_merge(model_filter_norm_data_scaled_2,c("VspMTF","VspP300","ATAC_Vsp",
                                                                    "loops_Vsp","merged_V"),
                                c("Vsp"))


#SSA_matrix_merge_A_2 <- SSA_merge(model_filter_norm_data_scaled_2,c("AspMTF","AspP300",
#                                                                    "AspGene","ATAC_Asp","loops_Asp",
#                                                                    "MTF_2_A","P300_2_A","ATAC_2_A","MTF_3_A","count_A_2","dis_A_2"),c("Asp"))
#SSA_matrix_merge_V_2 <- SSA_merge(model_filter_norm_data_scaled_2,c("VspMTF","VspP300","ATAC_Vsp","loops_Vsp",
#                                                                    "MTF_2_V","P300_2_V","ATAC_2_V","MTF_3_V","count_V_2","dis_V_2"),
#                                  c("Vsp"))


SSA_matrix_merged_2 <- rbind(SSA_matrix_merge_A_2,SSA_matrix_merge_V_2)
SSA_matrix_merge_A_matrix_2 <- SSA_merge_matrix(model_filter_norm_data_scaled_2,c("AspMTF","AspP300",
                                                                              "ATAC_Asp","loops_Asp",
                                                                              "merged_A"),
                                                c("Asp"))
SSA_matrix_merge_V_matrix_2 <- SSA_merge_matrix(model_filter_norm_data_scaled_2,c("VspMTF","VspP300","ATAC_Vsp","loops_Vsp",
                                                            "merged_V"),
                                              c("Vsp"))
write.table(SSA_matrix_merge_A_matrix_2, file="Asp_confmatrix.txt",sep = "\t",quote = F)
write.table(SSA_matrix_merge_V_matrix_2, file="Vsp_confmatrix.txt",sep = "\t",quote = F)

pdf("SSA_spMTF_spP300_spgene_spATAC_sploops.pdf",width = 7,height = 6)
A<-ggplot(SSA_matrix_merge_A_2,aes(y=factor(Factor,levels=rev(c("AspMTF","AspP300",
                                                              "ATAC_Asp","loops_Asp",
                                                              "merged_A"))),
                                 x=as.numeric(value),group=label,fill=factor(label)))+
  geom_bar(stat="identity")+ facet_wrap(vars(factor(label,levels=c("Specificity","Sensitivity","Accuracy"))),ncol = 3)+
  theme_classic()
V<-ggplot(SSA_matrix_merge_V_2,aes(y=factor(Factor,levels =rev(c("VspMTF","VspP300","ATAC_Vsp","loops_Vsp",
                                                                 "merged_V"))),
                                 x=as.numeric(value),group=label,fill=factor(label)))+
  geom_bar(stat="identity")+ facet_wrap(vars(factor(label,levels=c("Specificity","Sensitivity","Accuracy"))),ncol = 3)+
  theme_classic()

A/V
dev.off()
#################################################################
model_filter_norm_data_scaled_2$Gata4A <-model_filter_norm_data_scaled[rownames(model_filter_norm_data_scaled_2),"Gata4_A"]
model_filter_norm_data_scaled_2$Gata4V <-model_filter_norm_data_scaled[rownames(model_filter_norm_data_scaled_2),"Gata4_V"]
model_filter_norm_data_scaled_2$Mef2aA <-model_filter_norm_data_scaled[rownames(model_filter_norm_data_scaled_2),"Mef2a_A"]
model_filter_norm_data_scaled_2$Mef2aV <-model_filter_norm_data_scaled[rownames(model_filter_norm_data_scaled_2),"Mef2a_V"]
model_filter_norm_data_scaled_2$Mef2cA <-model_filter_norm_data_scaled[rownames(model_filter_norm_data_scaled_2),"Mef2c_A"]
model_filter_norm_data_scaled_2$Mef2cV <-model_filter_norm_data_scaled[rownames(model_filter_norm_data_scaled_2),"Mef2c_V"]
model_filter_norm_data_scaled_2$SrfA <-model_filter_norm_data_scaled[rownames(model_filter_norm_data_scaled_2),"Srf_A"]
model_filter_norm_data_scaled_2$SrfV <-model_filter_norm_data_scaled[rownames(model_filter_norm_data_scaled_2),"Srf_V"]
model_filter_norm_data_scaled_2$Tead1A <-model_filter_norm_data_scaled[rownames(model_filter_norm_data_scaled_2),"Tead1_A"]
model_filter_norm_data_scaled_2$Tead1V <-model_filter_norm_data_scaled[rownames(model_filter_norm_data_scaled_2),"Tead1_V"]
model_filter_norm_data_scaled_2$Nkx2.5A <-model_filter_norm_data_scaled[rownames(model_filter_norm_data_scaled_2),"Nkx2.5_A"]
model_filter_norm_data_scaled_2$Nkx2.5V <-model_filter_norm_data_scaled[rownames(model_filter_norm_data_scaled_2),"Nkx2.5_V"]
model_filter_norm_data_scaled_2$Tbx5A <-model_filter_norm_data_scaled[rownames(model_filter_norm_data_scaled_2),"Tbx5_A"]
model_filter_norm_data_scaled_2$Tbx5V <-model_filter_norm_data_scaled[rownames(model_filter_norm_data_scaled_2),"Tbx5_V"]

model_filter_norm_data_scaled_2$Gata4_A <-annotate_files_2[rownames(model_filter_norm_data_scaled_2),"Gata4_A"]
model_filter_norm_data_scaled_2$Gata4_V <-annotate_files_2[rownames(model_filter_norm_data_scaled_2),"Gata4_V"]
model_filter_norm_data_scaled_2$Mef2a_A <-annotate_files_2[rownames(model_filter_norm_data_scaled_2),"Mef2a_A"]
model_filter_norm_data_scaled_2$Mef2a_V <-annotate_files_2[rownames(model_filter_norm_data_scaled_2),"Mef2a_V"]
model_filter_norm_data_scaled_2$Mef2c_A <-annotate_files_2[rownames(model_filter_norm_data_scaled_2),"Mef2c_A"]
model_filter_norm_data_scaled_2$Mef2c_V <-annotate_files_2[rownames(model_filter_norm_data_scaled_2),"Mef2c_V"]
model_filter_norm_data_scaled_2$Srf_A <-annotate_files_2[rownames(model_filter_norm_data_scaled_2),"Srf_A"]
model_filter_norm_data_scaled_2$Srf_V <-annotate_files_2[rownames(model_filter_norm_data_scaled_2),"Srf_V"]
model_filter_norm_data_scaled_2$Tead1_A <-annotate_files_2[rownames(model_filter_norm_data_scaled_2),"Tead1_A"]
model_filter_norm_data_scaled_2$Tead1_V <-annotate_files_2[rownames(model_filter_norm_data_scaled_2),"Tead1_V"]
model_filter_norm_data_scaled_2$Nkx2.5_A <-annotate_files_2[rownames(model_filter_norm_data_scaled_2),"Nkx2.5_A"]
model_filter_norm_data_scaled_2$Nkx2.5_V <-annotate_files_2[rownames(model_filter_norm_data_scaled_2),"Nkx2.5_V"]
model_filter_norm_data_scaled_2$Tbx5_A <-annotate_files_2[rownames(model_filter_norm_data_scaled_2),"Tbx5_A"]
model_filter_norm_data_scaled_2$Tbx5_V <-annotate_files_2[rownames(model_filter_norm_data_scaled_2),"Tbx5_V"]
model_filter_norm_data_scaled_2$Gata4_2_V<-0
model_filter_norm_data_scaled_2[model_filter_norm_data_scaled_2$Gata4V==1 & (model_filter_norm_data_scaled_2$Gata4_V+1)/(model_filter_norm_data_scaled_2$Gata4_A+1)>=1.5,]$Gata4_2_V<-1
model_filter_norm_data_scaled_2$Gata4_2_A<-0
model_filter_norm_data_scaled_2[model_filter_norm_data_scaled_2$Gata4A==1 & (model_filter_norm_data_scaled_2$Gata4_A+1)/(model_filter_norm_data_scaled_2$Gata4_V+1)>=1.5,]$Gata4_2_A<-1

model_filter_norm_data_scaled_2$Mef2a_2_V<-0
model_filter_norm_data_scaled_2[model_filter_norm_data_scaled_2$Mef2aV==1 & (model_filter_norm_data_scaled_2$Mef2a_V+1)/(model_filter_norm_data_scaled_2$Mef2a_A+1)>=1.5,]$Mef2a_2_V<-1
model_filter_norm_data_scaled_2$Mef2a_2_A<-0
model_filter_norm_data_scaled_2[model_filter_norm_data_scaled_2$Mef2aA ==1 & (model_filter_norm_data_scaled_2$Mef2a_A+1)/(model_filter_norm_data_scaled_2$Mef2a_V+1)>=1.5,]$Mef2a_2_A<-1

model_filter_norm_data_scaled_2$Tead1_2_V<-0
model_filter_norm_data_scaled_2[model_filter_norm_data_scaled_2$Tead1V==1 & (model_filter_norm_data_scaled_2$Tead1_V+1)/(model_filter_norm_data_scaled_2$Tead1_A+1)>=1.5,]$Tead1_2_V<-1
model_filter_norm_data_scaled_2$Tead1_2_A<-0
model_filter_norm_data_scaled_2[model_filter_norm_data_scaled_2$Tead1A==1 & (model_filter_norm_data_scaled_2$Tead1_A+1)/(model_filter_norm_data_scaled_2$Tead1_V+1)>=1.5,]$Tead1_2_A<-1

model_filter_norm_data_scaled_2$Mef2c_2_V<-0
model_filter_norm_data_scaled_2[model_filter_norm_data_scaled_2$Mef2cV==1 & (model_filter_norm_data_scaled_2$Mef2c_V+1)/(model_filter_norm_data_scaled_2$Mef2c_A+1)>=1.5,]$Mef2c_2_V<-1
model_filter_norm_data_scaled_2$Mef2c_2_A<-0
model_filter_norm_data_scaled_2[model_filter_norm_data_scaled_2$Mef2cA==1 & (model_filter_norm_data_scaled_2$Mef2c_A+1)/(model_filter_norm_data_scaled_2$Mef2c_V+1)>=1.5,]$Mef2c_2_A<-1

model_filter_norm_data_scaled_2$Nkx2.5_2_V<-0
model_filter_norm_data_scaled_2[model_filter_norm_data_scaled_2$Nkx2.5V==1 & (model_filter_norm_data_scaled_2$Nkx2.5_V+1)/(model_filter_norm_data_scaled_2$Nkx2.5_A+1)>=1.5,]$Nkx2.5_2_V<-1
model_filter_norm_data_scaled_2$Nkx2.5_2_A<-0
model_filter_norm_data_scaled_2[model_filter_norm_data_scaled_2$Nkx2.5A==1 & (model_filter_norm_data_scaled_2$Nkx2.5_A+1)/(model_filter_norm_data_scaled_2$Nkx2.5_V+1)>=1.5,]$Nkx2.5_2_A<-1

model_filter_norm_data_scaled_2$Tbx5_2_V<-0
model_filter_norm_data_scaled_2[model_filter_norm_data_scaled_2$Tbx5V==1 & (model_filter_norm_data_scaled_2$Tbx5_V+1)/(model_filter_norm_data_scaled_2$Tbx5_A+1)>=1.5,]$Tbx5_2_V<-1
model_filter_norm_data_scaled_2$Tbx5_2_A<-0
model_filter_norm_data_scaled_2[model_filter_norm_data_scaled_2$Tbx5A==1 & (model_filter_norm_data_scaled_2$Tbx5_A+1)/(model_filter_norm_data_scaled_2$Tbx5_V+1)>=1.5,]$Tbx5_2_A<-1

model_filter_norm_data_scaled_2$Srf_2_V<-0
model_filter_norm_data_scaled_2[model_filter_norm_data_scaled_2$SrfV==1 & (model_filter_norm_data_scaled_2$Srf_V+1)/(model_filter_norm_data_scaled_2$Srf_A+1)>=1.5,]$Srf_2_V<-1
model_filter_norm_data_scaled_2$Srf_2_A<-0
model_filter_norm_data_scaled_2[model_filter_norm_data_scaled_2$SrfA==1 & (model_filter_norm_data_scaled_2$Srf_A+1)/(model_filter_norm_data_scaled_2$Srf_V+1)>=1.5,]$Srf_2_A<-1

SSA_matrix_merge_A_2 <- SSA_merge(model_filter_norm_data_scaled_2,c("Gata4_2_A","Mef2a_2_A","Mef2c_2_A",
                                                                    "Nkx2.5_2_A","Tbx5_2_A","Srf_2_A","Tead1_2_A"),c("Asp"))
SSA_matrix_merge_V_2 <- SSA_merge(model_filter_norm_data_scaled_2,c("Gata4_2_V","Mef2a_2_V","Mef2c_2_V",
                                                                    "Nkx2.5_2_V","Tbx5_2_V","Srf_2_V","Tead1_2_V"),
                                  c("Vsp"))


SSA_matrix_merged_2 <- rbind(SSA_matrix_merge_A_2,SSA_matrix_merge_V_2)
SSA_matrix_merge_A_matrix_2 <- SSA_merge_matrix(model_filter_norm_data_scaled_2,c("Gata4_2_A","Mef2a_2_A","Mef2c_2_A",
                                                                                  "Nkx2.5_2_A","Tbx5_2_A","Srf_2_A","Tead1_2_A"),
                                                c("Asp"))
SSA_matrix_merge_V_matrix_2 <- SSA_merge_matrix(model_filter_norm_data_scaled_2,c("Gata4_2_V","Mef2a_2_V","Mef2c_2_V",
                                                                                  "Nkx2.5_2_V","Tbx5_2_V","Srf_2_V","Tead1_2_V"  ),
                                                c("Vsp"))
write.table(SSA_matrix_merge_A_matrix_2, file="Asp_TFs_confmatrix.txt",sep = "\t",quote = F)
write.table(SSA_matrix_merge_V_matrix_2, file="Vsp_TFs_confmatrix.txt",sep = "\t",quote = F)

pdf("SSA_spMTF_spTFs.pdf",width = 7,height = 6)
A<-ggplot(SSA_matrix_merge_A_2,aes(y=factor(Factor,levels=rev(c("Gata4_2_A","Mef2a_2_A","Mef2c_2_A",
                                                                "Nkx2.5_2_A","Tbx5_2_A","Srf_2_A","Tead1_2_A"))),
                                   x=as.numeric(value),group=label,fill=factor(label)))+
  geom_bar(stat="identity")+ facet_wrap(vars(factor(label,levels=c("Specificity","Sensitivity","Accuracy"))),ncol = 3)+
  theme_classic()
V<-ggplot(SSA_matrix_merge_V_2,aes(y=factor(Factor,levels =rev(c("Gata4_2_V","Mef2a_2_V","Mef2c_2_V",
                                                                 "Nkx2.5_2_V","Tbx5_2_V","Srf_2_V","Tead1_2_V"))),
                                   x=as.numeric(value),group=label,fill=factor(label)))+
  geom_bar(stat="identity")+ facet_wrap(vars(factor(label,levels=c("Specificity","Sensitivity","Accuracy"))),ncol = 3)+
  theme_classic()

A/V
dev.off()
###################################################################################
model_filter_norm_data_scaled_2 <-inter_join_2[,c("AspMTF","VspMTF","Asp","Vsp","AspP300",
                                                  "VspP300","AspGene","VspGene")]

#ATAC_anntation<-read.table(file = "MPRA_7/new_annotation_sensitive/MPRA_with_A-V_ATAC_As_Vs_ATAC_full_list.txt")
#loops_anntation<-read.table(file = "MPRA_7/new_annotation_sensitive/MPRA_with_A-V_loops_As_Vs_loops_full_list.txt")
dim(ATAC_anntation)
dim(model_filter_norm_data_scaled_2)
model_filter_norm_data_scaled_2$ATAC_Asp <-ATAC_anntation[ATAC_anntation$V1 %in% rownames(model_filter_norm_data_scaled_2),"V4"]
model_filter_norm_data_scaled_2$ATAC_Vsp <-ATAC_anntation[ATAC_anntation$V1 %in% rownames(model_filter_norm_data_scaled_2),"V5"]
model_filter_norm_data_scaled_2$loops_Asp <-loops_anntation[loops_anntation$V1 %in% rownames(model_filter_norm_data_scaled_2),"V4"]
model_filter_norm_data_scaled_2$loops_Vsp <-loops_anntation[loops_anntation$V1 %in% rownames(model_filter_norm_data_scaled_2),"V5"]
SSA_matrix_merge_A_2 <- SSA_merge(model_filter_norm_data_scaled_2,c("AspMTF","AspP300",
                                                                    "AspGene","ATAC_Asp","loops_Asp"),c("Asp"))
SSA_matrix_merge_V_2 <- SSA_merge(model_filter_norm_data_scaled_2,c("VspMTF","VspP300","VspGene","ATAC_Vsp","loops_Vsp"),
                                  c("Vsp"))
SSA_matrix_merged_2 <- rbind(SSA_matrix_merge_A_2,SSA_matrix_merge_V_2)
SSA_matrix_merge_A_matrix_2 <- SSA_merge_matrix(model_filter_norm_data_scaled_2,c("AspMTF","AspP300",
                                                                                  "AspGene","ATAC_Asp","loops_Asp"),
                                                c("Asp"))
SSA_matrix_merge_V_matrix_2 <- SSA_merge_matrix(model_filter_norm_data_scaled_2,c("VspMTF","VspP300","VspGene","ATAC_Vsp","loops_Vsp"),
                                                c("Vsp"))
write.table(SSA_matrix_merge_A_matrix_2, file="Asp_confmatrix.txt",sep = "\t",quote = F)
write.table(SSA_matrix_merge_V_matrix_2, file="Vsp_confmatrix.txt",sep = "\t",quote = F)

pdf("SSA_spMTF_spP300_spgene_spATAC_sploops.pdf",width = 7,height = 6)
A<-ggplot(SSA_matrix_merge_A_2,aes(y=factor(Factor,levels=rev(c("AspMTF","AspP300",
                                                                "AspGene","ATAC_Asp","loops_Asp"))),
                                   x=as.numeric(value),group=label,fill=factor(label)))+
  geom_bar(stat="identity")+ facet_wrap(vars(factor(label,levels=c("Specificity","Sensitivity","Accuracy"))),ncol = 3)+
  theme_classic()
V<-ggplot(SSA_matrix_merge_V_2,aes(y=factor(Factor,levels =rev(c("VspMTF","VspP300","VspGene","ATAC_Vsp","loops_Vsp"))),
                                   x=as.numeric(value),group=label,fill=factor(label)))+
  geom_bar(stat="identity")+ facet_wrap(vars(factor(label,levels=c("Specificity","Sensitivity","Accuracy"))),ncol = 3)+
  theme_classic()

A/V
dev.off()
##################################################################################






#confmatrix<-table(model_filter_norm_data_scaled_2$AspGene_,model_filter_norm_data_scaled_2$A_active)
write.table(SSA_matrix_merged_2,"AVsplibrary_Sensitivity_specitivity_Accuracy_merged.txt",sep="\t",quote=F)
library(pROC)
roc(model_filter_norm_data_scaled_2$Asp, model_filter_norm_data_scaled_2$loops_Asp)

       