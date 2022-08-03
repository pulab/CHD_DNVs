setwd("~/Desktop/work/Bill/yangpo/MPRA_4")
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

#data<-read.table(file="MPRA_4_barcode_number.txt",header=T,row.names = 1)
#data<-read.table(file="MPRA_4_UMI_number_with_header.txt",header=T,row.names = 1)
data<-read.table(file="AAV_DNA-A-V_mutagenesis_readscount.txt",header=T,row.names = 1)
region_list<-read.table(file="new_mutagenesis_list.txt",header=T)
###################################DESeq2########################
region_list
head(data)
dim(data)
data_location <- data
reads_split<-str_split_fixed(rownames(data_location), pattern=":::",2)
data_location$group_name <-reads_split[,1]
data_location$rownames <-rownames(data_location)
data_sub_r <- as.data.frame(merge(data_location, region_list,by.x="group_name",by.y="region",no.dups=T))$rownames
data_sub_r <- data_sub_r[!is.na(data_sub_r)]
data_sub2 <-data[data_sub_r,]
dim(data_sub2)
Negative_C_raw <- t(dplyr::select(as.data.frame(t(data)), contains("ESC")))
Positive_C_raw <- t(dplyr::select(as.data.frame(t(data)), contains("postive")))
Posivie_list <- rownames(read.table("positive_region_list.txt",row.names = 1))
Positive_list_in<-intersect(rownames(data), Posivie_list)
Positive_C_raw <- data[Positive_list_in,]
dim(Positive_C_raw)
data_sub3<- rbind(data_sub2,Negative_C_raw,Positive_C_raw)
dim(data_sub3)

##########################################

##########################################
filter_for_degs<- data_sub3[,c(1:10,12:15)]
filter_for_degs<- filter_for_degs[rowSums(filter_for_degs)>0,]
#filter_for_degs<- data_sub[rowMins(as.matrix(norm_data[,13:18]), value = T) >=20,13:24]
coldata <-as.data.frame(colnames(filter_for_degs))

coldata$type <- coldata$`colnames(filter_for_degs)`
coldata$condition <-gsub("_\\d","",coldata$type)
genes <-as.matrix((filter_for_degs))
dds <- DESeqDataSetFromMatrix(countData = genes, colData=coldata, design=~condition)
dds <- DESeq(dds)
res1 <- results( dds, contrast = c("condition","A","DNA") )
res2 <- results( dds, contrast = c("condition","V","DNA") )
res1$qvalue_0_05 <- 0
res1$qvalue_0_05[res1$padj < 0.05 & res1$log2FoldChange >0] <- 1
res2$qvalue_0_05 <- 0
res2$qvalue_0_05[res2$padj < 0.05 & res2$log2FoldChange >0] <- 1
sum(res1$qvalue_0_05)
#res2$qvalue_0_05 <- 0
#res2$qvalue_0_05[res2$padj < 0.05 & res2$log2FoldChange >0] <- 1
sum(res2$qvalue_0_05)
######################to FPM#####################################
data <- data_sub3
data_sub <- as.matrix(data_sub3[,c(1:10,12:15)])
norm_data <- t(t(data_sub)/colSums(data_sub)*1000000)
#filter_norm_data <- cbind(norm_data, data[,6:7]) %>% dplyr::filter(DNA >= 5 &rowSums(.[,1:4]>0))
test <- dplyr::select(as.data.frame(norm_data), contains("DNA"))
filter_norm_data <- norm_data[rowMaxs(as.matrix(test), value = T) >=5,]
dim(filter_norm_data)
head(test)
new <- data.frame(matrix(ncol = 0, nrow = length(rownames(filter_norm_data))))
rownames(new) <- rownames(filter_norm_data)

for_plot <- as.data.frame(test) 
ggplot(for_plot, aes(x = DNA_1))+ geom_histogram(binwidth = 5,color="skyblue")+theme_bw()+theme_classic()+
  geom_vline(xintercept = 5,color="skyblue",size=0.4,linetype="dashed")+xlim(0,500)
ggplot(for_plot, aes(x = DNA_4))+ geom_histogram(binwidth = 5,color="skyblue")+theme_bw()+theme_classic()+
  geom_vline(xintercept = 5,color="skyblue",size=0.4,linetype="dashed")+xlim(0,500)
ggplot(for_plot, aes(x = DNA_2))+ geom_histogram(binwidth = 5,color="skyblue")+theme_bw()+theme_classic()+
  geom_vline(xintercept = 5,color="skyblue",size=0.4,linetype="dashed")+xlim(0,500)
ggplot(for_plot, aes(x = DNA_3))+ geom_histogram(binwidth = 5,color="skyblue")+theme_bw()+theme_classic()+
  geom_vline(xintercept = 5,color="skyblue",size=0.4,linetype="dashed")+xlim(0,500)
ggplot(for_plot, aes(x = DNA_5))+ geom_histogram(binwidth = 5,color="skyblue")+theme_bw()+theme_classic()+
  geom_vline(xintercept = 5,color="skyblue",size=0.4,linetype="dashed")+xlim(0,500)
#########################################################3
DNA_total <-as.data.frame(rowMaxs(norm_data[,c(1:5)],value=T))
colnames(DNA_total)<-c("DNA_max_value")
pdf("DNA_counts_distribution_Mutagenesis.pdf")
ggplot(DNA_total, aes(x = DNA_max_value))+ geom_histogram(binwidth = 5,color="skyblue")+theme_bw()+theme_classic()+
  geom_vline(xintercept = 5,color="skyblue",size=0.4,linetype="dashed")+xlim(0,500)
dev.off()

###########################################################
filter_norm_data <- as.data.frame(filter_norm_data)
filter_norm_data$DNA <-rowMeans(filter_norm_data[,1:5])
for (i in 6:14){
  ratio <- (filter_norm_data[,i]+1)/(filter_norm_data$DNA+1)
  ratio <- log2(ratio)
  new <- cbind(new,ratio)
}

new_colname <- paste0(colnames(filter_norm_data[,6:14]), "_vs_DNA_log")
colnames(new) <- new_colname
head(new)
new$log_A_ratio <- rowMeans(new[1:5])
new$log_V_ratio <- rowMeans(new[6:9])
as.data.frame(colnames(filter_norm_data))
as.data.frame(colnames(new))
####################################################################
total_number <- length(DNA_total[,1])
kept_number <- length(filter_norm_data[,1])
keptratio <- kept_number/total_number
DNA_counts <- data.frame(
  label=c("total_number","kept_number","keptratio"),
  value=c(as.character(total_number),as.character(kept_number),as.character(keptratio))
)

write.table(DNA_counts,"Mutagenesis_library_DNA_distribution_filter_log.txt",sep="\t",quote=F)
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

Negative_C <- t(dplyr::select(as.data.frame(t(new)), !contains("ESC_enhancers")))
#Positive_C <- t(dplyr::select(as.data.frame(t(new)), contains("postive")))
#Posivie_list <- rownames(read.table("positive_region_list.txt",row.names = 1))
Positive_list_in_2<-intersect(rownames(new), Posivie_list)
Positive_C2 <- new[Positive_list_in_2,]
both_C <- unique(c(rownames(Negative_C),rownames(Positive_C2)))
filter_norm_data_nc_normalize_exp <-MPRA_SizeFactors(new,both_C,c(1:9))
filter_norm_data_nc_normalize_exp <- as.data.frame(filter_norm_data_nc_normalize_exp)


A_V_cor <- cor(filter_norm_data_nc_normalize_exp)
A_V_cor2 <- cor(new[,1:9])
#A_V_cor2 <- cor(new[,c(1:5,7:10)])
corrplot(A_V_cor, method  = "number",tl.cex=0.5,hclust.method="median")
corrplot(A_V_cor2, method  = "number",tl.cex=0.5,hclust.method="median")
######################################pcc figure ##########################################
my_cols <- c("black") 
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

pdf("MPRA_mutagenesis_pcc.pdf",width = 11, height = 11)
A_V_cor <- cor(filter_norm_data_nc_normalize_exp)
p<-corrplot(A_V_cor2, method  = "number",order = "alphabet")
pairs(filter_norm_data_nc_normalize_exp, 
      lower.panel = panel.cor,
      upper.panel = upper.panel)
dev.off()

###########################################################################
###################barcode vs UMI#######################################
#UMI <- new[,c(1:5,7:10)]
#compare <- intersect(rownames(UMI),rownames(new[,c(1:5,7:10)]))
#data_test <- new[,c(1:5,7:10)]
#A_V_cor <- cbind(UMI[compare,],data_test[compare,])
#colnames(A_V_cor) <- c("UMI_A1","UMI_A2","UMI_A3","UMI_A4","UMI_A5","UMI_V2","UMI_V3","UMI_V4","UMI_V5","raw_A1","raw_A2","raw_A3","raw_A4","raw_A5","raw_V2","raw_V3","raw_V4","raw_V5")
#A_V_cor2 <- cor(A_V_cor)
#library(RColorBrewer)
#corrplot(A_V_cor2, method  = "number",order = "alphabet",number.font = 2,number.cex = 0.8)
#################################################################################
ggplot(filter_norm_data_nc_normalize_exp,aes(x=V_5_vs_DNA_log,y=V_4_vs_DNA_log))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(-7,7)+ylim(-7,7)
ggplot(filter_norm_data_nc_normalize_exp,aes(x=A_1_vs_DNA_log,y=A_2_vs_DNA_log))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(-10,10)+ylim(-10,10)
ggplot(filter_norm_data_nc_normalize_exp,aes(x=RNA1_vs_DNA_log,y=RNA4_vs_DNA_log))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(0,3)+ylim(0,3)
ggplot(filter_norm_data_nc_normalize_exp,aes(x=RNA2_vs_DNA_log,y=RNA3_vs_DNA_log))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(0,3)+ylim(0,3)

write.table(new,file="UMI_mutagenesis_log2_ratio_add1_new_region.txt",sep="\t")
write.table(filter_norm_data,file="UMI_mutagenesis_log2_FPM_add1_new_region.txt",sep="\t")
write.table(filter_norm_data_nc_normalize_exp,file="UMI_mutagenesis_log2_normalized_log2_ratio_add1_new_region.txt",sep="\t")
filter_norm_data_scaled<- filter_norm_data_nc_normalize_exp
filter_norm_data_scaled$log_A_ratio <- rowMeans(filter_norm_data_scaled[1:5])
filter_norm_data_scaled$log_V_ratio <- rowMeans(filter_norm_data_scaled[6:9])
filter_norm_data_scaled$A_active <- res1[rownames(filter_norm_data_scaled),]$qvalue_0_05
filter_norm_data_scaled$V_active <- res2[rownames(filter_norm_data_scaled),]$qvalue_0_05
filter_norm_data_scaled$order_log_A_mean_ratio <- rank(-as.numeric(as.character(filter_norm_data_scaled$log_A_ratio)),ties.method = "first")
filter_norm_data_scaled$order_log_V_mean_ratio <- rank(-as.numeric(as.character(filter_norm_data_scaled$log_V_ratio)),ties.method = "first")
#########################################FDR 0.05#######################

Negative_C <- t(dplyr::select(as.data.frame(t(filter_norm_data_scaled)), contains("ESC")))

Negative_C <-as.data.frame(Negative_C)
fitg<-fitdist(as.numeric(Negative_C$log_A_ratio) ,"norm")
FDR_A_0.95 <- qnorm(0.95, mean = fitg$estimate["mean"], sd = fitg$estimate["sd"])

fitg2<-fitdist(as.numeric(Negative_C$log_V_ratio) ,"norm")
FDR_V_0.95 <- qnorm(0.95, mean = fitg2$estimate["mean"], sd = fitg$estimate["sd"])

FDR_A_0.95
FDR_V_0.95
FDR_A_0.95<-0.6706277
FDR_V_0.95<-1.28182
select_REF <- t(filter_norm_data_scaled) %>% as.data.frame() %>% dplyr::select(contains(":::W:::")) %>% t()
select_ALT <- t(filter_norm_data_scaled) %>% as.data.frame() %>% dplyr::select(contains(":::M:::")) %>% t()
filter_norm_data_scaled$Group <-""
filter_norm_data_scaled[rownames(select_REF),]$Group <- "REF"
filter_norm_data_scaled[rownames(select_ALT),]$Group <- "ALT"
#filter_norm_data_scaled[rownames(Positive_C),]$Group <- "PositiveControl"
filter_norm_data_scaled[rownames(Negative_C),] $Group<- "NegativeControl"
filter_norm_data_scaled[rownames(Positive_C2),]$Group <- "PositiveControlmouse"
ggplot(filter_norm_data_scaled, aes(x = log_A_ratio, color=Group))+geom_density()+theme_bw()+theme_classic()
ggplot(Negative_C, aes(x = as.numeric(log_A_ratio)))+geom_density()+theme_bw()+theme_classic()
ggplot(Negative_C, aes(x = as.numeric(log_V_ratio)))+geom_density()+theme_bw()+theme_classic()
library(patchwork)
pdf("mutagenesis_library_figures_new_region.pdf", width = 11, height = 8)
A<-ggplot(filter_norm_data_scaled,aes(x=order_log_A_mean_ratio,y=log_A_ratio,color=factor(A_active)))+
  geom_point(size=2)+
  geom_hline(yintercept=FDR_A_0.95,color="red",linetype="dashed")+theme_bw()+
  scale_color_manual(breaks = c("1", "0"), values=c("red", "gray"))+
  theme_classic()+geom_text(aes(0,round(FDR_A_0.95,3),label = round(FDR_A_0.95,3), vjust = -1))
B<-ggplot(filter_norm_data_scaled, aes(x = order_log_A_mean_ratio, y = factor(Group,levels=c("NegativeControl",
                                                                                             "PositiveControlmouse",
                                                                                "ALT","REF"))))+
  geom_bin2d(aes(alpha = ..density..), fill = "slategray",binwidth = c(3,1)) +
  scale_fill_continuous(guide = "none")+
  theme_bw()
A/B
#stat_bin2d(aes(fill = after_stat(count)), binwidth = c(1,1),color="Black")+theme_bw()

#ggplot(filter_norm_data_scaled, aes(x = order_log_A_mean_ratio, y = Group))+stat_bin2d(aes(fill = after_stat(count)), binwidth = c(3,1))
#ggplot(filter_norm_data_scaled, aes(x = order_log_A_mean_ratio, y = MTF_2))+stat_bin2d(aes(fill = after_stat(count)), binwidth = c(3,1))
#ggplot(filter_norm_data_scaled, aes(x = order_log_A_mean_ratio, y = P300_2))+stat_bin2d(aes(fill = after_stat(count)), binwidth = c(3,1))
#ggplot(filter_norm_data_scaled, aes(x = order_log_A_mean_ratio, y = Gene))+stat_bin2d(aes(fill = after_stat(count)), binwidth = c(3,1))

C<-ggplot(filter_norm_data_scaled,aes(x=order_log_V_mean_ratio,y=log_V_ratio,color=factor(V_active)))+
  geom_point(size=2)+
  geom_hline(yintercept=FDR_V_0.95,color="red",linetype="dashed")+theme_bw()+
  scale_color_manual(breaks = c("1", "0"), values=c("red", "gray"))+
  theme_classic()+geom_text(aes(0,round(FDR_V_0.95,3),label = round(FDR_V_0.95,3), vjust = -1))
D<-ggplot(filter_norm_data_scaled, aes(x = order_log_V_mean_ratio, y = factor(Group,levels=c("NegativeControl",
                                                                                             "PositiveControlmouse",
                                                                                             "ALT","REF"))))+
  geom_bin2d(aes(alpha = ..density..), fill = "slategray",binwidth = c(3,1)) +
  scale_fill_continuous(guide = "none")+
  theme_bw()
C/D
dev.off()
write.table(filter_norm_data_scaled,file="mutagenesis_log2_normalized_log2_ratio_add1.txt",sep="\t")
filter_norm_data_scaled<-read.table(file="mutagenesis_log2_normalized_log2_ratio_add1.txt",header = T,row.names = 1)
dim(filter_norm_data_scaled[filter_norm_data_scaled$Group=="REF",]) 
dim(filter_norm_data_scaled[filter_norm_data_scaled$Group=="REF",])
#####################################################################
library(kmer)
colnames(filter_norm_data)
filter_norm_data$RNA_A <- rowMeans(filter_norm_data[6:10])
filter_norm_data$RNA_V <- rowMeans(filter_norm_data[11:15])
filter_norm_data_scaled$RNA_A <- filter_norm_data$RNA_A
filter_norm_data_scaled$RNA_V <- filter_norm_data$RNA_V
filter_norm_data_scaled$DNA <- filter_norm_data$DNA
ggplot(filter_norm_data_scaled,aes(x=RNA_A,y=RNA_V,color=factor(V_active)))+
  geom_point(size=2)+
  geom_hline(yintercept=FDR_V_0.95,color="red",linetype="dashed")+theme_bw()+
  scale_color_manual(breaks = c("1", "0"), values=c("red", "darkblue"))+
  theme_classic()+geom_text(aes(0,round(FDR_V_0.95,3),label = round(FDR_V_0.95,3), vjust = -1))

woodmouse <- filter_norm_data_scaled[,c("RNA_A","RNA_V","DNA")]
woodmouse <- woodmouse[, apply(woodmouse, 2, function(v) !any(v == 0xf0))]
woodmouse.kdist <- kdistance(woodmouse, k = 6)
print(as.matrix(woodmouse.kdist)[1:7, 1:7], digits = 2)
###################################################################

##############################A plot ###############################################33
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
################################A  plot ###########################
filter_norm_data_scaled_grouped <- filter_norm_data_scaled
filter_norm_data_scaled_grouped_test_A <-Enrich_score(filter_norm_data_scaled_grouped, "Group", "log_A_ratio", "order_log_A_mean_ratio")
dat_text_A <-filter_norm_data_scaled_grouped_test_A$dat_text
filter_norm_data_for_plot_A <- filter_norm_data_scaled_grouped_test_A$datasets
pdf("A_mutagenesis_enrich_score_new_with_rank.pdf")
p_A<-ggplot(filter_norm_data_for_plot_A, aes(order_log_A_mean_ratio, final_ES,color=factor(score_group)))+geom_line()+theme_bw()+theme_classic()+
  geom_hline(yintercept = 0,color="skyblue",size=0.4,linetype="dashed")+geom_text(
  data    = dat_text_A,
  mapping = aes(x = -Inf, y = -Inf, label = label),
  hjust   = -2,
  vjust   = -0.25) 
A<-ggplot(filter_norm_data_scaled,aes(x=order_log_A_mean_ratio,y=log_A_ratio,color=factor(A_active)))+
  geom_point(size=2)+
  geom_hline(yintercept=FDR_A_0.95,color="red",linetype="dashed")+theme_bw()+
  scale_color_manual(breaks = c("1", "0"), values=c("red", "gray"))+
  theme_classic()+geom_text(aes(0,round(FDR_A_0.95,3),label = round(FDR_A_0.95,3), vjust = -1))
B<-ggplot(filter_norm_data_scaled, aes(x = order_log_A_mean_ratio, y = factor(Group,levels=c("NegativeControl",
                                                                                             "PositiveControlmouse",
                                                                                             "ALT","REF"))))+
  geom_bin2d(aes(alpha = ..density..), fill = "slategray",binwidth = c(3,1)) +
  scale_fill_continuous(guide = "none")+
  theme_bw()
A/p_A/B

dev.off()
####################V plot####################################
filter_norm_data_scaled_grouped_test_V <-Enrich_score(filter_norm_data_scaled_grouped, "Group", "log_V_ratio", "order_log_V_mean_ratio")
dat_text_V <-filter_norm_data_scaled_grouped_test_V$dat_text
filter_norm_data_for_plot_V <- filter_norm_data_scaled_grouped_test_V$datasets
pdf("V_mutagenesis_enrich_score_new_with_rank.pdf")
p_V<-ggplot(filter_norm_data_for_plot_V, aes(order_log_V_mean_ratio, final_ES,color=factor(score_group)))+geom_line()+theme_bw()+theme_classic()+
  geom_hline(yintercept = 0,color="skyblue",size=0.4,linetype="dashed")+geom_text(
  data    = dat_text_V,
  mapping = aes(x = -Inf, y = -Inf, label = label),
  hjust   = -2,
  vjust   = -0.25) 
C<-ggplot(filter_norm_data_scaled,aes(x=order_log_V_mean_ratio,y=log_V_ratio,color=factor(V_active)))+
  geom_point(size=2)+
  geom_hline(yintercept=FDR_V_0.95,color="red",linetype="dashed")+theme_bw()+
  scale_color_manual(breaks = c("1", "0"), values=c("red", "gray"))+
  theme_classic()+geom_text(aes(0,round(FDR_V_0.95,3),label = round(FDR_V_0.95,3), vjust = -1))
D<-ggplot(filter_norm_data_scaled, aes(x = order_log_V_mean_ratio, y = factor(Group,levels=c("NegativeControl",
                                                                                             "PositiveControlmouse",
                                                                                             "ALT","REF"))))+
  geom_bin2d(aes(alpha = ..density..), fill = "slategray",binwidth = c(3,1)) +
  scale_fill_continuous(guide = "none")+
  theme_bw()
C/p_V/D
dev.off()
######################################################################################



################################diff enhancers for A####################3
#table <- read.table("mutagenesis_log2_normalized_log2_ratio_add1.txt", header = T, row.names = 1)
head(table)
table <- filter_norm_data_scaled_grouped
table_RNA <- table %>% dplyr::select(contains("DNA_log"))
table_RNA_A <- subset(table_RNA, select = startsWith(names(table_RNA), "A")) 
table_RNA_V <- subset(table_RNA, select = startsWith(names(table_RNA), "V")) 
table_active_A <- table %>% dplyr::select(contains("A_active"))
table_active_V <- table %>% dplyr::select(contains("V_active"))
select_REF_A <- t(table_RNA_A) %>% as.data.frame() %>% dplyr::select(contains(":::W:::")) %>% t()
select_REF_V <- t(table_RNA_V) %>% as.data.frame() %>% dplyr::select(contains(":::W:::")) %>% t()
select_ALT_A <- t(table_RNA_A) %>% as.data.frame() %>% dplyr::select(contains(":::M:::")) %>% t()
select_ALT_V <- t(table_RNA_V) %>% as.data.frame() %>% dplyr::select(contains(":::M:::")) %>% t()
select_REF_active_A <- t(table_active_A) %>% as.data.frame() %>% dplyr::select(contains(":::W:::")) %>% t()
select_REF_active_V <- t(table_active_V) %>% as.data.frame() %>% dplyr::select(contains(":::W:::")) %>% t()
sum(as.data.frame(select_REF_active_A)$A_active)
sum(as.data.frame(select_REF_active_V)$V_active)
select_ALT_active_A <- t(table_active_A) %>% as.data.frame() %>% dplyr::select(contains(":::M:::")) %>% t()
select_ALT_active_V <- t(table_active_V) %>% as.data.frame() %>% dplyr::select(contains(":::M:::")) %>% t()
sum(as.data.frame(select_ALT_active_A)$A_active)
sum(as.data.frame(select_ALT_active_V)$V_active)
#########################
head(table)
#table_RNA <- table %>% dplyr::select(contains("DNA_log"))
table_RNA <- table %>% dplyr::select(contains("DNA_log"))
table_RNA_A <- table %>% dplyr::select(contains( "DNA_log")) %>% dplyr::select(starts_with( "A_"))
table_RNA_V <- table %>% dplyr::select(contains( "DNA_log")) %>% dplyr::select(starts_with( "V_"))
select_REF <- as.data.frame(t(table_RNA) %>% as.data.frame() %>% dplyr::select(contains(":::W:::")) %>% t())
select_ALT <- as.data.frame(t(table_RNA) %>% as.data.frame() %>% dplyr::select(contains(":::M:::")) %>% t())
select_REF2 <-as.data.frame(t(table) %>% as.data.frame() %>% dplyr::select(contains(":::W:::")) %>% t())
select_ALT2 <- as.data.frame(t(table) %>% as.data.frame() %>% dplyr::select(contains(":::M:::")) %>% t())
rownames_REF<-str_split_fixed(rownames(select_REF), pattern=":::W",2)[,1]
rownames_ALT <- str_split_fixed(rownames(select_ALT), pattern=":::M",2)[,1]
intersect_loc <- intersect(rownames_REF, rownames_ALT)
#new_REF_name <- paste(intersect_loc, "REF", sep = ":::")
rownames(select_REF)<- rownames_REF
rownames(select_ALT)<- rownames_ALT
#new_ALT_name <- paste(intersect_loc, "ALT", sep = ":::")
#select_ALT$pairedname <- new_ALT_name
rownames_REF2 <- str_split_fixed(rownames(select_REF2), pattern=":::W",2)[,1]
rownames_ALT2 <- str_split_fixed(rownames(select_ALT2), pattern=":::M",2)[,1]
rownames(select_REF2)<- rownames_REF2
rownames(select_ALT2)<- rownames_ALT2
intersect_loc2 <- intersect(rownames_REF2, rownames_ALT2)
#new_REF_name2 <- paste(intersect_loc2, "REF", sep = ":::")
#new_ALT_name2 <- paste(intersect_loc2, "ALT", sep = ":::")
inter_REF <- select_REF[intersect_loc,]
inter_ALT <- select_ALT[intersect_loc,]
inter_REF2 <- select_REF2[intersect_loc2,]
inter_ALT2 <- select_ALT2[intersect_loc2,]
pvalue_A <- c()
pvalue_V <- c()
pvalue_WT_AV <- c()
pvalue_ALT_AV <- c()
for (i in 1:dim(inter_REF)[1]) {
  ref_value_A <- as.numeric(inter_REF[i,c(1:5)])
  alt_value_A <- as.numeric(inter_ALT[i,c(1:5)])
  ref_value_V <- as.numeric(inter_REF[i,c(6:9)])
  alt_value_V <- as.numeric(inter_ALT[i,c(6:9)])
  pvalue_i_A <- t.test(ref_value_A, alt_value_A)
  pvalue_A <- c(pvalue_A, pvalue_i_A$p.value)
  pvalue_i_V <- t.test(ref_value_V, alt_value_V)
  pvalue_V <- c(pvalue_V, pvalue_i_V$p.value)
  pvalue_i_WT_AV <- t.test(ref_value_A, ref_value_V)
  pvalue_i_ALT_AV <- t.test(alt_value_A, alt_value_V)
  pvalue_WT_AV <- c(pvalue_WT_AV, pvalue_i_WT_AV$p.value)
  pvalue_ALT_AV <- c(pvalue_ALT_AV, pvalue_i_ALT_AV$p.value)
}
qvalue_A <- pvalue_A * length(pvalue_A) /rank(pvalue_A)
qvalue_V <- pvalue_V * length(pvalue_V) /rank(pvalue_V)
qvalue_WT_AV <- pvalue_WT_AV * length(pvalue_WT_AV) /rank(pvalue_WT_AV)
qvalue_ALT_AV <- pvalue_ALT_AV * length(pvalue_ALT_AV) /rank(pvalue_ALT_AV)
inter_REF <- as.data.frame(inter_REF)
inter_ALT <- as.data.frame(inter_ALT)

inter_REF2 <- as.data.frame(inter_REF2)
inter_ALT2 <- as.data.frame(inter_ALT2)
head(inter_REF2)
head(inter_REF)

################################################################
inter_join <- cbind(inter_REF, inter_ALT)
colnames(inter_join) <- c("A1_REF","A2_REF","A3_REF","A4_REF","A5_REF","V2_REF","V3_REF","V4_REF","V5_REF",
                        "A1_ALT","A2_ALT","A3_ALT","A4_ALT","A5_ALT","V2_ALT","V3_ALT","V4_ALT","V5_ALT")
inter_join$REF_A <- rowMeans(inter_join[,1:5])
inter_join$ALT_A <- rowMeans(inter_join[,10:14])
inter_join$REF_V <- rowMeans(inter_join[,6:9])
inter_join$ALT_V <- rowMeans(inter_join[,15:19])
inter_join$pvalue_A <- as.numeric(as.character(pvalue_A))
inter_join$qvalue_A <- as.numeric(as.character(qvalue_A))
inter_join$pvalue_V <- as.numeric(as.character(pvalue_V))
inter_join$qvalue_V <- as.numeric(as.character(qvalue_V))
inter_join$pvalue_WT_AV <- as.numeric(as.character(pvalue_WT_AV))
inter_join$qvalue_WT_AV <- as.numeric(as.character(qvalue_WT_AV))
inter_join$pvalue_ALT_AV <- as.numeric(as.character(pvalue_ALT_AV))
inter_join$qvalue_ALT_AV <- as.numeric(as.character(qvalue_ALT_AV))
#inter_join$fold <- rowMeans(inter_join[,5:8])/rowMeans(inter_join[,1:4])
inter_join$fold_A <- inter_join$ALT_A-inter_join$REF_A
inter_join$fold_V <- inter_join$ALT_V-inter_join$REF_V
inter_join$fold_WT_AV <- inter_join$REF_V-inter_join$REF_A
inter_join$fold_ALT_AV <- inter_join$ALT_V-inter_join$ALT_A
inter_join$A_REF_active <- inter_REF2[rownames(inter_join),]$A_active
inter_join$V_REF_active <- inter_REF2[rownames(inter_join),]$V_active
inter_join$A_ALT_active <- inter_ALT2[rownames(inter_join),]$A_active
inter_join$V_ALT_active <- inter_ALT2[rownames(inter_join),]$V_active
inter_join$active_A <- 0
inter_join$active_A[inter_join$A_REF_active >0 | inter_join$A_ALT_active >0] <-1
inter_join$active_V <- 0
inter_join$active_V[inter_join$V_REF_active >0 | inter_join$V_ALT_active >0] <-1
inter_join$active_WT <- 0
inter_join$active_WT[inter_join$A_REF_active >0 | inter_join$V_REF_active >0] <-1
inter_join$active_ALT <- 0
inter_join$active_ALT[inter_join$A_ALT_active >0 | inter_join$V_ALT_active >0] <-1
inter_join$final_sig_A <-0
inter_join$final_sig_V <-0
inter_join$final_sig_A[inter_join$qvalue_A < 0.05 &  inter_join$fold_A >= 0.58 & inter_join$active_A >0] <-2
inter_join$final_sig_A[inter_join$qvalue_A <0.05  & inter_join$fold_A <= -0.58 & inter_join$active_A >0] <-1
inter_join$final_sig_V[inter_join$qvalue_V < 0.05 &  inter_join$fold_V >= 0.58 & inter_join$active_V >0] <-2
inter_join$final_sig_V[inter_join$qvalue_V < 0.05  & inter_join$fold_V <= -0.58 & inter_join$active_V >0] <-1

inter_join$final_sig_WT_AV <-0
inter_join$final_sig_ALT_AV <-0
inter_join$final_sig_WT_AV[inter_join$qvalue_WT_AV < 0.05 &  inter_join$fold_WT_AV >= 0.58 & inter_join$active_WT >0] <-2
inter_join$final_sig_WT_AV[inter_join$qvalue_WT_AV <0.05  & inter_join$fold_WT_AV <= -0.58 & inter_join$active_WT >0] <-1
inter_join$final_sig_ALT_AV[inter_join$qvalue_ALT_AV < 0.05 &  inter_join$fold_ALT_AV >= 0.58 & inter_join$active_ALT >0] <-2
inter_join$final_sig_ALT_AV[inter_join$qvalue_ALT_AV < 0.05  & inter_join$fold_ALT_AV <= -0.58 & inter_join$active_ALT >0] <-1
##################################################
write.table(inter_join,file="AV_mutagenesis_MPRA_4_normalized_log2_ratio_fold_pvalue_add1_with_AV_diff_new_region.txt",sep="\t")
inter_join<-read.table(file="AV_mutagenesis_MPRA_4_normalized_log2_ratio_fold_pvalue_add1_with_AV_diff_new_region.txt",header = T, row.names = 1)
#################################################
inter_join$qvalue_0_05_A <- 0
inter_join$qvalue_0_05_A[inter_join$qvalue_A < 0.05] <- 1
inter_join$qvalue_0_05_V <- 0
inter_join$qvalue_0_05_V[inter_join$qvalue_V < 0.05] <- 1

inter_join$qvalue_0_05_WT <- 0
inter_join$qvalue_0_05_WT[inter_join$qvalue_WT_AV < 0.05] <- 1
inter_join$qvalue_0_05_ALT <- 0
inter_join$qvalue_0_05_ALT[inter_join$qvalue_ALT_AV < 0.05] <- 1

######################################################################3

table(inter_join$final_sig_WT_AV)
inter_join_test<-(inter_join[! is.na(inter_join$final_sig_WT_AV),])
unique(rownames(inter_join_test))
####################################################################


ggplot(inter_join,aes(REF_A, ALT_A,color=factor(final_sig_A,levels=c('0','1','2')),alpha =factor(active_A,levels=c("0","1"))))+
  geom_point()+theme_bw()+theme_classic()+geom_abline(slope = 1,color="black",size=0.2,linetype="dashed")+
  geom_abline(slope = 1,intercept=0.58,color="black",size=0.2,linetype="dashed")+
  geom_abline(slope=1,intercept=-0.58,color="black",size=0.2,linetype="dashed")+
  scale_color_manual(breaks = c("0", "1", "2"), values=c("gray", "red", "blue"))+
  scale_alpha_manual(breaks = c("0", "1"), values=c(0.2, 1))+
  geom_vline(xintercept = FDR_A_0.95,color="skyblue",size=0.4,linetype="dashed")+
  geom_hline(yintercept = FDR_A_0.95,color="skyblue",size=0.4,linetype="dashed")+
  xlim(-6,6)+ylim(-6,6)
pdf("diffenhancers_mutagenesis_WTA-ALTA_WTV-ALTV_REFA-REFV_ALTA-ALTV_new_region.pdf",height = 9, width = 11)
ggplot(inter_join,aes(REF_A, ALT_A,color=factor(final_sig_A,levels=c('0','1','2')),alpha =factor(qvalue_0_05_A,levels=c("0","1"))))+
  geom_point()+theme_bw()+theme_classic()+geom_abline(slope = 1,color="black",size=0.2,linetype="dashed")+
  geom_abline(slope = 1,intercept=0.58,color="black",size=0.2,linetype="dashed")+
  geom_abline(slope=1,intercept=-0.58,color="black",size=0.2,linetype="dashed")+
  scale_color_manual(breaks = c("0", "1", "2"), values=c("gray", "red", "blue"))+
  scale_alpha_manual(breaks = c("0", "1"), values=c(0.2, 1))+
  geom_vline(xintercept = FDR_A_0.95,color="skyblue",size=0.4,linetype="dashed")+
  geom_hline(yintercept = FDR_A_0.95,color="skyblue",size=0.4,linetype="dashed")+
  xlim(-6,6)+ylim(-6,6)

ggplot(inter_join,aes(REF_V, ALT_V,color=factor(final_sig_V,levels=c('0','1','2')),alpha =factor(qvalue_0_05_V,levels=c("0","1"))))+
  geom_point()+theme_bw()+theme_classic()+geom_abline(slope = 1,color="black",size=0.2,linetype="dashed")+
  geom_abline(slope = 1,intercept=0.58,color="black",size=0.2,linetype="dashed")+
  geom_abline(slope=1,intercept=-0.58,color="black",size=0.2,linetype="dashed")+
  scale_color_manual(breaks = c("0", "1", "2"), values=c("gray", "red", "blue"))+
  scale_alpha_manual(breaks = c("0", "1"), values=c(0.2, 1))+
  geom_vline(xintercept = FDR_V_0.95,color="skyblue",size=0.4,linetype="dashed")+
  geom_hline(yintercept = FDR_V_0.95,color="skyblue",size=0.4,linetype="dashed")+
  xlim(-6,6)+ylim(-6,6)

ggplot(inter_join,aes(REF_V, REF_A,color=factor(final_sig_WT_AV,levels=c('0','1','2')),alpha =factor(qvalue_0_05_WT,levels=c("0","1"))))+
  geom_point()+theme_bw()+theme_classic()+geom_abline(slope = 1,color="black",size=0.2,linetype="dashed")+
  geom_abline(slope = 1,intercept=0.58,color="black",size=0.2,linetype="dashed")+
  geom_abline(slope=1,intercept=-0.58,color="black",size=0.2,linetype="dashed")+
  scale_color_manual(breaks = c("0", "1", "2"), values=c("gray", "red", "blue"))+
  scale_alpha_manual(breaks = c("0", "1"), values=c(0.2, 1))+
  geom_vline(xintercept = FDR_V_0.95,color="skyblue",size=0.4,linetype="dashed")+
  geom_hline(yintercept = FDR_A_0.95,color="skyblue",size=0.4,linetype="dashed")+
  xlim(-6,6)+ylim(-6,6)

ggplot(inter_join,aes(ALT_V, ALT_A,color=factor(final_sig_ALT_AV,levels=c('0','1','2')),alpha =factor(qvalue_0_05_ALT,levels=c("0","1"))))+
  geom_point()+theme_bw()+theme_classic()+geom_abline(slope = 1,color="black",size=0.2,linetype="dashed")+
  geom_abline(slope = 1,intercept=0.58,color="black",size=0.2,linetype="dashed")+
  geom_abline(slope=1,intercept=-0.58,color="black",size=0.2,linetype="dashed")+
  scale_color_manual(breaks = c("0", "1", "2"), values=c("gray", "red", "blue"))+
  scale_alpha_manual(breaks = c("0", "1"), values=c(0.2, 1))+
  geom_vline(xintercept = FDR_V_0.95,color="skyblue",size=0.4,linetype="dashed")+
  geom_hline(yintercept = FDR_A_0.95,color="skyblue",size=0.4,linetype="dashed")+
  xlim(-6,6)+ylim(-6,6)
dev.off()
number1 <- as.data.frame(table(inter_join$final_sig_A))
number2 <- as.data.frame(table(inter_join$final_sig_V))
number3 <- as.data.frame(table(inter_join$final_sig_WT_AV))
number4 <- as.data.frame(table(inter_join$final_sig_ALT_AV))

###########################################################################################
###########################################################################################
#number <- as.data.frame(table(inter_join$final_sig_A))

dat_text_2<- data.frame(
  label = c("nochange","REForA_specific","ALTorV_specific"),
  REF_ALT_A = number1$Freq
)
dat_text_2$REF_ALT_V <- number2$Freq
dat_text_2$REF_A_V <- number3$Freq
dat_text_2$ALT_A_V <- number4$Freq
write.table(dat_text_2,file = "Diff_enhancers_number_mutagenesis.txt",quote=F, sep = "\t")
W1<-as.data.frame(table(inter_join$A_REF_active))
W2<-as.data.frame(table(inter_join$V_REF_active))
A1 <-as.data.frame(table(inter_join$A_ALT_active))
A2 <-as.data.frame(table(inter_join$V_ALT_active))
dat_text_3 <- cbind(W1, W2,A1,A2) 
colnames(dat_text_3) <-c("WT_A_active","number","WT_V_active","number","ALT_A_active","number","ALT_V_active","number")
write.table(dat_text_3, "mutagenesis_active_regions_number.txt",quote=F, sep = "\t")
###########################################################################################
###########################################################################################
groups_split<-str_split_fixed(rownames(inter_join), pattern=":::",2)
groups_split2<-str_split_fixed(rownames(inter_join), pattern=":::chr",3)[,3]
groups_split3<-str_split_fixed(groups_split2,pattern=":::",2)
#paste0("chr",groups_split3[,1])
#groups_split <- str_split(rownames(inter_join),split = ":::")
#as.matrix(groups_split)
#inter_join_location <- inter_join[,21:53]
inter_join_location <- inter_join[,19:length(inter_join[1,])]
inter_join_location$group_name <-groups_split[,1]
inter_join_location$group_number <-groups_split3[,2]
inter_join_location$group_F <-paste0("chr",groups_split3[,1])
#region_motif <- read.table("region_changed_motifs.list",header = T)
#region_motif <- read.table("region_changed_motifs_gene_V.list",header = T)
region_motif <- read.table("region_changed_motifs_gene_V.txt",header = T)
head(region_motif$motif_name_2)
colnames(region_motif)<-c("region","motif_name_2")
#region_motif$motif_name_1<-str_split_fixed(region_motif$motif_name, pattern="_HUMAN",2)[,1]
#region_motif$motif_name_2<-str_split_fixed(region_motif$motif_name_1, pattern="_MOUSE",2)[,1]
rownames(region_motif) <- region_motif$region
inter_join_location$changed_motif <- region_motif[rownames(inter_join_location),]$motif_name_2
head(inter_join_location)
#select_bias <- inter_join_location %>% as.data.frame() %>% dplyr::select(contains(":::"))
#####################################################
common_400 <- c("chr1:77062496-77062896",
           "chr11:54879982-54880382",
           "chr2:155604926-155605326",
           "chr6:37398994-37399394",
           "chr7:100465696-100466096",
           "chr7:49176880-49177280",
           "chr8:35630234-35630634",
           "chr9:109730912-109731312")
A_400 <-c("chr1:130804279-130804679",
          "chr11:5889033-5889433",
          "chr11:86916688-86917088",
          "chr13:12133446-12133846",
          "chr13:35398106-35398506",
          "chr13:93992860-93993260",
          "chr15:12749980-12750380",
          "chr2:153609584-153609984",
          "chr6:30545609-30546009",
          "chr8:33992986-33993386",
          "chr8:77988334-77988734",
          "chr11:5889071-5889384",
          "chr6:87449513-87449967")
V_400 <-c("chr13:74098998-74099398",
          "chr14:101872884-101873284",
          "chr16:87263521-87263921",
          "chr17:9895269-9895669",
          "chr6:8732525-8732925",
          "chr7:112276741-112277141",
          "chr13:72023637-72023968",
          "chr9:110765610-110766148")
length(A_400)
length(V_400)
length(common_400)
A_mute_WT<-table(inter_join_location[inter_join_location$A_REF_active==1,]$group_name)
test<-group_by (inter_join_location,group_name,group_number)
A_A<-intersect(A_400,names(A_mute_WT))
V_mute_WT<-table(inter_join_location[inter_join_location$V_REF_active==1,]$group_name)
V_A<-intersect(V_400,names(V_mute_WT))
V_mute_WT[V_A]
AV_mute_WT<-table(inter_join_location[inter_join_location$V_REF_active==1 |inter_join_location$A_REF_active==1 ,]$group_name)
C_A<-intersect(common_400,names(AV_mute_WT))
unique(inter_join_location$group_name)
all_names<-c(A_400,common_400,V_400)
####################################################

length(unique(inter_join_location[inter_join_location$A_REF_active >0,"group_name"]))
length(unique(inter_join_location[inter_join_location$active_A >0,"group_name"]))
length(unique(inter_join_location[,"group_name"]))

E_names3_A_s<-c("chr6:30545609-30546009")
#E_names4_V_s<-c("chr7:145087775-145088175")
E_names4_V_s<-c("chr7:100465696-100466096")
E_names2 <- c("chr11:5889033-5889433")###Myl7
E_names <- c("chr6:30545609-30546009")
E_names <- c("chr7:80130433-80130833")
E_names <- c("chr9:110765610-110766148") ##Myl3
E_names <- E_names4_V_s
E_names2 <- E_names3_A_s
#active_paired <- unique(inter_join_location[inter_join_location$A_REF_active >0,c("group_name","group_number")])
active_paired <- unique(inter_join_location[,c("group_name","group_number","A_REF_active","REF_A")])
#ggplot(active_paired, aes(x=as.numeric(group_number),y=group_name))+geom_bin_2d()+
#  theme_bw()+theme_classic()

active_paired$A_REF_active=as.numeric(active_paired$A_REF_active)
#filter_active_paired <- active_paired
filter_active_paired <- active_paired[active_paired$A_REF_active == 1,]

rows_numner <- as.data.frame(seq(-16,97,1))

new_df <- NULL
for (i in 1:dim(filter_active_paired)[1]){
  start = as.numeric(as.character(filter_active_paired$group_number[i])) - 17
  end = as.numeric(as.character(filter_active_paired$group_number[i])) + 17
  new_data_frame <- data.frame(group_name = rep(filter_active_paired$group_name[i],35),
                               group_number = c(start:end))
  new_data_frame$A_REF_active <- filter_active_paired$A_REF_active[i]
  new_data_frame$REF_A <- filter_active_paired$REF_A[i]
  new_data_frame$center_ID <- filter_active_paired$group_number[i]
  new_df <- rbind(new_df, new_data_frame)
}


#################################################################
ggplot(new_df,aes(group_number, reorder(center_ID,-group_number), fill = REF_A))+geom_raster(interpolate = TRUE)+
  #scale_color_manual(breaks = c("0", "1", ""), values=c("gray", "red"))+
  geom_tile(colour = "grey50")+xlim(-17,98)+
  geom_vline(xintercept = 1,color="skyblue",size=0.4,linetype="dashed")+
  geom_vline(xintercept = 80,color="skyblue",size=0.4,linetype="dashed")+
  scale_fill_gradient2(limits =c(-5,5),low="gray", mid="white",high="red")+
  theme_classic2()+facet_wrap(~group_name,nrow =29,strip.position="left")+
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())
##########################################################################
one_enhancer_t <- new_df[new_df$group_name == E_names,]
one_enhancer_t_2 <- new_df[new_df$group_name == E_names2,]

#one_enhancer_t <- one_enhancer_t[!is.na(one_enhancer_t$REF_A),]
#my_palette <- colorRampPalette(c("blue","white", "orange"))(n = 299)
#pdf("myl7_bin_A.pdf")
pdf("myl3_bin_A.pdf")
pdf("top2_V_bin_A.pdf")
ggplot(one_enhancer_t,aes(group_number, reorder(center_ID,-group_number), fill = REF_A))+geom_raster()+
  #scale_color_manual(breaks = c("0", "1", ""), values=c("gray", "red"))+
  geom_tile(colour = "grey50")+xlim(-17,98)+
  geom_vline(xintercept = 1,color="skyblue",size=0.4,linetype="dashed")+
  geom_vline(xintercept = 80,color="skyblue",size=0.4,linetype="dashed")+
  scale_fill_gradient2(limits =c(-5,5),low="gray", mid="white",high="red")+
  theme_classic2()
dev.off()
pdf("myl7_bin_A.pdf")
pdf("top2_A_bin_A.pdf")
ggplot(one_enhancer_t_2,aes(group_number, reorder(center_ID,-group_number), fill = REF_A))+geom_raster()+
  #scale_color_manual(breaks = c("0", "1", ""), values=c("gray", "red"))+
  geom_tile(colour = "grey50")+xlim(-17,98)+
  geom_vline(xintercept = 1,color="skyblue",size=0.4,linetype="dashed")+
  geom_vline(xintercept = 80,color="skyblue",size=0.4,linetype="dashed")+
  scale_fill_gradient2(limits =c(-5,5),low="gray", mid="white",high="red")+
  theme_classic2()
dev.off()
###########################################################################333##########
  
sum_new_df <- new_df %>% group_by(group_name, group_number) %>% 
  dplyr::summarise(sum = sum(REF_A)) %>%
  as.data.frame()
new_df_matrix <- acast(data = sum_new_df, formula = group_name~group_number,value = sum, 
                       fill = 0)
#########################3

sum_new_df_active <- new_df %>% group_by(group_name, group_number) %>% 
  dplyr::summarise(sum = sum(A_REF_active)) %>%
  as.data.frame()
new_df_matrix_active_A <- acast(data = sum_new_df_active, formula = group_name~group_number,value = A_REF_active, 
                                fill = 0)
A_A
pdf("A_active_regions.pdf")
pheatmap(new_df_matrix_active_A[A_A,],cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("white", "orange", "firebrick3"))(50))
pheatmap(new_df_matrix_active_A[C_A,],cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("white", "orange", "firebrick3"))(50))
dev.off()
######################################################
#new_df_matrix <- acast(data = sum_new_df, formula = group_name~group_number,value = A_REF_active, 
#                       fill = 0)
whole_mat_A <- matrix(0,nrow=length(all_names),ncol=114)
rownames(whole_mat_A) <- as.character(all_names)
colnames(whole_mat_A) <- as.character(seq(-17,96,1))
whole_mat_A[as.character(rownames(new_df_matrix)),as.character(colnames(new_df_matrix))] <- new_df_matrix
#A_row_names_2 <- pheatmap(new_df_matrix,cluster_cols = F,cluster_rows = T,color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
#                        broder_color="Black",gaps_col="Yellow"  )
#A_row_names_2 <- pheatmap2(new_df_matrix,cluster_cols = F,cluster_rows = F,border_color = "Black"
#                          )
#E_names <- c("chr9:110765610-110766148")

rownames(rows_numner) <- rows_numner[,1]

library(gplots)
#my_palette <- colorRampPalette(c("blue","white", "orange"))(n = 299)
my_palette <- colorRampPalette(c("white", "orange","red"))(n = 299)
heatmap.2(new_df_matrix, breaks=seq(0, 40, length.out=300),dendrogram = "none",cexCol = 0.8,key = T,
          col=my_palette, Colv = F, Rowv = F, trace = "none", scale="none",sepcolor="black",colsep=1:ncol(new_df_matrix),
          rowsep=1:nrow(new_df_matrix),sepwidth=c(0.1,0.1))
new_df_matrix_A <-new_df_matrix
only_one<-as.data.frame(new_df_matrix[E_names,])
only_one_2<-as.data.frame(new_df_matrix[E_names2,])
rows_numner$value <- 0
rows_numner$value2<-0
rows_numner[rownames(only_one),]$value <- only_one[,1]
rows_numner[rownames(only_one_2),]$value2 <- only_one_2[,1]

pdf("Myl3_bottom_A.pdf",height = 12,width = 0.7)
pheatmap(rows_numner[,2],cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("gray", "white", "firebrick3"))(50),
                                broder_color="Black", gaps_row = c(17,97),show_colnames=F,breaks=seq(0,60,1.2))
row_names<-pheatmap(rows_numner[,2],cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("gray", "white", "firebrick3"))(50),
                    broder_color="Black", gaps_row = c(17,97),show_colnames=F,breaks=seq(0,60,1.2))
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

pdf("Myl3_bottom_A_2.pdf",height = 16,width = 0.7)

row_names<-pheatmap(rows_numner[c(18:97),2],cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("gray", "white", "firebrick3"))(50),
                    broder_color="Black", show_colnames=F,breaks=seq(0,60,1.2))
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
pdf("top2_V_bottom_A_2.pdf",height = 16,width = 0.7)

row_names<-pheatmap(rows_numner[c(18:97),2],cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("gray", "white", "firebrick3"))(50),
                    broder_color="Black", show_colnames=F,breaks=seq(0,60,1.2))
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
pdf("Myl7_bottom_A.pdf",height = 12,width = 0.7)
pheatmap(rows_numner[,3],cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("gray", "white", "firebrick3"))(50),
         broder_color="Black", gaps_row = c(17,97),show_colnames=F,breaks=seq(0,60,1.2))
row_names<-pheatmap(rows_numner[,3],cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("gray", "white", "firebrick3"))(50),
                    broder_color="Black", gaps_row = c(17,97),show_colnames=F,breaks=seq(0,60,1.2))
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
pdf("Myl7_bottom_A_2.pdf",height = 16,width = 0.7)
#pheatmap(rows_numner[c(17:97),3],cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("gray", "white", "firebrick3"))(50),
#         broder_color="Black",show_colnames=F,breaks=seq(0,60,1.2))
#row_names<-pheatmap(rows_numner[,3],cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("gray", "white", "firebrick3"))(50),
#                    broder_color="Black", gaps_row = c(17,97),show_colnames=F,breaks=seq(0,60,1.2))
row_names<-pheatmap(rows_numner[c(18:97),3],cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("gray", "white", "firebrick3"))(50),
                    broder_color="Black", show_colnames=F,breaks=seq(0,60,1.2))
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
pdf("top2_A_bottom_A_2.pdf",height = 16,width = 0.7)

row_names<-pheatmap(rows_numner[c(18:97),3],cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("gray", "white", "firebrick3"))(50),
                    broder_color="Black", show_colnames=F,breaks=seq(0,60,1.2))
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

active_paire_matrix <-acast(active_paired, group_name~group_number,value=A_REF_active, fill =0)
colnames(active_paire_matrix)<-as.numeric(colnames(active_paire_matrix))
active_paire_matrix_new_A <-active_paire_matrix[,order(as.numeric(colnames(active_paire_matrix)))]
##########################kmeans ##################################
#kmeansObj <- kmeans(active_paire_matrix_new_A, centers = 3)
#image(t(active_paire_matrix_new_A)[, order(kmeansObj$cluster)], yaxt = "n", main = "Clustered Data")
#####################################################################

A_row_names <- pheatmap(active_paire_matrix_new_A,cluster_cols = F,cluster_rows = T,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
A_row_names_sorted<-rownames(active_paire_matrix_new_A[A_row_names$tree_row[["order"]],])
#test <-as.data.frame(colSums(active_paire_matrix_new))

ggplot(active_paired, aes(x=as.numeric(group_number)))+geom_bar(color="gray")+
  theme_bw()+theme_classic()
####################################################
length(unique(inter_join_location[inter_join_location$V_REF_active >0,"group_name"]))
length(unique(inter_join_location[inter_join_location$active_V >0,"group_name"]))
length(unique(inter_join_location[,"group_name"]))

active_paired <- unique(inter_join_location[,c("group_name","group_number","V_REF_active","REF_V")])
#ggplot(active_paired, aes(x=as.numeric(group_number),y=group_name))+geom_bin_2d()+
#  theme_bw()+theme_classic()

active_paired$V_REF_active=as.numeric(active_paired$V_REF_active)
#filter_active_paired <- active_paired
filter_active_paired <- active_paired[active_paired$V_REF_active == 1,]
new_df <- NULL
for (i in 1:dim(filter_active_paired)[1]){
  start = as.numeric(as.character(filter_active_paired$group_number[i])) - 17
  end = as.numeric(as.character(filter_active_paired$group_number[i])) + 17
  new_data_frame <- data.frame(group_name = rep(filter_active_paired$group_name[i],35),
                               group_number = c(start:end))
  new_data_frame$V_REF_active <- filter_active_paired$V_REF_active[i]
  new_data_frame$REF_V <- filter_active_paired$REF_V[i]
  new_data_frame$center_ID <- filter_active_paired$group_number[i]
  new_df <- rbind(new_df, new_data_frame)
}
sum_new_df <- new_df %>% group_by(group_name, group_number) %>% 
  dplyr::summarise(sum = sum(REF_V)) %>%
  as.data.frame()
new_df_matrix <- acast(data = sum_new_df, formula = group_name~group_number,value = V_REF_active, 
                       fill = 0)
whole_mat_V <- matrix(0,nrow=length(all_names),ncol=114)
rownames(whole_mat_V) <- as.character(all_names)
colnames(whole_mat_V) <- as.character(seq(-17,96,1))
whole_mat_V[as.character(rownames(new_df_matrix)),as.character(colnames(new_df_matrix))] <- new_df_matrix
############################################
sum_new_df_active <- new_df %>% group_by(group_name, group_number) %>% 
  dplyr::summarise(sum = sum(V_REF_active)) %>%
  as.data.frame()
new_df_matrix_active_V <- acast(data = sum_new_df_active, formula = group_name~group_number,value = V_REF_active, 
                       fill = 0)
V_A
pdf("V_active_regions.pdf")
pheatmap(new_df_matrix_active_V[V_A,],cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("white", "orange", "firebrick3"))(50))
pheatmap(new_df_matrix_active_V[C_A,],cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("white", "orange", "firebrick3"))(50))
dev.off()
##########################################################################
one_enhancer_t <- new_df[new_df$group_name == E_names,]
one_enhancer_t_2 <- new_df[new_df$group_name == E_names2,]
#my_palette <- colorRampPalette(c("blue","white", "orange"))(n = 299)
pdf("myl3_bin_V.pdf")
pdf("top2_V_bin_V.pdf")
ggplot(one_enhancer_t,aes(group_number, reorder(center_ID,-group_number), fill = REF_V))+geom_raster()+
  #scale_color_manual(breaks = c("0", "1", ""), values=c("gray", "red"))+
  geom_tile(colour = "grey50")+xlim(-17,98)+
  geom_vline(xintercept = 0,color="skyblue",size=0.4,linetype="dashed")+
  geom_vline(xintercept = 80,color="skyblue",size=0.4,linetype="dashed")+
  scale_fill_gradient2(limits =c(-3.5,3.5),low="gray", mid="white",high="red")+
  theme_classic2()
dev.off()
pdf("myl7_bin_V.pdf")
pdf("top2_A_bin_V.pdf")
ggplot(one_enhancer_t_2,aes(group_number, reorder(center_ID,-group_number), fill = REF_V))+geom_raster()+
  #scale_color_manual(breaks = c("0", "1", ""), values=c("gray", "red"))+
  geom_tile(colour = "grey50")+xlim(-17,98)+
  geom_vline(xintercept = 0,color="skyblue",size=0.4,linetype="dashed")+
  geom_vline(xintercept = 80,color="skyblue",size=0.4,linetype="dashed")+
  scale_fill_gradient2(limits =c(-3.5,3.5),low="gray", mid="white",high="red")+
  theme_classic2()
dev.off()
###########################################################################333##########
V_row_names_2 <- pheatmap(whole_mat_V,cluster_cols = F,cluster_rows = T,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
#E_names <- c("chr9:110765610-110766148")

library(gplots)
my_palette <- colorRampPalette(c("lightgray","white", "orange"))(n = 299)
heatmap.2(whole_mat_V, breaks=seq(0, 40, length.out=300),dendrogram = "none",cexCol = 0.8,key = T,
          col=my_palette, Colv = F, Rowv = T, trace = "none", scale="none",sepcolor="black",colsep=1:ncol(whole_mat_V),
          rowsep=1:nrow(whole_mat_V),sepwidth=c(0.1,0.1))
rows_numner$value <- 0
rows_numner$value2 <- 0
rownames(rows_numner) <- rows_numner[,1]
########## without active
only_one<-as.data.frame(new_df_matrix[E_names,])
only_one_2<-as.data.frame(new_df_matrix[E_names2,])
rows_numner[rownames(only_one),]$value <- only_one[,1]
rows_numner[rownames(only_one_2),]$value2 <- only_one_2[,1]
pdf("myl3_bottom_V.pdf",height = 12,width = 0.7)

#pheatmap(rows_numner[,2],cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("gray", "white", "firebrick3"))(50),
#         broder_color="Black", gaps_row = c(17,97),show_colnames=F,breaks=seq(0,60,1.2))
row_names<-pheatmap(rows_numner[,2],cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("gray", "white", "firebrick3"))(50),
                    broder_color="Black", gaps_row = c(17,97),show_colnames=F,breaks=seq(0,60,1.2))
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
pdf("myl3_bottom_V_2.pdf",height = 16,width = 0.7)

#pheatmap(rows_numner[,2],cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("gray", "white", "firebrick3"))(50),
#         broder_color="Black", gaps_row = c(17,97),show_colnames=F,breaks=seq(0,60,1.2))
row_names<-pheatmap(rows_numner[c(18:97),2],cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("gray", "white", "firebrick3"))(50),
                    broder_color="Black", show_colnames=F,breaks=seq(0,60,1.2))
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
pdf("top2_V_bottom_V_2.pdf",height = 16,width = 0.7)

#pheatmap(rows_numner[,2],cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("gray", "white", "firebrick3"))(50),
#         broder_color="Black", gaps_row = c(17,97),show_colnames=F,breaks=seq(0,60,1.2))
row_names<-pheatmap(rows_numner[c(18:97),2],cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("gray", "white", "firebrick3"))(50),
                    broder_color="Black", show_colnames=F,breaks=seq(0,60,1.2))
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

pdf("myl7_bottom_V.pdf",height = 12,width = 0.7)
#pheatmap(rows_numner[,3],cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("gray", "white", "firebrick3"))(50),
#         broder_color="Black", gaps_row = c(17,97),show_colnames=F,breaks=seq(0,60,1.2))
#row_names<-pheatmap(rows_numner[,3],cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("gray", "white", "firebrick3"))(50),
#                    broder_color="Black", gaps_row = c(17,97),show_colnames=F,breaks=seq(0,60,1.2))
row_names<-pheatmap(rows_numner[c(18:97),3],cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("gray", "white", "firebrick3"))(50),
                    broder_color="Black", show_colnames=F,breaks=seq(0,60,1.2))
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
pdf("myl7_bottom_V_2.pdf",height = 16,width = 0.7)
#pheatmap(rows_numner[,3],cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("gray", "white", "firebrick3"))(50),
#         broder_color="Black", gaps_row = c(17,97),show_colnames=F,breaks=seq(0,60,1.2))
row_names<-pheatmap(rows_numner[c(18:97),3],cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("gray", "white", "firebrick3"))(50),
                    broder_color="Black",show_colnames=F,breaks=seq(0,60,1.2))
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
pdf("top2_A_bottom_V_2.pdf",height = 16,width = 0.7)
#pheatmap(rows_numner[,3],cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("gray", "white", "firebrick3"))(50),
#         broder_color="Black", gaps_row = c(17,97),show_colnames=F,breaks=seq(0,60,1.2))
row_names<-pheatmap(rows_numner[c(18:97),3],cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("gray", "white", "firebrick3"))(50),
                    broder_color="Black",show_colnames=F,breaks=seq(0,60,1.2))
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

active_paire_matrix <-acast(active_paired, group_name~group_number,value=V_REF_active, fill =0)
colnames(active_paire_matrix)<-as.numeric(colnames(active_paire_matrix))
active_paire_matrix_new_V <-active_paire_matrix[,order(as.numeric(colnames(active_paire_matrix)))]
V_row_names<- pheatmap(active_paire_matrix_new_V,cluster_cols = F,cluster_rows =T,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
#A_row_names <- pheatmap(active_paire_matrix_new_A,cluster_cols = F,cluster_rows = T,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
V_row_names_sorted<-rownames(active_paire_matrix_new_V[V_row_names$tree_row[["order"]],])
#test <-as.data.frame(colSums(active_paire_matrix_new))
ggplot(active_paired, aes(x=as.numeric(group_number)))+geom_bar(color="gray")+
  theme_bw()+theme_classic()
##########################################################
#merge_AV <- merge(new_df_matrix_A,new_df_matrix, by="row.names")
merge_AV <- merge(whole_mat_A,whole_mat_V, by="row.names")
merge_AV[is.na(merge_AV)]<-0
rownames(merge_AV) <-merge_AV$Row.names
merge_AV_new <- merge_AV[,2:length(merge_AV[1,])]
set.seed(1)
###################################
merge_AV_new$A<-rowMeans(merge_AV_new[,1:114])
merge_AV_new$V<-rowMeans(merge_AV_new[,115:228])
merge_AV_new$Amax<-rowMaxs(as.matrix(merge_AV_new[,1:114]),value = T)
merge_AV_new$Vmax<-rowMaxs(as.matrix(merge_AV_new[,115:228]),value = T)
########################################
#merge_AV_new$A<-rowMeans(merge_AV_new[,1:110])
#merge_AV_new$V<-rowMeans(merge_AV_new[,111:223])
#merge_AV_new$Amax<-rowMaxs(as.matrix(merge_AV_new[,1:110]),value = T)
#merge_AV_new$Vmax<-rowMaxs(as.matrix(merge_AV_new[,111:223]),value = T)
############################################
merge_AV_new$AV <-merge_AV_new$A-merge_AV_new$V
merge_AV_new$AV_2 <-0
#merge_AV_new[merge_AV_new$AV>0,]$AV_2 <-1
#merge_AV_new[merge_AV_new$AV<0,]$AV_2 <-2
merge_AV_new[A_400,]$AV_2 <-1
merge_AV_new[V_400,]$AV_2 <-2
#A_row_names <- pheatmap(active_paire_matrix_new_A,cluster_cols = F,cluster_rows = T,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
#row_names_sorted <- c("chr1:130804279-130804679", "chr11:5889033-5889433",    
#                      "chr11:5889071-5889384" ,  "chr6:30545609-30546009",    "chr13:35398106-35398506",   "chr9:109730912-109731312",
#                      "chr1:77062496-77062896",    "chr2:153609584-153609984",  "chr8:77988334-77988734", "chr12:79044668-79045068",   
#                      "chr8:33992986-33993386", 
#                      "chr11:54879982-54880382", "chr15:12749980-12750380", 
#                      "chr6:37398994-37399394",    "chr11:30257774-30258174",  
#                      "chr6:8732525-8732925",      
#                      "chr6:87449513-87449967",      
#                      "chr13:72023637-72023968",   "chr2:122652384-122652784",
#                      "chr2:154390099-154390499",  "chr15:77027673-77028073", 
#                      "chr9:110765610-110766148",  "chr7:100465696-100466096", "chr8:35630234-35630634" ,
#                      "chr2:155604926-155605326", "chr7:112276741-112277141", "chr11:4296432-4296832",    
#                      "chr7:145087775-145088175", "chr7:80130433-80130833", "chr13:74098998-74099398", "chr10:42463240-42463640","chr11:86916688-86917088" ,
#                      "chr16:87263521-87263921", "chr13:12133446-12133846", "chr14:101872884-101873284","chr2:167775302-167775702", 
#                      "chr7:49176880-49177280","chr18:61346368-61346768","chr13:93992860-93993260",   "chr17:9895269-9895669", "chr1:160825526-160825926", 
#                      "chr5:51591402-51591802"
                      
#                       )
#row_names_sorted<-rownames(merge_AV_new[row_names$tree_row[["order"]],])
#row_names_sorted<- rownames(merge_AV_new[order(merge_AV_new$AV_2, merge_AV_new$A,merge_AV_new$V),])
merge_AV_new_1 <-merge_AV_new[merge_AV_new$AV_2==1,]
row_names_sorted_1 <- rownames(merge_AV_new_1[order(-merge_AV_new_1$Amax),])

merge_AV_new_2 <-merge_AV_new[merge_AV_new$AV_2==2,]
row_names_sorted_2 <- rownames(merge_AV_new_2[order(merge_AV_new_2$Vmax),])
merge_AV_new_3 <-merge_AV_new[merge_AV_new$AV_2==0,]
row_names_sorted_3 <- rownames(merge_AV_new_3[order(-merge_AV_new_3$Amax),])
row_names_sorted <-c(row_names_sorted_1,row_names_sorted_3,row_names_sorted_2)
merge_AV_new<-merge_AV_new[row_names_sorted,]
#V_row_names_sorted <- row_names_sorted
#A_row_names_sorted <- row_names_sorted

library(viridis)
pdf("AV_merge_activity_new_region.pdf",width = 16, height = 6)
pdf("AV_merge_only_activity_new_region.pdf",width = 16, height = 6)
pheatmap(as.matrix(merge_AV_new[,1:228]),cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("gray", "white", "firebrick3"))(50),
         broder_color="blue", show_colnames=F,breaks=seq(0,40,0.8),gaps_col = c(17,97,114,131,211))

my_palette <- colorRampPalette(c("white",brewer.pal(9,"OrRd")))(n = 50)
library("RColorBrewer")
#heatmap.2(as.matrix(merge_AV_new[,1:228]), breaks=seq(0, 40, length.out=300),dendrogram = "none",cexCol = 0.5,key = T,
#          col=my_palette, Colv = F, Rowv = F, trace = "none", scale="none",sepcolor="black",colsep= c(17,98,115,132,213),
#          rowsep=1:nrow(merge_AV_new),sepwidth=c(0.01,0.01))
#row_names<-pheatmap(merge_AV_new[,1:228],cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("white", "pink","firebrick3","purple"))(50),
#                    broder_color="Black", show_colnames=T,breaks=seq(0,40,0.8),gaps_row = c(13,21),gaps_col = c(17,97,114,131,211))
row_names<-pheatmap(merge_AV_new[,1:228],cluster_cols = F,cluster_rows = F,color =my_palette,
                    broder_color="Black", show_colnames=T,breaks=seq(0,40,0.8),gaps_row = c(13,21),gaps_col = c(17,97,114,131,211))
row_names<-pheatmap(merge_AV_new[,1:228],cluster_cols = F,cluster_rows = F,color =viridis(50),
                    broder_color="Black", show_colnames=T,breaks=seq(0,40,0.8),gaps_row = c(13,21),gaps_col = c(17,97,114,131,211))
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
row_names<-pheatmap(merge_AV_new[,18:97],cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("gray", "white", "firebrick3"))(50),
                    broder_color="Black", show_colnames=T,breaks=seq(0,40,0.8),gaps_row = c(12))
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

row_names<-pheatmap(merge_AV_new[,c(132:211)],cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("gray", "white", "firebrick3"))(50),
                    broder_color="Black", show_colnames=T,breaks=seq(0,40,0.8),gaps_row = c(12))
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

data_AV <-read.table("../MPRA_5/MPRA_7/A_V_mouse_MPRA_7_UMI_normalized_count_fd_version_add_1.txt",header = T,row.names = 1)
data_AV$names <- str_split_fixed(rownames(data_AV), pattern="::",2)[,2]
intersect(data_AV$names,row_names_sorted)
rownames(data_AV) <- data_AV$names
data_AV_filter <-data_AV[intersect(data_AV$names,row_names_sorted) ,c("A","V")]
data_AV_filter_2 <-data_AV[intersect(data_AV$names,row_names_sorted) ,c("final_sig","A","V","fold","A_active","V_active")]
length(data_AV_filter_2$A)
dim(data_AV_filter_2)
data_AV_filter[row_names_sorted,]
pdf("mutagenes_to_400bp_new_region.pdf")
pheatmap(as.matrix(data_AV_filter[row_names_sorted,]),cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("gray", "white", "firebrick3"))(50),
         broder_color="blue", show_colnames=F,na_col = "lightblue")
dev.off()
data_AV_filter3<-data.frame(c(A_400,common_400,V_400))
rownames(data_AV_filter3)<-c(A_400,common_400,V_400)
A_new<- row_names_sorted[row_names_sorted %in% A_400]
common_new <-row_names_sorted[row_names_sorted %in% common_400]
V_new <-row_names_sorted[row_names_sorted %in% V_400]
row_names_sorted_re_sort <-c(A_new,common_new,V_new)
data_AV_filter3$type<-0
data_AV_filter3[A_400,]$type<-1
data_AV_filter3[V_400,]$type<--1
merge_AV_new<-merge_AV_new[row_names_sorted_re_sort,]

######################################
row_names<-pheatmap(merge_AV_new[,1:223],cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("white", "orange", "firebrick3"))(50),
                    broder_color="Black", show_colnames=T,breaks=seq(0,40,0.8),gaps_row = c(13,21),gaps_col = c(17,97,114,131,211))
#row_names<-pheatmap(merge_AV_new[,1:228],cluster_cols = F,cluster_rows = F,color = viridis(50),
#                    broder_color="Black", show_colnames=T,breaks=seq(0,40,0.8),gaps_row = c(12),gaps_col = c(17,97,114,131,211))
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
#########################



pheatmap(as.matrix(data_AV_filter3[row_names_sorted,"type"]),cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("blue","white","firebrick3"))(50),
         broder_color="blue", show_colnames=T,na_col = "lightblue",breaks=seq(-1,1,0.04))

pheatmap(as.matrix(data_AV_filter_2[row_names_sorted,]),cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("gray", "blue", "firebrick3"))(50),
         broder_color="blue", show_colnames=T,na_col = "lightblue",breaks=seq(0,2,0.04))
pheatmap(as.matrix(data_AV_filter_2[row_names_sorted,]$fold),cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("blue", "white", "firebrick3"))(50),
         broder_color="blue", show_colnames=T,na_col = "lightblue",breaks=seq(-2,2,0.08))
pheatmap(as.matrix(data_AV_filter_2[row_names_sorted,]),cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("gray", "firebrick3"))(50),
         broder_color="blue", show_colnames=T,na_col = "lightblue",breaks=seq(0,1,0.02))
write.table(data_AV_filter_2[row_names_sorted,], file = "mutagenesis_with_AV_specific.txt",sep = "\t",quote = F)
##########################################################
merge_AV <- cbind(active_paire_matrix_new_A,active_paire_matrix_new_V)
AV_row_names<- pheatmap(merge_AV,cluster_cols = F,cluster_rows =T,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
kmeansObj2 <- kmeans(merge_AV, centers = 5)
image(t(merge_AV)[, order(kmeansObj2$cluster)], yaxt = "n", main = "Clustered Data")
library(vegan)
library(jaccard)
jaccard_list <-c()
for (line in 1:dim(merge_AV)[1]){
  jaccard_tem <- jaccard(merge_AV[line,1:80], merge_AV[line, 81:160])
  jaccard_list <-c(jaccard_list,jaccard_tem)
}
merge_AV <- as.data.frame(merge_AV)
merge_AV$ja <- jaccard_list
image(t(merge_AV[,1:160])[, order(merge_AV$ja)], yaxt = "n", main = "Clustered Data")
########################################################
#E_names <- c("chr9:110765610-110766148")
active_paired <- unique(inter_join_location[,c("group_name","group_number","fold_A", "final_sig_A")])
active_paired$fold_A[active_paired$final_sig_A==0]="NA" 
active_paired$fold_A <-as.numeric(active_paired$fold_A)
active_paired <- unique(active_paired[,c("group_name","group_number","fold_A")])
active_paire_matrix <-acast(active_paired, group_name~group_number,value=fold_A, fill = NA)
#active_paire_matrix <-acast(active_paired, group_name~group_number,value=fold_A)
colnames(active_paire_matrix)<-as.numeric(colnames(active_paire_matrix))
active_paire_matrix_new_WT_ALT_A <-active_paire_matrix[,order(as.numeric(colnames(active_paire_matrix)))]
active_paire_matrix_new_WT_ALT_A <- active_paire_matrix_new_WT_ALT_A[A_row_names_sorted,]
library("RColorBrewer")
display.brewer.all()
pdf(file="A_WT_ALT_fold_sorted_with_active_score_new_region.pdf", width = 16, height = 6)
#pheatmap(active_paire_matrix_new_WT_ALT_A,cluster_cols = F,cluster_rows = F,breaks=seq(-4,4,0.08),color = colorRampPalette(c("Green", "white", "Gold"))(100))
#pheatmap(active_paire_matrix_new_WT_ALT_A,cluster_cols = F,cluster_rows = F,breaks=seq(-4,4,0.08),color = colorRampPalette(rev(brewer.pal(10, "BrBG")))(100))
row_names<-pheatmap(active_paire_matrix_new_WT_ALT_A,cluster_cols = F,cluster_rows = F,breaks=seq(-4,4,0.08),
                    gaps_row=c(12),color = colorRampPalette(rev(brewer.pal(10, "BrBG")))(100))

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
pheatmap(active_paire_matrix_new_WT_ALT_A[E_names,],cluster_cols = F,cluster_rows = F,breaks=seq(-4,4,0.08),color = colorRampPalette(c("Green", "white", "Gold"))(100))
#test <-as.data.frame(colSums(active_paire_matrix_new))
#ggplot(active_paired, aes(x=as.numeric(group_number)))+geom_bar(color="gray")+
#  theme_bw()+theme_classic()
###################################################

#############################################################
active_paired <- unique(inter_join_location[,c("group_name","group_number","fold_V", "final_sig_V")])
#active_paired$value=active_paired
#active_paire_matrix <-acast(active_paired, group_name~group_number,value=fold_V, fill =0)
active_paired$fold_V[active_paired$final_sig_V==0]="NA" 
active_paired$fold_V <-as.numeric(active_paired$fold_V)
active_paired <- unique(active_paired[,c("group_name","group_number","fold_V")])
active_paire_matrix <-acast(active_paired, group_name~group_number,value=fold_V, fill = NA)
colnames(active_paire_matrix)<-as.numeric(colnames(active_paire_matrix))
active_paire_matrix_new_WT_ALT_V <-active_paire_matrix[,order(as.numeric(colnames(active_paire_matrix)))]
active_paire_matrix_new_WT_ALT_V <- active_paire_matrix_new_WT_ALT_V[V_row_names_sorted,]
#pdf(file="V_WT_ALT_fold_sorted_with_active_score.pdf", width = 16, height = 6)
#pheatmap(active_paire_matrix_new_WT_ALT_V,cluster_cols = F,cluster_rows = F,breaks=seq(-4,4,0.16),color = colorRampPalette(c("Green", "white", "Gold"))(50))
#pheatmap(active_paire_matrix_new_WT_ALT_V,cluster_cols = F,cluster_rows = F,breaks=seq(-4,4,0.08),color = colorRampPalette(rev(brewer.pal(10, "BrBG")))(100))
#dev.off()
###################################################
pdf("V_WT_ALT_fold_sorted_with_active_score_new_region.pdf",width = 16, height = 6)
row_names<-pheatmap(active_paire_matrix_new_WT_ALT_V,cluster_cols = F,cluster_rows = F,breaks=seq(-4,4,0.08),
                    gaps_row=c(12),color = colorRampPalette(rev(brewer.pal(10, "BrBG")))(100))

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
#######################################################
pheatmap(active_paire_matrix_new_WT_ALT_V[E_names,],cluster_cols = F,cluster_rows = F,breaks=seq(-4,4,0.08),color = colorRampPalette(c("Green", "white", "Gold"))(100))
#######################################################



##################################################3
active_paired <- unique(inter_join_location[,c("group_name","group_number","fold_WT_AV", "final_sig_WT_AV")])
#active_paired$value=active_paired
#active_paire_matrix <-acast(active_paired, group_name~group_number,value=fold_WT_AV, fill =0)
active_paired$fold_WT_AV[active_paired$final_sig_WT_AV==0]="NA" 
active_paired$fold_WT_AV <-as.numeric(active_paired$fold_WT_AV)
active_paired <- unique(active_paired[,c("group_name","group_number","fold_WT_AV")])
active_paire_matrix <-acast(active_paired, group_name~group_number,value=fold_WT_AV, fill = NA)
colnames(active_paire_matrix)<-as.numeric(colnames(active_paire_matrix))
active_paire_matrix_new_WT_AV <-active_paire_matrix[,order(as.numeric(colnames(active_paire_matrix)))]
active_paire_matrix_new_WT_AV <- active_paire_matrix_new_WT_AV[V_row_names_sorted,]
#pdf(file="WT_A_V_fold_sorted_with_active_score.pdf", width = 16, height = 6)
#pheatmap(active_paire_matrix_new_WT_AV,cluster_cols = F,cluster_rows = F,breaks=seq(-4,4,0.16),color = colorRampPalette(c("Green", "white", "Gold"))(50))
#pheatmap(active_paire_matrix_new_WT_AV,cluster_cols = F,cluster_rows = F,breaks=seq(-4,4,0.08),color = colorRampPalette(rev(brewer.pal(10, "BrBG")))(100))
#dev.off()
###################################################
pdf("WT_A_V_fold_sorted_with_active_score_new_region.pdf",width = 16, height = 6)
row_names<-pheatmap(active_paire_matrix_new_WT_AV,cluster_cols = F,cluster_rows = F,breaks=seq(-4,4,0.08),
                    gaps_row=c(12),color = colorRampPalette(rev(brewer.pal(10, "BrBG")))(100))

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

#######################################################
pheatmap(active_paire_matrix_new_WT_AV[E_names,],cluster_cols = F,cluster_rows = F,breaks=seq(-4,4,0.08),color = colorRampPalette(c("Green", "white", "Gold"))(100))

active_paired <- unique(inter_join_location[,c("group_name","group_number","fold_ALT_AV","final_sig_ALT_AV")])
#active_paired$value=active_paired
#active_paire_matrix <-acast(active_paired, group_name~group_number,value=fold_ALT_AV, fill =0)
active_paired$fold_ALT_AV[active_paired$final_sig_ALT_AV==0]="NA" 
active_paired$fold_ALT_AV <-as.numeric(active_paired$fold_ALT_AV)
active_paired <- unique(active_paired[,c("group_name","group_number","fold_ALT_AV")])
active_paire_matrix <-acast(active_paired, group_name~group_number,value=fold_ALT_AV, fill = NA)
colnames(active_paire_matrix)<-as.numeric(colnames(active_paire_matrix))
active_paire_matrix_new_ALT_AV <-active_paire_matrix[,order(as.numeric(colnames(active_paire_matrix)))]
active_paire_matrix_new_ALT_AV <- active_paire_matrix_new_ALT_AV[V_row_names_sorted,]
#pdf(file="ALT_A_V_fold_sorted_with_active_score.pdf", width = 16, height = 6)
#pheatmap(active_paire_matrix_new_ALT_AV,cluster_cols = F,cluster_rows = F,breaks=seq(-4,4,0.16),color = colorRampPalette(c("Green", "white", "Gold"))(50))
#pheatmap(active_paire_matrix_new_ALT_AV,cluster_cols = F,cluster_rows = F,breaks=seq(-4,4,0.08),color = colorRampPalette(rev(brewer.pal(10, "BrBG")))(100))
#dev.off()
##########################################################
pdf("ALT_A_V_fold_sorted_with_active_score_new_region.pdf",width = 16, height = 6)
row_names<-pheatmap(active_paire_matrix_new_ALT_AV,cluster_cols = F,cluster_rows = F,breaks=seq(-4,4,0.08),
                    gaps_row=c(12),color = colorRampPalette(rev(brewer.pal(10, "BrBG")))(100))
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
#####################################################
pheatmap(active_paire_matrix_new_ALT_AV[E_names,],cluster_cols = F,cluster_rows = F,breaks=seq(-4,4,0.08),color = colorRampPalette(c("Green", "white", "Gold"))(100))
#merged_V <- cbind(active_paire_matrix_new_V,active_paire_matrix_new_WT_AV,active_paire_matrix_new_ALT_AV)
merged_V <- cbind(active_paire_matrix_new_WT_AV,active_paire_matrix_new_ALT_AV)
merged_V <- cbind(active_paire_matrix_new_WT_AV[E_names,],active_paire_matrix_new_ALT_AV[E_names,])
merged_V <- cbind(active_paire_matrix_new_ALT_AV[E_names,],active_paire_matrix_new_WT_AV[E_names,],active_paire_matrix_new_WT_ALT_V[E_names,],active_paire_matrix_new_WT_ALT_A[E_names,])
pdf("Myl3_figure_middle.pdf",height = 12,width = 1.5)
pdf("top2_V_figure_middle.pdf",height = 12,width = 1.5)
#pheatmap(merged_V,cluster_cols = F,cluster_rows = F,breaks=seq(-4,4,0.16),color = colorRampPalette(c("green", "white", "gold"))(50))
#pheatmap(merged_V,cluster_cols = F,cluster_rows = F,breaks=seq(-4,4,0.08),color = colorRampPalette(rev(brewer.pal(10, "BrBG")))(100))
row_names<-pheatmap(merged_V,cluster_cols = F,cluster_rows = F,breaks=seq(-4,4,0.08),color = colorRampPalette(rev(brewer.pal(10, "BrBG")))(100))
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
#E_names2 <- c("chr11:5889033-5889433")###Myl7
merged_V2 <- cbind(active_paire_matrix_new_WT_AV,active_paire_matrix_new_ALT_AV)
merged_V2 <- cbind(active_paire_matrix_new_WT_AV[E_names2,],active_paire_matrix_new_ALT_AV[E_names2,])
merged_V2 <- cbind(active_paire_matrix_new_ALT_AV[E_names2,],active_paire_matrix_new_WT_AV[E_names2,],active_paire_matrix_new_WT_ALT_V[E_names2,],active_paire_matrix_new_WT_ALT_A[E_names2,])
pdf("Myl7_figure_middle.pdf",height = 12,width = 1.5)
pdf("top2_A_figure_middle.pdf",height = 12,width = 1.5)
#pheatmap(merged_V2,cluster_cols = F,cluster_rows = F,breaks=seq(-4,4,0.16),color = colorRampPalette(c("green", "white", "gold"))(50))
#pheatmap(merged_V2,cluster_cols = F,cluster_rows = F,breaks=seq(-4,4,0.16),color = colorRampPalette(rev(brewer.pal(10, "BrBG")))(100))
row_names<-pheatmap(merged_V2,cluster_cols = F,cluster_rows = F,breaks=seq(-4,4,0.08),color = colorRampPalette(rev(brewer.pal(10, "BrBG")))(100))
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

#####################################################
###################################################
active_paired <- unique(inter_join_location[,c("group_name","group_number","fold_WT_AV")])
#active_paired$value=active_paired
active_paire_matrix <-acast(active_paired, group_name~group_number,value=fold_WT_AV, fill =0)
colnames(active_paire_matrix)<-as.numeric(colnames(active_paire_matrix))
active_paire_matrix_new_WT_AV <-active_paire_matrix[,order(as.numeric(colnames(active_paire_matrix)))]
active_paire_matrix_new_WT_AV <- active_paire_matrix_new_WT_AV[A_row_names_sorted,]
pheatmap(active_paire_matrix_new_WT_AV,cluster_cols = F,cluster_rows = F,breaks=seq(-4,4,0.16),color = colorRampPalette(c("navy", "white", "firebrick3"))(50))


active_paired <- unique(inter_join_location[,c("group_name","group_number","fold_ALT_AV")])
#active_paired$value=active_paired
active_paire_matrix <-acast(active_paired, group_name~group_number,value=fold_ALT_AV, fill =0)
colnames(active_paire_matrix)<-as.numeric(colnames(active_paire_matrix))
active_paire_matrix_new_ALT_AV <-active_paire_matrix[,order(as.numeric(colnames(active_paire_matrix)))]
active_paire_matrix_new_ALT_AV <- active_paire_matrix_new_ALT_AV[A_row_names_sorted,]
pheatmap(active_paire_matrix_new_ALT_AV,cluster_cols = F,cluster_rows = F,breaks=seq(-4,4,0.16),color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

merged_V <- cbind(active_paire_matrix_new_A,active_paire_matrix_new_WT_AV,active_paire_matrix_new_ALT_AV)
pheatmap(merged_V,cluster_cols = F,cluster_rows = T,breaks=seq(-4,4,0.16),color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
#####################################################
#inter_join_location$final_sig_A <- inter_join[rownames(inter_join_location),]$final_sig
length(unique(inter_join_location[inter_join_location$final_sig_A==1,]$group_name))
length(unique(inter_join_location[inter_join_location$final_sig_A==2,]$group_name))
#inter_join_location_1<-inter_join_location[inter_join_location$group_name=="chr20:49039671-49039842",]
#ggplot(inter_join_location_1,aes(as.numeric (group_number), log2(fold),color=factor(fold_1_5,levels=c('0','1','2')),alpha =factor(qvalue_0_05,levels=c("0","1"))))+
#  geom_point()+theme_bw()+theme_classic()+
#  geom_hline(yintercept = 0,color="black",size=0.4,linetype="dashed")+
#  scale_color_manual(breaks = c("0", "1", "2"), values=c("gray", "red", "blue"))+
#  scale_alpha_manual(breaks = c("0", "1"), values=c(0.2, 1))+
#  geom_hline(yintercept = log2(1.5),color="skyblue",size=0.4,linetype="dashed")+
#  geom_hline(yintercept = log2(0.677),color="skyblue",size=0.4,linetype="dashed")+
#  facet_grid(rows = vars(group_name),cols = vars(group_F),scales="free")
####################################################################################
#head(inter_join_location)
#KO_median<- as.data.frame(base::tapply( inter_join_location$mean_KO,inter_join_location$WT, FUN =median))
#colnames(KO_median)<-"KO_median"
#head(KO_median)
#inter_join_location$KO_median <- KO_median[inter_join_location$WT,]
###################################################################################

inter_join_location_1<-inter_join_location[inter_join_location$full_names=="chr20:49530503-49530903",]
inter_join_location_1<-inter_join_location[inter_join_location$full_names=="chr10:35047312-35047712",]
inter_join_location_1<-inter_join_location[inter_join_location$full_names=="chr6:39745032-39745432",]
inter_join_location_1<-inter_join_location[inter_join_location$full_names=="chr3:34053580-34053980",]
full_names_1 <- inter_join_location[inter_join_location$final_sig_A>0,][1:5]

table(unique(full_names_1)$group_name)[as.character(full_names_1$group_name)]

#inter_join_location_1<-inter_join_location[inter_join_location$full_names %in% full_names_1,]
inter_join_location_names <- head(inter_join_location[inter_join_location$final_sig_A>0,]$group_name)
inter_join_location_names <- c("chr7:93870775-93871175","chr14:59536249-59536649","chr4:148231410-148231810","chr1:245748592-245748992")
inter_join_location_names <- c("chr15:96808175-96808575","chr2:47925801-47926201","chr18:22819712-22820112","chr13:39332360-39332760")
inter_join_location_names <- c("chr9:110765610-110766148",
                               "chr7:80130433-80130833",
                               "chr1:130804279-130804679",
                               "chr1:77062496-77062896")
inter_join_location_1<-inter_join_location[inter_join_location$group_name %in% inter_join_location_names,]
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
#melt_inter_join_location_1<- melt(inter_join_location_1[,c(1:2,5:6)])
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
###############################################################################################
inter_join_location_1$mean <- (inter_join_location_1$REF_A+inter_join_location_1$ALT_A)/2
ggplot(inter_join_location_1)+
  geom_point(aes(as.numeric (group_number), as.numeric(as.character(ALT_A)),color=factor(final_sig_A,levels=c('0','1','2')), alpha =factor(qvalue_0_05_A,levels=c("0","1")),size=factor(A_ALT_active,levels=c('0','1'))))+
  theme_bw()+theme_classic()+
  geom_hline(yintercept = 0,color="black",size=0.4,linetype="dashed")+
  scale_color_manual(breaks = c("0", "1", "2"), values=c("gray", "red", "blue"))+
  scale_alpha_manual(breaks = c("0", "1"), values=c(0.2, 1))+
  scale_size_manual(breaks = c("0", "1"), values=c(2, 4))+
  geom_point(aes(as.numeric (group_number), as.numeric(as.character(REF_A)),  
                 alpha =factor(qvalue_0_05_A,levels=c("0","1")),size=factor(A_REF_active,levels=c('0','1'))),
             shape=17, color="darkblue")+
  geom_line(aes(as.numeric(as.character(group_number)), as.numeric(as.character(mean))),
                color="darkgreen",linetype="dashed")+
  #geom_line(aes(as.numeric (group_number), as.numeric(as.character(ALT_A)),alpha =factor(A_ALT_active,levels=c("0","1"))),color="coral",linetype="dotdash")+
  facet_grid(rows = vars(group_name),scales="free")
#############################################################################################
inter_join_location_names <- c("chr9:110765610-110766148") ###myl3 V use it
#inter_join_location_names <- c("chr7:80130433-80130833") ##########ERR
#inter_join_location_names <- c("chr1:130804279-130804679")
#inter_join_location_names <- c("chr1:77062496-77062896")
#inter_join_location_names <- c("chr1:77062496-77062896")
#inter_join_location_names <- c("chr11:5889028-5889427") ########Myl7 1 A??? from yangpo
#inter_join_location_names <- c("chr11:5888940-5889515") ########Myl7 2 A???? from yangpo
inter_join_location_names2 <- c("chr11:5889033-5889433") ########Myl7 2 A use it
#inter_join_location_names <- c("chr11:5889071-5889384") ########Myl7 1 A
#inter_join_location_names <- c("chr13:72023637-72023968")##############Irx1
#inter_join_location_names <- c("chr6:87449513-87449967") ################BMP10
#inter_join_location_names <- c("chr6:30545609-30546009") ################Tbx5
#inter_join_location_names <- c("chr18:61346368-61346768")
inter_join_location_names4_top2_A<-c("chr6:30545609-30546009")
#inter_join_location_names3_top2_V<-c("chr7:145087775-145088175")
inter_join_location_names3_top2_V<-c("chr7:100465696-100466096")
#########
inter_join_location_names<-inter_join_location_names3_top2_V
inter_join_location_names2<- inter_join_location_names4_top2_A
###########
inter_join_location_2<-inter_join_location[inter_join_location$group_name %in% inter_join_location_names,]
inter_join_location_2_A <- inter_join_location_2[,c("REF_A" ,"ALT_A","final_sig_A","qvalue_0_05_A" ,"A_ALT_active","group_number","A_REF_active","changed_motif")]
inter_join_location_2_A$mean <- (inter_join_location_2_A$REF_A+inter_join_location_2_A$ALT_A)/2
inter_join_location_2_A$sample <- "A"
colnames(inter_join_location_2_A) <- c("REF","ALT","final_sig","qvalue_0_05","ALT_active","group_number","REF_active","changed_motif","mean","sample")
inter_join_location_2_V <- inter_join_location_2[,c("REF_V" ,"ALT_V","final_sig_V","qvalue_0_05_V" ,"V_ALT_active","group_number","V_REF_active","changed_motif")]
inter_join_location_2_V$mean <- (inter_join_location_2_V$REF_V+inter_join_location_2_V$ALT_V)/2
inter_join_location_2_V$sample <- "V"
colnames(inter_join_location_2_V) <- c("REF","ALT","final_sig","qvalue_0_05","ALT_active","group_number","REF_active","changed_motif","mean","sample")
inter_join_location_2_AV <- rbind(inter_join_location_2_A,inter_join_location_2_V)
pdf("Myl3_full_plot.pdf",width = 26, height = 16)
pdf("top2_V_full_plot.pdf",width = 26, height = 16)
ggplot(inter_join_location_2_AV)+
  geom_point(aes(as.numeric (group_number), as.numeric(as.character(ALT)),
                 color=factor(final_sig,levels=c('0','1','2')), 
                 alpha =factor(qvalue_0_05,levels=c("0","1")),
                 size=factor(ALT_active,levels=c('0','1'))))+
  theme_bw()+theme_classic()+
  geom_hline(yintercept = 0,color="black",size=0.4,linetype="dashed")+
  scale_color_manual(breaks = c("0", "1", "2"), values=c("gray", "red", "blue"))+
  scale_alpha_manual(breaks = c("0", "1"), values=c(0.2, 1))+
  scale_size_manual(breaks = c("0", "1"), values=c(2, 4))+
  geom_point(aes(as.numeric (group_number), as.numeric(as.character(REF)),  
                 alpha =factor(qvalue_0_05,levels=c("0","1")),size=factor(REF_active,levels=c('0','1'))),
             shape=17, color="darkblue")+
  geom_text(aes(as.numeric (group_number), as.numeric(as.character(REF)),label=changed_motif),
            size=3,angle = 90,check_overlap = T,nudge_y=3)+
  geom_line(aes(as.numeric(as.character(group_number)), as.numeric(as.character(mean))),
            color="darkgreen",linetype="dashed")+
  #geom_line(aes(as.numeric (group_number), as.numeric(as.character(ALT_A)),alpha =factor(A_ALT_active,levels=c("0","1"))),color="coral",linetype="dotdash")+
  facet_grid(rows = vars(sample),scales="free")
dev.off()
colnames(inter_join_location_2_AV)
for_heatmap_plot <- inter_join_location_2_A <- inter_join_location_2[,c("group_number","REF_A" ,"ALT_A","REF_V" ,"ALT_V","changed_motif")]
rownames(for_heatmap_plot) <- for_heatmap_plot$group_number
cols <-as.data.frame(seq(1,80,1))
rownames(cols) <- as.character( cols$`seq(1, 80, 1)`)
colnames(cols)<-c("group_number")
for_heatmap_plot <- merge(for_heatmap_plot,cols,by = "group_number",all.y=T)

for_heatmap_plot <- for_heatmap_plot[order(as.numeric(for_heatmap_plot$group_number)),]
for_heatmap_plot_M <- as.matrix(for_heatmap_plot)
rownames(for_heatmap_plot_M) <- for_heatmap_plot_M[,"changed_motif"]
for_heatmap_plot_M <- for_heatmap_plot_M[,c("ALT_V","ALT_A","REF_V","REF_A")]
rn <- rownames(for_heatmap_plot_M)
for_heatmap_plot_M <- apply(for_heatmap_plot_M, 2, as.numeric)
rownames(for_heatmap_plot_M) <- rn
pdf("myl3_top_figure.pdf",height = 12,width = 5)
pheatmap(for_heatmap_plot_M,cluster_cols = F,cluster_rows = F,breaks=seq(-4,4,0.16), border_color="black",fontsize=8,
         color = colorRampPalette(c("gold", "pink","purple"))(50))
pheatmap(for_heatmap_plot_M,cluster_cols = F,cluster_rows = F,breaks=seq(-4,4,0.16), border_color="black",fontsize=8,na_col="darkgray",
         color = colorRampPalette(c("gray", "white", "firebrick3"))(50))
dev.off()
pdf("myl3_top_figure_2.pdf",width = 5,height = 12)
pdf("top2_V_top_figure_2.pdf",width = 7,height = 12)
row_names2<-pheatmap(for_heatmap_plot_M,cluster_cols = F,cluster_rows = F,breaks=seq(-4,4,0.16), border_color="black",fontsize=8,
                     color = colorRampPalette(c("gold", "pink","purple"))(50))
## Extract the right grob
grob_classes <- purrr::map(row_names2$gtable$grobs, class)
idx_grob <- which(purrr::map_lgl(grob_classes, function(cl) 'gTree' %in% cl))[1]
grob_names <- names(row_names2$gtable$grobs[[idx_grob]]$children)
idx_rect <- grob_names[grep('rect', grob_names)][1]

## Remove borders around cells
row_names2$gtable$grobs[[idx_grob]]$children[[idx_rect]]$gp$col <- "black"
row_names2$gtable$grobs[[idx_grob]]$children[[idx_rect]]$gp$lwd <- 0.5

## Plot result
graphics::plot.new()
print(row_names2)
row_names3<-pheatmap(for_heatmap_plot_M,cluster_cols = F,cluster_rows = F,breaks=seq(-4,4,0.16), border_color="black",fontsize=8,na_col="darkgray",
                     color = colorRampPalette(c("gray", "white", "firebrick3"))(50))
## Extract the right grob
grob_classes <- purrr::map(row_names3$gtable$grobs, class)
idx_grob <- which(purrr::map_lgl(grob_classes, function(cl) 'gTree' %in% cl))[1]
grob_names <- names(row_names3$gtable$grobs[[idx_grob]]$children)
idx_rect <- grob_names[grep('rect', grob_names)][1]

## Remove borders around cells
row_names3$gtable$grobs[[idx_grob]]$children[[idx_rect]]$gp$col <- "black"
row_names3$gtable$grobs[[idx_grob]]$children[[idx_rect]]$gp$lwd <- 0.5

## Plot result
graphics::plot.new()
print(row_names3)
dev.off()
######################################################################################3
inter_join_location_3<-inter_join_location[inter_join_location$group_name %in% inter_join_location_names2,]
inter_join_location_3_A <- inter_join_location_3[,c("REF_A" ,"ALT_A","final_sig_A","qvalue_0_05_A" ,"A_ALT_active","group_number","A_REF_active","changed_motif")]
inter_join_location_3_A$mean <- (inter_join_location_3_A$REF_A+inter_join_location_3_A$ALT_A)/2
inter_join_location_3_A$sample <- "A"
colnames(inter_join_location_3_A) <- c("REF","ALT","final_sig","qvalue_0_05","ALT_active","group_number","REF_active","changed_motif","mean","sample")
inter_join_location_3_V <- inter_join_location_3[,c("REF_V" ,"ALT_V","final_sig_V","qvalue_0_05_V" ,"V_ALT_active","group_number","V_REF_active","changed_motif")]
inter_join_location_3_V$mean <- (inter_join_location_3_V$REF_V+inter_join_location_3_V$ALT_V)/2
inter_join_location_3_V$sample <- "V"
colnames(inter_join_location_3_V) <- c("REF","ALT","final_sig","qvalue_0_05","ALT_active","group_number","REF_active","changed_motif","mean","sample")
inter_join_location_3_AV <- rbind(inter_join_location_3_A,inter_join_location_3_V)
pdf("Myl7_full_plot.pdf",width = 26, height = 16)
pdf("top2_A_full_plot.pdf",width = 26, height = 16)
ggplot(inter_join_location_3_AV)+
  geom_point(aes(as.numeric (group_number), as.numeric(as.character(ALT)),
                 color=factor(final_sig,levels=c('0','1','2')), 
                 alpha =factor(qvalue_0_05,levels=c("0","1")),
                 size=factor(ALT_active,levels=c('0','1'))))+
  theme_bw()+theme_classic()+
  geom_hline(yintercept = 0,color="black",size=0.4,linetype="dashed")+
  scale_color_manual(breaks = c("0", "1", "2"), values=c("gray", "red", "blue"))+
  scale_alpha_manual(breaks = c("0", "1"), values=c(0.2, 1))+
  scale_size_manual(breaks = c("0", "1"), values=c(2, 4))+
  geom_point(aes(as.numeric (group_number), as.numeric(as.character(REF)),  
                 alpha =factor(qvalue_0_05,levels=c("0","1")),size=factor(REF_active,levels=c('0','1'))),
             shape=17, color="darkblue")+
  geom_text(aes(as.numeric (group_number), as.numeric(as.character(REF)),label=changed_motif),
            size=3,angle = 90,check_overlap = T,nudge_y=3)+
  geom_line(aes(as.numeric(as.character(group_number)), as.numeric(as.character(mean))),
            color="darkgreen",linetype="dashed")+
  #geom_line(aes(as.numeric (group_number), as.numeric(as.character(ALT_A)),alpha =factor(A_ALT_active,levels=c("0","1"))),color="coral",linetype="dotdash")+
  facet_grid(rows = vars(sample),scales="free")
dev.off()
colnames(inter_join_location_3_AV)
for_heatmap_plot <- inter_join_location_3_A <- inter_join_location_3[,c("group_number","REF_A" ,"ALT_A","REF_V" ,"ALT_V","changed_motif")]
rownames(for_heatmap_plot) <- for_heatmap_plot$group_number
cols <-as.data.frame(seq(1,80,1))
rownames(cols) <- as.character( cols$`seq(1, 80, 1)`)
colnames(cols)<-c("group_number")
for_heatmap_plot <- merge(for_heatmap_plot,cols,by = "group_number",all.y=T)

for_heatmap_plot <- for_heatmap_plot[order(as.numeric(for_heatmap_plot$group_number)),]
for_heatmap_plot_M <- as.matrix(for_heatmap_plot)
rownames(for_heatmap_plot_M) <- for_heatmap_plot_M[,"changed_motif"]
for_heatmap_plot_M <- for_heatmap_plot_M[,c("ALT_V","ALT_A","REF_V","REF_A")]
rn <- rownames(for_heatmap_plot_M)
for_heatmap_plot_M <- apply(for_heatmap_plot_M, 2, as.numeric)
rownames(for_heatmap_plot_M) <- rn
pdf("myl7_top_figure.pdf",width = 15,height = 12)
pheatmap(for_heatmap_plot_M,cluster_cols = F,cluster_rows = F,breaks=seq(-4,4,0.16), border_color="black",fontsize=8,
         color = colorRampPalette(c("gold", "pink","purple"))(50))
pheatmap(for_heatmap_plot_M,cluster_cols = F,cluster_rows = F,breaks=seq(-4,4,0.16), border_color="black",fontsize=8,na_col="darkgray",
         color = colorRampPalette(c("gray", "white", "firebrick3"))(50))

dev.off()
pdf("myl7_top_figure_2.pdf",width = 15,height = 12)
pdf("top2_A_top_figure_2.pdf",width = 15,height = 12)
row_names2<-pheatmap(for_heatmap_plot_M,cluster_cols = F,cluster_rows = F,breaks=seq(-4,4,0.16), border_color="black",fontsize=8,
                    color = colorRampPalette(c("gold", "pink","purple"))(50))
## Extract the right grob
grob_classes <- purrr::map(row_names2$gtable$grobs, class)
idx_grob <- which(purrr::map_lgl(grob_classes, function(cl) 'gTree' %in% cl))[1]
grob_names <- names(row_names2$gtable$grobs[[idx_grob]]$children)
idx_rect <- grob_names[grep('rect', grob_names)][1]

## Remove borders around cells
row_names2$gtable$grobs[[idx_grob]]$children[[idx_rect]]$gp$col <- "black"
row_names2$gtable$grobs[[idx_grob]]$children[[idx_rect]]$gp$lwd <- 0.5

## Plot result
graphics::plot.new()
print(row_names2)
row_names3<-pheatmap(for_heatmap_plot_M,cluster_cols = F,cluster_rows = F,breaks=seq(-4,4,0.16), border_color="black",fontsize=8,na_col="darkgray",
                     color = colorRampPalette(c("gray", "white", "firebrick3"))(50))
## Extract the right grob
grob_classes <- purrr::map(row_names3$gtable$grobs, class)
idx_grob <- which(purrr::map_lgl(grob_classes, function(cl) 'gTree' %in% cl))[1]
grob_names <- names(row_names3$gtable$grobs[[idx_grob]]$children)
idx_rect <- grob_names[grep('rect', grob_names)][1]

## Remove borders around cells
row_names3$gtable$grobs[[idx_grob]]$children[[idx_rect]]$gp$col <- "black"
row_names3$gtable$grobs[[idx_grob]]$children[[idx_rect]]$gp$lwd <- 0.5

## Plot result
graphics::plot.new()
print(row_names3)
dev.off()

########################################################################################3
ggplot(filter_norm_data_scaled_grouped_2[filter_norm_data_scaled_grouped_2$active_WT>0,], aes(x = group_F))+ geom_bar(color="skyblue")+theme_bw()+theme_classic()
geom_bar(aes(log2(fold),color=factor(fold_1_5,levels=c('0','1','2')),alpha =factor(qvalue_0_05,levels=c("0","1"))))



write.table(inter_join_location,file="AV_mutategenesis_diff_Enhancers_paired_add1_added_AV_new_region.txt",sep="\t")  








inter_join_location[inter_join_location$mean_WT == 1.82784864351185,]

head(inter_join)
WT_selct<- unique(inter_join[,c(2,7:10,12,15)])
row.names(WT_selct)<-WT_selct$WT
WT_selct <- merge(filter_norm_data,WT_selct,by=0)
write.table(WT_selct,file="10bp_Mutagenesis_WT_final.txt",sep="\t",quote = T)
###########################################################################################


