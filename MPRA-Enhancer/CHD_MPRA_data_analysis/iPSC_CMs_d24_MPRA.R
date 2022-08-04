setwd("~/Desktop/work/Bill/Fengxiao/batch1-3")
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
library(patchwork)
#data<-read.table(file="format_table.txt",header=T,row.names = 1)
data<-read.table(file="batch1-2_rawcounts_d24_with_merge.txt",header=T,row.names = 1)
######################to FPM#####################################
dim(data)
data$RNA.5 <- data$d24_RNA.5_1
data_sub <- as.matrix(data[,2:33])
norm_data <- t(t(data_sub)/colSums(data_sub)*1000000)
total_sum <-as.data.frame(colSums(data_sub))
total_sum
ggplot(total_sum, aes(x=row.names(total_sum),y=colSums(data_sub)))+ geom_bar(color="skyblue",stat="identity")+theme_bw()+theme_classic()
#filter_norm_data <- cbind(norm_data, data[,6:7]) %>% dplyr::filter(DNA >= 5 &rowSums(.[,1:4]>0))
test <- dplyr::select(as.data.frame(norm_data), contains("DNA"))
#filter_norm_data <- norm_data[rowMins(as.matrix(test), value = T) >=5,]
filter_norm_data <- test
dim(filter_norm_data)
head(test)

for_plot <- as.data.frame(test) 
ggplot(for_plot, aes(x = DNA.1))+ geom_histogram(binwidth = 5,color="skyblue")+theme_bw()+theme_classic()+
  geom_vline(xintercept = 20,color="skyblue",size=0.4,linetype="dashed")+xlim(0,500)
ggplot(for_plot, aes(x = for_plot$d24_DNA.1_2))+ geom_histogram(binwidth = 5,color="skyblue")+theme_bw()+theme_classic()+
  geom_vline(xintercept = 20,color="skyblue",size=0.4,linetype="dashed")+xlim(0,500)

##############################distribution of DNA####################333
dim(test)
for_plot2 <- melt(test[,c(7:12,19:24)]) 
#ggplot(for_plot2, aes(x = value))+ geom_density(color=for_plot2$variable)+theme_bw()+theme_classic()+
#  geom_vline(xintercept = 10,color="skyblue",size=0.4,linetype="dashed")+xlim(0,500)
ggplot(for_plot2, aes(x = value, color=variable))+geom_density()+theme_bw()+theme_classic()+
  geom_vline(xintercept = 20,color="skyblue",size=0.4,linetype="dashed")+xlim(0,500)
for_plot3 <- melt(test[,c(1:12,19:24)]) 
#ggplot(for_plot2, aes(x = value))+ geom_density(color=for_plot2$variable)+theme_bw()+theme_classic()+
#  geom_vline(xintercept = 10,color="skyblue",size=0.4,linetype="dashed")+xlim(0,500)
ggplot(for_plot3, aes(x = value, color=variable))+geom_density()+theme_bw()+theme_classic()+
  geom_vline(xintercept = 20,color="skyblue",size=0.4,linetype="dashed")+xlim(0,500)
for_plot4 <- melt(test) 
#ggplot(for_plot2, aes(x = value))+ geom_density(color=for_plot2$variable)+theme_bw()+theme_classic()+
#  geom_vline(xintercept = 10,color="skyblue",size=0.4,linetype="dashed")+xlim(0,500)
ggplot(for_plot4, aes(x = value, color=variable))+geom_density()+theme_bw()+theme_classic()+
  geom_vline(xintercept = 20,color="skyblue",size=0.4,linetype="dashed")+xlim(0,500)
########################new DNA figures###############################
DNA_total <-as.data.frame(rowMin(as.matrix(test[,c(11:16)],value=T)))
colnames(DNA_total)<-c("DNA_min_value")
pdf("DNA_counts_distribution_d24_new.pdf")
ggplot(DNA_total, aes(x = DNA_min_value))+ geom_histogram(binwidth = 5,color="skyblue")+theme_bw()+theme_classic()+
  geom_vline(xintercept = 20,color="skyblue",size=0.4,linetype="dashed")+xlim(0,500)
dev.off()

#ggplot(for_plot2, aes(x = value))+ geom_histogram(binwidth = 5,color="skyblue")+theme_bw()+theme_classic()+
#  geom_vline(xintercept = 20,color="skyblue",size=0.4,linetype="dashed")+xlim(0,300)
#head(for_plot)


####################################################################3

filter_norm_data <- norm_data[rowMins(as.matrix(test[,11:16]), value = T) >=20,21:32]

#######################################filter ratio#########################
total_number <- length(DNA_total[,1])
kept_number <- length(filter_norm_data[,1])
keptratio <- kept_number/total_number
DNA_counts <- data.frame(
  label=c("total_number","kept_number","keptratio"),
  value=c(as.character(total_number),as.character(kept_number),as.character(keptratio))
)

write.table(DNA_counts,"D24_DNA_distribution_filter_log.txt",sep="\t",quote=F)
###########################################################################


filter_for_degs<- data_sub[rowMins(as.matrix(test[,11:16]), value = T) >=20,21:32]
#filter_for_degs<- data_sub[rowMins(as.matrix(norm_data[,13:18]), value = T) >=20,13:24]
coldata <-as.data.frame(colnames(filter_for_degs))

coldata$type <- coldata$`colnames(filter_for_degs)`
coldata$condition <-gsub("\\.\\d","",coldata$type)
genes=as.matrix(filter_for_degs)
dds <- DESeqDataSetFromMatrix(countData = genes, colData=coldata, design=~condition)
dds <- DESeq(dds)
res <- results( dds, contrast = c("condition","RNA","DNA") )
#res2 <- results( dds, contrast = c("condition","d24_DNA_2","d24_RNA_2") )
res$qvalue_0_05 <- 0
res$qvalue_0_05[res$padj < 0.05 & res$log2FoldChange >0] <- 1
sum(res$qvalue_0_05)
dim(res)

write.table(res, file="D24_DES_result.txt",sep = "\t",quote = F)
#res2$qvalue_0_05 <- 0
#res2$qvalue_0_05[res2$padj < 0.05 & res2$log2FoldChange >0] <- 1
############################################################################

new <- data.frame(matrix(ncol = 0, nrow = length(rownames(filter_norm_data))))
rownames(new) <- rownames(filter_norm_data)
dim(filter_norm_data)
for (i in 1:6){
  j <- i +6
  ratio <- (filter_norm_data[,j]+1)/(filter_norm_data[,i]+1)
  ratio <- log2(ratio)
  new <- cbind(new,ratio)
}

new_colname <- paste0(colnames(filter_norm_data)[grep("RNA", colnames(filter_norm_data))], "_vs_DNA_log")
colnames(new) <- new_colname
head(new)
new$group <- data[as.character(rownames(new)),1]
head(new)
as.data.frame(colnames(new))
######################size factor#####################################
filter_norm_data <- as.data.frame(filter_norm_data)
d24_ratio <- as.data.frame(as.numeric(as.character(rowMeans(new[,1:6]))))

#filter_norm_data <- as.data.frame(cbind(filter_norm_data,d24_ratio))
filter_norm_data$log_d24_ratio <- d24_ratio[,1]
filter_norm_data$order_d24_ratio <- rank(-as.numeric(as.character(filter_norm_data$log_d24_ratio)))
#filter_norm_data$order_log_V_ratio <- rank(-as.numeric(as.character(filter_norm_data$d24_ratio)))
head(new)
new <-cbind(new,filter_norm_data)
new2 <- cbind(new,res$qvalue_0_05)
dim (new)
dim(res)
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
Negative_C <- filter_norm_data[which(new$group == "d0_"),]
both_C <- rownames(new)
head(new)
filter_norm_data_nc_normalize_exp <-MPRA_SizeFactors(new,both_C,c(1:6))
filter_norm_data_nc_normalize_exp <- as.data.frame(filter_norm_data_nc_normalize_exp)

A_V_cor <- cor(filter_norm_data_nc_normalize_exp)
A_V_cor2 <- cor(new[1:6])
corrplot(A_V_cor, method  = "number",cl.length = 21, tl.cex=0.5,hclust.method="median")
corrplot(A_V_cor2, method  = "number",tl.cex=0.5,hclust.method="median")
ggplot(filter_norm_data_nc_normalize_exp,aes(x=DNA.1,y=RNA.1))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(0,3000)+ylim(0,3000)
ggplot(new,aes(x=DNA.1,y=RNA.1))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(0,3000)+ylim(0,3000)
ggplot(filter_norm_data_nc_normalize_exp,aes(x=RNA1_vs_DNA_log,y=RNA3_vs_DNA_log))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(0,3)+ylim(0,3)
ggplot(filter_norm_data_nc_normalize_exp,aes(x=RNA1_vs_DNA_log,y=RNA4_vs_DNA_log))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(0,3)+ylim(0,3)
ggplot(filter_norm_data_nc_normalize_exp,aes(x=RNA2_vs_DNA_log,y=RNA3_vs_DNA_log))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(0,3)+ylim(0,3)
##############new pcc figures################################
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

pdf("D24_pcc_new5.pdf",width = 14, height = 14)
A_V_cor <- cor(filter_norm_data_nc_normalize_exp)
p<-corrplot(A_V_cor, method  = "number",order = "alphabet",type = 'lower')
pairs(filter_norm_data_nc_normalize_exp, 
      lower.panel = panel.cor,
      upper.panel = upper.panel)
dev.off()

#############################d17 and d24######################################
d24_filter_norm_data_nc_normalize_exp<-filter_norm_data_nc_normalize_exp
d17_filter_norm_data_nc_normalize_exp <- read.table(file="d17_batch1-2__normalized_log2_ratio_add1.txt",header=T,row.names = 1)
d17_d24_filter_norm_data_nc_normalize_exp<- merge(d17_filter_norm_data_nc_normalize_exp,d24_filter_norm_data_nc_normalize_exp,
                                                  by="row.names")
rownames(d17_d24_filter_norm_data_nc_normalize_exp)<-d17_d24_filter_norm_data_nc_normalize_exp$Row.names
d17_d24_filter_norm_data_nc_normalize_exp<-d17_d24_filter_norm_data_nc_normalize_exp[,c(2:length(d17_d24_filter_norm_data_nc_normalize_exp[1,]))]
colnames(d17_d24_filter_norm_data_nc_normalize_exp)<-c("D17_rep1","D17_rep2","D17_rep3","D17_rep4","D17_rep5","D17_rep6",
                                                       "D24_rep1","D24_rep2","D24_rep3","D24_rep4","D24_rep5","D24_rep6")
pdf("D17_D24_pcc_new5.pdf",width = 14, height = 14)

A_V_cor_merge <- cor(d17_d24_filter_norm_data_nc_normalize_exp)
p<-corrplot(A_V_cor_merge, method  = "number",order = "alphabet",type = 'lower')
pheatmap(p$corr,cluster_cols=F, cluster_rows=F,col=colorRampPalette(c("white","firebrick3"))(100),breaks=seq(0.6,1,0.004))
#corrplot(A_V_cor_merge,method="color",addgrid.col ="White",
#         col.lim=c(0,1), addCoef.col = 'black',col=colorRampPalette(c("blue","white","red"))(200))
pairs(d17_d24_filter_norm_data_nc_normalize_exp, 
      lower.panel = panel.cor,
      upper.panel = upper.panel)
dev.off()
######################## merged keep all ###############################

d17_filter_norm_data_nc_normalize_exp_new<-read.table(file="d17_batch1-2_log2_ratio_add1.txt",header=T,row.names = 1)
d24_filter_norm_data_nc_normalize_exp_new <- read.table(file="d24_batch1-2_log2_ratio_add1.txt",header=T,row.names = 1)
d17_d24_filter_norm_data_nc_normalize_exp_new<- merge(d17_filter_norm_data_nc_normalize_exp_new,d24_filter_norm_data_nc_normalize_exp_new,
                                                  by="row.names",all=T)
#d17_d24_filter_norm_data_nc_normalize_exp_new<- merge(d17_filter_norm_data_nc_normalize_exp_new,d24_filter_norm_data_nc_normalize_exp_new,
#                                                      by="row.names")
d17_d24_filter_norm_data_nc_normalize_exp_new<-d17_d24_filter_norm_data_nc_normalize_exp_new[!is.na(d17_d24_filter_norm_data_nc_normalize_exp_new$Row.names),]
d17_d24_filter_norm_data_nc_normalize_exp_new$group<-0
row.names(d17_d24_filter_norm_data_nc_normalize_exp_new)<-d17_d24_filter_norm_data_nc_normalize_exp_new$Row.names
d17_d24_filter_norm_data_nc_normalize_exp_new[rownames(d17_filter_norm_data_nc_normalize_exp_new),]$group<-d17_d24_filter_norm_data_nc_normalize_exp_new[rownames(d17_filter_norm_data_nc_normalize_exp_new),]$group.x
d17_d24_filter_norm_data_nc_normalize_exp_new[rownames(d24_filter_norm_data_nc_normalize_exp_new),]$group<-d17_d24_filter_norm_data_nc_normalize_exp_new[rownames(d24_filter_norm_data_nc_normalize_exp_new),]$group.y

head(d17_d24_filter_norm_data_nc_normalize_exp_new)
table(d17_d24_filter_norm_data_nc_normalize_exp_new[(d17_d24_filter_norm_data_nc_normalize_exp_new$res.qvalue_0_05.x==1 | d17_d24_filter_norm_data_nc_normalize_exp_new$res.qvalue_0_05.y==1),c("group")])

table(d17_d24_filter_norm_data_nc_normalize_exp_new[d17_d24_filter_norm_data_nc_normalize_exp_new$res.qvalue_0_05.x==1,c("group")])
table(d17_d24_filter_norm_data_nc_normalize_exp_new[,c("group.y")])
table(d17_filter_norm_data_nc_normalize_exp_new[,c("group")])

###################################################################

write.table(new2,file="d24_batch1-2_log2_ratio_add1.txt",sep="\t")
write.table(filter_norm_data,file="d24_batch1-2_log2_FPM_add1.txt",sep="\t")
write.table(filter_norm_data_nc_normalize_exp,file="d24_batch1-2__normalized_log2_ratio_add1.txt",sep="\t")
write.table(d17_d24_filter_norm_data_nc_normalize_exp,file="d17-d24_batch1-2__normalized_log2_ratio_add1.txt",sep="\t")
#################################rank figure#####################################
filter_norm_data_nc_normalize_exp <- read.table(file="d24_batch1-2__normalized_log2_ratio_add1.txt",header=T,row.names = 1)
new2<-read.table(file="d24_batch1-2_log2_ratio_add1.txt",header=T,row.names = 1)
#filter_norm_data_scaled <- dplyr::rename(new2, d24_active=`res$qvalue_0_05`)
filter_norm_data_scaled <- dplyr::rename(new2, d24_active=`res.qvalue_0_05`)
#filter_norm_data_scaled <- cbind(new[,c(1:7,20:22)],filter_norm_data_nc_normalize_exp)
#filter_norm_data_scaled <- cbind(new2[,c(1:7,20:23)],filter_norm_data_nc_normalize_exp)
filter_norm_data_scaled<-as.data.frame(filter_norm_data_scaled)
filter_norm_data_scaled$rank_value <- filter_norm_data_scaled$order_d24_ratio

filter_norm_data_scaled_grouped <- as.data.frame(filter_norm_data_scaled)
select_NC <- filter_norm_data_scaled_grouped[(filter_norm_data_scaled_grouped$group=="d0_noncoding"| filter_norm_data_scaled_grouped$group=="d0_exon_nonHHE"),]
filter_norm_data_scaled_grouped$group2<-"HHE_ATAC+"
filter_norm_data_scaled_grouped[rownames(select_NC),]$group2<-"Neg"
##########################FDR 0.05################################################
fitg<-fitdist(select_NC$log_d24_ratio ,"norm")
head(select_NC)
FDR_0.95 <- qnorm(0.95, mean = fitg$estimate["mean"], sd = fitg$estimate["sd"])
FDR_0.95
ggplot(filter_norm_data_scaled_grouped, aes(x = log_d24_ratio, color=group))+geom_density()+theme_bw()+theme_classic()
ggplot(filter_norm_data_scaled_grouped, aes(x = order_d24_ratio, color=group))+geom_density()+theme_bw()+theme_classic()
library(patchwork)
pdf("D24_library_figures.pdf", width = 11, height = 8)
A<-ggplot(filter_norm_data_scaled,aes(x=rank_value,y=log_d24_ratio,color=factor(d24_active)))+
  geom_point(size=2)+
  theme_bw()+
  scale_color_manual(breaks = c("1", "0"), values=c("red", "gray"))+
  theme_classic()

B<-ggplot(filter_norm_data_scaled_grouped, aes(x = order_d24_ratio, y = factor(group,levels=c("d0_exon_nonHHE","d0_noncoding",
                                                                                                           "d8_d17_noncoding_HHE "))))+
  stat_bin2d(aes(fill = after_stat(count)), binwidth = c(3,1))+theme_bw()
A/B
dev.off()
pdf("D24_library_figures_2.pdf", width = 11, height = 8)
A<-ggplot(filter_norm_data_scaled,aes(x=rank_value,y=log_d24_ratio,color=factor(d24_active)))+
  geom_point(size=2)+
  theme_bw()+
  scale_color_manual(breaks = c("1", "0"), values=c("red", "gray"))+
  theme_classic()

B<-ggplot(filter_norm_data_scaled_grouped, aes(x = order_d24_ratio, y = factor(group2,levels=c("Neg",
                                                                                              "HHE_ATAC+"))))+
  stat_bin2d(aes(fill = after_stat(count)), binwidth = c(3,1))+theme_bw()
A/B
dev.off()
##############################################
first_col<-c("GBE1",
             "PTPN1",
             "COL5A1",
             "SALL1",
             "PKP2",
             "COL1A2",
             "TGFBR1",
             "NR2F2",
             "NOTCH2",
             "SHROOM3")
second_col<-c("chr3:82149507-82149907",
"chr20:49157407-49157807",
"chr9:137494462-137494862",
"chr16:51046302-51046702",
"chr12:33136295-33136695",
"chr7:93870775-93871175",
"chr9:101837589-101837989",
"chr15:96808175-96808575",
"chr1:120557155-120557555",
"chr4:77586893-77587293")
highlight <- data.frame(first_col, second_col)
rownames(highlight)<-highlight$second_col
head(filter_norm_data_scaled)
filter_norm_data_scaled$label<-NA
filter_norm_data_scaled[rownames(highlight),]$label<-highlight$first_col
filter_norm_data_scaled$have<-1
filter_norm_data_scaled[rownames(highlight),]$have<-2
library(ggrepel)
pdf("D24_library_figures_10_top.pdf", width = 8, height = 8,useDingbats = F)
ggplot(filter_norm_data_scaled,aes(x=rank_value,y=log_d24_ratio,color=factor(d24_active),label=label,size=factor(have)))+
  geom_point()+
  theme_bw()+
  scale_color_manual(breaks = c("1", "0"), values=c("red", "gray"))+
  scale_size_manual(breaks = c("1", "2"), values=c(1, 3))+
  geom_text_repel(aes(label=label,color=factor(d24_active)), nudge_x=20,direction="x",color="red",size=3)+ 
  theme_classic()
ggplot(filter_norm_data_scaled,aes(x=rank_value,y=log_d24_ratio,label=label,size=factor(have)))+
  geom_point(color="darkblue")+
  theme_bw()+
  #scale_color_manual(breaks = c("1", "0"), values=c("red", "gray"))+
  scale_size_manual(breaks = c("1", "2"), values=c(1, 3))+
  geom_text_repel(aes(label=label),nudge_x=20,direction="x",color="red",size=3)+ 
  theme_classic()
dev.off()
############################################3
ggplot(filter_norm_data_scaled_grouped,aes(x=rank_value,y=log_d24_ratio,color=factor(d24_active)))+geom_point(size=2)+
  geom_hline(yintercept=FDR_0.95,color="red",linetype="dashed")+theme_bw()+
  scale_color_manual(breaks = c("1", "0"), values=c("red", "darkblue"))+
  theme_classic()+geom_text(aes(0,round(FDR_0.95,3),label = round(FDR_0.95,3), vjust = -1))
ggplot(filter_norm_data_scaled_grouped,aes(x=order_d24_ratio,y=log_d24_ratio,color=factor(d24_active)))+geom_point(size=2)+
  theme_bw()+
  scale_color_manual(breaks = c("1", "0"), values=c("red", "darkblue"))+
  theme_classic()
ggplot(filter_norm_data_scaled_grouped, aes(x = order_d24_ratio, y = group))+stat_bin2d(aes(fill = after_stat(count)), binwidth = c(3,1))
ggplot(filter_norm_data_scaled_grouped, aes(x = order_d24_ratio, y = group))+stat_bin2d(aes(fill = after_stat(count)), binwidth = c(1,1))
dev.off()
sum(filter_norm_data_scaled_grouped$log_d24_ratio>=FDR_0.95)
W1<-as.data.frame(table(filter_norm_data_scaled_grouped$d24_active))

colnames(W1) <-c("D24_active","number")
write.table(W1, "D24_active_regions_number.txt",quote=F, sep = "\t")
##############################################################################
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
################################3A###########################
#filter_norm_data_scaled_grouped_test_A <-Enrich_score(filter_norm_data_scaled_grouped, "group", "log_A_ratio", "order_log_A_mean_ratio")
filter_norm_data_scaled_grouped_test <-Enrich_score(filter_norm_data_scaled_grouped, "group2", "log_d24_ratio", "order_d24_ratio")
dat_text_A <-filter_norm_data_scaled_grouped_test$dat_text
filter_norm_data_for_plot_A <- filter_norm_data_scaled_grouped_test$datasets
#pdf("A_MTF_2_ES.pdf")

ggplot(filter_norm_data_scaled_grouped, aes(x = order_d24_ratio, y = factor(group2,levels=c("Neg",
                                                                                            "HHE_ATAC+"))))+
  #stat_bin2d(aes(fill = ..density..),bins = 200)+
  #scale_fill_gradient(low = "light blue", high = "dark red")+
  geom_bin2d(aes(alpha = ..count..), fill = "slategray",binwidth = c(5,1)) +
  scale_fill_continuous(guide = "none")+
  #scale_fill_gradient(low = "light blue", high = "dark red")+
  #stat_bin2d(aes(fill = after_stat(count)), binwidth = c(3,1))+
  theme_bw()
ggplot(filter_norm_data_scaled_grouped, aes(x = order_d24_ratio, y = factor(group2,levels=c("Neg",
                                                                                            "HHE_ATAC+"))))+
  #stat_bin2d(aes(fill = ..density..),bins = 200)+
  #scale_fill_gradient(low = "light blue", high = "dark red")+
  geom_bin2d(aes(alpha = ..density..), fill = "slategray",binwidth = c(5,1)) +
  scale_fill_continuous(guide = "none")+
  #scale_fill_gradient(low = "light blue", high = "dark red")+
  #stat_bin2d(aes(fill = after_stat(count)), binwidth = c(3,1))+
  theme_bw()
pdf("d24_ES.pdf")
p_A<-ggplot(filter_norm_data_for_plot_A, aes(order_d24_ratio, final_ES,color=factor(score_group)))+geom_line()+theme_bw()+theme_classic()+
  #facet_wrap(~score_group,scales = "free",nrow =8)+ 
  theme(panel.grid.major.y = element_blank(),
                                                          panel.grid.minor.y = element_blank(),
                                                          panel.grid.major.x = element_blank(),
                                                          panel.grid.minor.x = element_blank())+
  scale_color_manual(breaks = c("HHE_ATAC+", "Neg"), values=c("darkred", "darkblue"))+
  geom_hline(yintercept = 0,color="skyblue",size=0.4,linetype="dashed")+
  geom_text(
  data    = dat_text_A,
  mapping = aes(x = -Inf, y = -Inf, label = label),
  hjust   = -2,
  vjust   = -0.25) 
A<-ggplot(filter_norm_data_scaled,aes(x=rank_value,y=log_d24_ratio,color=factor(d24_active)))+
  geom_point(size=2)+
  theme_bw()+
  scale_color_manual(breaks = c("1", "0"), values=c("red", "gray"))+
  theme_classic()

B<-ggplot(filter_norm_data_scaled_grouped, aes(x = order_d24_ratio, y = factor(group2,levels=c("Neg",
                                                                                               "HHE_ATAC+"))))+
  
  geom_bin2d(aes(alpha = ..density..), fill = "slategray",binwidth = c(5,1)) +
  scale_fill_continuous(guide = "none")+

  theme_bw()
A/p_A/B
dev.off()
