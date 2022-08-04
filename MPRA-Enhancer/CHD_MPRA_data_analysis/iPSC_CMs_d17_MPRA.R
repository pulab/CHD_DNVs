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
#data<-read.table(file="format_table.txt",header=T,row.names = 1)
data<-read.table(file="batch1-2_rawcounts_d17_with_merge.txt",header=T,row.names = 1)
######################to FPM#####################################
dim(data)
data_sub <- as.matrix(data[,2:37])
norm_data <- t(t(data_sub)/colSums(data_sub)*1000000)
total_sum <-as.data.frame(colSums(data_sub))

ggplot(total_sum, aes(x=row.names(total_sum),y=colSums(data_sub)))+ geom_bar(color="skyblue",stat="identity")+theme_bw()+theme_classic()
#filter_norm_data <- cbind(norm_data, data[,6:7]) %>% dplyr::filter(DNA >= 5 &rowSums(.[,1:4]>0))
test <- dplyr::select(as.data.frame(norm_data), contains("DNA"))
#filter_norm_data <- norm_data[rowMins(as.matrix(test), value = T) >=5,]
filter_norm_data <- test
dim(filter_norm_data)
head(test)


for_plot <- as.data.frame(test) 
ggplot(for_plot, aes(x = DNA.4))+ geom_histogram(binwidth = 5,color="skyblue")+theme_bw()+theme_classic()+
  geom_vline(xintercept = 20,color="skyblue",size=0.4,linetype="dashed")+xlim(0,500)
ggplot(for_plot, aes(x = d17_DNA.4_2))+ geom_histogram(binwidth = 5,color="skyblue")+theme_bw()+theme_classic()+
  geom_vline(xintercept = 20,color="skyblue",size=0.4,linetype="dashed")+xlim(0,500)
##############################distribution of DNA####################333
#dim(test)
#for_plot2 <- melt(test[,c(7:12,19:24)]) 
##ggplot(for_plot2, aes(x = value))+ geom_density(color=for_plot2$variable)+theme_bw()+theme_classic()+
##  geom_vline(xintercept = 10,color="skyblue",size=0.4,linetype="dashed")+xlim(0,500)
#ggplot(for_plot2, aes(x = value, color=variable))+geom_density()+theme_bw()+theme_classic()+
#  geom_vline(xintercept = 20,color="skyblue",size=0.4,linetype="dashed")+xlim(0,500)
#for_plot3 <- melt(test[,c(1:12,19:24)]) 
##ggplot(for_plot2, aes(x = value))+ geom_density(color=for_plot2$variable)+theme_bw()+theme_classic()+
##  geom_vline(xintercept = 10,color="skyblue",size=0.4,linetype="dashed")+xlim(0,500)
#ggplot(for_plot3, aes(x = value, color=variable))+geom_density()+theme_bw()+theme_classic()+
#  geom_vline(xintercept = 20,color="skyblue",size=0.4,linetype="dashed")+xlim(0,500)
#for_plot4 <- melt(test) 
#ggplot(for_plot2, aes(x = value))+ geom_density(color=for_plot2$variable)+theme_bw()+theme_classic()+
#  geom_vline(xintercept = 10,color="skyblue",size=0.4,linetype="dashed")+xlim(0,500)
#ggplot(for_plot4, aes(x = value, color=variable))+geom_density()+theme_bw()+theme_classic()+
#  geom_vline(xintercept = 20,color="skyblue",size=0.4,linetype="dashed")+xlim(0,500)
####################################################################
########################new DNA figures###############################
DNA_total <-as.data.frame(rowMin(as.matrix(test[,c(13:18)],value=T)))
colnames(DNA_total)<-c("DNA_min_value")
pdf("DNA_counts_distribution_D17.pdf")
ggplot(DNA_total, aes(x = DNA_min_value))+ geom_histogram(binwidth = 5,color="skyblue")+theme_bw()+theme_classic()+
  geom_vline(xintercept = 20,color="skyblue",size=0.4,linetype="dashed")+xlim(0,500)
dev.off()

####################################################################
filter_norm_data <- norm_data[rowMins(as.matrix(test[,13:18]), value = T) >=20,25:36]
#######################################filter ratio#########################
total_number <- length(DNA_total[,1])
kept_number <- length(filter_norm_data[,1])
keptratio <- kept_number/total_number
DNA_counts <- data.frame(
  label=c("total_number","kept_number","keptratio"),
  value=c(as.character(total_number),as.character(kept_number),as.character(keptratio))
)

write.table(DNA_counts,"D17_DNA_distribution_filter_log.txt",sep="\t",quote=F)
##############################################################################3
filter_for_degs<- data_sub[rowMins(as.matrix(test[,13:18]), value = T) >=20,25:36]
#filter_for_degs<- data_sub[rowMins(as.matrix(norm_data[,13:18]), value = T) >=20,13:24]
coldata <-as.data.frame(colnames(filter_for_degs))

coldata$type <- coldata$`colnames(filter_for_degs)`
coldata$condition <-gsub("\\.\\d","",coldata$type)
genes=as.matrix(filter_for_degs)
dds <- DESeqDataSetFromMatrix(countData = genes, colData=coldata, design=~condition)
dds <- DESeq(dds)
res <- results( dds, contrast = c("condition","RNA","DNA") )
#res2 <- results( dds, contrast = c("condition","d17_DNA_2","d17_RNA_2") )
res$qvalue_0_05 <- 0
res$qvalue_0_05[res$padj < 0.05 & res$log2FoldChange >0] <- 1
sum(res$qvalue_0_05)
dim(res)

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
d17_ratio <- as.data.frame(as.numeric(as.character(rowMeans(new[,1:6]))))

#filter_norm_data <- as.data.frame(cbind(filter_norm_data,d24_ratio))
filter_norm_data$log_d17_ratio <- d17_ratio[,1]
filter_norm_data$order_d17_ratio <- rank(-as.numeric(as.character(filter_norm_data$log_d17_ratio)))
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
#both_C <- rownames(rbind(Negative_C,Positive_C))

#d17_ratio <- as.numeric(as.character(rowMeans(filter_norm_data[,7:12])/rowMeans(filter_norm_data[,1:6])))

#filter_norm_data <- as.data.frame(cbind(filter_norm_data,d17_ratio))
#filter_norm_data$log_d17_ratio <- log2(filter_norm_data$d17_ratio+1)
#filter_norm_data$order_d17_ratio <- rank(-as.numeric(as.character(filter_norm_data$d17_ratio)))
#filter_norm_data$order_log_V_ratio <- rank(-as.numeric(as.character(filter_norm_data$d17_ratio)))
head(new)
new <-cbind(new,filter_norm_data)
new2 <- cbind(new,res$qvalue_0_05)
dim (new)
dim(res)
#Negative_C <- filter_norm_data[which(new$group == "d0_noncoding"),]
#Positive_C <- t(dplyr::select(as.data.frame(t(new)), contains("::F")))
#Positive_C <-new
#filter_norm_data_nc <- Positive_C[,c(8:19)] %>% as.matrix()
#filter_norm_data_full_nc <- new[,c(8:19)] %>% as.matrix()
#filter_norm_data_nc <- filter_norm_data_nc[rowMins(filter_norm_data_nc,value = T)>0,]

filter_norm_data_nc_normalize_exp <-MPRA_SizeFactors(new,both_C,c(1:6))
filter_norm_data_nc_normalize_exp <- as.data.frame(filter_norm_data_nc_normalize_exp)

#filter_norm_data_nc_ratio_pseudo_exp <- geometricmeanRow(filter_norm_data_nc)
#filter_norm_data_nc_ratio_scaled <- t(t(filter_norm_data_nc)/filter_norm_data_nc_ratio_pseudo_exp)

#filter_norm_data_nc_normalize_factor <- colMedians(filter_norm_data_nc_ratio_scaled)
#filter_norm_data_nc_normalize_exp <- t(t(new[,8:19])/filter_norm_data_nc_normalize_factor)
#def sizeFactor()
#filter_norm_data_scaled <- cbind(filter_norm_data[,c(5:7)], filter_norm_data_nc_normalize_exp)
#filter_norm_data_scaled <- cbind(new[,c(1:4)], filter_norm_data_nc_normalize_exp)
#filter_norm_data_scaled$A_ratio <- (filter_norm_data_scaled$A1 + filter_norm_data_scaled$A2)/2/filter_norm_data_scaled$DNA
#filter_norm_data_scaled$V_ratio <- (filter_norm_data_scaled$V1 + filter_norm_data_scaled$V2)/2/filter_norm_data_scaled$DNA
#filter_norm_data_nc_normalize_exp <- as.data.frame(filter_norm_data_nc_normalize_exp)
A_V_cor <- cor(filter_norm_data_nc_normalize_exp)
A_V_cor2 <- cor(new[1:6])
corrplot(A_V_cor2, method  = "number",cl.length = 21, tl.cex=0.5,hclust.method="median")
pheatmap(A_V_cor2,cluster_rows=FALSE,cluster_cols=FALSE,color = colorRampPalette(c("white","white","firebrick3"))(200),breaks=seq(0,1,0.005),fontsize=9,display_numbers=T)
ggplot(filter_norm_data_nc_normalize_exp,aes(x=DNA.1,y=RNA.1))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(0,3000)+ylim(0,3000)
ggplot(new,aes(x=DNA.1,y=RNA.1))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(0,3000)+ylim(0,3000)
ggplot(filter_norm_data_nc_normalize_exp,aes(x=RNA1_vs_DNA_log,y=RNA3_vs_DNA_log))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(0,3)+ylim(0,3)
ggplot(filter_norm_data_nc_normalize_exp,aes(x=RNA1_vs_DNA_log,y=RNA4_vs_DNA_log))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(0,3)+ylim(0,3)
ggplot(filter_norm_data_nc_normalize_exp,aes(x=RNA2_vs_DNA_log,y=RNA3_vs_DNA_log))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(0,3)+ylim(0,3)
#######################
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

pdf("D17_pcc.pdf",width = 14, height = 14)
A_V_cor <- cor(filter_norm_data_nc_normalize_exp)
p<-corrplot(A_V_cor, method  = "number",order = "alphabet",type = 'lower')
pairs(filter_norm_data_nc_normalize_exp, 
      lower.panel = panel.cor,
      upper.panel = upper.panel)
dev.off()



d17_filter_norm_data_nc_normalize_exp <-filter_norm_data_nc_normalize_exp

#write.table(new,file="d17_batch1-2_log2_ratio.txt",sep="\t")
write.table(new2,file="d17_batch1-2_log2_ratio_add1.txt",sep="\t")
write.table(filter_norm_data,file="d17_batch1-2_log2_FPM_add1.txt",sep="\t")
write.table(filter_norm_data_nc_normalize_exp,file="d17_batch1-2__normalized_log2_ratio_add1.txt",sep="\t")

#################################rank figure#####################################
filter_norm_data_nc_normalize_exp <- read.table(file="d17_batch1-2__normalized_log2_ratio_add1.txt",header=T,row.names = 1)
new2<-read.table("d17_batch1-2_log2_ratio_add1.txt",header = T,row.names = 1)
#new <- rename(new2, d17_active=`res$qvalue_0_05`)
filter_norm_data_scaled <- dplyr::rename(new2, d17_active=`res$qvalue_0_05`)
filter_norm_data_scaled <- dplyr::rename(new2, d17_active=`res.qvalue_0_05`)
#filter_norm_data_scaled <- cbind(new[,c(1:7,20:23)],filter_norm_data_nc_normalize_exp)
filter_norm_data_scaled<-as.data.frame(filter_norm_data_scaled)
filter_norm_data_scaled$rank_value <- filter_norm_data_scaled$order_d17_ratio

filter_norm_data_scaled_grouped <- as.data.frame(filter_norm_data_scaled)
#select_NC <- filter_norm_data_scaled_grouped[filter_norm_data_scaled_grouped$group=="d0_noncoding",]
select_NC <- filter_norm_data_scaled_grouped[(filter_norm_data_scaled_grouped$group=="d0_noncoding"| filter_norm_data_scaled_grouped$group=="d0_exon_nonHHE"),]
filter_norm_data_scaled_grouped$group2<-"HHE_ATAC+"
filter_norm_data_scaled_grouped[rownames(select_NC),]$group2<-"Neg"
##########################FDR 0.05################################################
fitg<-fitdist(select_NC$log_d17_ratio ,"norm")
head(select_NC)
FDR_0.95 <- qnorm(0.95, mean = fitg$estimate["mean"], sd = fitg$estimate["sd"])
FDR_0.95
ggplot(filter_norm_data_scaled_grouped, aes(x = log_d17_ratio, color=group))+geom_density()+theme_bw()+theme_classic()

pdf("D17_library_figures.pdf", width = 11, height = 8)
A<-ggplot(filter_norm_data_scaled,aes(x=rank_value,y=log_d17_ratio,color=factor(d17_active)))+
  geom_point(size=2)+
  theme_bw()+
  scale_color_manual(breaks = c("1", "0"), values=c("red", "gray"))+
  theme_classic()

B<-ggplot(filter_norm_data_scaled_grouped, aes(x = order_d17_ratio, y = factor(group,levels=c("d0_exon_nonHHE","d0_noncoding",
                                                                                              "d8_d17_noncoding_HHE "))))+
  stat_bin2d(aes(fill = after_stat(count)), binwidth = c(3,1))+theme_bw()
A/B
dev.off()
pdf("D17_library_figures_2.pdf", width = 11, height = 8)
A<-ggplot(filter_norm_data_scaled,aes(x=rank_value,y=log_d17_ratio,color=factor(d17_active)))+
  geom_point(size=2)+
  theme_bw()+
  scale_color_manual(breaks = c("1", "0"), values=c("red", "gray"))+
  theme_classic()

B<-ggplot(filter_norm_data_scaled_grouped, aes(x = order_d17_ratio, y = factor(group2,levels=c("Neg",
                                                                                               "HHE_ATAC+"))))+
  stat_bin2d(aes(fill = after_stat(count)), binwidth = c(3,1))+theme_bw()
A/B
dev.off()

ggplot(filter_norm_data_scaled_grouped,aes(x=order_d17_ratio,y=log_d17_ratio,color=factor(d17_active)))+geom_point(size=2)+
  geom_hline(yintercept=FDR_0.95,color="red",linetype="dashed")+theme_bw()+
  scale_color_manual(breaks = c("1", "0"), values=c("red", "darkblue"))+
  theme_classic()+geom_text(aes(0,round(FDR_0.95,3),label = round(FDR_0.95,3), vjust = -1))
ggplot(filter_norm_data_scaled_grouped,aes(x=order_d17_ratio,y=log_d17_ratio,color=factor(d17_active)))+geom_point(size=2)+
  theme_bw()+
  scale_color_manual(breaks = c("1", "0"), values=c("red", "darkblue"))+
  theme_classic()
ggplot(filter_norm_data_scaled_grouped, aes(x = order_d17_ratio, y = group))+stat_bin2d(aes(fill = after_stat(count)), binwidth = c(3,1))
ggplot(filter_norm_data_scaled_grouped, aes(x = order_d17_ratio, y = group))+stat_bin2d(aes(fill = after_stat(count)), binwidth = c(1,1))
dev.off()
sum(filter_norm_data_scaled_grouped$log_d17_ratio>=FDR_0.95)

ggplot(filter_norm_data_scaled_grouped,aes(y= log_d17_ratio, x= log2(DNA.1)))+geom_point()
W1<-as.data.frame(table(filter_norm_data_scaled_grouped$d17_active))

colnames(W1) <-c("D17_active","number")
write.table(W1, "D17_active_regions_number.txt",quote=F, sep = "\t")
##########################################################################
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
filter_norm_data_scaled_grouped_test <-Enrich_score(filter_norm_data_scaled_grouped, "group2", "log_d17_ratio", "order_d17_ratio")
dat_text_A <-filter_norm_data_scaled_grouped_test$dat_text
filter_norm_data_for_plot_A <- filter_norm_data_scaled_grouped_test$datasets
#pdf("A_MTF_2_ES.pdf")

ggplot(filter_norm_data_scaled_grouped, aes(x = order_d17_ratio, y = factor(group2,levels=c("Neg",
                                                                                            "HHE_ATAC+"))))+
  #stat_bin2d(aes(fill = ..density..),bins = 200)+
  #scale_fill_gradient(low = "light blue", high = "dark red")+
  geom_bin2d(aes(alpha = ..count..), fill = "slategray",binwidth = c(5,1)) +
  scale_fill_continuous(guide = "none")+
  #scale_fill_gradient(low = "light blue", high = "dark red")+
  #stat_bin2d(aes(fill = after_stat(count)), binwidth = c(3,1))+
  theme_bw()
ggplot(filter_norm_data_scaled_grouped, aes(x = order_d17_ratio, y = factor(group2,levels=c("Neg",
                                                                                            "HHE_ATAC+"))))+
  #stat_bin2d(aes(fill = ..density..),bins = 200)+
  #scale_fill_gradient(low = "light blue", high = "dark red")+
  geom_bin2d(aes(alpha = ..density..), fill = "slategray",binwidth = c(5,1)) +
  scale_fill_continuous(guide = "none")+
  #scale_fill_gradient(low = "light blue", high = "dark red")+
  #stat_bin2d(aes(fill = after_stat(count)), binwidth = c(3,1))+
  theme_bw()
pdf("A_MTF_2_ES_d17.pdf")
p_A<-ggplot(filter_norm_data_for_plot_A, aes(order_d17_ratio, final_ES,color=factor(score_group)))+geom_line()+theme_bw()+theme_classic()+
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
A<-ggplot(filter_norm_data_scaled,aes(x=rank_value,y=log_d17_ratio,color=factor(d17_active)))+
  geom_point(size=2)+
  theme_bw()+
  scale_color_manual(breaks = c("1", "0"), values=c("red", "gray"))+
  theme_classic()

B<-ggplot(filter_norm_data_scaled_grouped, aes(x = order_d17_ratio, y = factor(group2,levels=c("Neg",
                                                                                               "HHE_ATAC+"))))+
  
  geom_bin2d(aes(alpha = ..density..), fill = "slategray",binwidth = c(5,1)) +
  scale_fill_continuous(guide = "none")+
  
  theme_bw()
A/p_A/B
dev.off()
