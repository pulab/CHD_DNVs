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
library(stringr)
#library(Rbase)
data<-read.table(file="CHD_variane_raw_counts_DNA1-4_RNA1-4_with_footprint.txt",header=T,row.names = 1)
#data<-read.table(file="Mutagenesis_raw_counts_DNA1-4_RNA1-4.txt",header=T,row.names = 1)
######################to FPM#####################################
data_sub <- as.matrix(data[,1:8])
head(data_sub)
#setClass("MPRA_data",slots=list(readcount="numberic",age="numeric"))
norm_data <- t(t(data_sub)/colSums(data_sub)*1000000)
#filter_norm_data <- cbind(norm_data, data[,6:7]) %>% dplyr::filter(DNA >= 5 &rowSums(.[,1:4]>0))
test <- dplyr::select(as.data.frame(norm_data), contains("DNA"))
#filter_norm_data <- norm_data[rowMins(as.matrix(test), value = T) >=15,]
filter_norm_data <- norm_data[rowMins(as.matrix(test), value = T) >=0,]
head(test)


for_plot <- as.data.frame(test) 
#head(for_plot)
for_plot2 <- melt(for_plot)

DNA_total <-as.data.frame(rowMin(norm_data[,c(1:4)]))
colnames(DNA_total)<-c("DNA_min_value")
pdf("DNA_counts_distribution_CHD.pdf")
ggplot(DNA_total, aes(x = DNA_min_value))+ geom_histogram(binwidth = 5,color="skyblue")+theme_bw()+theme_classic()+
  geom_vline(xintercept = 20,color="skyblue",size=0.4,linetype="dashed")+xlim(0,500)
dev.off()

ggplot(for_plot2, aes(x = value))+ geom_histogram(binwidth = 5,color="skyblue")+theme_bw()+theme_classic()+
  geom_vline(xintercept = 20,color="skyblue",size=0.4,linetype="dashed")+xlim(0,300)
head(for_plot)
#######################################################################
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

filter_for_degs<- data_sub[rownames(filter_norm_data),]

#filter_for_degs<-filter_norm_data
#filter_for_degs<- data_sub[rowMins(as.matrix(norm_data[,13:18]), value = T) >=20,13:24]
coldata <-as.data.frame(colnames(filter_for_degs))

coldata$type <- coldata$`colnames(filter_for_degs)`
coldata$condition <-gsub("\\d","",coldata$type)
res<-Active_enhancer(filter_for_degs,coldata,"condition")
#genes=as.matrix(filter_for_degs)
#dds <- DESeqDataSetFromMatrix(countData = genes, colData=coldata, design=~condition)
#dds <- DESeq(dds)
#res <- results( dds, contrast = c("condition","RNA","DNA") )
##res2 <- results( dds, contrast = c("condition","d17_DNA_2","d17_RNA_2") )
#res$qvalue_0_05 <- 0
#res$qvalue_0_05[res$padj < 0.05 & res$log2FoldChange >0] <- 1
sum(res$qvalue_0_05)
dim(res)
write.table(res, file="CHD_DES_result.txt",sep = "\t",quote = F)
#######################################################################
filter_norm_data <- norm_data[rowMins(as.matrix(test), value = T) >=20,]
total_number <- length(DNA_total[,1])
kept_number <- length(filter_norm_data[,1])
keptratio <- kept_number/total_number
DNA_counts <- data.frame(
  label=c("total_number","kept_number","keptratio"),
  value=c(as.character(total_number),as.character(kept_number),as.character(keptratio))
)

write.table(DNA_counts,"CHD_library_DNA_distribution_filter_log.txt",sep="\t",quote=F)
###########################################################################
new <- data.frame(matrix(ncol = 0, nrow = length(rownames(filter_norm_data))))
rownames(new) <- rownames(filter_norm_data)
#for (i in 1:4){
#  j <- i +4
#  ratio <- filter_norm_data[,j]/filter_norm_data[,i]
#  ratio <- log2(ratio+1)
#  new <- cbind(new,ratio)
#}
for (i in 1:4){
  j <- i +4
  ratio <- (filter_norm_data[,j]+1)/(filter_norm_data[,i]+1)
  ratio <- log2(ratio)
  new <- cbind(new,ratio)
}
new_colname <- paste0(colnames(filter_norm_data)[grep("RNA", colnames(filter_norm_data))], "_vs_DNA_log")
colnames(new) <- new_colname
head(res)
new$qvalue_0_05 <- res[rownames(new),]$qvalue_0_05
head(new)
head(new2)
###########################DESeq2#############################

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

#corrplot(A_V_cor, method  = "number",tl.cex=0.5,hclust.method="median")
######################################################################################

data_CHD<-filter_norm_data_nc_normalize_exp[both_C,]
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

pdf("MPRA_CHD_pcc.pdf",width = 11, height = 11)
A_V_cor <- cor(filter_norm_data_nc_normalize_exp)
p<-corrplot(A_V_cor, method  = "number",order = "alphabet",type = 'lower')
pairs(filter_norm_data_nc_normalize_exp, 
      lower.panel = panel.cor,
      upper.panel = upper.panel)
dev.off()
#########################################################################
ggplot(filter_norm_data_nc_normalize_exp,aes(x=RNA1_vs_DNA_log,y=RNA2_vs_DNA_log))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(0,3)+ylim(0,3)
ggplot(new,aes(x=RNA1_vs_DNA_log,y=RNA2_vs_DNA_log))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(0,3)+ylim(0,3)
ggplot(filter_norm_data_nc_normalize_exp,aes(x=RNA1_vs_DNA_log,y=RNA3_vs_DNA_log))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(0,3)+ylim(0,3)
ggplot(filter_norm_data_nc_normalize_exp,aes(x=RNA1_vs_DNA_log,y=RNA4_vs_DNA_log))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(0,3)+ylim(0,3)
ggplot(filter_norm_data_nc_normalize_exp,aes(x=RNA2_vs_DNA_log,y=RNA3_vs_DNA_log))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(0,3)+ylim(0,3)

head(new)
head(filter_norm_data_nc_normalize_exp)
write.table(new,file="CHD_variane_log2_ratio_FPM20_new_add_1.txt",sep="\t")
write.table(filter_norm_data,file="CHD_variane_FPM20_new_add_1.txt",sep="\t")
write.table(filter_norm_data_nc_normalize_exp,file="CHD_variane_normalized_log2_ratio_FPM20_add_1.txt",sep="\t")
##################################residual correction##########################


#################################rank figure#####################################
filter_norm_data_nc_normalize_exp <- read.table(file="CHD_variane_normalized_log2_ratio_FPM20_add_1.txt",header=T,row.names = 1)

mean_log_ratio <- rowMeans(filter_norm_data_nc_normalize_exp)
filter_norm_data_scaled <- cbind(filter_norm_data_nc_normalize_exp,mean_log_ratio)
filter_norm_data_scaled<-as.data.frame(filter_norm_data_scaled)
filter_norm_data_scaled$rank_value <- rank(-as.numeric(as.character(filter_norm_data_scaled$mean_log_ratio)))

select_REF <- t(filter_norm_data_scaled) %>% as.data.frame() %>% dplyr::select(contains("REF")) %>% t()
select_REF <- as.data.frame(select_REF)
select_REF$group <- "REF"
select_ALT <- t(filter_norm_data_scaled) %>% as.data.frame() %>% dplyr::select(contains("ALT")) %>% t()
select_ALT <- as.data.frame(select_ALT)
select_ALT$group <- "ALT"
select_NC <- t(filter_norm_data_scaled) %>% as.data.frame() %>% dplyr::select(!contains("::")) %>% t()
select_NC <- as.data.frame(select_NC)
select_NC$group <- "NegativeControl"
select_PC <- t(filter_norm_data_scaled) %>% as.data.frame() %>% dplyr::select(contains("::F")) %>% t()
select_PC <- as.data.frame(select_PC)
select_PC$group <- "PostiveControl"
filter_norm_data_scaled_grouped <-rbind(select_REF,select_ALT,select_NC,select_PC)
filter_norm_data_scaled_grouped <- as.data.frame(filter_norm_data_scaled_grouped)
##########################FDR 0.05################################################
FDR_negative <- function(Negative_C,colname,distribution="norm"){
  fitg<-fitdist(Negative_C[,colname] ,distribution)
  
  FDR_0.95 <- qnorm(0.95, mean = fitg$estimate["mean"], sd = fitg$estimate["sd"])
  return(FDR_0.95)
}
#fitg<-fitdist(select_NC$mean_log_ratio ,"norm")

#FDR_0.95 <- qnorm(0.95, mean = fitg$estimate["mean"], sd = fitg$estimate["sd"])
#head(filter_norm_data_scaled_grouped)
ggplot(filter_norm_data_scaled_grouped, aes(x = mean_log_ratio, color=group))+geom_density()+theme_bw()+theme_classic()
ggplot(filter_norm_data_scaled_grouped, aes(x = mean_log_ratio))+geom_density()+theme_bw()+theme_classic()
FDR_0.95<-FDR_negative(select_NC,"mean_log_ratio")

############################## heatmap #############################

filter_norm_data_scaled_grouped$order <- rank(-as.numeric(as.character(filter_norm_data_scaled_grouped$mean_log_ratio)))
filter_norm_data_scaled_grouped$active <- res[rownames(filter_norm_data_scaled_grouped),]$qvalue_0_05
##############################################################################
filter_norm_data_scaled_grouped_test <-Enrich_score(filter_norm_data_scaled_grouped, "group2", "log_d24_ratio", "order_d24_ratio")
dat_text_A <-filter_norm_data_scaled_grouped_test$dat_text
filter_norm_data_for_plot_A <- filter_norm_data_scaled_grouped_test$datasets
####################################rank figure####################################
library(patchwork)
pdf("CHD_library_figures.pdf", width = 11, height = 8)
A<-ggplot(filter_norm_data_scaled_grouped,aes(x=order,y=mean_log_ratio,color=factor(active)))+geom_point(size=2)+
  scale_color_manual(breaks = c("1", "0"), values=c("red", "darkblue"))+
  geom_hline(yintercept=FDR_0.95,color="red",linetype="dashed")+theme_bw()+theme_classic()+geom_text(aes(0,round(FDR_0.95,3),label = round(FDR_0.95,3), vjust = -1))

#ggplot(filter_norm_data_scaled_grouped, aes(x = order, y = group))+stat_bin2d(aes(fill = after_stat(count)), binwidth = c(3,1))
B<-ggplot(filter_norm_data_scaled_grouped, aes(x = order, y = group))+stat_bin2d(aes(fill = after_stat(count)), 
                                                                              binwidth = c(3,1))
A/B
B
dev.off()

#################################################################################
library(jmuOutlier) 
#library(distr6)
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
filter_norm_data_scaled_grouped_test <-Enrich_score(filter_norm_data_scaled_grouped, "group", "mean_log_ratio", "order")
dat_text <-filter_norm_data_scaled_grouped_test$dat_text
filter_norm_data_for_plot <- filter_norm_data_scaled_grouped_test$datasets
pdf("Enrichscore_of_CHD_library.pdf")

A<-ggplot(filter_norm_data_scaled_grouped,aes(x=order,y=mean_log_ratio,color=factor(active)))+
  geom_point(size=2)+
  theme_bw()+
  scale_color_manual(breaks = c("1", "0"), values=c("red", "gray"))+
  theme_classic()

B<-ggplot(filter_norm_data_scaled_grouped, aes(x = order, y = factor(group,levels=c("NegativeControl","PostiveControl",
                                                                                               "ALT","REF"))))+
  
  geom_bin2d(aes(alpha = ..density..), fill = "slategray",binwidth = c(5,1)) +
  scale_fill_continuous(guide = "none")+
  
  theme_bw()

p_E<- ggplot(filter_norm_data_for_plot, aes(rank_value, final_ES,color=factor(score_group)))+geom_line()+theme_bw()+theme_classic()+
  #facet_wrap(~score_group,scales = "free",nrow =8)+ 
  theme(panel.grid.major.y = element_blank(),
                                                          panel.grid.minor.y = element_blank(),
                                                          panel.grid.major.x = element_blank(),
                                                          panel.grid.minor.x = element_blank())+
  scale_color_manual(breaks = c("NegativeControl","PostiveControl",
                                "ALT","REF"), values=c("darkblue", "darkred","blue","red"))+
  geom_hline(yintercept = 0,color="skyblue",size=0.4,linetype="dashed")+
  geom_text(
  data    = dat_text,
  mapping = aes(x = -Inf, y = -Inf, label = label),
  hjust   = -2,
  vjust   = -0.25) 

A/p_E/B
dev.off()
#################################### pcc ##################################3
new <- filter_norm_data_scaled_grouped
A_V_cor <- cor(new[1:4])
corrplot(A_V_cor, method  = "number",tl.cex=0.5,hclust.method="median")
ggplot(new,aes(x=RNA1_vs_DNA_log,y=RNA2_vs_DNA_log))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(0,2)+ylim(0,2)
ggplot(new,aes(x=RNA1_vs_DNA_log,y=RNA3_vs_DNA_log))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(0,2)+ylim(0,2)
ggplot(new,aes(x=RNA1_vs_DNA_log,y=RNA4_vs_DNA_log))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(0,2)+ylim(0,2)
ggplot(new,aes(x=RNA2_vs_DNA_log,y=RNA3_vs_DNA_log))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(0,2)+ylim(0,2)

#write.table(new,file="CHD_variant_FPM20_log_rank_new.txt",sep="\t")
write.table(filter_norm_data,file="CHD_variant_FPM20_add_1.txt",sep="\t")
Anewdata <-cbind(filter_norm_data[rownames(new),],new)

write.table(Anewdata,file="CHD_variant_FPM20_log_rank_add_1.txt",sep="\t",quote = F)
Anewdata <-read.table(file="CHD_variant_FPM20_log_rank_add_1.txt",header = T,row.names = 1)
#pheatmap(new,cluster_rows=T,cluster_cols=F,scale="row",color = colorRampPalette(c("navy", "white","gold"))(50))
#pheatmap(new,cluster_rows=T,cluster_cols=F,color = colorRampPalette(c("navy", "white","gold"))(50))
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
write.table(filter_norm_data_scaled_grouped,"CHD_FPM20_new_grouped",sep="\t",quote=F)
table <- read.table("CHD_FPM20_new_grouped", header = T, row.names = 1)
#table <- read.table("CHD_variane_log2_ratio_FPM15.txt", header = T, row.names = 1)

table <-filter_norm_data_scaled_grouped
head(table)
table_RNA <- table %>% dplyr::select(contains("DNA_log"))

select_REF <- t(table_RNA) %>% as.data.frame() %>% dplyr::select(contains("REF")) %>% t()
select_ALT <- t(table_RNA) %>% as.data.frame() %>% dplyr::select(contains("ALT")) %>% t()
select_REF2 <- t(table) %>% as.data.frame() %>% dplyr::select(contains("REF")) %>% t()
select_ALT2 <- t(table) %>% as.data.frame() %>% dplyr::select(contains("ALT")) %>% t()
rownames_REF <- gsub("::REF", "", rownames(select_REF))
rownames_ALT <- gsub("::ALT", "", rownames(select_ALT))
intersect_loc <- intersect(rownames_REF, rownames_ALT)
new_REF_name <- paste(intersect_loc, "REF", sep = "::")
new_ALT_name <- paste(intersect_loc, "ALT", sep = "::")

rownames_REF2 <- gsub("::REF", "", rownames(select_REF2))
rownames_ALT2 <- gsub("::ALT", "", rownames(select_ALT2))
intersect_loc2 <- intersect(rownames_REF2, rownames_ALT2)
new_REF_name2 <- paste(intersect_loc2, "REF", sep = "::")
new_ALT_name2 <- paste(intersect_loc2, "ALT", sep = "::")

inter_REF <- select_REF[new_REF_name,]
inter_ALT <- select_ALT[new_ALT_name,]
inter_REF2 <- select_REF2[new_REF_name2,]
inter_ALT2 <- select_ALT2[new_ALT_name2,]

pvalue <- c()
for (i in 1:dim(inter_REF)[1]) {
  ref_value <- as.numeric(inter_REF[i,])
  alt_value <- as.numeric(inter_ALT[i,])
  pvalue_i <- t.test(ref_value, alt_value)
  pvalue <- c(pvalue, pvalue_i$p.value)
}
qvalue <- pvalue * length(pvalue) /rank(pvalue)
inter_REF <- as.data.frame(inter_REF)
inter_ALT <- as.data.frame(inter_ALT)

inter_REF2 <- as.data.frame(inter_REF2)
inter_ALT2 <- as.data.frame(inter_ALT2)
head(inter_REF2)

##################################3
rownames(inter_REF) <- intersect_loc
rownames(inter_ALT) <- intersect_loc
rownames(inter_REF2) <- intersect_loc2
rownames(inter_ALT2) <- intersect_loc2
inter_join <- cbind(inter_REF, inter_ALT)
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
inter_join2 <- Diff_enhancers(inter_join,c(1:4),c(5:8))

inter_join$REF <- rowMeans(inter_join[,1:4])
inter_join$ALT <- rowMeans(inter_join[,5:8])
inter_join$qvalue <- as.numeric(as.character(qvalue))
inter_join$pvalue <- as.numeric(as.character(pvalue))
#inter_join$fold <- rowMeans(inter_join[,5:8])/rowMeans(inter_join[,1:4])
inter_join$fold <- rowMeans(inter_join[,5:8])-rowMeans(inter_join[,1:4])
inter_join$qvalue_0_05 <- 0
inter_join$qvalue_0_05[inter_join$qvalue < 0.05] <- 1
inter_join$fold_1_5 <- 0
inter_join$fold_1_5 [inter_join$fold <= -0.58] <- 1
inter_join$fold_1_5 [inter_join$fold >= 0.58] <- 2
inter_join$active <- 0
inter_join[rownames(inter_ALT2[inter_ALT2$active==1,]),]$active <-1
inter_join[rownames(inter_REF2[inter_REF2$active==1,]),]$active <-1
inter_join$active_WT <- 0
inter_join$active_ALT <- 0
inter_join[rownames(inter_ALT2[inter_ALT2$active==1,]),]$active_ALT <-1
inter_join[rownames(inter_REF2[inter_REF2$active==1,]),]$active_WT <-1
#inter_join$active <- 0
#inter_join$active[rowMeans(inter_join[,1:4])>FDR_0.95] <-1
#inter_join$active[rowMeans(inter_join[,5:8])>FDR_0.95] <-1
inter_join$final_sig <-0
inter_join$final_sig[inter_join$qvalue_0_05 >0 &  inter_join$fold >= 0.58 & inter_join$active >0] <-2
inter_join$final_sig[inter_join$qvalue_0_05 >0 & inter_join$fold <= -0.58 & inter_join$active >0] <-1
inter_join$final_sig_2 <-0
inter_join$final_sig_2[inter_join$qvalue_0_05 >0 & inter_join$fold >= 1  & inter_join$active >0] <-2
inter_join$final_sig_2[inter_join$qvalue_0_05 >0 & inter_join$fold <= -1  & inter_join$active >0] <-1
inter_join$footprint <- inter_ALT2[rownames(inter_join),]$footprint
inter_join$WTfootprint <- inter_REF2[rownames(inter_join),]$footprint
#inter_join$ratiofootprint <- (as.numeric(inter_join$footprint))/(as.numeric(inter_join$WTfootprint))
head(inter_join)
#dim(inter_join[inter_join$REF > FDR_0.95,])
#dim(inter_join[inter_join$ALT > FDR_0.95,])
#length(inter_join$active[rowMeans(inter_join[,5:8])>FDR_0.95])
#dim(inter_join[inter_join$final_sig_2 ==1 & inter_join$active >0,])
#dim(inter_join[inter_join$final_sig_2 ==2 & inter_join$active >0,])
write.table(inter_join,file="CHD_variant_fold_pvalue_new_FPM20_add_1_newest.txt",sep="\t",quote = F)
inter_join<-read.table(file="CHD_variant_fold_pvalue_new_FPM20_add_1_newest.txt",header = T, row.names = 1)
colnames(inter_join)[5]<-"RNA1_vs_DNA_log_ALT"
colnames(inter_join)[6]<-"RNA2_vs_DNA_log_ALT"
colnames(inter_join)[7]<-"RNA3_vs_DNA_log_ALT"
colnames(inter_join)[8]<-"RNA4_vs_DNA_log_ALT"
ggplot(inter_join,aes(log2(fold), -log10(qvalue),color=factor(fold_1_5,levels=c('0','1','2')),alpha =factor(active,levels=c("0","1"))))+
  geom_point()+theme_bw()+theme_classic()+geom_vline(xintercept = 0.58,color="black",size=0.2,linetype="dashed")+
  geom_vline(xintercept = -0.58,color="black",size=0.2,linetype="dashed")+geom_hline(yintercept =-log10(0.05),color="black",size=0.2,linetype="dashed")+scale_color_manual(breaks = c("0", "1", "2"), values=c("gray", "red", "blue"))

ggplot(inter_join,aes(REF, ALT,color=factor(fold_1_5,levels=c('0','1','2')),alpha =factor(qvalue_0_05,levels=c("0","1"))))+
  geom_point()+theme_bw()+theme_classic()+geom_abline(slope = 1,color="black",size=0.2,linetype="dashed")+
  geom_abline(slope = 1.5,color="black",size=0.2,linetype="dashed")+
  geom_abline(slope=0.667,color="black",size=0.2,linetype="dashed")+
  scale_color_manual(breaks = c("0", "1", "2"), values=c("gray", "red", "blue"))+
  scale_alpha_manual(breaks = c("0", "1"), values=c(0.2, 1))
#+
 # geom_vline(xintercept = FDR_0.95,color="skyblue",size=0.4,linetype="dashed")+
  #geom_hline(yintercept = FDR_0.95,color="skyblue",size=0.4,linetype="dashed")

ggplot(inter_join,aes(REF, ALT,color=factor(final_sig,levels=c('0','1','2')),alpha =factor(qvalue_0_05,levels=c("0","1"))))+
  geom_point()+theme_bw()+theme_classic()+geom_abline(slope = 1,color="black",size=0.2,linetype="dashed")+
  geom_abline(slope = 1, intercept = 0.58,color="black",size=0.2,linetype="dashed")+
  geom_abline(slope=1,intercept = -0.58,color="black",size=0.2,linetype="dashed")+
  scale_color_manual(breaks = c("0", "1", "2"), values=c("gray", "red", "blue"))+
  scale_alpha_manual(breaks = c("0", "1"), values=c(0.2, 1))

#+
 # geom_vline(xintercept = FDR_0.95,color="skyblue",size=0.4,linetype="dashed")+
  #geom_hline(yintercept = FDR_0.95,color="skyblue",size=0.4,linetype="dashed")
pdf("Diff_enhancers_CHD.pdf",width = 9,height = 7)
ggplot(inter_join,aes(REF, ALT,color=factor(final_sig,levels=c('0','1','2')), alpha =factor(qvalue_0_05,levels=c("0","1"))))+
  geom_point()+theme_bw()+theme_classic()+geom_abline(slope = 1,color="black",size=0.2,linetype="dashed")+
  geom_abline(slope = 1, intercept = 0.58,color="black",size=0.2,linetype="dashed")+
  geom_abline(slope=1,intercept = -0.58,color="black",size=0.2,linetype="dashed")+
  scale_color_manual(breaks = c("0", "1", "2"), values=c("gray", "red", "blue"))+
  scale_alpha_manual(breaks = c("0", "1"), values=c(0.2, 1))+
  xlim(-4,4)+ylim(-4,4)

dev.off()
ggplot(inter_join,aes(REF, ALT,color=factor(final_sig,levels=c('0','1','2')),size=factor(footprint), alpha =factor(qvalue_0_05,levels=c("0","1"))))+
  geom_point()+theme_bw()+theme_classic()+geom_abline(slope = 1,color="black",size=0.2,linetype="dashed")+
  geom_abline(slope = 1.5,color="black",size=0.2,linetype="dashed")+
  geom_abline(slope=0.667,color="black",size=0.2,linetype="dashed")+
  scale_color_manual(breaks = c("0", "1", "2"), values=c("gray", "red", "blue"))+
  scale_alpha_manual(breaks = c("0", "1"), values=c(0.2, 1))+
  

ggplot(inter_join,aes(log2(fold), footprint,color=factor(final_sig,levels=c('0','1','2')), alpha =factor(qvalue_0_05,levels=c("0","1"))))+
  geom_point()+theme_bw()+theme_classic()+
  geom_vline(xintercept = 0,color="black",size=0.2,linetype="dashed")+
  geom_hline(yintercept =1,color="black",size=0.2,linetype="dashed")+
  scale_color_manual(breaks = c("0", "1", "2"), values=c("gray", "red", "blue"))+
  scale_alpha_manual(breaks = c("0", "1"), values=c(0.2, 1))

ggplot(inter_join[inter_join$active==1,],aes(fold, footprint,color=factor(final_sig,levels=c('0','1','2')), alpha =factor(qvalue_0_05,levels=c("0","1"))))+
  geom_point()+theme_bw()+theme_classic()+
  geom_vline(xintercept = 0,color="black",size=0.2,linetype="dashed")+
  geom_hline(yintercept =1,color="black",size=0.2,linetype="dashed")+
  scale_color_manual(breaks = c("0", "1", "2"), values=c("gray", "red", "blue"))+
  scale_alpha_manual(breaks = c("0", "1"), values=c(0.2, 1))

ggplot(inter_join[inter_join$active==1,],aes(fold, footprint/(WTfootprint),color=factor(final_sig,levels=c('0','1','2')), alpha =factor(qvalue_0_05,levels=c("0","1"))))+
  geom_point()+theme_bw()+theme_classic()+
  geom_vline(xintercept = 0,color="black",size=0.2,linetype="dashed")+
  scale_color_manual(breaks = c("0", "1", "2"), values=c("gray", "red", "blue"))+
  scale_alpha_manual(breaks = c("0", "1"), values=c(0.2, 1))

#ggplot(inter_join,aes(fold,ratiofootprint,color=factor(final_sig,levels=c('0','1','2')), alpha =factor(qvalue_0_05,levels=c("0","1"))))+
    geom_point()+theme_bw()+theme_classic()+
    geom_vline(xintercept = 0,color="black",size=0.2,linetype="dashed")+
  scale_color_manual(breaks = c("0", "1", "2"), values=c("gray", "red", "blue"))+
    scale_alpha_manual(breaks = c("0", "1"), values=c(0.2, 1))


ggplot(inter_join,aes(x=footprint))+
  geom_histogram(stat="count")+theme_bw()+theme_classic()+
  facet_wrap(~final_sig,nrow =3)+ theme(panel.grid.major.y = element_blank(),
                                  panel.grid.minor.y = element_blank(),
                                  panel.grid.major.x = element_blank(),
                                  panel.grid.minor.x = element_blank())

sub <-inter_join[inter_join$final_sig==1,1:8]      
pheatmap(sub,cluster_rows=T,cluster_cols=F,scale="row",color = colorRampPalette(c("navy", "white","gold"))(50))

sub2 <-inter_join[inter_join$final_sig !=0,c("REF","ALT","final_sig")]   
sub2$final_sig3 <-paste0("x",sub2$final_sig)
sub2[sub2$final_sig==1,]$final_sig3 <-"LoF"
sub2[sub2$final_sig==2,]$final_sig3 <-"GoF"
melt_sub2 <-melt(sub2[,c(1:2,4)])
pdf("REF_ALT_GoF_LoF_boxplot.pdf")
ggplot(melt_sub2, aes(y=factor(final_sig3,levels=c("LoF","GoF")),x=value)) +
  geom_boxplot(aes(colour = variable)) + theme_bw()+
  facet_wrap(~variable,nrow =2)+ theme(panel.grid.major.y = element_blank(),
                                       panel.grid.minor.y = element_blank(),
                                       panel.grid.major.x = element_blank(),
                                       panel.grid.minor.x = element_blank())
dev.off()
#sub3<- as.matrix(table[(table$group=="NegativeControl" |table$group=="PostiveControl") ,c("group","mean_log_ratio")])
#melt_sub2$group<-paste0(melt_sub2$final_sig3,"_",melt_sub2$variable)
sub3<- as.matrix(table[(table$group=="NegativeControl") ,c("group","mean_log_ratio")])
melt_sub2$group<-paste0(melt_sub2$final_sig3,"_",melt_sub2$variable)
sub4<-melt_sub2[,c("group","value")]
colnames(sub4)<-c("group","mean_log_ratio")
sub4_for_plot<-rbind(sub3,sub4)
#library(Rmisc)
#CI(Anewdata[Anewdata$active==1,"mean_log_ratio"],ci=0.95)
#CI(Anewdata[(Anewdata$active==1),"mean_log_ratio"],ci=0.95)
Ql<-quantile(Anewdata[(Anewdata$active==1),"mean_log_ratio"], 0.05)
Qh<-quantile(Anewdata[(Anewdata$active==1),"mean_log_ratio"], 0.95)
pdf("REF_ALT_GoF_LoF_boxplot_2.pdf")
ggplot(sub4_for_plot, aes(y=factor(group,levels=rev(c("GoF_REF","LoF_REF","GoF_ALT","LoF_ALT","NegativeControl",
                                                  "PostiveControl"))),x=as.numeric(mean_log_ratio))) +
  geom_violin(aes(fill = group))+
  geom_boxplot(aes(fill = group),color="white",width = 0.2) + theme_bw()+
  geom_vline(xintercept = Ql, color="darkblue",size=0.2,linetype="dashed")+
  geom_vline(xintercept = Qh, color="darkblue",size=0.2,linetype="dashed")

dev.off()  
library(ggpubr)
library(tidyverse)
library(broom)
library(AICcmodavg)  
summary(sub4_for_plot)
sub4_for_plot$value <- as.numeric(sub4_for_plot$mean_log_ratio)
one.way <- aov(value ~ group, data = sub4_for_plot)
tukey.one.way<-TukeyHSD(one.way)

tukey.one.way
#library(DescTools)
#DunnettTest(sub4_for_plot$value, sub4_for_plot$group)
length(Anewdata[Anewdata$active==1,"mean_log_ratio"])
sub5<-Anewdata[Anewdata$active==1,]
sub5$order2<-rank(Anewdata[Anewdata$active==1,"mean_log_ratio"])
sub5[sub5$order2==133,]$mean_log_ratio
ggplot(inter_join[diffnames,],aes(REF, ALT,color=factor(final_sig,levels=c('0','1','2')), alpha =factor(qvalue_0_05,levels=c("0","1"))))+
  geom_point()+theme_bw()+theme_classic()+geom_abline(slope = 1,color="black",size=0.2,linetype="dashed")+
  geom_abline(slope = 1, intercept = 0.58,color="black",size=0.2,linetype="dashed")+
  geom_abline(slope=1,intercept = -0.58,color="black",size=0.2,linetype="dashed")+
  scale_color_manual(breaks = c("0", "1", "2"), values=c("gray", "red", "blue"))+
  scale_alpha_manual(breaks = c("0", "1"), values=c(0.2, 1))+
  xlim(-4,4)+ylim(-4,4)

ggplot(inter_join[inter_join$active==1,], aes(x=fold,color=factor(footprint)))+geom_density(size=1)+theme_bw()+theme_classic()
#  facet_wrap(~inter_join$footprint,nrow =4)+ theme(panel.grid.major.y = element_blank(),
#                                    panel.grid.minor.y = element_blank(),
#                                    panel.grid.major.x = element_blank(),
#                                    panel.grid.minor.x = element_blank())
inter_join$foot <-0
inter_join[inter_join$footprint>0,]$foot <- 1
ggplot(inter_join[inter_join$WTfootprint>=1 & inter_join$active==1 ,],aes(final_sig3,color=factor(foot)))+
  geom_density()+theme_bw()+theme_classic()
#############################################################
data3 <- read.table("mutation_located_at_CDS_region.bed",header = F)
WT_FPM_more20 <- intersect(data3[,2],rownames_REF2)
ALT_FPM_more20 <- intersect(data3[,2],rownames_ALT2)

rownames(data3) <- data3[,2]
inter_join$CDS_region <- 0
inter_join[intersect(rownames(inter_join),rownames(data3)),]$CDS_region <- 1 
inter_join_no_CDS <- inter_join[inter_join$CDS_region==0,]
write.table(inter_join,file="CHD_variant_fold_pvalue_new_FPM20_add_1_newest_with_CDS.txt",sep="\t",quote = F)
number <- as.data.frame(table(inter_join_no_CDS$final_sig))
dat_text_2 <- data.frame(
  label = c("nochange","GoF","LoF"),
  score_group = number$Freq
)
pdf("Diff_enhancers_CHD.pdf",width = 9,height = 7)
  ggplot(inter_join_no_CDS,aes(REF, ALT,color=factor(final_sig,levels=c('0','1','2')), alpha =factor(qvalue_0_05,levels=c("0","1"))))+
  geom_point()+theme_bw()+theme_classic()+geom_abline(slope = 1,color="black",size=0.2,linetype="dashed")+
  geom_abline(slope = 1, intercept = 0.58,color="black",size=0.2,linetype="dashed")+
  geom_abline(slope=1,intercept = -0.58,color="black",size=0.2,linetype="dashed")+
  scale_color_manual(breaks = c("0", "1", "2"), values=c("gray", "red", "blue"))+
  scale_alpha_manual(breaks = c("0", "1"), values=c(0.2, 1))+
  xlim(-4,4)+ylim(-4,4)
dev.off()
write.table(dat_text_2,file = "Diff_enhancers_number.txt",quote=F, sep = "\t")
W1<-as.data.frame(table(inter_join_no_CDS$active_WT))
A1 <-as.data.frame(table(inter_join_no_CDS$active_ALT))
dat_text_3 <- cbind(W1, A1) 
colnames(dat_text_3) <-c("WT_active","number","ALT_active","number")
write.table(dat_text_3, "CHD_active_regions_number.txt",quote=F, sep = "\t")
