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
####make new inputfiles ##
data<-read.table(file="../data/Mutagenesis_raw_counts_DNA1-4_RNA1-4.txt",header=T,row.names = 1)
head(data)
select_ALT <- t(data) %>% as.data.frame() %>% dplyr::select(contains("::") & !contains("_17")) %>% t()
select_ALT <- as.data.frame(select_ALT)
select_ALT$group <- "ALT"
select_REF <- t(data) %>% as.data.frame() %>% dplyr::select(contains("_17")) %>% t()
select_REF <- as.data.frame(select_REF)
select_REF$group <- "WT"

select_NC <- t(data) %>% as.data.frame() %>% dplyr::select(!contains("::")) %>% t()
select_NC <- as.data.frame(select_NC)
select_NC$group <- "NegativeControl"
data_grouped <-rbind(select_REF,select_ALT,select_NC)
groups_split<-str_split_fixed(rownames(data_grouped), pattern="::|_",3)
data_grouped$group_name <-groups_split[,1]
data_grouped$group_number <-groups_split[,3]
data_grouped$group_F <-groups_split[,2]
data_grouped[data_grouped$group_F=="",]$group_F<-data_grouped[data_grouped$group_F=="",]$group_name
data_grouped[data_grouped$group_number=="",]$group_F<-data_grouped[data_grouped$group_number=="",]$group_name

data2<-read.table("../data/library3_171_full_names.txt",header= F, row.names = 1)

data_grouped$full_names <- data2[as.character(data_grouped$group_name),]
data_grouped[is.na(data$full_names),]$full_names<-rownames(data_grouped[is.na(data$full_names),])
write.table(data_grouped,file="Mutagenesis_raw_data_with_group.txt",sep="\t",quote=F)
###################
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
sum(res$qvalue_0_05)#1898
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
DNA_total <-as.data.frame(rowMin(norm_data[,c(1:4)]))
colnames(DNA_total)<-c("DNA_min_value")
pdf("DNA_counts_distribution_Mutagenesis.pdf")
ggplot(DNA_total, aes(x = DNA_min_value))+ geom_histogram(binwidth = 5,color="skyblue")+theme_bw()+theme_classic()+
  geom_vline(xintercept = 20,color="skyblue",size=0.4,linetype="dashed")+xlim(0,500)
dev.off()
############################################################################
total_number <- length(DNA_total[,1])
kept_number <- length(filter_norm_data[,1])
keptratio <- kept_number/total_number
DNA_counts <- data.frame(
  label=c("total_number","kept_number","keptratio"),
  value=c(as.character(total_number),as.character(kept_number),as.character(keptratio))
)

write.table(DNA_counts,"Mutagenesis_library_DNA_distribution_filter_log.txt",sep="\t",quote=F)
###########################################################################


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

A_V_cor <- cor(filter_norm_data_nc_normalize_exp)
A_V_cor2 <- cor(new)
corrplot(A_V_cor, method  = "number",tl.cex=0.5,hclust.method="median")
corrplot(A_V_cor2, method  = "number",tl.cex=0.5,hclust.method="median")
ggplot(filter_norm_data_nc_normalize_exp,aes(x=RNA1_vs_DNA_log,y=RNA2_vs_DNA_log))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(0,3)+ylim(0,3)
ggplot(filter_norm_data_nc_normalize_exp,aes(x=RNA1_vs_DNA_log,y=RNA3_vs_DNA_log))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(0,3)+ylim(0,3)
ggplot(filter_norm_data_nc_normalize_exp,aes(x=RNA1_vs_DNA_log,y=RNA4_vs_DNA_log))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(0,3)+ylim(0,3)
ggplot(filter_norm_data_nc_normalize_exp,aes(x=RNA2_vs_DNA_log,y=RNA3_vs_DNA_log))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(0,3)+ylim(0,3)
#####################################################################################3
data_mutagenesis <-filter_norm_data_nc_normalize_exp[both_C,]
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

combin <-merge(data_mutagenesis,data_CHD,by="row.names")

rownames(combin)<-combin$Row.names
combin<-combin[,c(2:9)]
combin_p <-t(dplyr::select(as.data.frame(t(combin)), contains("::F")))
combin_n <-t(dplyr::select(as.data.frame(t(combin)), !contains("::F")))
pairs(combin_p, 
      lower.panel = panel.cor,
      upper.panel = upper.panel)
combin_cor <- cor(combin)
pdf("CHD_vs_mutegenesis_x_muta_y_CHD_pcc.pdf",width = 11, height = 11)
p<-corrplot(combin_cor, method  = "number",order = "alphabet")
pairs(combin, 
      lower.panel = panel.cor,
      upper.panel = upper.panel)
dev.off()
write.table(combin,file="Mutagenesis_and_CHD_overlapped.txt",sep="\t",quote = F)
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
p<-corrplot(A_V_cor, method  = "number",order = "alphabet",type = 'lower')
pairs(filter_norm_data_nc_normalize_exp, 
      lower.panel = panel.cor,
      upper.panel = upper.panel)
dev.off()

###########################################################################
write.table(new,file="10bp_tiled_log2_ratio_add1.txt",sep="\t")
write.table(filter_norm_data,file="10bp_tiled_log2_FPM_add1.txt",sep="\t")
filter_norm_data<-read.table(file="10bp_tiled_log2_FPM_add1.txt",header=T,row.names = 1)
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

#################################################################
library(jmuOutlier) 
filter_norm_data_scaled_grouped_test <-Enrich_score(filter_norm_data_scaled_grouped, "group", "mean_log_ratio", "order")
dat_text <-filter_norm_data_scaled_grouped_test$dat_text
filter_norm_data_for_plot <- filter_norm_data_scaled_grouped_test$datasets
##################################################################

library(patchwork)
pdf("mutagenesis_library_figures.pdf", width = 11, height = 8)
#ggplot(filter_norm_data_scaled_grouped,aes(x=order,y=mean_log_ratio,color=factor(active)))+geom_point(size=2)+
#  scale_color_manual(breaks = c("1", "0"), values=c("red", "darkblue"))+
#  geom_hline(yintercept=FDR_0.95,color="red",linetype="dashed")+theme_bw()+theme_classic()+geom_text(aes(0,round(FDR_0.95,3),label = round(FDR_0.95,3), vjust = -1))

A<-ggplot(filter_norm_data_scaled_grouped,aes(x=order,y=mean_log_ratio,color=factor(active)))+
  geom_point(size=2)+
  theme_bw()+
  scale_color_manual(breaks = c("1", "0"), values=c("red", "gray"))+
  theme_classic()

B<-ggplot(filter_norm_data_scaled_grouped, aes(x = order, y = factor(group,levels=c("NegativeControl",
                                                                                    "ALT","WT"))))+
  
  geom_bin2d(aes(alpha = ..density..), fill = "slategray",binwidth = c(5,1)) +
  scale_fill_continuous(guide = "none")+
  
  theme_bw()
p_E<- ggplot(filter_norm_data_for_plot, aes(rank_value, final_ES,color=factor(score_group)))+geom_line()+theme_bw()+theme_classic()+
  #facet_wrap(~score_group,scales = "free",nrow =8)+ 
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  scale_color_manual(breaks = c("NegativeControl",
                                "ALT","WT"), values=c("darkblue", "blue","red"))+
  geom_hline(yintercept = 0,color="skyblue",size=0.4,linetype="dashed")+
  geom_text(
    data    = dat_text,
    mapping = aes(x = -Inf, y = -Inf, label = label),
    hjust   = -2,
    vjust   = -0.25) 
A/p_E/B
dev.off()

write.table(filter_norm_data_scaled_grouped,file="10bp_tiled_log2_normalized_log2_ratio_with_active_add1.txt",sep="\t",quote = F)
f

filter_norm_data_scaled_grouped<-read.table(file="10bp_tiled_log2_normalized_log2_ratio_with_active_add1.txt",header = T,row.names = 1)
mutagenesis_merged<- merge(filter_norm_data,filter_norm_data_scaled_grouped,by="row.names")
write.table(mutagenesis_merged,file="10bp_tiled_log2_normalized_log2_ratio_with_active_add1_merged.txt",sep="\t",quote = F)
f
head(filter_norm_data_scaled_grouped)
table(filter_norm_data_scaled_grouped$group)
#ggplot(filter_norm_data_scaled_grouped, aes(x = order, y = group))+stat_bin2d(aes(fill = after_stat(count)), binwidth = c(1,1))
#g.iris <- ggplotGrob(g1)
#g2 <-ggplot(filter_norm_data_scaled_grouped, aes(x = order, y = group))+stat_bin2d(aes(fill = after_stat(count)), binwidth = c(3,1))
#g.mpg <- ggplotGrob(g2)
#iris.widths <- g.iris$widths[1:3]
#mpg.widths <- g.mpg$widths[1:3]
#max.widths <- unit.pmax(iris.widths, mpg.widths)
#g.iris$widths[1:3] <- max.widths # assign max. widths to iris gtable
#g.mpg$widths[1:3] <- max.widths # assign max widths to mpg gtable

# plot_grid() can work directly with gtables, so this works
#plot_grid(g.iris, g.mpg, align="v", ncol = 1,rel_widths=c(2,2))
#plot_grid(g1,g2, align="v",nrow=2)
#melt_filter_norm_data_2 <- melt(filter_norm_data[,c(1:5,7)])
#ggplot(melt_filter_norm_data_2, aes(x = value, color=variable))+geom_density()+scale_x_log10()+theme_bw()+theme_classic()


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
inter_join$active_WT <- filter_norm_data_nc@res[inter_join$WT,]$Active
inter_join$active_KO <- filter_norm_data_nc@res[inter_join$KO,]$Active
inter_join$active_WT <- res[inter_join$WT,]$qvalue_0_05
inter_join$active_KO <- res[inter_join$KO,]$qvalue_0_05
inter_join$fold <- as.numeric(as.character(inter_join$mean_KO))-as.numeric(as.character(inter_join$mean_WT))
inter_join$qvalue_0_05 <- 0
inter_join$qvalue_0_05[inter_join$qvalue < 0.05] <- 1
inter_join$fold_1_5 <- 0
inter_join$fold_1_5 [inter_join$fold <= -0.58] <- 1
inter_join$fold_1_5 [inter_join$fold >= 0.58] <- 2
inter_join$active <- 0
inter_join$active[as.numeric(as.character(inter_join$active_WT))>0] <-1
inter_join$active[as.numeric(as.character(inter_join$active_KO))>0] <-1
inter_join$final_sig <-0
inter_join$final_sig[inter_join$qvalue_0_05 >0 & inter_join$fold <= -0.58 & inter_join$active >0] <-1
inter_join$final_sig[inter_join$qvalue_0_05 >0 & inter_join$fold >= 0.58 & inter_join$active >0] <-2
inter_join$final_sig_2 <-0
inter_join$final_sig_2[inter_join$qvalue_0_05 >0 & inter_join$fold >= 0.58] <-2
inter_join$final_sig_2[inter_join$qvalue_0_05 >0 & inter_join$fold <= -0.58] <-1
write.table(inter_join,file="10bp_tiled_log2_normalized_log2_ratio_fold_pvalue_add1.txt",sep="\t")

####################################################################
inter_join<-read.table(file="10bp_tiled_log2_normalized_log2_ratio_fold_pvalue.txt",header = T, row.names = 1)

inter_join$mean_KO <- as.numeric(as.character(inter_join$mean_KO))
inter_join$mean_WT <- as.numeric(as.character(inter_join$mean_WT))

ggplot(inter_join,aes(log2(fold), -log10(qvalue),color=factor(fold_1_5,levels=c('0','1','2')),alpha =factor(active,levels=c("0","1"))))+
  geom_point()+theme_bw()+theme_classic()+geom_vline(xintercept = 0.58,color="black",size=0.2,linetype="dashed")+
  geom_vline(xintercept = -0.58,color="black",size=0.2,linetype="dashed")+geom_hline(yintercept =-log10(0.05),color="black",size=0.2,linetype="dashed")+scale_color_manual(breaks = c("0", "1", "2"), values=c("gray", "red", "blue"))

ggplot(inter_join,aes(mean_WT, mean_KO,color=factor(final_sig_2,levels=c('0','1','2')),alpha =factor(active,levels=c("0","1"))))+
  geom_point()+theme_bw()+theme_classic()+geom_abline(slope = 1,color="black",size=0.2,linetype="dashed")+
  geom_abline(slope = 1,intercept=0.58,color="black",size=0.2,linetype="dashed")+
  geom_abline(slope=0.667,intercept=-0.58,color="black",size=0.2,linetype="dashed")+
  scale_color_manual(breaks = c("0", "1", "2"), values=c("gray", "red", "blue"))+
  scale_alpha_manual(breaks = c("0", "1"), values=c(0.2, 1))+
  geom_vline(xintercept = FDR_0.95,color="skyblue",size=0.4,linetype="dashed")+
  geom_hline(yintercept = FDR_0.95,color="skyblue",size=0.4,linetype="dashed")+
  xlim(-5,5)+ylim(-5,5)
pdf("Diff_enhancers_mutagenesis.pdf",width = 9,height = 7)
ggplot(inter_join,aes(mean_WT, mean_KO,color=factor(final_sig,levels=c('0','1','2')),alpha =factor(qvalue_0_05,levels=c("0","1"))))+
  geom_point()+theme_bw()+theme_classic()+geom_abline(slope = 1,color="black",size=0.2,linetype="dashed")+
  geom_abline(slope = 1,intercept=0.58,color="black",size=0.2,linetype="dashed")+
  geom_abline(slope=1,intercept=-0.58,color="black",size=0.2,linetype="dashed")+
  scale_color_manual(breaks = c("0", "1", "2"), values=c("gray", "red", "blue"))+
  scale_alpha_manual(breaks = c("0", "1"), values=c(0.2, 1))+
  geom_vline(xintercept = FDR_0.95,color="skyblue",size=0.4,linetype="dashed")+
  geom_hline(yintercept = FDR_0.95,color="skyblue",size=0.4,linetype="dashed")+
  xlim(-5,5)+ylim(-5,5)
dev.off()
ggplot(inter_join,aes(mean_WT, mean_KO,color=factor(final_sig,levels=c('0','1','2')), alpha =factor(qvalue_0_05,levels=c("0","1"))))+
  geom_point()+theme_bw()+theme_classic()+geom_abline(slope = 1,color="black",size=0.2,linetype="dashed")+
  geom_abline(slope = 1.5,color="black",size=0.2,linetype="dashed")+
  geom_abline(slope=0.667,color="black",size=0.2,linetype="dashed")+
  scale_color_manual(breaks = c("0", "1", "2"), values=c("gray", "red", "blue"))+
  scale_alpha_manual(breaks = c("0", "1"), values=c(0.2, 1))
#sub <-inter_join[inter_join$final_sig==1,1:8]      
#pheatmap(sub,cluster_rows=T,cluster_cols=F,scale="row",color = colorRampPalette(c("navy", "white","gold"))(50))
sum(inter_join$final_sig==1)
sum(inter_join$final_sig==2)


###########################################################################################
number <- as.data.frame(table(inter_join$final_sig))
dat_text_2 <- data.frame(
  label = c("nochange","LoF","GoF"),
  score_group = number$Freq
)
write.table(dat_text_2,file = "Diff_enhancers_number_mutagenesis.txt",quote=F, sep = "\t")
W1<-as.data.frame(table(inter_join$active_WT))
A1 <-as.data.frame(table(inter_join$active_KO))
dat_text_3 <- cbind(W1, A1) 
colnames(dat_text_3) <-c("WT_active","number","ALT_active","number")
write.table(dat_text_3, "mutagenesis_active_regions_number.txt",quote=F, sep = "\t")
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
pdf("active_regions.pdf")
pheatmap(active_paire_matrix,cluster_cols = F,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
ggplot(active_paired, aes(group_F))+geom_histogram(stat = "count")+theme_classic()
dev.off()
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

#inter_join_location_1<-inter_join_location[inter_join_location$full_names=="chr20:49530503-49530903",]
#inter_join_location_1<-inter_join_location[inter_join_location$full_names=="chr10:35047312-35047712",]
#inter_join_location_1<-inter_join_location[inter_join_location$full_names=="chr6:39745032-39745432",]
#inter_join_location_1<-inter_join_location[inter_join_location$full_names=="chr3:34053580-34053980",]
full_names_1 <- inter_join_location[inter_join_location$final_sig>0,]

table(unique(full_names_1)$full_names)[as.character(full_names_1$full_names)]

inter_join_location_1<-inter_join_location[inter_join_location$full_names %in% full_names_1,]
#inter_join_location_names <- head(inter_join_location[inter_join_location$final_sig>0,]$full_names)
#inter_join_location_names <- c("chr7:93870775-93871175","chr14:59536249-59536649","chr4:148231410-148231810","chr1:245748592-245748992")
#inter_join_location_names <- c("chr15:96808175-96808575","chr2:47925801-47926201","chr18:22819712-22820112","chr13:39332360-39332760")
inter_join_location_names <- c("chr3:82149507-82149907",
                               "chr20:49157407-49157807",
                               "chr9:137494462-137494862",
                               "chr9:101837589-101837989",
                               "chr4:77586893-77587293")
#inter_join_location_names<-c("chr9:11135936-11136336",
#                             "chr15:30133105-30133505",
#                             "chr15:96666221-96666621",
#                             "chr3:34053580-34053980",
#                             "chr4:17273756-17274156"
#)
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

##############################################################################################
head(inter_join_location)
sig_motif <-read.table(file = "mutagenesis_WT-ALT_fimo_new_2_expressed_new_new.txt",header=F)
head(sig_motif)
sig_motif$color<-"0"
sig_motif$color<-"1"
sig_motif[(sig_motif$V7-sig_motif$V5)<=-2,]$color<-"2"

sig_motif_for_merge <-sig_motif[abs(sig_motif$V7-sig_motif$V5)>=4,c("V2","V8","color")]

inter_join_location_names<-c(#"chr9:11135936-11136336",
                             #"chr15:30133105-30133505",
                             #"chr15:96666221-96666621",
                             "chr3:82149507-82149907",
                             "chr20:49157407-49157807",
                             "chr9:137494462-137494862",
                             "chr9:101837589-101837989")
                             #"chr4:77586893-77587293"
#)

inter_join_location_motif <- merge(inter_join_location,sig_motif_for_merge,by.x="KO",by.y="V2",all.x=T)
inter_join_location_1<-inter_join_location_motif[inter_join_location_motif$full_names %in% inter_join_location_names,]
inter_join_location_1$motif_show<-inter_join_location_1$motif_show
motif_split<-str_split_fixed(inter_join_location_1$V8, pattern=",",2)
inter_join_location_1$motif_show<-motif_split[,1]
inter_join_location_1[inter_join_location_1$final_sig==0,]$motif_show<-NA

pdf(file="mutagenesis_example_2.pdf",width = 10,height = 6)
ggplot(inter_join_location_1,aes(as.numeric (group_number), as.numeric(as.character(mean_KO)),label=motif_show))+
  geom_point(aes(color=factor(final_sig,levels=c('0','1','2')), alpha =factor(qvalue_0_05,levels=c("0","1")),size=factor(active,levels=c('0','1'))))+
  theme_bw()+theme_classic()+
  geom_hline(yintercept = 0,color="black",size=0.4,linetype="dashed")+
  scale_color_manual(breaks = c("0", "1", "2"), values=c("gray", "red", "blue"))+
  scale_alpha_manual(breaks = c("0", "1"), values=c(0.2, 1))+
  scale_size_manual(breaks = c("0", "1"), values=c(2, 4))+
  geom_line(aes(as.numeric (group_number), as.numeric(as.character(mean_WT)),alpha =factor(active_WT,levels=c("0","1"))),color="blue",linetype="dashed")+
  geom_line(aes(as.numeric (group_number), as.numeric(as.character(KO_median)),alpha =factor(active_KO,levels=c("0","1"))),color="coral",linetype="dotdash")+
  facet_grid(rows = vars(full_names),cols = vars(group_F),scales="free")+
  geom_line(aes(x=as.numeric (group_number), y=as.numeric(as.character(mean_KO))), alpha=0.2)+
  geom_text(angle = 90, aes(colour =factor(color)),position = position_dodge(1),
            vjust = 2)
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
dev.off()
write.table(inter_join_location_motif,file="mutagenesis_library_with_changed_motif.txt",sep="\t",quote=F)
write.table(inter_join_location_1,file="mutagenesis_library_with_changed_example_2.txt",sep="\t",quote=F)
################################################################################################


head(inter_join_location_1)
ggplot(inter_join_location_1, aes(y=full_names,group_number, fill=fold, alpha=factor(final_sig)),)+ 
  geom_tile(colour = "grey50")+
  facet_grid(cols = vars(group_F),scales="free")+
  scale_alpha_manual(breaks = c("0", "1","2"), values=c(0.2, 1,1))+
  scale_fill_gradient2(
    low = "navy",
    mid = "white",
    high = "red",
    midpoint = 0)+
  theme_classic()
inter_join_location_1$group_number_new<-inter_join_location_1$group_number
inter_join_location_1[inter_join_location_1$group_F=="F2",]$group_number_new<-as.numeric(inter_join_location_1[inter_join_location_1$group_F=="F2",]$group_number)+17
inter_join_location_1[inter_join_location_1$group_F=="F3",]$group_number_new<-as.numeric(inter_join_location_1[inter_join_location_1$group_F=="F3",]$group_number)+34
melt_interjoin_location_1 <- melt(inter_join_location_1[,c("full_names","group_number_new","fold")])
inter_join_location_acast <-acast(melt_interjoin_location_1[,c(1,2,4)],formula =full_names~group_number_new, fill = 0)
pheatmap(inter_join_location_acast,cluster_cols=FALSE,color = colorRampPalette(c("navy", "white","gold"))(50),breaks=seq(-3,3,0.12))
order(colnames(inter_join_location_acast))
inter_join_location_acast_sorted<-rownames(inter_join_location_acast[inter_join_location_acast$tree_row[["order"]],])
ggplot(inter_join_location_1, aes(y=full_names,group_number, fill=fold, alpha=factor(final_sig)),)+ 
  geom_tile(colour = "grey50")+
  facet_grid(cols = vars(group_F),scales="free")+
  scale_alpha_manual(breaks = c("0", "1","2"), values=c(0.2, 1,1))+
  scale_fill_gradient2(
    low = "navy",
    mid = "white",
    high = "red",
    midpoint = 0)+
  theme_classic()



ggplot(filter_norm_data_scaled_grouped_2[filter_norm_data_scaled_grouped_2$active_WT>0,], aes(x = group_F))+ geom_bar(color="skyblue")+theme_bw()+theme_classic()
geom_bar(aes(log2(fold),color=factor(fold_1_5,levels=c('0','1','2')),alpha =factor(qvalue_0_05,levels=c("0","1"))))



write.table(inter_join_location,file="10bp_tiled_log2_normalized_log2_ratio_fold_pvalue_final_add1.txt",sep="\t")  








inter_join_location[inter_join_location$mean_WT == 1.82784864351185,]

head(inter_join)
WT_selct<- unique(inter_join[,c(2,7:10,12,15)])
row.names(WT_selct)<-WT_selct$WT
WT_selct <- merge(filter_norm_data,WT_selct,by=0)
write.table(WT_selct,file="10bp_Mutagenesis_WT_final.txt",sep="\t",quote = T)
###########################################################################################

