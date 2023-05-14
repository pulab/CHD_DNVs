setwd("~/Desktop/work/Bill/Fengxiao/muatiaon_MPRA/raw_counts")
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
library(stringr)


#########################################

RNA=data[,c(5:8)]
DNA=data[,c(1:4)]

data_norm<-create_MPRA_class(RNA=RNA,DNA=DNA)

##############

DNA_total <- as.data.frame(matrixStats::rowMins(as.matrix(data_norm@DNA)))
colnames(DNA_total)<-"minDNA"
dim(data_norm@DNA)
#for_plot <- as.data.frame(data_norm@DNA) %>% melt()
pdf("DNA_counts_distribution_CHD.pdf")
ggplot(DNA_total, aes(x = minDNA))+ geom_histogram(binwidth = 5,color="skyblue")+theme_bw()+theme_classic()+
  geom_vline(xintercept = 20,color="skyblue",size=0.4,linetype="dashed")+xlim(0,300)
dev.off()

#######################################################################
filter_norm_data <- filter_data(data_norm, dna_filter=0)
dim(filter_norm_data@DNA)

filter_norm_data <-DEGseq2_active(filter_norm_data=filter_norm_data)
sum(filter_norm_data@res$Active)
#write.table(filter_norm_data@res$Active, file="CHD_DES_result.txt",sep = "\t",quote = F)
#######################################################################
filter_norm_data <- filter_data(filter_norm_data, dna_filter=20)
total_number <- length(DNA_total[,1])
kept_number <- length(filter_norm_data@DNA[,1])
keptratio <- kept_number/total_number
DNA_counts <- data.frame(
  label=c("total_number","kept_number","keptratio"),
  value=c(as.character(total_number),as.character(kept_number),as.character(keptratio))
)
DNA_counts
#write.table(DNA_counts,"CHD_library_DNA_distribution_filter_log.txt",sep="\t",quote=F)
###########################################################################
filter_norm_data_new <-get_ratio(filter_norm_data)

head(filter_norm_data_new@res)

######################size factor#####################################

Negative_C <- t(dplyr::select(as.data.frame(t(filter_norm_data_new@now_value)), !contains("::")))
Positive_C <- t(dplyr::select(as.data.frame(t(filter_norm_data_new@now_value)), contains("::F")))
both_C <- rownames(rbind(Negative_C,Positive_C))

filter_norm_data_nc <-MPRA_SizeFactors(filter_norm_data_new,both_C,c(1:4))
#filter_norm_data_nc_normalize_exp <- as.data.frame(filter_norm_data_nc@res[,c(6:9)])
filter_norm_data_nc_normalize_exp<-filter_norm_data_nc@now_value

######################################pcc figure ##########################################

# Create the PCC plots

#pdf("MPRA_CHD_pcc.pdf",width = 11, height = 11)
A_V_cor <- cor(filter_norm_data_nc_normalize_exp, method = "spearman")
p<-corrplot(A_V_cor, method  = "number",order = "alphabet")
scatter_plot_cor(filter_norm_data_nc_normalize_exp,cor_method="spearman")
#dev.off()
#########################################################################
ggplot(filter_norm_data_nc_normalize_exp,aes(x=RNA1_ratio_scaled,y=RNA2_ratio_scaled))+geom_point()+theme_bw()+theme_classic()+geom_abline(intercept = 0, slope = 1,color="black",size=0.5)+xlim(0,3)+ylim(0,3)

head(filter_norm_data_nc_normalize_exp)

#################################rank figure#####################################
#filter_norm_data_nc_normalize_exp <- read.table(file="CHD_variane_normalized_log2_ratio_FPM20_add_1.txt",header=T,row.names = 1)
mean_log_ratio <- rowMeans(filter_norm_data_nc_normalize_exp)
filter_norm_data_scaled <- cbind(filter_norm_data_nc_normalize_exp,mean_log_ratio)
filter_norm_data_scaled<-as.data.frame(filter_norm_data_scaled)
filter_norm_data_scaled$rank_value <- rank(-as.numeric(as.character(filter_norm_data_scaled$mean_log_ratio)))
filter_norm_data_nc@res$mean_log_ratio <-mean_log_ratio
filter_norm_data_nc@res$order <-filter_norm_data_scaled$rank_value
########add a Group factor to the annotation matrix
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
filter_norm_data_nc@annot$Group<- filter_norm_data_scaled_grouped[rownames(filter_norm_data_nc@annot),]$group
#colnames(filter_norm_data_nc@annot)[ncol(filter_norm_data_nc@annot)]<-"Group"
########3
filter_norm_data_scaled_grouped <- as.data.frame(filter_norm_data_scaled_grouped)
##########################FDR 0.05################################################
FDR_negative <- function(Negative_C,colname,distribution="norm",fdr_value=0.05){
  th_value=1-fdr_value
  fitg<-fitdist(Negative_C[,colname] ,distribution)
  
  FDR_0.95 <- qnorm(0.95, mean = fitg$estimate["mean"], sd = fitg$estimate["sd"])
  return(FDR_0.95)
}

ggplot(filter_norm_data_scaled_grouped, aes(x = mean_log_ratio, color=group))+geom_density()+theme_bw()+theme_classic()
#ggplot(filter_norm_data_scaled_grouped, aes(x = mean_log_ratio))+geom_density()+theme_bw()+theme_classic()
FDR_0.95<-FDR_negative(select_NC,"mean_log_ratio")

filter_norm_data_nc@res$FDR_active<-0
filter_norm_data_nc@res[filter_norm_data_nc@res$mean_log_ratio>FDR_0.95,]$FDR_active <- 1
group_data<-as.data.frame(filter_norm_data_nc@annot$Group)
names(group_data)<-"group"
filter_norm_data_scaled_grouped <- cbind(filter_norm_data_nc@res[,c("mean_log_ratio","order","Active")],
                                         group_data)


################################3A###########################
filter_norm_data_scaled_grouped$rank_value <- rank(-as.numeric(as.character(filter_norm_data_scaled_grouped$mean_log_ratio)))
filter_norm_data_scaled_grouped_test <-Enrich_score(filter_norm_data_scaled_grouped, "group", "mean_log_ratio", "order")

dat_text <-filter_norm_data_scaled_grouped_test$dat_text
filter_norm_data_for_plot <- filter_norm_data_scaled_grouped_test$datasets
#pdf("Enrichscore_of_CHD_library.pdf")

A<-ggplot(filter_norm_data_scaled_grouped,aes(x=order,y=mean_log_ratio,color=factor(Active)))+
  geom_point(size=2)+
  theme_bw()+
  scale_color_manual(breaks = c("1", "0"), values=c("red", "gray"))+
  theme_classic()

B<-ggplot(filter_norm_data_scaled_grouped, aes(x = order, y = factor(group,levels=c("NegativeControl","PostiveControl",
                                                                                    "ALT","REF"))))+
  
  geom_bin2d(aes(alpha = ..density..), fill = "slategray",binwidth = c(5,1)) +
  scale_fill_continuous(guide = "none")+
  
  theme_bw()

p_E<- ggplot(filter_norm_data_for_plot, aes(order, final_ES,color=factor(score_group)))+geom_line()+theme_bw()+theme_classic()+
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
#dev.off()

################################diff enhancers####################3
#####make the diff table, skip it if you read it as annotated table
#write.table(filter_norm_data_scaled_grouped,"CHD_FPM20_new_grouped",sep="\t",quote=F)
#table <- read.table("CHD_FPM20_new_grouped", header = T, row.names = 1)

#head(table)
table_RNA<-filter_norm_data_nc@now_value
#table_RNA <- table %>% dplyr::select(contains("DNA_log"))
select_REF <- t(table_RNA) %>% as.data.frame() %>% dplyr::select(contains("REF")) %>% t()
select_ALT <- t(table_RNA) %>% as.data.frame() %>% dplyr::select(contains("ALT")) %>% t()


filter_norm_data_nc@annot$type <-NA
filter_norm_data_nc@annot[rownames(select_REF),]$type<-"REF"
filter_norm_data_nc@annot[rownames(select_ALT),]$type<-"ALT"

###########filter data out
intersect_loc <- intersect(filter_norm_data_nc@annot[filter_norm_data_nc@annot$type=="REF","region_name"], 
                           filter_norm_data_nc@annot[filter_norm_data_nc@annot$type=="ALT","region_name"])
new_name<-as.data.frame(filter_norm_data_nc@annot[(filter_norm_data_nc@annot$region_name %in% intersect_loc),
                                                  c("type","region_name")])

#new_REF_name <- rownames(filter_norm_data_nc@annot[(filter_norm_data_nc@annot$region_name %in% intersect_loc & 
#                                                     filter_norm_data_nc@annot$type=="REF"),])
#new_ALT_name <- rownames(filter_norm_data_nc@annot[(filter_norm_data_nc@annot$region_name %in% intersect_loc & 
#                                                      filter_norm_data_nc@annot$type=="ALT"),])
filter_data<-merge(table_RNA,new_name,by=0,all.y=T)
filter_data <- filter_data[!is.na(filter_data$region_name),]

inter_REF<-filter_data[filter_data$type=="REF",]
rownames(inter_REF)<-inter_REF$region_name
inter_REF<-inter_REF[,c(2:5)]
inter_ALT<-filter_data[filter_data$type=="ALT",]
rownames(inter_ALT)<-inter_ALT$region_name
inter_ALT<-inter_ALT[,c(2:5)]


inter_join <- cbind(inter_REF, inter_ALT)
head(inter_join)
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
