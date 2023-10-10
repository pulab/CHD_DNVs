library(DESeq2)
library(GenomicRanges)
library(BiocGenerics)
library(MASS)
create_MPRA_class <- function(RNA, DNA, annot = NULL,  norm = TRUE) {
  # check if RNA and DNA have the same number of rows
  stopifnot(nrow(RNA) == nrow(DNA))
  stopifnot(identical(row.names(RNA), row.names(DNA)))
  
  # check if annot is a data frame if not null
  if (!is.null(annot) && !is.data.frame(annot)) {
    stop("'annot' must be a data frame")
  }
  
  # check if annot has all row names in DNA if not null
  if (!is.null(annot)) {
    if (!all(rownames(DNA) %in% rownames(annot))) {
      message("Warning: not all row names in DNA are present in annot")
    }
  }else{
    annot <- data.frame(row.names = rownames(DNA))
  }
  # normalize RNA and DNA to FPKM if norm is TRUE
  if (norm) {
    
    rna_fpkms <- as.data.frame(t(t(RNA)/colSums(RNA)*1000000))
    dna_fpkms <- as.data.frame(t(t(DNA)/colSums(DNA)*1000000))
  } else {
    rna_fpkms <- RNA
    dna_fpkms <- DNA
  }
  
  setClass("RNA_DNA",
           slots = list(RNA = "data.frame",
                        DNA = "data.frame",
                        raw_DNA ="data.frame",
                        raw_RNA ="data.frame",
                        annot = "data.frame",
                        res = "data.frame",
                        now_value="data.frame"))
  #res <- data.frame(matrix(NA, nrow=nrow(RNA), ncol=ncol(RNA)), row.names=row.names(RNA))
  #res <- data.frame(matrix(NA, nrow=nrow(RNA)), row.names=row.names(RNA))
  res <- data.frame(row.names = rownames(DNA))
  now_value <- data.frame(row.names = rownames(DNA))
  #history_data <- data.frame(row.names = rownames(DNA))
  rna_dna <- new("RNA_DNA", RNA = rna_fpkms, DNA = dna_fpkms, raw_DNA=DNA,raw_RNA=RNA,
                 annot = annot, res=res, now_value=now_value)
  
  # return filtered DNA
  #filtered_dna <- subset(rna_dna@DNA, count > dna_filter)
  return(rna_dna)
}
filter_data <- function(my_data, dna_filter = 0, filter_by="Min") {
  dna_norm <- my_data@DNA
  rna_norm <- my_data@RNA
  annot_norm <- my_data@annot
  res_norm <- my_data@res
  raw_dna <- my_data@raw_DNA
  raw_rna <- my_data@raw_RNA
  now_value <- my_data@now_value
  #print (dna_filter)
  # Filter DNA data
  if (filter_by =="Min"){
    #dna_filtered<- dna_norm[rowMins(as.matrix(my_data@DNA),value=T) >= dna_filter,]
    #dna_filtered <- dna_norm[rowMins(as.matrix(dna_norm), value = T) >=20,]
    dna_filtered <- dna_norm[rowMins(as.matrix(dna_norm), value = T) >=dna_filter,]
    print(dim(dna_filtered))
  }
  if (filter_by =="Mean"){
    dna_filtered <- dna_norm[rowMeans(as.matrix(dna_norm), value = T) >=dna_filter,]
  }
  if (filter_by =="Median"){
    dna_filtered <- dna_norm[rowMedians(as.matrix(dna_norm), value = T) >=dna_filter,]
  }
  if (filter_by =="Max"){
    dna_filtered <- dna_norm[rowMaxs(as.matrix(dna_norm), value = T) >=dna_filter,]
  }
  if (filter_by =="Sum"){
    dna_filtered <- dna_norm[rowSums(as.matrix(dna_norm), value = T) >=dna_filter,]
  }
  
  # Filter RNA data
  rna_filtered <- subset(rna_norm, rownames(rna_norm) %in% rownames(dna_filtered))
  res_filtered <- subset(res_norm, rownames(res_norm) %in% rownames(dna_filtered))
  annot_filtered <- subset(annot_norm, rownames(annot_norm) %in% rownames(dna_filtered))
  raw_rna_filtered <- subset(raw_rna, rownames(raw_rna) %in% rownames(dna_filtered))
  raw_dna_filtered <- subset(raw_dna, rownames(raw_dna) %in% rownames(dna_filtered))
  now_value_filtered <- subset(now_value, rownames(now_value) %in% rownames(dna_filtered))
  
  # Create a new instance of the RNA_DNA class with the filtered data
  filtered_data <- new("RNA_DNA", RNA = rna_filtered, DNA = dna_filtered, raw_DNA=raw_dna_filtered,
                       raw_RNA=raw_rna_filtered,
                       annot = annot_filtered,res=res_filtered, now_value=now_value_filtered)
  
  # Return the filtered data
  return(filtered_data)
}
####old function for active enhancers
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
######new function for DESeq2 active enhancers
DEGseq2_active <- function(filter_norm_data, fd = 0, qval = 0.05, addname="Active"){
  # the fold is RNA:DNA fold change,and the threshold is fd change > the value, and qvalue < the set value.
  filter_for_degs<-cbind(filter_norm_data@raw_DNA,filter_norm_data@raw_RNA)
  coldata <-as.data.frame(colnames(filter_for_degs))
  coldata$type <- coldata$`colnames(filter_for_degs)`
  #coldata$condition <-gsub("\\d","",coldata$type)
  coldata$condition <- c(rep("DNA", ncol(filter_norm_data@raw_DNA)), 
                         rep("RNA", ncol(filter_norm_data@raw_RNA)))
  genes=as.matrix(filter_for_degs)
  dds <- DESeqDataSetFromMatrix(countData = genes, colData=coldata, design=~condition)
  dds <- DESeq(dds)
  res <- results( dds, contrast = c("condition", "RNA", "DNA") )
  #res2 <- results( dds, contrast = c("condition","d17_DNA_2","d17_RNA_2") )
  res$qvalue_0_05 <- 0
  res$qvalue_0_05[res$padj < qval & res$log2FoldChange >fd] <- 1
  filter_norm_data@res[,addname]<-res$qvalue_0_05
  return (filter_norm_data)
  
}

get_ratio <-function(data, DNA_merge="Paired"){
  ncols <- length(data@RNA[1,])
  if (DNA_merge=="Paired"){
    stopifnot(ncol(data@RNA) == ncol(data@DNA))
    for (i in 1:ncols){
      ratio <- (data@RNA[,i]+1)/(data@DNA[,i]+1)
      ratio <- log2(ratio)
      colname_i<-paste0(colnames(data@RNA)[i], "_", "ratio")
      data@res <- cbind(data@res,ratio)
      names(data@res)[ncol(data@res)] <- colname_i
      data@now_value <-cbind(data@now_value,ratio)
      names(data@now_value)[ncol(data@now_value)] <- colname_i
      }
  }
  if (DNA_merge=="Median"){
    for (i in 1:ncols){
      data@DNA$median <-rowMedians(as.matrix(data@DNA))
      ratio <- (data@RNA[,i]+1)/(data@DNA$median+1)
      ratio <- log2(ratio)
      colname_i<-paste0(colnames(data@RNA)[i], "_", "ratio")
      data@res <- cbind(data@res,ratio)
      names(data@res)[ncol(data@res)] <- colname_i
      data@now_value <-cbind(data@now_value,ratio)
      names(data@now_value)[ncol(data@now_value)] <- colname_i

    }
  }
  if (DNA_merge=="Mean"){
    for (i in 1:ncols){
      data@DNA$mean <-rowMedians(as.matrix(data@DNA))
      ratio <- (data@RNA[,i]+1)/(data@DNA$mean+1)
      ratio <- log2(ratio)
      colname_i<-paste0(colnames(data@RNA)[i], "_", "ratio")
      data@res <- cbind(data@res,ratio)
      names(data@res)[ncol(data@res)] <- colname_i
      data@now_value <-cbind(data@now_value,ratio)
      names(data@now_value)[ncol(data@now_value)] <- colname_i
    }
  }
  if (DNA_merge=="Sum"){
    for (i in 1:ncols){
      data@DNA$sum <-rowSums(as.matrix(data@DNA))
      ratio <- (data@RNA[,i]+1)/(data@DNA$sum+1)
      ratio <- log2(ratio)
      colname_i<-paste0(colnames(data@RNA)[i], "_", "ratio")
      data@res<- cbind(data@res,ratio)
      names(data@res)[ncol(data@res)] <- colname_i
      data@now_value <-cbind(data@now_value,ratio)
      names(data@now_value)[ncol(data@now_value)] <- colname_i
    }
  }
  return(data)
}
#############this function is used to do a size factor scaled
MPRA_SizeFactors <- function(datas,groups,cols){
  datasets <- datas@now_value
  filter_norm_data_nc <- datasets[groups,cols] %>% as.matrix()
  filter_norm_data_full_nc <- datasets[,cols] %>% as.matrix()
  filter_norm_data_nc <- filter_norm_data_nc[rowMins(filter_norm_data_nc,value = T)>0,]
  filter_norm_data_nc_ratio_pseudo_exp <- geometricmeanRow(filter_norm_data_nc)
  filter_norm_data_nc_ratio_scaled <- t(t(filter_norm_data_nc)/filter_norm_data_nc_ratio_pseudo_exp)
  filter_norm_data_nc_normalize_factor <- colMedians(filter_norm_data_nc_ratio_scaled)
  print (filter_norm_data_nc_normalize_factor)
  filter_norm_data_nc_normalize_exp <- t(t(filter_norm_data_full_nc)/filter_norm_data_nc_normalize_factor)
  colnames(filter_norm_data_nc_normalize_exp)<-paste0(colnames(filter_norm_data_nc_normalize_exp),"_","scaled")
  datas@res<-cbind(datas@res,filter_norm_data_nc_normalize_exp)  
  datas@now_value <- as.data.frame(filter_norm_data_nc_normalize_exp)
  return (datas)
}

###########this function is used to get correlation bewteen samples
##method should be chosen from c("pearson", "kendall", "spearman")
scatter_plot_cor <- function(data, cor_method="pearson"){
  my_cols <- c("black") 
  panel.cor <- function(x, y,corr_method){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- round(cor(x, y,method=corr_method), digits=2)
    txt <- paste0("R = ", r)
    cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor*1 )
  }
  
  # Customize upper panel
  upper.panel<-function(x, y){
    points(x,y, pch = 19, col = my_cols[iris$Species],cex=0.1)
  }
  
  pairs(filter_norm_data_nc_normalize_exp, 
        panel = function(x, y){
          panel.cor(x, y, corr_method = cor_method)
        },
        upper.panel = upper.panel)
}

########FDR method
FDR_negative <- function(Negative_C,colname,distribution="norm",fdr_value=0.05){
  th_value=1-fdr_value
  fitg<-fitdist(Negative_C[,colname] ,distribution)
  
  FDR_0.95 <- qnorm(0.95, mean = fitg$estimate["mean"], sd = fitg$estimate["sd"])
  return(FDR_0.95)
}
library(jmuOutlier) 
#library(distr6)
Enrich_score <- function(datasets, groups_col , value_col , order_col,random=1000){
  #set.seed(11)
  ALL <- datasets[,value_col]
  Max <- length(ALL)
  datasets$Enrich_score <-0
  datasets$group_order <-0
  pvalue <- c()
  final_ES_merged <-c()
  #random_final_E <- data.frame()
  #random_final_E <- data.frame(matrix(nrow = random, ncol = 0))
  groups_names <-c()
  rows_number <-length(unique(datasets[,groups_col]))
  datasets2 <- datasets[FALSE,]
  datasets3 <- datasets
  all_random_ES <- data.frame(matrix(nrow = random, ncol = 0))
  for (cols_name in unique(datasets[,groups_col])) {
    datasets3 <- datasets
    PC1 <- datasets3[datasets3[,groups_col]==cols_name,]
    #p_PC1 <- perm.test(PC1[,value_col],ALL)$p.value
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
    random_final_E <-c()
    for (num in c(1:random)){
      #random_set<-datasets3
      selected_regions <- sort(sample(datasets3[,order_col], size = Max_pc1, replace = FALSE))
      selected_E<-as.data.frame(selected_regions)
      colnames(selected_E)<-"order_col"
      selected_E$group_order <-0
      selected_E$group_order <- rank(as.numeric(as.character(selected_E[,"order_col"])),ties.method = "first")
      Max_order_r<-max(selected_E[,"group_order"])
      Max_pc1_r <- length(selected_E[,"order_col"])

      selected_E$Enrich_score <- selected_E[,"group_order"]/Max_pc1_r - selected_E[,"order_col"]/Max
      random_final_E_num<-mean(selected_E$Enrich_score)
      #random_final_E_num<-sum(selected_E$Enrich_score)
      #random_final_E_num<- selected_E$Enrich_score[which.max(abs(selected_E$Enrich_score))]
      random_final_E <- c(random_final_E,random_final_E_num)
    }
    plotPdata<-as.data.frame(random_final_E)
    
    #ggplot(plotPdata, aes(random_final_E))+geom_histogram(bins = 0.05)+
    #  geom_vline(xintercept = final_ES)
    #  theme_classic()
    all_random_ES[,cols_name]<-random_final_E
    
    groups_names <- c(groups_names,cols_name)
    datasets3$score_group <- cols_name
    #datasets3$final_ES <- datasets3$Enrich_score - datasets3$negative_Enrich_score
    datasets3$final_ES <- (datasets3$Enrich_score) 
    final_ES <-mean(datasets3[datasets3[,groups_col]==cols_name,]$Enrich_score)
    
    if (final_ES <0 ){
      p_test <- (sum(final_ES >= random_final_E))/random
    }
    if (final_ES >=0 ){
      p_test <- (sum(final_ES <= random_final_E))/random
    }
    
    p_PC1 <- p_test
    pvalue <- c(pvalue, p_PC1)
    final_ES_merged<-c(final_ES_merged,final_ES)

    
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
    score_group = groups_names,
    final_ES =final_ES_merged
  )
  
  newlist <-list("dat_text"=dat_text, "datasets"=datasets2,"random_ES"=all_random_ES)
  return(newlist)
}

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
split_for_pair <- function(df,condition,vals=null,group_by,kept_annot=T,kept_res=T) {
  ####### this function is used to split the data by condition
  ###### then merged the sub datasets by "group_by", only paired data will be kept 
  ###### kept_annot/kept_res could be a factor or "TRUE/FALSE"
  library(dplyr)
  library(purrr)
  if (is.null(vals)){
    used_values<-levels(df@annot$condition)
  }
  else{
    used_values <- vals
  }
  all_Merge<-df@now_value
  all_Merge$group_by<-df@annot[,group_by]
  all_Merge$condition<-df@annot[,condition]
  value_columns <- colnames(df@now_value)
  if (kept_annot !=FALSE ){
    if (kept_annot==TRUE){
      values_to_remove<-c(group_by,condition)
      other_annot<-subset(colnames(df@annot), !colnames(df@annot) %in% values_to_remove)
      if (! is.null(length(other_annot))){
        all_Merge[,other_annot] <- df@annot[,other_annot]
        value_columns<-c(value_columns,other_annot)
      }else{
        message("No other annot values found.")
      }
    }
    else {
      if (! is.null(length(kept_annot))) {
      all_Merge[,kept_annot] <- df@annot[,kept_annot]
      value_columns<-c(value_columns,kept_annot)
      if (! kept_annot %in% (colnames(df@annot))){
        message("Please check res values! ")
      }
    }
    }
  }

  if (kept_res !=FALSE ){
    if (kept_res==TRUE){
      #values_to_remove<-c(group_by,condition)
      other_res<-colnames(df@res)
      if (! is.null(length(other_res))){
        all_Merge[,other_res] <- df@res[,other_res]
        value_columns<-c(value_columns,other_res)
      } 
      else{
        message("No res values found.")
      }
    }
    else{
      if (! is.null(length(kept_res))) {
        all_Merge[,kept_res] <- df@res[,kept_res]
        value_columns<-c(value_columns,kept_res)
        if (! kept_res %in% (colnames(df@res))){
          message("Please check res values! ")
        }
      }
    }
  }
   

  #new_colnames <- paste(value_columns, unique(all_Merge$condition), sep = "_")
  #names(all_Merge)[startsWith(names(all_Merge), "value")] <- new_colnames
  new_colnames <- expand.grid(value_columns, levels(factor(all_Merge$"condition")))
  concatenated <- with(new_colnames, paste(Var1, Var2, sep = "_"))
  # Split the filtered data.frame by 'condition' and rename only the value columns
  
  split_dfs <- all_Merge %>%
    filter(condition %in% used_values) %>%
    group_split(condition) %>%
    map(~ {
      condition_val <- unique(.x$condition)
      .x %>%
        rename_with(~ if_else(.x %in% value_columns, paste0(.x, "_", condition_val), .x))
    }) 
  #print (split_dfs)
  # Find the shared 'group_by' values across all sub data.frames
  shared_group_values <- Reduce(intersect, map(split_dfs, function(sub_df) unique(sub_df$group_by)))
  
  # Check if there are any shared group values
  if (length(shared_group_values) > 0) {
    # Filter the sub data.frames to keep only the shared 'group_by' values
    filtered_dfs <- split_dfs %>%
      map(~ .x %>%
            filter(group_by %in% shared_group_values))
    
    
    # Merge the filtered sub data.frames by 'group_by' column
    merged_df <- bind_cols(filtered_dfs, .id = "group_by")
    group_by_names <- colnames(merged_df)[grep ("group_by",colnames(merged_df))]
    # Set 'group_by' as new row names
    merged_df<- as.data.frame(merged_df)
    rownames(merged_df) <- merged_df[,group_by_names[1]]
    merged_df <- merged_df[,colnames(merged_df) %in% concatenated]
    # View the resulting merged data.frame
    return (merged_df)
  } else {
    message("No shared 'group_by' values found.")
    NULL
  }
}

library(MPRAnalyze)
trans_to_MPRAnalyze<-function(data_norm,condition_DNA=NULL,batch_DNA=NULL,barcode_DNA=NULL,
                              condition_RNA=NULL,batch_RNA=NULL,barcode_RNA=NULL, control=NULL,...)
{
  nDNA<-length(data_norm@raw_DNA)
  region<-dim(data_norm@raw_DNA)[1]
  nRNA<-length(data_norm@raw_RNA)
  controls<-c(rep(FALSE,region))
  colAnnot <- data.frame(batch=c(rep(1,nDNA)),
                         condition=c(rep("samples",nDNA)),
                         barcode=c(rep(1,nDNA))
  )
  colAnnot2 <- data.frame(batch=c(rep(1,nRNA)),
                          conditon=c(rep("samples",nRNA)),
                          barcode=c(rep(1,nRNA))
  )
  rownames(colAnnot)<-colnames(data_norm@raw_DNA)
  rownames(colAnnot2)<-colnames(data_norm@raw_RNA)
  
  if (!is.null (condition_DNA)){
    colAnnot$condition<-condition_DNA
  }
  if (!is.null (batch_DNA)){
    colAnnot$batch<-batch_DNA
  }
  if (!is.null (barcode_DNA)){
    colAnnot$barcode<-barcode_DNA
  }
  if (!is.null (condition_RNA)){
    colAnnot2$condition<-condition_RNA
  }
  if (!is.null (batch_RNA)){
    colAnnot2$batch<-batch_RNA
  }
  if (!is.null (barcode_RNA)){
    colAnnot2$barcode<-barcode_RNA
  }
  if (!is.null (control)){
    controls<-control
  }
  obj <- MpraObject(dnaCounts = as.matrix(data_norm@raw_DNA), rnaCounts = as.matrix(data_norm@raw_RNA), 
                    dnaAnnot = colAnnot, rnaAnnot = colAnnot2, 
                    controls = controls,...)
  return(obj)
}

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

Get_odds_ratio <-function(data, guess_cols, background_cols,min_value=0) {
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
      res_list_odds[rownames(Odd_input[Odd_input[,cols_name]<= min_value,]),1]<-NA
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








