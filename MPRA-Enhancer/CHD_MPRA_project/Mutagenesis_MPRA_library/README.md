# Mutagenesis MPRA Analysis

This is all the analysis notebooks and input files for Paper:
Functional dissection of human cardiac enhancers and non-coding de novo variants in congenital heart disease
Feng Xiao, Xiaoran Zhang, Sarah U. Morton, Seong Won Kim, Youfei Fan, Joshua M. Gorham, Huan Zhang, Paul J. Berkson, Neil Mazumdar, Yangpo Cao, Jian Chen, Jacob Hagen, Xujie Liu, Pingzhu Zhou, Felix Richter, Yufeng Shen, Tarsha Ward, Bruce Gelb, Christine E. Seidman, Jonathan G. Seidman, William T. Pu

This is the a part of the MADAP (MPRA Design and Analysis Package).

## Package/tools needed

R >=4.1, 
packages in the codes folder: MPRA_analyzer.R
packages: DESeq2,GenomicRanges,BiocGenerics,dplyr,purrr etc. 

## A schematic overview of MPRA libraries
----------------------------------------------------------------------
### Mutagenesis MPRA Library
The MPRA library analysis for a systematic tiling mutagenesis library.
Mutagenesis_MPRA_library/codes/Mutagenesis_MPRA_analysis.Rmd

<img width="372" alt="image" src="https://github.com/pulab/CHD_DNVs/assets/66787411/a4f14c1a-198b-4da8-b668-76bfe877aaaa">


### Inputfiles
Inputfiles for Mutagenesis_MPRA_analysis.Rmd could be found in the data folder.
The raw motif file as an input for Mutagenesis_motif_analysis.Rmd, it could be downloaded from:
https://zenodo.org/record/8162058
Other input files could be found in data folder.