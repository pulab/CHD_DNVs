# CHD MPRA Analysis

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
### CHD MPRA Library
The MPRA library analysis for human cardiac enhancers and non-coding de novo variants in congenital heart disease.
CHD_MPRA_library/codes/CHD_MPRA.Rmd

<img width="490" alt="image" src="https://github.com/pulab/CHD_DNVs/assets/66787411/9f992c5e-de4f-4932-8adb-dbb226e7416c">

### Inputfiles
Inputfiles for CHD_MPRA.Rmd could be found in the data folder.
The raw motif file as an input for CHD_motif_analysis.Rmd, it could be downloaded from:
https://zenodo.org/record/8162058
Other input files could be found in data folder.


