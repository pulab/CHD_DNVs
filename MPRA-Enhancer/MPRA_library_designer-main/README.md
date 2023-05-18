# MPRA_library_designer

This is the first part of the MADAP (MPRA Design and Analysis Package).

############ package/tools needed##############

python2.7+, 
packages: Math,string,sys,getopt,os

bedtools v >=2.28.0

###############################################
A schematic overview of MPRA_library_designer
----------------------------------------------------------------------
![MPRA_library_designer](https://github.com/pulab/CHD_DNVs/assets/66787411/7939882a-bf7e-4722-9b48-542bba6533c7)
########################

## Installation

Download the code and add the bedtools to your PATH.
For MPRA common library, Tiling library, Mutation library, Tiling-mutation library oligo sequence

## Usage and example

```
Usage: python design_library_v1.0.py -i Inputfile -o outputfile -g chromesizefile -a fasta [options]
-h --help useange()

-i: inputfile

-o: outputfile

-g: chromsize file

-a --fasta: genome fasta file

-l: length, the length of the whole peak, default=0,

-f: the float percentange to extend of the whole peak, default:use -l

-s: the col number of the summit, default:none

-m: the col number of the summit is the length from the start end, default:not

-e --enzyme: enzyme sites to motified, E1,E2:newE1,newE2, default:no

-t: the input file have a header in the frist line, default:not

-k: keep the input file have a header in the frist line, default:not

-b --binsize int : the bin for Tiling window

--sublength int : the sublength for one Tiling oligo,it must be less than -l

--tiling_mutation int: the length of KO region in each tiling region,it must be less than --sublength

--mutation vcf format file: mutation mode

--tiling: tiling mode,default Fasle

--only_middle: only for tiling mutation mode,keep the ko region in the oligo middle

--keeplength: default not
```
examples:
```bash
python ./design_library_v1.0.py -i input_peaks_example.txt -o test --enzyme ATGC -l 400 -s 4 -g ~/database/hg38/hg38.chrom.sizes --fasta ~/database/hg38/hg38.fasta --tiling -b 5 --sublength 180 --only_middle --tiling_mutation 5
```
python ./design_library_v1.0.py --mutation mutation_file_example.txt -o test2 --enzyme ATGC -l 170 --fasta hg38.fasta -t
```
###############################################################
##1. All type of library:
###Input format: (“\t” delimited files)
----------------------------------------------------------------------
![image](https://user-images.githubusercontent.com/66787411/127383828-93042f84-9ea8-40b5-8e17-2d13d6e9d247.png)

1.	Iterms in red are needed!
2.	Iterms in black are chosen. Any features can be added. 
3.	Each feature is a column. If you add a column, NA/0/Unknown is OK, but it can't be left blank.
4.	Do not include ":::" and Spaces in the feature description (For example:  “Negative_control” ✔️，“Negative control” ❌).
5.	The name column is not allowed to be duplicated in the same library.

###Annotation file:
![image](https://user-images.githubusercontent.com/66787411/127383926-a80aa9b8-90c5-4c68-86c2-407b60249930.png)
----------------------------------------------------------------------
1.	All described features need to be commented.
2.	Note the value if it has NA.
3.	You can add any description of features you want. There are no restrictions on the format of this table.
##2. Mutation library:
The first four columns of the file must be in Stand VCF file format,
The description of the format can be found in the link below:
https://www.internationalgenome.org/wiki/Analysis/vcf4.0/
----------------------------------------------------------------------
![image](https://user-images.githubusercontent.com/66787411/127384089-76b9df55-dfcb-4ed7-bad6-c36c3c7d56ad.png)

Tips: the input file should be 0 base and the mutation file is 1 base. 

##3. Tiling library:

###Information needed: 
----------------------------------------------------------------------
![image](https://user-images.githubusercontent.com/66787411/127384122-cf43215f-1c77-4e8e-b4f4-68423707073b.png)

###Type I:
----------------------------------------------------------------------
![image](https://user-images.githubusercontent.com/66787411/127384161-ab9711c6-449b-47f4-951d-3c17378ce312.png)
###Type II :
----------------------------------------------------------------------
![image](https://user-images.githubusercontent.com/66787411/127384204-f92bfd5a-9af4-4575-b2cc-02d1b7b8962f.png)
###Type III:

![image](https://user-images.githubusercontent.com/66787411/127384257-de32dfdc-c247-499c-ac35-087de23f7bec.png)


Contact: xrzhang0525@gmail.com 




