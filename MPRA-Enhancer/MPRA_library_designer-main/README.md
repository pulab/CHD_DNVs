# MPRA_library_designer

This is the first part of the MADAP (MPRA Design and Analysis Package).

############ package/tools needed##############

python2.7+, 
packages: Math,string,sys,getopt,os

bedtools v >=2.28.0

###############################################

Download the code and add the bedtools to your PATH.
For MPRA common library, Tiling library, Mutation library, Tiling-mutation library oligo sequence

examples:

python ./design_library_v1.0.py -i input_peaks_example.txt -o test --enzyme ATGC -l 400 -s 4 -g ~/database/hg38/hg38.chrom.sizes --fasta ~/database/hg38/hg38.fasta --tiling -b 5 --sublength 180 --only_middle --tiling_mutation 5

python ./design_library_v1.0.py --mutation mutation_file_example.txt -o test2 --enzyme ATGC -l 170 --fasta hg38.fasta -t

###############################################################
1. All type of library:
Input format: (“\t” delimited files)

![image](https://user-images.githubusercontent.com/66787411/127383828-93042f84-9ea8-40b5-8e17-2d13d6e9d247.png)

1.	Iterms in red are needed!
2.	Iterms in black are chosen. Any features can be added. 
3.	Each feature is a column. If you add a column, NA/0/Unknown is OK, but it can't be left blank.
4.	Do not include ":::" and Spaces in the feature description (For example:  “Negative_control” ✔️，“Negative control” ❌).
5.	The name column is not allowed to be duplicated in the same library.

Annotation file:
![image](https://user-images.githubusercontent.com/66787411/127383926-a80aa9b8-90c5-4c68-86c2-407b60249930.png)

1.	All described features need to be commented.
2.	Note the value if it has NA.
3.	You can add any description of features you want. There are no restrictions on the format of this table.
2. Mutation library:
The first four columns of the file must be in Stand VCF file format,
The description of the format can be found in the link below:
https://www.internationalgenome.org/wiki/Analysis/vcf4.0/

![image](https://user-images.githubusercontent.com/66787411/127384089-76b9df55-dfcb-4ed7-bad6-c36c3c7d56ad.png)

Tips: the input file should be 0 base and the mutation file is 1 base. 

3. Tiling library:

Information needed: 

![image](https://user-images.githubusercontent.com/66787411/127384122-cf43215f-1c77-4e8e-b4f4-68423707073b.png)

Type I:

![image](https://user-images.githubusercontent.com/66787411/127384161-ab9711c6-449b-47f4-951d-3c17378ce312.png)
Type II :

![image](https://user-images.githubusercontent.com/66787411/127384204-f92bfd5a-9af4-4575-b2cc-02d1b7b8962f.png)
Type III:

![image](https://user-images.githubusercontent.com/66787411/127384257-de32dfdc-c247-499c-ac35-087de23f7bec.png)


Contact: xrzhang0525@gmail.com 




