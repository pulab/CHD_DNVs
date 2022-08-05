#!/usr/bin/bash
#for i in `ls *paired_sorted.txt`; do cut -f2-3 $i |sort -u >${i/sorted/uniq}; done
for i in `ls *paired_sorted.txt`; do cut -f2-3 $i |sort |uniq -c >${i/sorted/uniq_number}; done
