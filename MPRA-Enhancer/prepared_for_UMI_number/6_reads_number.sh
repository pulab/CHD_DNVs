#!/usr/bin/bash
#for i in `ls *paired_sorted.txt`;do cut -f2 $i|sort |uniq -c |perl -alne'print "$F[1]\t$F[0]"' >${i/paired_sorted/barcode_number}; done
for i in `ls *_paired_uniq.txt`;do cut -f1 $i|sort |uniq -c |perl -alne'print "$F[1]\t$F[0]"' >${i/paired_uniq/UMI_number}; done
