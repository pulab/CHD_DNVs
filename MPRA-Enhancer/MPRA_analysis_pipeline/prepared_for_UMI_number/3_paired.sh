#!/usr/bin/bash
for i in `ls *paired.txt`;do bedForPair.py -i $i -c 1 -o ${i/paired.txt/paired_true.txt};done
