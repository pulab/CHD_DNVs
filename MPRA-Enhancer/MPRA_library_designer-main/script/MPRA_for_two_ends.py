#!/usr/bin/env python
#---coding:utf-8---

####################################################################
# File Name: MPRA_for_two_ends.py
# Author: Xiaoran Zhang
# mail: xrzhang0525@gmail.com
# Created Time: Mon 23 Nov 2020 12:12:19 AM EST
####################################################################

import math, string
import sys, getopt

def usage():
    print("-h --help useange()", "\n")
    print("-i: inputfile", "\n")
    print("-o: outputfile prefiex, _f and _r files will be the outfiles", "\n")
    print("-1: 5 end length","\n")
    print("-2: 3 end length","\n")
def reverse_seq(site):
    re_seq=''
    seq_len=len(site)
    for i in range(0,seq_len):
        if (site[i]=="A" or site[i] =="a"):
            re_seq+="T"
        elif (site[i]=="T" or site[i] =="t"):
            re_seq+="A"
        elif (site[i]=="G" or site[i] =="g"):
            re_seq+="C"
        elif (site[i]=="C" or site[i] =="c"):
            re_seq+="G"
    re_seq=re_seq[::-1]
    return re_seq
if __name__=="__main__":
    opts, args=getopt.getopt(sys.argv[1:], "hi:o:1:2:", ["help"])
    input=""
    output=""
    l1=0
    l2=0
    for op, value in opts:
        if op in ("-h", "--help"):
            usage()
            sys.exit()
        if op=="-i":
            input=value
        if op=="-o":
            output=value
        if op=="-1":
            l1=int(value)
        if op=="-2":
            l2=int(value)
    out1=output+"_f.fasta"
    out2=output+"_r.fasta"
    input=open(input,'r')
    output1=open(out1,'w')
    output2=open(out2,'w')
    for line in input:
        line=line.rstrip().split()
        line1_new=line[1][0:l1]
        l2_r=-1*l2
        line2_new_reverse=line[1][l2_r:]
        line2_new=reverse_seq(line2_new_reverse)
        output1.write(line[0]+"_f"+"\t"+line1_new+"\n")
        output2.write(line[0]+"_r"+"\t"+line2_new+"\n")
    input.close()
    output1.close()
    output2.close()
        


