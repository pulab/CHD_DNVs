#!/usr/bin/env python
#---coding:utf-8---

####################################################################
# File Name: /home/xz191/Bin/program/barcode.py
# Author: Xiaoran Zhang
# mail: xrzhang0525@gmail.com
# Created Time: Tue 17 Nov 2020 03:36:44 PM EST
####################################################################

import math, string
import sys, getopt
import random

def usage():
    print("-h --help useange()", "\n")
    print("-i: inputfile", "\n")
    print("-o: outputfile", "\n")
    print("-l: length of barcode","\n")
    print("-m: mismatch of barcodes,default:0","\n")
    print("--gc: GC percentage cut off,default:0,0~0.5,the gc percentage will be gc%~0.5" ,'\n')

def reverse_seq(seq):
    re_seq=''
    seq_len=len(seq)
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
def f(n,x):

#n为待转换的十进制数，x为机制，取值为2-16
    a=[0,1,2,3,4,5,6,7,8,9,'A','b','C','D','E','F']
    b=[]
    while True:
        s=n//x  # 商
        y=n%x  # 余数
        b=b+[y]
        if s==0:
            break
        n=s
    b.reverse()
    return b


def make_all(length,gc_percentage):
    seq_database=[]
    total_number=int(math.pow(4,length))
    GC_number_A=int(length*gc_percentage)
    GC_number_B=length-GC_number_A

    for i in range(0,total_number):
        a=["A","T","C","G"]
        b=f(i,4)  
        GC_number=0
        candidate=["A"]*length
        for j in range(1,len(b)+1):
            jj=j*-1
            if b[jj] >=2:
                GC_number+=1
            candidate[j-1]=str(a[b[jj]])
        if GC_number >= GC_number_A and GC_number <= GC_number_B:
            seq_database+=["".join(candidate)]
    return seq_database

if __name__=="__main__":
    opts, args=getopt.getopt(sys.argv[1:], "hi:o:m:n:l:", ["help","gc="])
    input=""
    output=""
    barcode_len=0
    mismatch_num=0
    total_num=0
    gc_percentage=0.0
    for op, value in opts:
        if op in ("-h", "--help"):
            usage()
            sys.exit()
        if op=="-i":
            input=value
        if op=="-o":
            output=value
        if op=="-l":
            barcode_len=int(value)
        if op=="-m":
            mismatch_num=int(value)
        if op=="-n":
            total_num=int(value)
        if op=="--gc":
            gc_percentage=float(value)
            if gc_percentage >0.5 or gc_percentage <0:
                sys.stderr.write("The GC% should be :0~0.5!",'\n')
                raise SystemExit(1)
    barcode_num=1+mismatch_num
    barcode_all=[]
    barcode_random=[]
    len1=int(barcode_len/barcode_num)
    len2=barcode_num-int(barcode_len%barcode_num)
    out=open(output,'w')
    for i in range(0,barcode_num):
        barcode_all=[]
        if i<len2:
            barcode_all=make_all(len1,gc_percentage)
            new_list=[]
            if len(barcode_all)<total_num:
                sys.stderr.write("The total number is out of range, make a longer length or a smaller GC%!",'\n')
                #print("The total number is out of range, make a longer length or a smaller GC%!",'\n')
                raise SystemExit(1)
            else:
                barcode_random.append(random.sample(barcode_all,total_num))
        else:
            len3=len1+1
            barcode_all=make_all(len3,gc_percentage)
            barcode_random.append(random.sample(barcode_all,total_num))
    for i in range(0,total_num):
        line=[]
        for j in range(0,barcode_num):
            line.append(str(barcode_random[j][i]))
        line="\t".join(line)
        out.write(line+"\n")
    out.close()
            




