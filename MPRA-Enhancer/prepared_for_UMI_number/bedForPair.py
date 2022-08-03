#!/usr/bin/env python
#---coding:utf-8---

####################################################################
# File Name: bedForPair.py
# Author: zxr
# mail: zhangxiaoran@picb.ac.cn 
# Created Time: Sun 31 Aug 2014 02:20:30 PM CST
####################################################################

import math, string
import sys, getopt
def usage():
    print("-h useange()", '\n')
    print("-i inputfile, merged filename,sorted by the name", '\n')
    print("-o outputfile, file name", '\n')
    print("-c columon, number,default=4","\n")
def forpair2(inputfile, output,c):
    colnum=c-1
    input=open(inputfile,'r')
    a1=0
    output=open(output, 'w')
    line2=[]
    id2=""
    for line in input:
        line1=line.rstrip().split()
        if a1==0:
            id2=line1[colnum]
            line2=line1
            a1+=1
        else:
            id1=line1[colnum]
            if id1==id2:
                str1="\t".join(line1)
                str2="\t".join(line2)
                if len(line2)>len(line1):
                    output.write(str2+"\t"+str1+"\n")
                    a1=0
                    line1=[]
                    line2=[]
                else:    
                    output.write(str1+"\t"+str2+"\n")
                    a1=0
                    line1=[]
                    line2=[]
            else:
                a1=1
                id2=line1[colnum]
                line2=line1
        continue
    input.close()
    output.close()

if __name__=="__main__":
    opts, args=getopt.getopt(sys.argv[1:], "hi:o:c:")
    input=""
    output=""
    c=3
    for op, value in opts:
        if op=="-h":
            usage()
            sys.exit()
        if op=="-i":
            input=value
        if op=="-o":
            output=value
        if op=="-c":
            c=int(value)
    forpair2(input,output,c)


