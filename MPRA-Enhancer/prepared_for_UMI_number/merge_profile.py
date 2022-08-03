#!/usr/bin/env python
#---coding:utf-8---

####################################################################
# File Name: merge_profile.py
# Author: zxr
# mail: zhangxiaoran@picb.ac.cn 
# Created Time: Sun 27 Apr 2014 05:16:38 PM CST
####################################################################

import math, string
import sys, getopt
import os

def sort_profile(output, n, cut, all, keep):
    print ("now sorting", '\n')
    sortsys=''
    for i in n:
        newi=i.split("r")
        if len(newi)>1:
            sortsys+="-k"+newi[0]+","+newi[0]+"nr "
        else:
            sortsys+="-k"+newi[0]+","+newi[0]+"n"+" "
    temp=str(output)+"_temp.txt"
    temp2=str(output)+"_temp2.txt"
    temp3=str(output)+"_temp3.txt"
    temp4=str(output)+"_temp4.txt"
    syst="sort "+sortsys+temp+" >"+temp2
    cutsys=[]
    for i in range(1,all+1):
        if str(i) in cut:
            continue
        else:
            cutsys.append(str(i))
    cutsyss=','.join(cutsys)
    sys="cut -f"+cutsyss+" "+temp2+" >"+temp4
    print (syst, '\n')
    os.system(syst)
    if keep==1 :
        header=['gene_name']
        if cut==[]:
            for i in range(1, all):
                header.append(str(i))
        else:
            for i in range(1, all-len(cut)):
                header.append(str(i))
                
        outfile=open(temp3, 'w')
        outfile.write('\t'.join(header)+'\n')
        outfile.close()
        if cut==[]:
            print ("cat "+temp3+" "+temp2+" >", str(output), '\n')
            os.system("cat "+temp3+" "+temp2+" >"+str(output) )
            os.system("rm -f "+temp2)
            os.system("rm -f "+temp)
            os.system("rm -f "+temp3)
        else:
            print (sys, '\n')
            os.system(sys)
            print ("cat "+temp3+" "+temp4+" >"+str(output), '\n')
            syst3="cat "+temp3+" "+temp4+" >"+str(output) 
            os.system(syst3)
            os.system("rm -f "+temp2)
            os.system("rm -f "+temp)
            os.system("rm -f "+temp3)
            os.system("rm -f "+temp4)
    else:
        if cut==[]:
            os.system("mv "+temp2+" "+str(output) )
            os.system("rm -f "+temp)
    
        else:
            print (sys, '\n')
            os.system(sys)
            os.system("mv "+temp4+" "+str(output) )
            os.system("rm -f "+temp2)
            os.system("rm -f "+temp)
            os.system("rm -f "+temp3)

def usage():
    print("usage:merge_profile.py -i inputfile -p plusfile1, plusfile2, ..., plusfilen -o outputfile", '\n') 
    print("-n ( if you want sort you can choose this parameter,  by the n column, default :not sort, and if you sort , you'd better choose -b 1 )", '\n') 
    print("   if you choose -n, do only one merge_profile.py in one folder", '\n')
    print("-b (1 or 2, 1 rm the frist row, 2 not, and default :2)", '\n')
    print("-k if you choose it, means you choose -b 1 and keep a new sample header like gene_name 1 2 3... in the result, default:not", '\n')
    print("-c cut this column in the result, default: none", '\n')
    print("-1 the column of the key in inputfile, default: 1", '\n')
    print("-2 the column of the key in plusfile, default: 1", '\n')
if __name__=="__main__":
    opts, args=getopt.getopt(sys.argv[1:], "hi:o:p:n:b:c:1:2:k")
    input=""
    output=""
    plus=[]
    n=['1']
    kind=0
    a=1
    b=2
    cut=[]
    c1=0
    c2=0
    keep=0
    for op, value in opts:
        if op=="-h":
            usage()
            sys.exit()
        if op=="-i":
            input=value
        if op=="-o":
            output=value
        if op=="-p":
            plus=value.split(',')
        if op=="-n":
            n=value.split(',')
            kind=1
        if op=="-b":
            b=int(value)
        if op=="-c":
            cut=value.split(",")
        if op=="-1":
            c1=int(value)-1
        if op=="-2":
            c2=int(value)-1
        if op=="-k":
            keep=1
    input=open(input, 'r')
    data={}
    outdata=[]
    all=0
    for i in range(0, len(plus)):
        plusfile=plus[i]
        plusfile=open(plusfile, 'r')
        for line in plusfile:
            line=line.rstrip()
            line1=line.split()
            c=len(line1)-1
            if c<0:
                continue
            data[line1[c2]]=line1[0:c2]
            for x in range(c2+1, len(line1)):
                data[line1[c2]].append(line1[x])
        if a==1:
            for line in input:
                if b==1:
                    b+=1
                    continue
                else:
                    line=line.rstrip()
                    line1=line.split()
                    line='\t'.join(line1)
                    cin=len(line1)
                    if len(line1)==0:
                        continue
                    if line1[c1] in data.keys():
                        for j in range(0, len(data[line1[c1]])):
                            line+='\t'
                            line+=data[line1[c1]][j]
                    else :
                        for j in range(0, c):
                            line+='\t'
                            line+=("0")
                    outdata.append(line)
            all+=cin
        else :
            for k  in range(0, len(outdata)):
                line1=outdata[k].split()
                if line1[c1] in data.keys():
                    for j in range(0, len(data[line1[c1]])):
                        outdata[k]+='\t'
                        outdata[k]+=data[line1[c1]][j]
                else :
                    for j in range(0, c):
                        outdata[k]+='\t'
                        outdata[k]+=("0")
        a+=1
        all+=c
        data={}
        plusfile.close()
    input.close()
    temp=str(output)+"_temp.txt"
    out=open(temp, 'w')
    for i in range (0, len(outdata)):
        out.write(str(outdata[i])+'\n')
    out.close()
    if kind==1:
        sort_profile(output, n, cut, all, keep)
    else:
        if cut==[]:
            os.system("mv "+temp+" "+str(output))
        else:
            cutsys=[]
            for i in range(1,all+1):
                if str(i) in cut:
                    continue
                else:
                    cutsys.append(str(i))
            cutsyss=','.join(cutsys)
            sys="cut -f"+cutsyss+" "+temp+" >"+str(output)
            print (sys, '\n')
            os.system(sys)
            os.system("rm -f "+temp)    



