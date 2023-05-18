#!/usr/bin/env python
#---coding:utf-8---

####################################################################
# File Name: extend_form_summit.py
# Author: zxr
# mail: zhangxiaoran@picb.ac.cn 
# Created Time: Tue 05 Jan 2020 03:22:12 PM CST
####################################################################

import math, string
import sys, getopt
import os
import re

def usage():
    print("-h --help useange()", "\n")
    print("-i: inputfile", "\n")
    print("-o: outputfile", "\n")
    print("-g: chromsize file", '\n')
    print("-a --fasta: genome fasta file", '\n')
    print("-l: length, the length of the whole peak, default=0, ", '\n')
    print("-f: the float percentange to extend of the whole peak, default:use -l", '\n')
    print("-s: the col number of the summit, default:none", '\n')
    print("-m: the col number of the summit is the length from the start end, default:not", '\n')
    print("-e --enzyme: enzyme sites to motified, E1,E2:newE1,newE2, default:no", '\n')
    print("-t: the input file have a header in the frist line, default:not", '\n')
    print("-k: keep the input file have a header in the frist line, default:not", '\n')
    print("-b --binsize int : the bin for Tiling window", '\n')
    print("--sublength int : the sublength for one Tiling oligo,it must be less than -l",'\n')
    print("--tiling_mutation int: the length of KO region in each tiling region,it must be less than --sublength","\n")
    print("--mutation vcf format file: mutation mode",'\n')
    print("--tiling: tiling mode",'\n')
    print("--only_middle: only for tiling mutation mode,keep the ko region in the oligo middle",'\n')
    print("--keeplength: default not","\n")
def extend_from_summit(i,o,g,l,s,m,f,t,k):
    input=i
    output=o
    genome=g
    length=0
    length=int(l)
    summit=0
    sumcol=0
    sumcol=int(s)
    if sumcol>0:
        summit=1
    mode=m
    middle=0
    extendto=0
    fp=0.0
    fp=float(f)
    if fp==0.0:
        usef=0
    else:
        usef=1
    title=0
    title=int(t)
    keep=0
    keep=int(k)
    A=0
    odd=0
    chromsize={}
    genomefile=open(genome)
    if usef==0:
        half=length/2
        if length%2==1:
            odd=1
    for line in genomefile:
        line=line.rstrip().split()
        chromsize[line[0]]=int(line[1])
    input=open(input,'r', buffering=-1)
    output=open(output, 'w')
    for line in input:
        line=line.rstrip().split()
        if title==1 and A==0:
            A+=1
            if keep==0:
                continue
            else:
                rawregion=':'.join(line)
                output.write("chr\tstar\tend\t"+rawregion+"\n")
                continue
        max=chromsize[line[0]]
        wholelen=int(line[2])-int(line[1])
        if extendto==0:
            if summit==1:
                if mode==1:
                    middle=int(line[sumcol])+int(line[1])
                else:
                    middle=int(line[sumcol])
            else:
                middle=int((int(line[1])+int(line[2]))/2)
            if usef==0:
                left=int(middle-half)
                right=int(middle+half+odd)
            else:
                half=int(wholelen*fp)
                left=int(middle-half)
                right=int(middle+half)
        else:
            if usef==0:
                left=int(line[1])-half
                right=int(line[2])+half+odd
            else:
                half=int(wholelen*fp)
                left=int(line[1])-half
                right=int(line[2])+half
        if left<0:
            ref2=ref[0].split(":")
            left=0
            right=length
        elif right>max:
            right=max
            if usef==0:
                left=max-length
                if left<0:
                    left=0
        rawregion=line[0]+":"+line[1]+"-"+line[2]+":::"+":::".join(line[3:])
        output.write(line[0]+"\t"+str(left)+"\t"+str(right)+"\t"+rawregion+"\n")
    input.close()
    output.close()
def keep_length(input,output,title):
    input=open(input,'r')
    output=open(output,'w')
    t=title
    for line in input:
        if t==1:
            t+=1
            continue
        else:
            line=line.rstrip().split()
            rawregion=line[0]+":"+line[1]+"-"+":".join(line[2:])
            output.write(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+rawregion+"\n")
    input.close()
    output.close()
def getfasta(bed,fasta,output):
    print (fasta,'\n')
    command_line="bedtools getfasta -fi "+fasta+" -fo "+output+" -tab -name -bed "+bed
    print(command_line,'\n')
    os.system(command_line)
def check_enzyme_sites(fasta,output,enzymesites,new_enzymesites):
    new_e=0
    fastafile=open(fasta)
    outfile1=output+"/"+"without_enzymesites.fasta"
    outfile2=output+"/"+"with_enzymesites.fasta"
    if len(new_enzymesites) ==len(enzymesites) :
        outfile3=output+"/"+"with_new_enzymesites.fasta"
        new_e=1
    out1=open(outfile1,'w')
    out2=open(outfile2,'w')
    out3=open(outfile3,'w')
    for line in fastafile:
        with_enzyme=0
        line=line.rstrip().split()
        for site in range(0,len(enzymesites)):
            line[1]=line[1].upper()
            if site==0:
                line_new=line[1]
            x =line[1].find(enzymesites[site])
            if x >=0:
                with_enzyme+=1
                if new_e==1:
                    line_new=line_new.replace(enzymesites[site],new_enzymesites[site])
                #out3.write(str(line[0])+'\t'+str(line_new)+'\n')
                #print (enzymesites[site]+'\t'+new_enzymesites[site]+'\n')
        if with_enzyme>0:
            out2.write(str(line[0])+'\t'+str(line[1])+'\n')
            if new_e==1:
                out3.write(str(line[0])+'\t'+str(line_new)+'\n')
        else:
            out1.write(str(line[0])+'\t'+str(line[1])+'\n')
            if new_e==1:
                out3.write(str(line[0])+'\t'+str(line_new)+'\n')
            
    out1.close()
    out2.close()
    out3.close()

def tiling(bed,windowsize,length,outfile):
    input=open(bed,'r')
    output=open(outfile,'w')
    out=[]
    for line in input:
        line=line.rstrip().split()
        peaklen=int(int(line[2])-int(line[1]))
        if peaklen <= length:
            print ("The region length is less than the tiling region length!!",'\n')
        else:
            i= int((peaklen-length)/windowsize)
            for j in range(0, i+1):
                start=int(line[1])+j*windowsize
                end=start+length
                if end >int(line[2]):
                    end=int(line[2])
                    break
                newline3=line[3]+":::"+line[0]+":"+line[1]+"-"+line[2]+":::"+line[0]+":"+str(start)+"-"+str(end)+":::"+str(j+1)
                output.write(line[0]+'\t'+str(start)+'\t'+str(end)+'\t'+newline3+'\n')

    output.close()
def tiling_mutated(fastafile,output,only_middle,tiling_mutation):
    out1=open(output,'w')
    input=open(fastafile,'r')
    for line in input:
        line=line.rstrip().split()
        peaklen=len(line[1])-tiling_mutation
        halfred=int(peaklen/2)
        bias=peaklen%2
        if peaklen <= 0:
            print ("The tiling mutation length is more than the tiling region length!!",'\n')
        elif only_middle==1:
            #seq1=line[1][0:peaklen].upper()
            #seq2=line[1][peaklen:peaklen+tiling_mutation].upper()
            #seq3=line[1][peaklen+tiling_mutation:].upper()
            seq1=line[1][0:halfred].upper()
            seq2=line[1][halfred:halfred+tiling_mutation].upper()
            seq3=line[1][halfred+tiling_mutation:].upper()
            seq_in="NO"
            out1.write(line[0]+":::"+"W"+":::"+seq2+'\t'+line[1].upper()+"\n")
            out1.write(line[0]+":::"+"M"+":::"+seq_in+'\t'+seq1+seq3+"\n")
       # else :

    input.close()
    out1.close()



    
def mutation(mutationfile,outputbed,outfastafile,length,t,fastafile,outfile):
 
    input=open(mutationfile,'r')
    out1=open(outputbed,'w')
    c=0
    for line in input:
        if t==1 and c==0:
            c+=1
            continue
        else:
            line=line.rstrip().split()
            M=line[3]
            R=line[2]
            if len(M)>len(R):
                mutation_type="I"
                readlen=length-len(M)+len(R)
                halfreadlen=int(readlen/2)
                bias=readlen%2
                start=int(line[1])-halfreadlen
                end=int(line[1])+halfreadlen+bias
                out1.write(line[0]+"\t"+str(start)+'\t'+str(end)+"\t"+line[0]+":"+line[1]+":"+line[2]+":"+":".join(line[3:])+":I"+":::"+line[0]+":"+str(start)+"-"+str(end)+":"+str(halfreadlen)+"\n")
            elif len(M) == len(R):
                mutation_type="R"
                halfreadlen=int(length/2)
                bias=length%2
                start=int(line[1])-halfreadlen
                end=int(line[1])+halfreadlen+bias
                out1.write(line[0]+"\t"+str(start)+'\t'+str(end)+"\t"+line[0]+":"+line[1]+":"+line[2]+":"+":".join(line[3:])+":R"+":::"+line[0]+":"+str(start)+"-"+str(end)+":"+str(halfreadlen)+"\n")
            else:
                mutation_type="D"
                halfreadlen=int(length/2)
                bias=length%2
                start=int(line[1])-halfreadlen
                end=int(line[1])+halfreadlen+bias
                out1.write(line[0]+"\t"+str(start)+'\t'+str(end)+"\t"+line[0]+":"+line[1]+":"+line[2]+":"+":".join(line[3:])+":D"+":::"+line[0]+":"+str(start)+"-"+str(end)+":"+str(halfreadlen)+"\n")
    input.close()
    out1.close()
    getfasta(outbed,fastafile,outfastafile)

    input2=open(outfastafile,'r')
    out2=open(outfile,'w')
    for line in input2:
        line=line.rstrip().split()
        ref=line[0].split(":::")
        ref2=ref[0].split(":")
        ref3=ref[1].split(":")

        M=ref2[3]
        R=ref2[2]
        startpoint=int(ref3[-1])
        endpoint=startpoint+len(R)
        seq1=line[1][0:startpoint].upper()
        seq2=line[1][startpoint:endpoint].upper()
        seq3=line[1][endpoint:].upper()

        WTseq=line[1].upper()
        if seq2.upper() != R.upper() :
            print ("The location of can not be find in the refgenome!",'\n')
            print(line[0],'\t',seq2,'\n')
        else:
            out2.write(line[0]+"\t"+WTseq+"\t"+seq1+M+seq3+"\n")
    out2.close()
    input2.close()
    
if __name__=="__main__":
    opts, args=getopt.getopt(sys.argv[1:], "hi:o:g:l:s:f:me:tka:b:", ["help","fasta=","enzyme=","tiling","binsize=","sublength=","mutation=","tiling_mutation=","only_middle","keeplength"])
    length=0
    fp=0.0
    extendto=0
    title=0
    keep=0
    raw=0
    mode=0
    enzyme=[]
    sum=0
    summit=0
    sumcol=0
    binsize=0
    tilingmode=0
    sublength=0
    mutationfile=''
    tiling_mutation=0
    only_middle=0
    keeplength=0
    for op, value in opts:
        if op in ("-h", "--help"):
            usage()
            sys.exit()
        if op=="-i":
            input=value
        if op=="-o":
            output=value
            if not os.path.isdir(output) :
                os.makedirs(output)
        if op=="-g":
            genome=value
        if op=="-l":
            length=int(value)
        if op=="-s":
            summit=1
            sumcol=int(value)-1
        if op=="-m":
            mode=1
        if op in ("-e","--enzyme"):
            enzyme1=value.split(':')[0]
            enzyme2=value.split(':')[1]
            enzyme=enzyme1.split(',')
            new_enzyme=enzyme2.split(',')
        if op=="-f":
            fp=float(value)
            usef=1
        if op=="-t":
            title=1
        if op=="-k":
            keep=1
        if op in ("-a","--fasta"):
            fastafile=value
        if op=="--tiling":
            tilingmode=1
        if op in ("-b","--binsize"):
            binsize=int(value)
        if op == "--sublength":
            sublength=int(value)
        if op == "--mutation":
            mutationfile=value
        if op == "--tiling_mutation":
            tiling_mutation=int(value)
        if op == "--only_middle":
            only_middle=1
        if op == "--keeplength":
            keeplength=1
            print ("keeplength mode on!","\n")
    outbed=""
    output1=output+"/extended.bed"
    print ("Extend bed is starting...",'\n')
    outbed=output1
    output2=output+"/get_fasta.txt"
    outfasta=output2
    if mutationfile !='':
        #if keeplength==1:
        #    keep_length(input,output1,title)
        #else:
        
        #extend_from_summit(input,output1,genome,length,sumcol,mode,fp,title,keep)
        #print ("Get fasta...",'\n')
        
        output4=output+"/WT_muataion_fasta.txt"
        print("Mutation modeing……",'\n')
        mutation(mutationfile,output1,output2,length,title,fastafile,output4)
        if enzyme != []:
            check_enzyme_sites(output4,output,enzyme)
    else:       
        if keeplength==1:
            keep_length(input,output1,title)
        else:
            extend_from_summit(input,output1,genome,length,sumcol,mode,fp,title,keep)
        print ("Get fasta...",'\n')
        if tilingmode == 1:
            output3=output+"/tiling.bed"
            tiling(output1,binsize,sublength,output3)
            outbed=output3
        getfasta(outbed,fastafile,output2)
        if tiling_mutation>0:
            outtilingfasta=output+"/tiling_fasta.txt"
            tiling_mutated(output2,outtilingfasta,only_middle,tiling_mutation)
            outfasta=outtilingfasta
        if enzyme != []:
            check_enzyme_sites(outfasta,output,enzyme,new_enzyme)

