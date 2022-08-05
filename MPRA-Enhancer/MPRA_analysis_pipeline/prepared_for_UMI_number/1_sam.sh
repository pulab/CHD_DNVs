#!/usr/bin/bash
#for i in `ls *.sam`;do perl -alne' next if $F[0]=~/^@/; print "$F[0]\t$F[2]" if $F[2] ne "*" ' $i|sort |uniq -c |perl -alne 'print "$F[1]\t$F[2]" if $F[0]==2' >${i/.sam/_readsname.txt}; done
#for i in `ls *.sam`;do perl -alne' next if $F[0]=~/^@/; print "$F[0]\t$F[2]" if ($F[2] ne "*" and $F[5] ne "*")' $i >${i/.sam/_sam.txt}; done
sort MPRA_400bp_A1_sam.txt |uniq -c |perl -alne'print "$F[1]\t$F[2]"if $F[0]==2' >MPRA_400bp_A1_readsname.txt &
sort MPRA_400bp_A2_sam.txt |uniq -c |perl -alne'print "$F[1]\t$F[2]" if $F[0]==2' >MPRA_400bp_A2_readsname.txt &
sort MPRA_400bp_A3_sam.txt |uniq -c |perl -alne'print "$F[1]\t$F[2]" if $F[0]==2' >MPRA_400bp_A3_readsname.txt &
sort MPRA_400bp_A4_sam.txt |uniq -c |perl -alne'print "$F[1]\t$F[2]" if $F[0]==2' >MPRA_400bp_A4_readsname.txt &
sort MPRA_400bp_A5_sam.txt |uniq -c |perl -alne'print "$F[1]\t$F[2]" if $F[0]==2' >MPRA_400bp_A5_readsname.txt &
sort MPRA_400bp_V1_sam.txt |uniq -c |perl -alne'print "$F[1]\t$F[2]" if $F[0]==2' >MPRA_400bp_V1_readsname.txt &
sort MPRA_400bp_V2_sam.txt |uniq -c |perl -alne'print "$F[1]\t$F[2]" if $F[0]==2' >MPRA_400bp_V2_readsname.txt &
sort MPRA_400bp_V3_sam.txt |uniq -c |perl -alne'print "$F[1]\t$F[2]" if $F[0]==2' >MPRA_400bp_V3_readsname.txt &
sort MPRA_400bp_V4_sam.txt |uniq -c |perl -alne'print "$F[1]\t$F[2]" if $F[0]==2' >MPRA_400bp_V4_readsname.txt &
sort MPRA_400bp_V5_sam.txt |uniq -c |perl -alne'print "$F[1]\t$F[2]" if $F[0]==2' >MPRA_400bp_V5_readsname.txt &
sort MPRA_400bp_DNA1_sam.txt |uniq -c |perl -alne'print "$F[1]\t$F[2]"if $F[0]==2' >MPRA_400bp_DNA1_readsname.txt &
sort MPRA_400bp_DNA2_sam.txt |uniq -c |perl -alne'print "$F[1]\t$F[2]" if $F[0]==2' >MPRA_400bp_DNA2_readsname.txt &
sort MPRA_400bp_DNA3_sam.txt |uniq -c |perl -alne'print "$F[1]\t$F[2]" if $F[0]==2' >MPRA_400bp_DNA3_readsname.txt &
sort MPRA_400bp_DNA4_sam.txt |uniq -c |perl -alne'print "$F[1]\t$F[2]" if $F[0]==2' >MPRA_400bp_DNA4_readsname.txt &
sort MPRA_400bp_DNA5_sam.txt |uniq -c |perl -alne'print "$F[1]\t$F[2]" if $F[0]==2' >MPRA_400bp_DNA5_readsname.txt ;wait
