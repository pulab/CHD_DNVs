#!/usr/bin/bash
perl -alne'print "$F[0]\t$F[1]\t$F[3]" if $F[1]=~/chr/; print "$F[0]\t$F[3]\t$F[1]" if $F[3]=~/chr/' MPRA_400bp_A1_R2_UMI_readsname_paired_true.txt >MPRA_400bp_A1_R2_UMI_readsname_paired_sorted.txt &
perl -alne'print "$F[0]\t$F[1]\t$F[3]" if $F[1]=~/chr/; print "$F[0]\t$F[3]\t$F[1]" if $F[3]=~/chr/' MPRA_400bp_A2_R2_UMI_readsname_paired_true.txt >MPRA_400bp_A2_R2_UMI_readsname_paired_sorted.txt &
perl -alne'print "$F[0]\t$F[1]\t$F[3]" if $F[1]=~/chr/; print "$F[0]\t$F[3]\t$F[1]" if $F[3]=~/chr/' MPRA_400bp_A3_R2_UMI_readsname_paired_true.txt >MPRA_400bp_A3_R2_UMI_readsname_paired_sorted.txt &
perl -alne'print "$F[0]\t$F[1]\t$F[3]" if $F[1]=~/chr/; print "$F[0]\t$F[3]\t$F[1]" if $F[3]=~/chr/' MPRA_400bp_A4_R2_UMI_readsname_paired_true.txt >MPRA_400bp_A4_R2_UMI_readsname_paired_sorted.txt &
perl -alne'print "$F[0]\t$F[1]\t$F[3]" if $F[1]=~/chr/; print "$F[0]\t$F[3]\t$F[1]" if $F[3]=~/chr/' MPRA_400bp_A5_R2_UMI_readsname_paired_true.txt >MPRA_400bp_A5_R2_UMI_readsname_paired_sorted.txt &
perl -alne'print "$F[0]\t$F[1]\t$F[3]" if $F[1]=~/chr/; print "$F[0]\t$F[3]\t$F[1]" if $F[3]=~/chr/' MPRA_400bp_V1_R2_UMI_readsname_paired_true.txt >MPRA_400bp_V1_R2_UMI_readsname_paired_sorted.txt &
perl -alne'print "$F[0]\t$F[1]\t$F[3]" if $F[1]=~/chr/; print "$F[0]\t$F[3]\t$F[1]" if $F[3]=~/chr/' MPRA_400bp_V2_R2_UMI_readsname_paired_true.txt >MPRA_400bp_V2_R2_UMI_readsname_paired_sorted.txt &
perl -alne'print "$F[0]\t$F[1]\t$F[3]" if $F[1]=~/chr/; print "$F[0]\t$F[3]\t$F[1]" if $F[3]=~/chr/' MPRA_400bp_V3_R2_UMI_readsname_paired_true.txt >MPRA_400bp_V3_R2_UMI_readsname_paired_sorted.txt &
perl -alne'print "$F[0]\t$F[1]\t$F[3]" if $F[1]=~/chr/; print "$F[0]\t$F[3]\t$F[1]" if $F[3]=~/chr/' MPRA_400bp_V4_R2_UMI_readsname_paired_true.txt >MPRA_400bp_V4_R2_UMI_readsname_paired_sorted.txt &
perl -alne'print "$F[0]\t$F[1]\t$F[3]" if $F[1]=~/chr/; print "$F[0]\t$F[3]\t$F[1]" if $F[3]=~/chr/' MPRA_400bp_V5_R2_UMI_readsname_paired_true.txt >MPRA_400bp_V5_R2_UMI_readsname_paired_sorted.txt &
perl -alne'print "$F[0]\t$F[1]\t$F[3]" if $F[1]=~/chr/; print "$F[0]\t$F[3]\t$F[1]" if $F[3]=~/chr/' MPRA_400bp_DNA1_R2_UMI_readsname_paired_true.txt >MPRA_400bp_DNA1_R2_UMI_readsname_paired_sorted.txt &
perl -alne'print "$F[0]\t$F[1]\t$F[3]" if $F[1]=~/chr/; print "$F[0]\t$F[3]\t$F[1]" if $F[3]=~/chr/' MPRA_400bp_DNA2_R2_UMI_readsname_paired_true.txt >MPRA_400bp_DNA2_R2_UMI_readsname_paired_sorted.txt &
perl -alne'print "$F[0]\t$F[1]\t$F[3]" if $F[1]=~/chr/; print "$F[0]\t$F[3]\t$F[1]" if $F[3]=~/chr/' MPRA_400bp_DNA3_R2_UMI_readsname_paired_true.txt >MPRA_400bp_DNA3_R2_UMI_readsname_paired_sorted.txt &
perl -alne'print "$F[0]\t$F[1]\t$F[3]" if $F[1]=~/chr/; print "$F[0]\t$F[3]\t$F[1]" if $F[3]=~/chr/' MPRA_400bp_DNA4_R2_UMI_readsname_paired_true.txt >MPRA_400bp_DNA4_R2_UMI_readsname_paired_sorted.txt &
perl -alne'print "$F[0]\t$F[1]\t$F[3]" if $F[1]=~/chr/; print "$F[0]\t$F[3]\t$F[1]" if $F[3]=~/chr/' MPRA_400bp_DNA5_R2_UMI_readsname_paired_true.txt >MPRA_400bp_DNA5_R2_UMI_readsname_paired_sorted.txt ;wait
