#!/usr/bin/bash
cat MPRA_400bp_A1_R2_UMI_10bp.txt MPRA_400bp_A1_readsname.txt |sort -k1,1 >MPRA_400bp_A1_R2_UMI_readsname_paired.txt &
cat MPRA_400bp_A2_R2_UMI_10bp.txt MPRA_400bp_A2_readsname.txt |sort -k1,1 >MPRA_400bp_A2_R2_UMI_readsname_paired.txt &
cat MPRA_400bp_A3_R2_UMI_10bp.txt MPRA_400bp_A3_readsname.txt |sort -k1,1 >MPRA_400bp_A3_R2_UMI_readsname_paired.txt &
cat MPRA_400bp_A4_R2_UMI_10bp.txt MPRA_400bp_A4_readsname.txt |sort -k1,1 >MPRA_400bp_A4_R2_UMI_readsname_paired.txt &
cat MPRA_400bp_A5_R2_UMI_10bp.txt MPRA_400bp_A5_readsname.txt |sort -k1,1 >MPRA_400bp_A5_R2_UMI_readsname_paired.txt &
cat MPRA_400bp_V1_R2_UMI_10bp.txt MPRA_400bp_V1_readsname.txt |sort -k1,1 >MPRA_400bp_V1_R2_UMI_readsname_paired.txt &
cat MPRA_400bp_V2_R2_UMI_10bp.txt MPRA_400bp_V2_readsname.txt |sort -k1,1 >MPRA_400bp_V2_R2_UMI_readsname_paired.txt &
cat MPRA_400bp_V3_R2_UMI_10bp.txt MPRA_400bp_V3_readsname.txt |sort -k1,1 >MPRA_400bp_V3_R2_UMI_readsname_paired.txt &
cat MPRA_400bp_V4_R2_UMI_10bp.txt MPRA_400bp_V4_readsname.txt |sort -k1,1 >MPRA_400bp_V4_R2_UMI_readsname_paired.txt &
cat MPRA_400bp_V5_R2_UMI_10bp.txt MPRA_400bp_V5_readsname.txt |sort -k1,1 >MPRA_400bp_V5_R2_UMI_readsname_paired.txt &
cat MPRA_400bp_DNA1_R2_UMI_10bp.txt MPRA_400bp_DNA1_readsname.txt |sort -k1,1 >MPRA_400bp_DNA1_R2_UMI_readsname_paired.txt &
cat MPRA_400bp_DNA2_R2_UMI_10bp.txt MPRA_400bp_DNA2_readsname.txt |sort -k1,1 >MPRA_400bp_DNA2_R2_UMI_readsname_paired.txt &
cat MPRA_400bp_DNA3_R2_UMI_10bp.txt MPRA_400bp_DNA3_readsname.txt |sort -k1,1 >MPRA_400bp_DNA3_R2_UMI_readsname_paired.txt &
cat MPRA_400bp_DNA4_R2_UMI_10bp.txt MPRA_400bp_DNA4_readsname.txt |sort -k1,1 >MPRA_400bp_DNA4_R2_UMI_readsname_paired.txt &
cat MPRA_400bp_DNA5_R2_UMI_10bp.txt MPRA_400bp_DNA5_readsname.txt |sort -k1,1 >MPRA_400bp_DNA5_R2_UMI_readsname_paired.txt ;wait
