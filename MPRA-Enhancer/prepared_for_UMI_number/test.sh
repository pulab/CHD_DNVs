#!/usr/bin/bash
cut -f1 MPRA_400bp_A1_R2_UMI_readsname_paired_sorted.txt |sort -u |wc -l >test.txt
