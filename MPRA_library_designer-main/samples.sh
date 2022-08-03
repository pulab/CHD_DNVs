#!/usr/bin/bash
python ./design_library_v1.0.py -i input_peaks_example.txt -o test --enzyme ATGC -l 400 -s 4 -g ~/database/hg38/hg38.chrom.sizes --fasta ~/database/hg38/hg38.fasta --tiling -b 5 --sublength 180 --only_middle --tiling_mutation 5
python ./design_library_v1.0.py --mutation mutation_file_example.txt -o test2 --enzyme ATGC -l 170 --fasta hg38.fasta -t
