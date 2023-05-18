#!/usr/bin/bash
perl -alne'print if $F[0]=~/:::S:::/' output_of_designer/with_new_enzymesites.fasta |random_selection.pl --file - --num 200|perl -alne 'print ">$F[0]\n$F[1]"' >S_random_200.fasta
## fasta-dinucleotide-shuffle-py3 https://meme-suite.org/meme/doc/fasta-dinucleotide-shuffle.html
fasta-dinucleotide-shuffle-py3 -f S_random_200.fasta -t _shuffled -s 1 -c 1 > S_random_200_shuffled.fasta
