#!/bin/bash

out="TesasSemble_contigs"
in="merged_fastas"
python /Users/cameronmartino/bin/TesasSemble/scripts/run_optim.py -i fraction -e 5 --optim-type simulated_annealing --alpha=0.5 -r 250 --k 3 --files $in/*.fasta --output-dir $out/$kmerlen

