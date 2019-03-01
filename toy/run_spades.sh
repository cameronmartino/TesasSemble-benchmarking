#!/bin/bash

cat fastqs/*_R1_001.fastq > fastqs/pooled_L001_R1_001.fastq 

cat fastqs/*_R2_001.fastq > fastqs/pooled_L001_R2_001.fastq 

/Users/cameronmartino/bin/SPAdes-3.13.0-Darwin/bin/spades.py --meta -k 21,33,55,77 -o /Users/cameronmartino/bin/TesasSemble-benchmarking/toy/assembly/assemble_s1 -1 /Users/cameronmartino/bin/TesasSemble-benchmarking/toy/fastqs/S0_L001_R1_001.fastq -2 /Users/cameronmartino/bin/TesasSemble-benchmarking/toy/fastqs/S0_L001_R2_001.fastq --phred-offset 33

/Users/cameronmartino/bin/SPAdes-3.13.0-Darwin/bin/spades.py --meta -k 21,33,55,77 -o /Users/cameronmartino/bin/TesasSemble-benchmarking/toy/assembly/assemble_s2 -1 /Users/cameronmartino/bin/TesasSemble-benchmarking/toy/fastqs/S1_L001_R1_001.fastq -2 /Users/cameronmartino/bin/TesasSemble-benchmarking/toy/fastqs/S1_L001_R2_001.fastq --phred-offset 33

/Users/cameronmartino/bin/SPAdes-3.13.0-Darwin/bin/spades.py --meta -k 21,33,55,77 -o /Users/cameronmartino/bin/TesasSemble-benchmarking/toy/assembly/assemble_pool -1 /Users/cameronmartino/bin/TesasSemble-benchmarking/toy/fastqs/pooled_L001_R1_001.fastq -2 /Users/cameronmartino/bin/TesasSemble-benchmarking/toy/fastqs/pooled_L001_R2_001.fastq --phred-offset 33