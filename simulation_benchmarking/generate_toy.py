import biom
import skbio
import itertools
import numpy as np
import pandas as pd
from toy_sims import subsample_reads,generate_strain_reads

## import contigs ## 
biom_path = 'time_species_resolution.biom'
total_reads = 5000

## import contigs ## 
g1 = 'genomes/GCA_000013785.1_ASM1378v1.fa' #Pseudomonas stutzeri
g2 = 'genomes/GCA_000816805.1_ASM81680v1.fa' #Bacillus subtilis
g3 = 'genomes/GCA_001598635.1_ASM159863v1.fa' #Enterococcus faecalis
g4 = 'genomes/GCA_000271365.1_ASM27136v1.fa' #Pseudomonas aeruginosa
genomes = [g1,g2,g3,g4]

## import contigs ## 
contigs = []
for gpath in genomes:
    contig = ''.join([line.rstrip('\n') for line in \
                    open(gpath)][1:])
    contigs.append(contig)
    
## import biom ##    
table = biom.load_table(biom_path)
table = table.to_dataframe().values
    
## reads by sample ##  
sample_reads = []
for i in range(table.shape[1]):
    read_dist = np.around(total_reads*table[:,i])
    sample_reads_tmp = []
    for c,r in zip(contigs,read_dist):
        sample_reads_tmp.append(subsample_reads(c,read_depth=int(r)))
    sample_reads.append(sample_reads_tmp)
        
## merge reads ##         
merged = []
for smp in sample_reads:
    merged_ = []
    for gr in smp:
        merged_.append([read[0][-20:]+read[1] 
                        for read in gr])
    merged.append(merged_)

## write files out ## 
for i,preads in enumerate(merged):
    # pool species
    preads = [j for i in preads for j in i]
    
    #fastq
    R1 = ['@sample_'+str(i)+'_read_'+str(j)+'\n'+str(pread_i)+'\n+\n'+''.join(['C'*len(pread_i)])
          for j,pread_i in enumerate(preads)]
    R1 = '\n'.join(R1)
    # write R1
    with open('fastqs/S'+str(i)+'_L001_001.fastq', 'w') as f:
        f.write(R1)
    f.close() 
    
    #fastq
    R1 = ['>sample_'+str(i)+'_read_'+str(j)+'\n'+str(pread_i)
          for j,pread_i in enumerate(preads)]
    R1 = '\n'.join(R1)
    # write R1
    with open('fastas/S'+str(i)+'_L001_001.fasta', 'w') as f:
        f.write(R1)
    f.close() 
