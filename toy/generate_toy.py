import numpy as np
import pandas as pd
from toy_sims import subsample_reads,generate_strain_reads

## import contigs ## 
gpath = 'genomes/GCF_000143685.1_ASM14368v1_genomic.fna'
contig = ''.join([line.rstrip('\n') for line in \
                  open(gpath)][1:])
## generate toy reads ## 
sample_1 = subsample_reads(contig)
new_reads,new_time = generate_strain_reads(sample_1)
## write files out ## 
for i,preads in enumerate(new_reads):
    R1 = ['@sample_'+str(i)+'_read_'+str(j)+'\n'+str(pread_i[0])+'\n+\n'+''.join(['C'*len(pread_i[0])])
          for j,pread_i in enumerate(preads)]
    R2 = ['@sample_'+str(i)+'_read_'+str(j)+'\n'+str(pread_i[1])+'\n+\n'+''.join(['C'*len(pread_i[1])])
          for j,pread_i in enumerate(preads)]
    R1 = '\n'.join(R1)
    R2 = '\n'.join(R2)
    # write R1
    with open('fastqs/S'+str(i)+'_L001_R1_001.fastq', 'w') as f:
        f.write(R1)
    f.close() 
    # write R2
    with open('fastqs/S'+str(i)+'_L001_R2_001.fastq', 'w') as f:
        f.write(R2)
    f.close()    
# write time
new_timedf = pd.DataFrame(new_time).T
new_timedf.index.name = 'read_number'
new_timedf.to_csv('time_read.tsv',sep='\t')
