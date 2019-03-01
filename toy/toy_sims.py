import copy
import numpy as np
import pandas as pd
from numpy.random import randint

def generate_strain_reads(reads_in, n_samples = 1, 
                          n_read_chnaged = 100, 
                          n_nuc_chnaged = .2, 
                          nuc_allowed = ['A','G','C','T'],
                         paired=True):
    
    generate_chnages = {i:[j for j in nuc_allowed if i!=j] 
                        for i in nuc_allowed}
    random_change = np.min([len(l) \
                            for l in \
                            generate_chnages.values()])
    samples = [reads_in.copy()]
    time = [[0] * len(samples[0])]
    for k in range(1,n_samples + 1):
        samples.append(copy.deepcopy(samples[k-1]))
        time.append(time[k-1].copy())
        r_start = randint(len(samples[k])-n_read_chnaged)
        for i in range(r_start,r_start+n_read_chnaged):
            time[k][i] = k
            # get read
            if paired==True:
                read_i_R1 = list(samples[k-1][i][0]).copy()
                read_i_R2 = list(samples[k-1][i][1]).copy()
            else:
                read_i_R1 = list(samples[k-1][i]).copy()
            # randomly chnage nuc.
            for j in randint(len(read_i_R1),size = int(len(read_i_R1) * n_nuc_chnaged)):
                if paired==True:
                    #R1
                    read_i_R1[j] = generate_chnages[read_i_R1[j]]\
                                    [np.random.randint(random_change)] 
                    #R2
                    read_i_R2[j] = generate_chnages[read_i_R2[j]]\
                                    [np.random.randint(random_change)] 
                else:
                    #R1
                    read_i_R1[j] = generate_chnages[read_i_R1[j]]\
                                [np.random.randint(random_change)]
            if paired==True:
                samples[k][i][0] = ''.join(read_i_R1.copy())
                samples[k][i][1] = ''.join(read_i_R2.copy())  
            else:                                       
                samples[k][i] = ''.join(read_i_R1)

    return samples,time

def subsample_reads(contig, read_length = 250,
                    read_depth = 3000, 
                    distance = 20, k = 55):
    
    # generate k-mer cov of a string ar read length
    reads = []
    i = 0 # overalp
    while i<=len(contig)-read_length:
        reads.append([contig[i:i+read_length],
                      contig[i+(read_length-distance):\
                             (i+(read_length-distance))+read_length]])
        i+=k
    # sub-sample randomly to read depth
    r_start = randint(len(reads)-read_depth)
    sub_reads = [[reads[i][0],reads[i][1]] 
                 for i in range(r_start,r_start+read_depth)]
    return sub_reads

