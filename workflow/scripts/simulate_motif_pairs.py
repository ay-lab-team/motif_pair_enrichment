import pandas as pd
import itertools as it
from collections import defaultdict,Counter
import numpy as np
import os
import csv
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--inputfile', type=str, required=True)
parser.add_argument('--dirfile',type=str, required=True)
parser.add_argument('--sims',type=int,required=True)
args = parser.parse_args()

columns = ['loop_anchor1_chr','loop_anchor1_start','loop_anchor1_end','loop_anchor2_chr','loop_anchor2_start','loop_anchor2_end','motif_chr','motif_start','motif_end','motif_id','motif_name']
df = pd.read_csv(args.inputfile,sep="\t",header =None,names=columns)

# Get rid of all duplicate motif-pair
# (doesn't remove dupliate motifs in a single anchort)
df = df.drop_duplicates().reset_index(drop=True)

# Determine which anchor the motif falls in
df["anchor_1"] = (df['loop_anchor1_start'] <= df['motif_start']) & (df['motif_end'] <= df['loop_anchor1_end']) 
df["anchor_2"] = (df['loop_anchor2_start'] <= df['motif_start']) & (df['motif_end'] <= df['loop_anchor2_end'])
df['anchor_id'] = 1

# Get unique loop ID
df = df.drop_duplicates().reset_index(drop=True)
df['anchor1_id'] = df['loop_anchor1_chr'] + ':' + df['loop_anchor1_start'].astype(str)
df['anchor2_id'] = df['loop_anchor2_chr'] + ':' + df['loop_anchor2_start'].astype(str)
df['loop_id'] = df['loop_anchor1_chr'] + ':' + df['loop_anchor1_start'].astype(str) + ':' + df['loop_anchor2_start'].astype(str)

# Drop duplicates motifs if they popu up
df.drop_duplicates(subset=['loop_id', 'anchor1_id', 'motif_id'], inplace=True)
df.drop_duplicates(subset=['loop_id', 'anchor2_id', 'motif_id'], inplace=True)

# Count up motif pairs
motif_pair_counter = Counter()
for loop_id, loop_df in df.groupby('loop_id'):
    
    anchor1_df = loop_df.loc[loop_df.anchor_1 == True]
    anchor2_df = loop_df.loc[loop_df.anchor_2 == True]
    
    if anchor1_df.shape[0] == 0 or anchor2_df.shape[0] == 0:
        continue
        
    motifs1 = list(anchor1_df.motif_name)
    motifs2 = list(anchor2_df.motif_name)
    perms = list(it.product(motifs1, motifs2))
    
    for p in perms:
        motif_pair_counter[p] += 1


##Inputs needed, dataframe or text file of dataframe, number of sims you want, and directory path
bootstrap_df = df
sims = argparse.sims
dirpath = argparse.dirfile
##Permutate All motifs across the Chromosome
sim = 1
i = 0
results = []
all_loop_pairs = []
##Create Simulations
while sim < 3:
    # Get length of bootstrap
    n = len(bootstrap_df["motif_name"])
    #Randomly sort out permutation column with motifs
    bootstrap_df["Random_Motif"] = np.random.choice(df['motif_name'],size=n,replace=True)

    #Call in counter function
    sim_motif_pair_counter = Counter()
    
    #Groupby loop id, determine anchors based on positoin
    for loop_id, loop_df in df.groupby('loop_id'):
        anchor1_df = loop_df.loc[loop_df.anchor_1 == True]
        anchor2_df = loop_df.loc[loop_df.anchor_2 == True]

        #If missing anchors, skip loop
        if anchor1_df.shape[0] == 0 or anchor2_df.shape[0] == 0:
            continue

        #Get motifs1
        motifs1 = list(anchor1_df.Random_Motif)

        #Get motifs2
        motifs2 = list(anchor2_df.Random_Motif)

        # get motif-pair combinations
        combs = list(it.product(motifs1, motifs2))

        # count each combination
        for c in combs:
            sim_motif_pair_counter[c] += 1
    results.append(sim_motif_pair_counter)
    print(sim)
    sim += 1

t = 1
for motif_pair in motif_pair_counter.items():
    sim = 0
    while sim < sims:
        # Clean out characters for file neames
        folder_name = str(motif_pair[0]).replace("'","").replace("(","").replace(")","").replace(" ","_").replace(",","").replace(":",".")
        with open(dirpath+'\\'+str(folder_name)+'\\'+'simulations.txt', 'r') as infile:
            reader = csv.reader(infile)
            data = list(reader)
            

        # Add the new header to the first row of the list of lists
        data[0].append("sim"+str(sim))

        # Add the new data to the remaining rows of the list of lists
        data[1].append(results[sim][motif_pair[0]])
        dataentry = [data[1]]
        # Append the updated list of lists to an existing CSV file
        with open(dirpath+'\\'+str(folder_name)+'\\'+'simulations.txt', 'w', newline='') as outfile:
            writer = csv.writer(outfile)
            writer.writerow(data[0])
            writer.writerows(dataentry)
        sim+=1
    t+=1
    print(t)
    if t == 100:
        break
            