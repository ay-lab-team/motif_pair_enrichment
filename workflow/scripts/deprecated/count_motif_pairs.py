import pandas as pd
import itertools as it
from collections import defaultdict,Counter
import numpy as np
import os
import csv
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--inputpath', type=str, required=True)
parser.add_argument('--outputpath', type=str, required=True)
args = parser.parse_args()

columns = ['loop_anchor1_chr','loop_anchor1_start','loop_anchor1_end','loop_anchor2_chr','loop_anchor2_start','loop_anchor2_end','motif_chr','motif_start','motif_end','motif_id','motif_name']
df = pd.read_csv(args.inputpath,sep="\t",header =None,names=columns)

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

# Drop duplicates if they occur (none should be dropped, just as a precaution))
df.drop_duplicates(subset=['loop_id', 'anchor1_id', 'motif_id'], inplace=True)
df.drop_duplicates(subset=['loop_id', 'anchor2_id', 'motif_id'], inplace=True)

motif_pair_counter = Counter()

# Test file only looks at first 100 to save on time
# PLEASE DELETE "t" WHEN DOING ACTUAL SIMULATION

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
t = 1
for motif_pair in motif_pair_counter.items():
    # Clean out characters for file neames
    folder_name = str(motif_pair[0]).replace("'","").replace("(","").replace(")","").replace(" ","_").replace(",","").replace(":",".")
    # Create File of paired observations
    if os.path.exists(args.outputpath+"\\"+str(folder_name)):
        pass
    else:
        os.makedirs(args.outputpath+"\\"+str(folder_name))
    # Save observed file in text file
    with open(args.outputpath+'\\'+str(folder_name)+'\\'+'simulations.txt', 'w', newline='') as csv_file:
        data = [[motif_pair[1]]]
        writer = csv.writer(csv_file)
        writer.writerow(['Observed'])
        writer.writerows(data)
    t+=1
    if t == 100:
        break
    