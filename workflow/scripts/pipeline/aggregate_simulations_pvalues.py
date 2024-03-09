import os
import json
import time
import itertools as it
from collections import defaultdict, Counter
import argparse

import numpy as np
import pandas as pd


###############################################################################
# Setup the commandline interface
###############################################################################
parser = argparse.ArgumentParser()
parser.add_argument('--inputfilepath', type=str, required=True)
parser.add_argument('--simpath',type=str,required=True)
parser.add_argument('--simnums', type=int, required=True)
parser.add_argument('--skip-missing-sims', dest='skip_sim', action='store_true', required=False)
args = parser.parse_args()

# inputfilepath should be the location to a fimo based output
#inputfilepath: "/mnt/bioadhoc-temp/Groups/vd-ay/kfetter/hichip-db-loop-calling/results/motif_analysis/meme/fimo/GM12878.GSE101498.Homo_Sapiens.H3K27ac.b1/new_summarize_results/new_summary.txt"
columns = ['loop_chr1', 'loop_start1', 'loop_end1', 'loop_chr2',
           'loop_start2', 'loop_end2', 'motif_ID_1', 'motif_name_1',
           'motif_ID_2', 'motif_name_2']
df = pd.read_csv(args.inputfilepath, sep='\t', header=0)

# Drop about half of the loops
# Drop any anchors with no motifs in either anchor
filter_df = df[(df['motif_ID_1'] != 'None' ) & (df['motif_ID_2'] != 'None')].reset_index(drop=True)

###############################################################################
# Count up unique motif pairs
###############################################################################
print("# Count up unique motif pairs")

motif_pair_counter = Counter()
uniq_motif_pairs = []
for num in range(len(filter_df)):

    # Are we suppoed to use set()???
    # Get motifs 1
    motifs1 = list(set(filter_df.motif_name_1[num].split(',')))
    
    # Get motifs 2
    motifs2 = list(set(filter_df.motif_name_2[num].split(',')))
    
    # Get the cross product of all 
    combos = list(it.product(motifs1, motifs2))
    
    # Filter data to take into account (X,Y) vs (Y,X)
    combos = [(x[0], x[1]) if x[0] < x[1] else (x[1], x[0]) for x in combos]
    
    # Get chromsome attributes for Chromsome 1
    chr1_name = str(filter_df['chr1'][num])
    chr1_start = str(filter_df['start1'][num])
    chr1_end = str(filter_df['end1'][num])
    
    # Get chromosome attributes for Chromsome 2
    chr2_name = str(filter_df['chr2'][num])
    chr2_start = str(filter_df['start2'][num])
    chr2_end = str(filter_df['end2'][num])
        
    # Record them in counter
    for p in combos:
        motif_pair_counter[p] += 1
        uniq_motif_pairs.append(p)
uniq_motif_pairs = set(uniq_motif_pairs)
        
###############################################################################
# Aggregate simulations
###############################################################################
print("# Aggregate simulations")

results = {}
simcount = args.simnums + 1
for i in range(1, simcount):

    # get file path of the current batch
    fn = str(args.simpath) + '/batch' + str(i) + '.results.txt'

    # if the current file doesn't exists and skip has been activated
    # then this batch will be skipped
    if not os.path.exists(fn) and args.skip_sim:
        print(f'Skipped: {fn}')
        continue

    # load the current batch
    curr_sim = json.loads(open(fn).read())

    for motif_pair in uniq_motif_pairs:

        obs = motif_pair_counter[motif_pair]

        # get the counts for sims if present or else return a negative number
        # which will ensure thie sim will not be higher than observer
        sims = curr_sim.get(str(motif_pair), [-1000])

        # Go through list to see if it is above the observered count
        listnum = [True if x >= obs else False for x in sims]

        # Count total number of sims above observed
        above_obs = sum(listnum)

        # If first batch, create result
        if i == 1:
            results[motif_pair] = above_obs

        # If any other batch, add to previous batch
        else:
            results[motif_pair] += above_obs
    print(i)

###############################################################################
# Calculate p-value
###############################################################################
print("# Calculate p-value")

p_values = []
for motif_pair, count in results.items():

    # Get p-value by dividing by total number of simulations
    totalsims = 1000 * args.simnums
    p_value = count / totalsims

    # Get motif-pair, p-value, and observed count of p-value
    entry = {'Motif1': motif_pair[0], 'Motif2': motif_pair[1], 'Sim_Count': count,
             'P_value': p_value, 'Obs_Count': motif_pair_counter[motif_pair]}
    p_values.append(entry)

# Create dataframe for all counts
mp_data = pd.DataFrame(p_values)

# Save to dataframe
outfn = str(args.simpath) + '/P_values_agg.tsv'
mp_data.to_csv(outfn, sep='\t', index=False)
