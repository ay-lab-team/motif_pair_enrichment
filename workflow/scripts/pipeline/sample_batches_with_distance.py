import itertools as it
import pandas as pd
from collections import defaultdict,Counter
import json
import numpy as np
import time
import argparse

###############################################################################
# Setup the commandline interface
###############################################################################
parser = argparse.ArgumentParser()
parser.add_argument('--inputpath', type=str, required=True)
parser.add_argument('--savepath', type=str, required=True)
parser.add_argument('--sims', type=int, required=True)
parser.add_argument('--account_dist', type=str, required=True, default="Yes")
parser.add_argument('--distance', type=int, required=False, default=20000000)
args = parser.parse_args()

inputpath = args.inputpath
savepath = args.savepath
sims = args.sims
account_dist = args.account_dist
within_dist = args.distance

# load the input data
columns = ['loop_chr1', 'loop_start1', 'loop_end1',
           'loop_chr2', 'loop_start2', 'loop_end2',
           'motif_ID_1', 'motif_name_1', 'motif_ID_2',
           'motif_name_2']
df = pd.read_csv(inputpath, sep = '\t', header=0)

# Drop any anchors with no motifs in either anchor
filter_df = df[(df['motif_ID_1'] != "None" ) & (df['motif_ID_2'] != "None")].reset_index(drop=True)

###############################################################################
# Counting Observed
###############################################################################
print("# Counting Observed")

# Set up counter
motif_pair_counter = Counter()

# Loop through all enteries in dataframe
uniq_motif_pairs = []
for num in range(len(filter_df)):
    
    # Get motifs 1
    motifs1 = list(filter_df.motif_name_1[num].split(','))
    
    # Get motifs 2
    motifs2 = list(filter_df.motif_name_2[num].split(','))
    
    # Get the cross product of all 
    combos = list(it.product(motifs1, motifs2))
    
    # Filter data to take into account (X,Y) vs (Y,X)
    combos = [(x[0], x[1]) if x[0] < x[1] else (x[1], x[0]) for x in combos]
    
    # Get chromsome attributes for Chromsome 1
    chr1_name = str(filter_df["chr1"][num])
    chr1_start = str(filter_df["start1"][num])
    chr1_end = str(filter_df["end1"][num])
    
    # Get chromosome attributes for Chromsome 2
    chr2_name = str(filter_df["chr2"][num])
    chr2_start = str(filter_df["start2"][num])
    chr2_end = str(filter_df["end2"][num])
        
    # Record them in counter
    for p in combos:
        motif_pair_counter[p] += 1
        uniq_motif_pairs.append(p)
        
uniq_motif_pairs = list(set(uniq_motif_pairs))

###############################################################################
# Construct Anchor Slots
###############################################################################
print("# Construct Anchor Slots")

# Initialize the anchor slots data structure
anchor_slots = {}
anchor_shuffle_list = {}
anchor_ids = []
motifs_count_sims = []

# Loop through all enteries in dataframe
for i, sr in filter_df.iterrows():

    # Create anchor ids
    anchor1_id = '{}-{}-{}'.format(sr["chr1"], sr["start1"], sr["end1"])
    anchor2_id = '{}-{}-{}'.format(sr["chr2"], sr["start2"], sr["end2"])

    # Create blank list of anchors for shuffling based on distance
    entry_shuffle = []
    
    # Get chromosme
    chr1 = sr["chr1"]
    chr2 = sr["chr2"]

    # Put anchor ids and slots into list
    if anchor1_id not in anchor_slots:
        
        motifs1 = list(sr.motif_name_1.split(','))
        anchor_slots[anchor1_id] = list(set(motifs1))
        anchor_ids.append(anchor1_id)

        # Create column in dataframe to get distance measure for anchor 1
        filter_df["tmp_dist"] = abs(filter_df["start1"] - sr["start1"])
        if account_dist == "Yes":

            # Only keep anchors within certain distance and within certain chromosome
            tmp_dist_df = filter_df[(filter_df["tmp_dist"]<=within_dist) & (filter_df["chr1"]==chr1) & (filter_df["chr2"]==chr2)].copy()

        else:
            # Only keep motifs within specific chromosome
            tmp_dist_df = filter_df[(filter_df["chr1"]==chr1) & (filter_df["chr2"]==chr2)].copy()
        
        # Create anchor id for anchors within distance of anchor 1
        tmp_dist_df["anchor1_id_name"] = tmp_dist_df["chr1"] + "-" + \
                                            tmp_dist_df["start1"].astype("str") + "-" + \
                                            tmp_dist_df["end1"].astype("str")
        entry_shuffle += list(tmp_dist_df["anchor1_id_name"])

        # Create list of anchors wihin distance
        anchor_shuffle_list[anchor1_id] = entry_shuffle

        # Get summary stats of amount of motifs in anchors
        motifs_count_sims.append(len(entry_shuffle))
    
    if anchor2_id not in anchor_slots:
        
        motifs2 = list(sr.motif_name_2.split(','))
        anchor_slots[anchor2_id] = list(set(motifs2))
        anchor_ids.append(anchor2_id)
        
        # Create column in dataframe to get distance measure for anchor 1
        filter_df["tmp_dist"] = abs(filter_df["start2"] - sr["start2"])
        if account_dist == "Yes":
            # Only keep anchors within certain distance and within certain chromosome
            tmp_dist_df = filter_df[(filter_df["tmp_dist"]<=within_dist) & (filter_df["chr1"]==chr1) & (filter_df["chr2"]==chr2)].copy()
        else:
            # Only keep motifs within specific chromosome
            tmp_dist_df = filter_df[(filter_df["chr1"]==chr1) & (filter_df["chr2"]==chr2)].copy()

        # Create anchor id for anchors within distance of anchor 2
        tmp_dist_df["anchor2_id_name"] = tmp_dist_df["chr2"] + "-" + \
                                            tmp_dist_df["start2"].astype("str") + "-" + \
                                            tmp_dist_df["end2"].astype("str")
        entry_shuffle+=list(tmp_dist_df["anchor2_id_name"])

        # Create list of anchors within distance
        anchor_shuffle_list[anchor2_id] = entry_shuffle

        # Get summary stats of amount of motifs in anchors
        motifs_count_sims.append(len(entry_shuffle))

###############################################################################
# Conduct Simulations
###############################################################################
print("# Conduct Simulations")

sims = args.sims
savepath = args.savepath
sim_results = []
sim_idx = 0
num_loops = len(filter_df)

while sim_idx < sims:

    sims_motif_pair_counter = Counter()
    for num, entry in filter_df.iterrows():

        # Get name of anchor
        anch_name = '{}-{}-{}'.format(entry["chr1"], entry["start1"], entry["end1"])

        ## Get length of how many anchors
        
        # Simulate anchors based on chromsome and distance
        sim_anchors1 = np.random.choice(anchor_shuffle_list[anch_name], size=1, replace=True) 
        sim_anchors2 = np.random.choice(anchor_shuffle_list[anch_name], size=1, replace=True) 
        
        # Get motifs from simulated anchors
        sim_motifs1 = [anchor_slots[x] for x in sim_anchors1]
        sim_motifs2 = [anchor_slots[x] for x in sim_anchors2]
    
        # Get the cross product of all motifs 
        shuffle_combos = list(it.product(sim_motifs1[0], sim_motifs2[0]))

        # Filter data to take into account (X,Y) vs (Y,X)
        shuffle_combos = [(x[0], x[1]) if x[0] < x[1] else (x[1], x[0]) for x in shuffle_combos]

        # Record them in counter
        for p in shuffle_combos:
            sims_motif_pair_counter[p] += 1
    sim_results.append(sims_motif_pair_counter)
    sim_idx += 1
    
###############################################################################
# Save json file to output
###############################################################################
print("# Save json file to output"

# Keys are from observed
# Aggregate simulations into one list only for observed motif pairs   
mp_sim_aggs = defaultdict(list)
for sim_idx in range(sims):
    for mp in uniq_motif_pairs:
        mp_sim_aggs[str(mp)].append(sim_results[sim_idx][mp])

# Save json
json_d = json.dumps(mp_sim_aggs)
with open(savepath, 'w') as fw:
    fw.write(json_d)

