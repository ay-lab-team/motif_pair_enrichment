import itertools as it
import pandas as pd
from collections import defaultdict, Counter
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
parser.add_argument('--distance', type=int, required=False, default=2000000)
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

# drop any anchors with no motifs in either anchor (contains NaN)
filter_df =  df.loc[(~df['motif_name_1'].isna()) & (~df['motif_name_2'].isna())].reset_index(drop=True)

# add the anchor names
def create_anchor_name(sr, anchor_no):
    anchor_name = '{}-{}-{}'.format(sr[f'chr{anchor_no}'], sr[f'start{anchor_no}'], sr[f'end{anchor_no}'])
    return(anchor_name)
filter_df.loc[:, 'anchor_name1'] = filter_df.apply(create_anchor_name, args=[1], axis=1)
filter_df.loc[:, 'anchor_name2'] = filter_df.apply(create_anchor_name, args=[2], axis=1)


###############################################################################
# Counting Observed Motif Pairs
###############################################################################
print('# Counting Observed')

def sort_motifs(motifs):
    sorted_motifs = tuple(sorted(motifs))
    return sorted_motifs

# Set up counter
motif_pair_counter = Counter()

# Loop through all entries in dataframe
uniq_motif_pairs = []
for num in range(len(filter_df)):

    # Get motifs 1
    motifs1 = filter_df.motif_name_1[num].split(',')
    
    # JR-Fix: I tried to remove duplicates but removed self motifs in the process.
    # need to fix this.
    # Get motifs 2; to avoid duplicates motif pairs, motifs in anchor1 are removed
    # from motifs in anchor 2
    motifs2 = filter_df.motif_name_2[num].split(',')
    
    # Form the combination of all motif1 and motif2 pairs by
    # using the cross product
    combos = it.product(motifs1, motifs2)
    
    # Organize motifs into alphabetical order and drop duplicates
    combos = [sort_motifs(tmotifs) for tmotifs in combos]
    combos = sorted(set(combos))
    
    # Record them in counter
    for p in combos:
        motif_pair_counter[p] += 1
        uniq_motif_pairs.append(p)
        
uniq_motif_pairs = list(set(uniq_motif_pairs))

###############################################################################
# Construct Anchor Slots and Dictionary
###############################################################################
# An anchor slot is composed of a genetic region corresponding to a loop anchor
# and its constitutent motifs as a list. These are stored in dictionarys as:
# {'chr1-10000-15000': ['ACT', 'GGG', 'CGT']}. These anchor slots are used for
# simulation steps where the anchor slots are randomized.
print('# Construct Anchor Slots')

# Initialize the anchor slots data structure
anchor_slots = {}
anchor_shuffle_dict = {}
anchor_ids = []
motifs_count_sims = []

# Loop through all entries in dataframe
for i, sr in filter_df.iterrows():

    # Create anchor ids
    anchor1_id = sr['anchor_name1']
    anchor2_id = sr['anchor_name2']

    # Create blank list of anchors for shuffling based on distance
    entry_shuffle = []
    
    # Get chromosome
    chr1 = sr['chr1']
    chr2 = sr['chr2']

    # Put anchor ids and slots into list
    if anchor1_id not in anchor_slots:
        
        # Create the current slot
        motifs1 = list(sr.motif_name_1.split(','))
        anchor_slots[anchor1_id] = list(set(motifs1))
        anchor_ids.append(anchor1_id)

        # For the given slot, determine what anchors can be shuffled with it. If distance has been
        # provided then this will be applied.
        # Create column in dataframe to get distance measure for anchor 1
        # JR note: how is this calculating distance? Why are we using start1 of filter_df with start1 of sr?
        # JR note: I think I see how the distance is being calculated but you would have to check the end as well!
        filter_df['tmp_dist'] = abs(filter_df['start1'] - sr['start1'])
        if account_dist == 'Yes':

            # Only keep anchors within certain distance and within certain chromosome
            tmp_dist_df = filter_df[(filter_df['tmp_dist'] <= within_dist) & \
                                    (filter_df['chr1'] == chr1) & \
                                    (filter_df['chr2'] == chr2)].copy()

        else:
            # Only keep motifs within specific chromosome
            tmp_dist_df = filter_df[(filter_df['chr1'] == chr1) & \
                                    (filter_df['chr2'] == chr2)].copy()
        
        # Create anchor id for anchors within distance of anchor 1
        tmp_dist_df['anchor1_id_name'] = tmp_dist_df['chr1'] + '-' + \
                                            tmp_dist_df['start1'].astype('str') + '-' + \
                                            tmp_dist_df['end1'].astype('str')
        entry_shuffle += list(tmp_dist_df['anchor1_id_name'])

        # Create list of anchors wihin distance
        anchor_shuffle_dict[anchor1_id] = entry_shuffle

        # Get summary stats of amount of motifs in anchors
        motifs_count_sims.append(len(entry_shuffle))
    
    if anchor2_id not in anchor_slots:
        
        motifs2 = list(sr.motif_name_2.split(','))
        anchor_slots[anchor2_id] = list(set(motifs2))
        anchor_ids.append(anchor2_id)
        
        # Create column in dataframe to get distance measure for anchor 1
        filter_df['tmp_dist'] = abs(filter_df['start2'] - sr['start2'])
        if account_dist == 'Yes':

            # Only keep anchors within certain distance and within certain chromosome
            tmp_dist_df = filter_df[(filter_df['tmp_dist'] <= within_dist) & \
                                    (filter_df['chr1'] == chr1) & \
                                    (filter_df['chr2'] == chr2)].copy()

        else:
            # Only keep motifs within specific chromosome
            tmp_dist_df = filter_df[(filter_df['chr1'] == chr1) & \
                                    (filter_df['chr2'] == chr2)].copy()

        # Create anchor id for anchors within distance of anchor 2
        tmp_dist_df['anchor2_id_name'] = tmp_dist_df['chr2'] + '-' + \
                                            tmp_dist_df['start2'].astype('str') + '-' + \
                                            tmp_dist_df['end2'].astype('str')
        entry_shuffle += list(tmp_dist_df['anchor2_id_name'])

        # Create list of anchors within distance
        anchor_shuffle_dict[anchor2_id] = entry_shuffle

        # Get summary stats of amount of motifs in anchors
        motifs_count_sims.append(len(entry_shuffle))


###############################################################################
# Conduct Simulations
###############################################################################
print('# Conduct Simulations')

sim_results = []
sim_idx = 0
num_loops = len(filter_df)

while sim_idx < sims:

    sims_motif_pair_counter = Counter()
    for (anchor_name1, anchor_name2) in filter_df[['anchor_name1', 'anchor_name2']].values:

        # Simulate anchors based on chromsome and distance
        sim_anchors1 = np.random.choice(anchor_shuffle_dict[anchor_name1], size=1, replace=True) 
        sim_anchors2 = np.random.choice(anchor_shuffle_dict[anchor_name2], size=1, replace=True) 
        
        # Get motifs from simulated anchors
        # To remove duplicates motif pairs, motifs in anchor1 are removed from motifs in anchor 2
        sim_motifs1 = [anchor_slots[x] for x in sim_anchors1]
        sim_motifs2 = [anchor_slots[x] for x in sim_anchors2]
    
        # Form the combination of all motif1 and motif2 pairs by
        # using the cross product
        # JR Note: Why do we have use index [0]?
        shuffle_combos = it.product(sim_motifs1[0], sim_motifs2[0])

        # Organize motifs into alphabetical order
        shuffle_combos = [sort_motifs(tmotifs) for tmotifs in shuffle_combos]
        shuffle_combos = sorted(set(shuffle_combos))

        # Record them in counter
        for p in shuffle_combos:
            sims_motif_pair_counter[p] += 1

    sim_results.append(sims_motif_pair_counter)
    sim_idx += 1
    
###############################################################################
# Save json file to output
###############################################################################
print('# Save json file to output')

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

