import pandas as pd
import itertools as it
from collections import defaultdict,Counter
import numpy as np
import os
import csv
import argparse



parser = argparse.ArgumentParser()
parser.add_argument('--savepath', type=str, required=True)
parser.add_argument('--dirfile',type=str, required=True)
parser.add_argument('--siglevel',type=float,required=True)
args = parser.parse_args()

## Loop through all files

sig_level = argparse.siglevel
dirpath = argparse.dirfile
savepath = argparse.savepath

# Determine if file exists.  If does, delete

outputpath = savepath+"\\"+"p_values.txt"
if os.path.exists(outputpath):
        os.remove(outputpath)
else:
    pass

# Put in names of columns of dataframe
with open(outputpath, 'w', newline='') as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow(['Motif_Pair','Observed','P-value','Significant'])

for folder_name in os.listdir(dirpath):
    print(folder_name)
    outputpath = savepath+"\\"+"p_values.txt"
    with open(dirpath+'\\'+str(folder_name)+'\\'+'simulations.txt', 'r') as infile:
        reader = csv.reader(infile)
        data = list(reader)

    

    with open(outputpath, 'a', newline='') as csv_file:
        listnum = [int(num) for num in data[1][1:]]
        sims = len(listnum[1:])
        p_value = len([num for num in listnum[1:] if num <= listnum[0]])/sims
        sig = p_value >= sig_level

        # If significant, output message, "Yes".  "No" otherwise
        if sig == 1:
            message = "Yes"
        else:
            message = "No"
        writer = csv.writer(csv_file)
        writer.writerows([[data[1][0],data[1][1],p_value,message]])
