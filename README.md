# Enrichment Analysis for Motif Pairs

## Description
The main idea behind this analysis is that 3D contacts may be happening where pairs of 
motifs meet. To analyze the significance of motif pairs appearing we extract all pairs
of motifs that happen within a given set of 3D contacts. The original set of pairs are
counted, simulations are performed to build a distribution and the original counts are
compared to their distribution to assess signficance. 

The most important script is: sample_batches_with_distance.py which run a set of
simulations.

Note: We initially developed this pipeline for the case of motif pairs joined by 3D
contacts but this pipeline could be applied to motif pairs generated through any form
of analysis.

## Installation
Install Python 3.10 and seaborn

## Usage
1) Run the different simulations with: `workflow/scripts/pipeline/simulate_motif_pairs.within_chroms.commander.sh`
2) Aggregate the simulations and calculate p-values: `workflow/scripts/pipeline/aggregate_simulations_pvalues.commander.sh`
3) Analyze the results with `workflow/scripts/summary/QQPlot_Analysis_For_All_Groups.ipynb`
4) Analyze the results with `workflow/scripts/summary/Summarize_Paired_Motif_Analysis.ipynb`

## Usage (fast mode)
1) Run the different simulations with: `bash workflow/scripts/pipeline/simulate_motif_pairs.fast.within_chroms.commander.sh`
2) Aggregate the simulations and calculate p-values: `bash workflow/scripts/pipeline/aggregate_simulations_pvalues.fast.commander.sh`
3) Analyze the results with `workflow/scripts/summary/QQPlot_Analysis.ipynb`
4) Analyze the results with `workflow/scripts/summary/Summarize_Paired_Motif_Analysis.ipynb`
Test) Run: `/mnt/BioAdHoc/Groups/vd-ay/jreyna/software/Herman-Cluster-v3/mambaforge-pypy3/envs/lc-pipelines/bin/python workflow/scripts/pipeline/simulate_motif_pairs.fast.py --sims 1 --account_dist No --inputpath "results/fimo_single_sample/Aortic-VIC.GSE154513.Homo_Sapiens.H3K27ac.b1/summarize_results/summary.txt" --savepath "results/motif_pairs/chrom_correction_only.fast/Aortic-VIC.GSE154513.Homo_Sapiens.H3K27ac.b1/batch1.results.txt"`

ATTENTION: Using the blacklist version of these scripts under: `workflow/scripts/pipeline/simulate_mps_fast_with_blacklist`

## Project status
In progress. Romeo lead this work and Joaquin is re-running it for new samples.
