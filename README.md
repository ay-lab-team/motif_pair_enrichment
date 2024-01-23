# motif_pair_enrichment

## Description
The main idea behind this analysis is that 3D contacts may be happening where pairs of 
motifs meet. To analyze the significance of motif pairs appearing we extract all pairs
of motifs that happen within a given set of 3D contacts. The original set of pairs are
counted, simulations are performed to build a distribution and the original counts are
compared to their distribution to assess signficance. 

The most important script is: sample_batches_with_distance.py which run a set of
simulations.

## Installation
Install Python 3.10 and seaborn

## Usage
1) Run the different simulations with: `batch-process-all-only-chr-a_plus.qsh`
2) Aggregate the simulations and calculate p-values: `aggregate-all-sims-only-chr-a_plus.qsh`

## Project status
In progress. Romeo lead this work and Joaquin is re-running it for new samples.
