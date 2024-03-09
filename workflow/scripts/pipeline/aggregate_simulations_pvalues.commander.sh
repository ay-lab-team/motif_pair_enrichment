
# run with bash
samplesheet="results/samplesheet/motif_pairs/samplesheet.txt"
#bash workflow/scripts/pipeline/aggregate_simulations_pvalues.qarray.qsh $samplesheet 1

# run with sbatch
samplesheet="results/samplesheet/motif_pairs/samplesheet.txt"
export_vars="samplesheet=$samplesheet"
sbatch --array=1-53 --export="$export_vars" workflow/scripts/pipeline/aggregate_simulations_pvalues.qarray.qsh
