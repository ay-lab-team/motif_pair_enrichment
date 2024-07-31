
###############################################################################
# run with bash
###############################################################################
samplesheet="results/samplesheet/motif_pairs/samplesheet.txt"
#bash workflow/scripts/pipeline/aggregate_simulations_pvalues.qarray.qsh $samplesheet 1

###############################################################################
# run with SLURM
###############################################################################
samplesheet="results/samplesheet/motif_pairs/samplesheet.txt"
export_vars="samplesheet=$samplesheet"

# set output directoreis 
output_log="results/motif_pairs/chrom_correction_only/logs/agg_simulations/agg_simulations.%A.%a.out"
error_log="results/motif_pairs/chrom_correction_only/logs/agg_simulations/agg_simulations.%A.%a.out"
mkdir -p $(dirname $output_log)

# submit jobs 
# need to run 4,34 and 53
#sbatch --array=1-52 \
#sbatch --array=4,34 \
sbatch --array=53 \
        --export="$export_vars" \
        --output "$output_log" \
        --error "$error_log" \
        workflow/scripts/pipeline/aggregate_simulations_pvalues.qarray.qsh
