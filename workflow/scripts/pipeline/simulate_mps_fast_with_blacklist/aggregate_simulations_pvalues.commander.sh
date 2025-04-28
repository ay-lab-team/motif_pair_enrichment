
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
output_log="results/motif_pairs/chrom_correction_only.fast.with_blacklist/logs/agg_simulations/agg_simulations.%A.%a.out"
error_log="results/motif_pairs/chrom_correction_only.fast.with_blacklist/logs/agg_simulations/agg_simulations.%A.%a.out"
mkdir -p $(dirname $output_log)

# submit jobs 
#sbatch --array=1 \
sbatch --array=1-65 \
        --export="$export_vars" \
        --output "$output_log" \
        --error "$error_log" \
        workflow/scripts/pipeline/simulate_mps_fast_with_blacklist/aggregate_simulations_pvalues.qarray.qsh
