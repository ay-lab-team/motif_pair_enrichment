###############################################################################
# Details 
###############################################################################
# SLURM_ARRAY_TASK_IDs are used for each batch rather than samples


###############################################################################
# Running with bash
###############################################################################

# run with bash
#samplesheet="results/samplesheet/motif_pairs/samplesheet.txt"
#sample_name=$(sed -n '1p' $samplesheet)
#bash workflow/scripts/pipeline/simulate_motif_pairs.within_chroms.qarray.qsh $samplesheet $sample_name 1 


###############################################################################
# Running with SLURM
###############################################################################


# run with qsub
distance=2000000
samplesheet="results/samplesheet/motif_pairs/samplesheet.txt"
for i in $(seq 3 53);
do
    sample_name=$(sed -n "${i}p" $samplesheet)
    echo $sample_name

    # setting the variables for exporting
    export_vars="samplesheet=$samplesheet,sample_name=$sample_name,distance=$distance"

    # setting the logs
    output_log="results/motif_pairs/${distance}_distance_correction/logs/${sample_name}/run_simulations.%A.%a.out"
    error_log="results/motif_pairs/${distance}_distance_correction/logs/${sample_name}/run_simulations.%A.%a.err"
    mkdir -p $(dirname $output_log)

    # submit the job
    sbatch --array=1-100 \
        --export=$export_vars \
        --output="${output_log}" \
        --error="${error_log}" \
        workflow/scripts/pipeline/simulate_motif_pairs.within_2mb_dist.qarray.qsh
done
