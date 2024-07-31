# SLURM_ARRAY_TASK_IDs are used for each batch rather than samples

# run with bash
#samplesheet="results/samplesheet/motif_pairs/samplesheet.txt"
#sample_name=$(sed -n '1p' $samplesheet)
#bash workflow/scripts/pipeline/simulate_motif_pairs.within_chroms.qarray.qsh $samplesheet $sample_name 1 


###############################################################################
# run with qsub
###############################################################################
samplesheet="results/samplesheet/motif_pairs/samplesheet.txt"
#for i in $(seq 1);
#for i in $(seq 2 53);
do
    sample_name=$(sed -n "${i}p" $samplesheet)
    echo $sample_name

    # setting the variables for exporting
    export_vars="samplesheet=$samplesheet,sample_name=$sample_name"

    # setting the logs
    output_log="results/motif_pairs/chrom_correction_only/logs/${sample_name}/run_simulations.%A.%a.out"
    error_log="results/motif_pairs/chrom_correction_only/logs/${sample_name}/run_simulations.%A.%a.err"
    mkdir -p $(dirname $output_log)

    # submit the job
    #sbatch --array=1-100 \
    sbatch --array=1-100 \
        --export=$export_vars \
        --output="${output_log}" \
        --error="${error_log}" \
        workflow/scripts/pipeline/simulate_motif_pairs.within_chroms.qarray.qsh
done
