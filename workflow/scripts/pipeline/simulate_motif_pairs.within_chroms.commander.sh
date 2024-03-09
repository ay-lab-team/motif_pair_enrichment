
# SLURM_ARRAY_TASK_IDs are used for each batch rather than samples

# run with bash
#samplesheet="results/samplesheet/motif_pairs/samplesheet.txt"
#sample_name=$(sed -n '1p' $samplesheet)
#bash workflow/scripts/pipeline/simulate_motif_pairs.within_chroms.qarray.qsh $samplesheet $sample_name 1 

# run with qsub
samplesheet="results/samplesheet/motif_pairs/samplesheet.txt"
for i in $(seq 1 53);
do
    sample_name=$(sed -n "${i}p" $samplesheet)
    echo $sample_name
    export_vars="samplesheet=$samplesheet,sample_name=$sample_name"
    sbatch --array=1-100 --export=$export_vars workflow/scripts/pipeline/simulate_motif_pairs.within_chroms.qarray.qsh
done
