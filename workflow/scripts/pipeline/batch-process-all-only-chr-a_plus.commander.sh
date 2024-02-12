
SLURM_TASK_ARRAY_ID=1

## bash execution
#bash workflow/scripts/pipeline/batch-process-all-only-chr-a_plus.qsh \
#    CD4_Naive_1800-RH-1.phs001703v3p1.Homo_Sapiens.H3K27ac.biorep_merged $SLURM_TASK_ARRAY_ID

# slurm execution
#sbatch --array=1 workflow/scripts/pipeline/batch-process-all-only-chr-a_plus.qsh \
#    CD4_Naive_1800-RH-1.phs001703v3p1.Homo_Sapiens.H3K27ac.biorep_merged


echo -e "View Samplesheet: vim results/samplesheet/motif_pairs/a_plus.txt"
#for sample_name in $(cat results/samplesheet/motif_pairs/a_plus.txt | cut -f 1 | sed '1d');
for sample_name in $(cat results/samplesheet/motif_pairs/a_plus.txt | cut -f 1);
do
    echo $sample_name
    sbatch --array=1-100 workflow/scripts/pipeline/batch-process-all-only-chr-a_plus.qsh $sample_name
    break
done

