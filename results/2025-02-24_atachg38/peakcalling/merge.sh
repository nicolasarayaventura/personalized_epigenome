set -x -e

scratch="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-02-24_atachg38"
work="/sc/arion/work/arayan01/project/personalized_epigenome/results/2025-02-24_atachg38"

merged_output="merged_hg38_SEM-7-8.bam"
sorted_output="merged_sorted_hg38_SEM-7-8.bam"

bam_files=(*_SC_subset_dedup.bam)

if [ ${#bam_files[@]} -gt 0 ]; then
    bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=16000]" -o "${work}/merged_job.txt" \
    "samtools merge -@ 4 ${merged_output} ${bam_files[*]} && \
     samtools sort -@ 4 -o ${sorted_output} ${merged_output} && \
     samtools index ${sorted_output}"
fi