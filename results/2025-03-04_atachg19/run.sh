set -x -e

adapter_path="/sc/arion/work/arayan01/test-env/envs/atacseq/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa"
scratch="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-03-04_atachg19"
work="/sc/arion/work/arayan01/project/personalized_epigenome/results/2025-03-04_atachg19"
gen_dir="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/hg19_rfgen"
samplelist="/sc/arion/work/arayan01/project/personalized_epigenome/results/2025-02-24_atachg38/sample_ids.txt"

function filteredref {
    samtools faidx "${gen_dir}/hg19.fa"
    grep -vE "Un|random|alt" "${gen_dir}/hg19.fa.fai" > "${gen_dir}/hg19_filtered.fa.fai"
    cut -f1 "${gen_dir}/hg19_filtered.fa.fai" > "${gen_dir}/filtered_chromosomes.txt"

    rm -f "${gen_dir}/hg19_filtered.fa"
    while read chr; do
        samtools faidx "${gen_dir}/hg19.fa" "$chr" >> "${gen_dir}/hg19_filtered.fa"
    done < "${gen_dir}/filtered_chromosomes.txt"

    bsub -P acc_oscarlr -q premium -W 24:00 -R "rusage[mem=16000]" \
    -o "${work}/job_hg19filterd_indexing.txt" -eo "${work}/job_hg19filterd_indexing_err.txt" \
    "bwa index -p ${gen_dir}/hg19_filtered ${gen_dir}/hg19_filtered.fa"
}

function mapping {
    rm -rf ${scratch}/mapping
    mkdir -p ${scratch}/mapping
	ref_gen="${gen_dir}/hg19_filtered"
	mapped_dir="${scratch}/mapping"
	sample_dir="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-02-24_atachg38/trimming"

	while IFS=$'\t' read -r sample R1 R2; do
        basename="${sample}"
        sample1="${sample}_tr_1P.fastq.gz"
        sample2="${sample}_tr_2P.fastq.gz"

        # Submit job to bsub
        bsub -P acc_oscarlr -q premium -n 4 -W 24:00 -R "rusage[mem=8000]" -o "${work}/job_atacmapping_${basename}.txt" \
        "bwa mem -t 4 ${ref_gen} ${sample_dir}/${sample1} ${sample_dir}/${sample2} \
        | samtools sort -@4 -o ${mapped_dir}/${basename}.bam - \
        && samtools index ${mapped_dir}/${basename}.bam \
        && samtools flagstat ${mapped_dir}/${basename}.bam > ${mapped_dir}/${basename}_map_stats.txt"
    done < "${samplelist}"
}
function pre_peakcalling_processing {
    mapped_dir="${scratch}/mapping"

    while IFS=$'\t' read -r sample R1 R2; do
        basename="${sample}"
        rm -f \
            "${mapped_dir}/${basename}_SC_subset.bam" \
            "${mapped_dir}/${basename}_SC_subset_dedup.bam" \
            "${mapped_dir}/${basename}_SC_subset_dedup.bam.bai" \
            "${mapped_dir}/${basename}_markdup_metrics.txt" \
            "${mapped_dir}/${basename}_SC_subset_dedup_map_stats.txt" \
            "job2_prepeakcalling_${basename}.txt"
        bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "${work}/prepeakcalling_${basename}.txt" \
            "samtools view -bh ${mapped_dir}/${basename}.bam > ${mapped_dir}/${basename}_SC_subset.bam && \
            picard MarkDuplicates I=${mapped_dir}/${basename}_SC_subset.bam O=${mapped_dir}/${basename}_SC_subset_dedup.bam M=${mapped_dir}/${basename}_markdup_metrics.txt && \
            samtools index ${mapped_dir}/${basename}_SC_subset_dedup.bam && \
            samtools flagstat ${mapped_dir}/${basename}_SC_subset_dedup.bam > ${mapped_dir}/${basename}_SC_subset_dedup_map_stats.txt"
    done < "${samplelist}"
}

function peakcalling {
    mapped_dir="${scratch}/mapping"
    rm -rf "${work}/peakcalling"
    mkdir -p "${work}/peakcalling"
    peakcall="${work}/peakcalling"

    while IFS=$'\t' read -r sample R1 R2; do
        basename="${sample}"
        bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "${work}/peakcalling_${basename}.txt" \
            "macs2 callpeak -t ${mapped_dir}/${basename}_SC_subset_dedup.bam -f BAMPE -n ${peakcall}/${basename}_peak -g hs --keep-dup all"
    done < "${samplelist}"
}

#filteredref
mapping
#pre_peakcalling_processing
#peakcalling
