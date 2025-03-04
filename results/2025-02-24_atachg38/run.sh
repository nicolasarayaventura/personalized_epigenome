set -x -e

adapter_path="/sc/arion/work/arayan01/test-env/envs/atacseq/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa"
scratch="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-02-24_atac"
work="/sc/arion/work/arayan01/project/personalized_epigenome/results/2025-02-24_atachg38"

function samplelist {
	rm -rf "${work}/sample_ids.txt"
	sampledir="${scratch}/2025-02-19_ranjan_data/2022-11-22_data_from_Ranjan_Xiang/atac"
	output="${work}/sample_ids.txt"
	
	samplepath=$(realpath "${sampledir}")
    find "${sampledir}" -name "*_R1_001.fastq.gz" | while read R1; do
        sample=$(basename "${R1}" | sed 's/_R1_001.fastq.gz//')
        R2="${R1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
        echo -e "${sample}\t${R1}.fasta.gz\t${R2}.fasta.gz" >> "${output}"
    done

}

function fastqc_initial {
	samplelist="${work}/sample_ids.txt"
	sampledir="${scratch}/2025-02-19_ranjan_data/2022-11-22_data_from_Ranjan_Xiang/atac"
	qc_out="${work}/qc/init_fastqc" #once job is finished run multiqc in directory!

	while IFS=$'\t' read -r sample fastq_file; do
		bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "${scratch}/job_atacfastqc_${sample}.txt" \
		fastqc -t 2 "${sampledir}/${fastq_file}" -o "${qc_out}"
	done < "${samplelist}"
}

function trimming {
	rm -rf "${scratch}/trimmed"
	mkdir -p "${scratch}/trimmed"
	
	samplelist="${work}/sample_ids.txt"
	adapter_path="/sc/arion/work/arayan01/test-env/envs/atacseq/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa"

	while IFS=$'\t' read -r sample R1 R2; do
        bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "${scratch}/job_atactrimming_${sample}.txt" \
        "trimmomatic PE -threads 2 \
        ${R1} ${R2} \
        ${scratch}/trimmed/${sample}_tr_1P.fastq.gz \
        ${scratch}/trimmed/${sample}_tr_1U.fastq.gz \
        ${scratch}/trimmed/${sample}_tr_2P.fastq.gz \
        ${scratch}/trimmed/${sample}_tr_2U.fastq.gz \
        ILLUMINACLIP:${adapter_path}:2:30:10 MINLEN:87" #2:30:10 tell Trimmomatic to cut the adapters allowing a maximum of 2 mismatches, a score of 30 for paired reads, and score of 10 for single reads 
														#MINLEN option tells Trimmomatic to drop reads shorter than 30 base pairs.
	done < "${samplelist}"
}

function fastqc_trimming {
	rm -rf "${work}/qc/trimming"
	mkdir -p "${work}/qc/trimming"

	samplelist="${work}/sample_ids.txt"
	qc_output="${work}/qc/trimming"
	trimmed_dir="${scratch}/trimming"

	while IFS=$'\t' read -r sample R1 R2; do
		fastqc -t 2 "${trimmed_dir}/${sample}_tr_1P.fastq.gz" "${trimmed_dir}/${sample}_tr_2P.fastq.gz" -o "${qc_output}"
	done < "${samplelist}"

}

function indexing {
	gen_dir="${scratch}/refgenome"
	bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "${scratch}/job_atacgenomeindexing.txt" \
	"bwa index -p ${gen_dir}/hg38  ${gen_dir}/hg38.fasta"

}

function mapping {
	ref_gen="${scratch}/refgenome/hg38"
	mapped_dir="${scratch}/mapping"
	sample_dir="${scratch}/trimming"

	while IFS=$'\t' read -r sample R1 R2; do
        basename="${sample}"
        sample1="${sample}_tr_1P.fastq.gz"
        sample2="${sample}_tr_2P.fastq.gz"

        # Submit job to bsub
        bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "${scratch}/job_atacmapping_${basename}.txt" \
        "bwa mem -t 2 ${ref_gen} ${sample_dir}/${sample1} ${sample_dir}/${sample2} \
        | samtools sort -@2 -o ${mapped_dir}/${basename}.bam - \
        && samtools index ${mapped_dir}/${basename}.bam \
        && samtools flagstat ${mapped_dir}/${basename}.bam > ${mapped_dir}/${basename}_map_stats.txt"
    done < "${samplelist}"
}

function pre_peakcalling_processing {
    samplelist="${work}/sample_ids.txt"
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
        bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "job2_prepeakcalling_${basename}.txt" \
            "samtools view -bh ${mapped_dir}/${basename}.bam > ${mapped_dir}/${basename}_SC_subset.bam && \
            picard MarkDuplicates I=${mapped_dir}/${basename}_SC_subset.bam O=${mapped_dir}/${basename}_SC_subset_dedup.bam M=${mapped_dir}/${basename}_markdup_metrics.txt && \
            samtools index ${mapped_dir}/${basename}_SC_subset_dedup.bam && \
            samtools flagstat ${mapped_dir}/${basename}_SC_subset_dedup.bam > ${mapped_dir}/${basename}_SC_subset_dedup_map_stats.txt"
    done < "${samplelist}"
}

function peakcalling {
    samplelist="${work}/sample_ids.txt"
    mapped_dir="${scratch}/mapping"
    rm -rf "${work}/peakcalling"
    mkdir -p "${work}/peakcalling"
    peakcall="${work}/peakcalling"

    while IFS=$'\t' read -r sample R1 R2; do
        basename="${sample}"
        bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "job3_${basename}.txt" \
            "macs2 callpeak -t ${mapped_dir}/${basename}_SC_subset_dedup.bam -f BAMPE -n ${peakcall}/${basename}_peak -g hs --keep-dup all"
    done < "${samplelist}"
}

#samplelist
#fastqc_initial
#trimming
#trimming2
#fastqc_trimming
#refgen_indexing
#mapping1
#mapping2
#pre_peakcalling_processing
#peakcalling
