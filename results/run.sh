set -x -e

adapter_path="/sc/arion/work/arayan01/test-env/envs/atacseq/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa"
scratch="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-02-24_atac"
work="/sc/arion/work/arayan01/project/personalized_epigenome/results/2025-02-24_atac"

function samplelist {
	sampledir="${scratch}/2025-02-19_ranjan_data/2022-11-22_data_from_Ranjan_Xiang/atac"
	output="/sc/arion/work/arayan01/project/personalized_epigenome/results/2025-02-24_atac/sample_ids.txt"
		ls ${sampledir} | cut -f 1 -d "." > sample_ids.txt
}

function fastqc_initial {
samplelist="/sc/arion/work/arayan01/project/personalized_epigenome/results/sample_ids.txt"
#once job is finished run multiqc in directory!
	sampledir="${scratch}/2025-02-19_ranjan_data/2022-11-22_data_from_Ranjan_Xiang/atac"
	qc_out="${work}/qc/init_fastqc"

	while read -r sample ; do
		bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "${scratch}/job_atacfastqc_${sample}.txt" \
    				fastqc -t 2 ${sampledir}/${sample}.fastq.gz -o ${qc_out}

done < "${samplelist}"
}

function trimming {
rm -rf "${scratch}/trimming"
mkdir -p "${scratch}/trimming"
samplelist="/sc/arion/work/arayan01/project/personalized_epigenome/results/sample_ids.txt"
sampledir="${scratch}/2025-02-19_ranjan_data/2022-11-22_data_from_Ranjan_Xiang/atac"

adapter_path="/sc/arion/work/arayan01/test-env/envs/atacseq/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa"

	sample1="7_SEM-1_S53_L004"
	output1="${scratch}/trimming/${sample1}_R1_001"
	output2="${scratch}/trimming/${sample1}_R2_001"

	bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "${scratch}/job_atactrimming_${sample1}.txt" \
		trimmomatic PE -threads 2 \
		${sampledir}/${sample1}_R1_001.fastq.gz ${sampledir}/${sample1}_R2_001.fastq.gz \
		${output1}_tr_1P.fastq.gz ${output1}_tr_1U.fastq.gz \
		${output2}_tr_2P.fastq.gz ${output2}_tr_2U.fastq.gz \
		ILLUMINACLIP:${adapter_path}:2:30:10 MINLEN:87
}

function trimming2 {
samplelist="/sc/arion/work/arayan01/project/personalized_epigenome/results/sample_ids.txt"
sampledir="${scratch}/2025-02-19_ranjan_data/2022-11-22_data_from_Ranjan_Xiang/atac"

        adapter_path="/sc/arion/work/arayan01/test-env/envs/atacseq/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa"
	sample2="8_SEM-2_S54_L004"
	output1="${scratch}/trimming/${sample2}_R1_001"
	output2="${scratch}/trimming/${sample2}_R2_001"

        bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "${scratch}/job_atactrimming_${sample2}.txt" \
		trimmomatic PE -threads 2 \
                ${sampledir}/${sample2}_R1_001.fastq.gz ${sampledir}/${sample2}_R2_001.fastq.gz \
                ${output1}_tr_1P.fastq.gz ${output1}_tr_1U.fastq.gz \
                ${output2}_tr_2P.fastq.gz ${output2}_tr_2U.fastq.gz \
                ILLUMINACLIP:${adapter_path}:2:30:10 MINLEN:87
}

function trimming {
rm -rf "${work}/qc/trimming"
mkdir -p "${work}/qc/trimming"

samplelist="/sc/arion/work/arayan01/project/personalized_epigenome/results/sample_ids.txt"

	qc_output="${work}/qc/trimming"
	trimmed_dir="${scratch}/trimming"

	while read -r sample; do
		fastqc -t 2 ${trimmed_dir}/*_SEM-*_S**_L004_R*_001_tr_*P.fastq.gz -o ${qc_output}

done < "${samplelist}"
}

function indexing {
	gen_dir="${scratch}/refgenome"

		bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "${scratch}/job_atacgenomeindexing.txt" \
			bwa index -p ${gen_dir}/hg38  ${gen_dir}/hg38.fasta
}

function mapping1 {
ref_gen="${scratch}/refgenome/hg38"
mapped_dir="${scratch}/mapping"
sample_dir="${scratch}/trimming"

sample1="7_SEM-1_S53_L004_R1_001"
sample2="7_SEM-1_S53_L004_R2_001"
basename1="7_SEM-1_S53_L004"

		bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "job_atacmapping_${basename1}.txt" \
			"bwa mem -t 2 ${ref_gen} ${sample_dir}/${sample1}_tr_1P.fastq.gz ${sample_dir}/${sample2}_tr_2P.fastq.gz \
			| samtools sort -@2 -o ${mapped_dir}/${basename1}.bam - \
			&& samtools index ${mapped_dir}/${basename1}.bam \
			&& samtools flagstat ${mapped_dir}/${basename1}.bam > ${mapped_dir}/${basename1}_map_stats.txt"

}

function mapping2 {
ref_gen="${scratch}/refgenome/hg38"
mapped_dir="${scratch}/mapping"
sample_dir="${scratch}/trimming"

sample3="8_SEM-2_S54_L004_R1_001"
sample4="8_SEM-2_S54_L004_R2_001"
basename2="8_SEM-2_S54_L004"

                bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "job_atacmapping_${basename2}.txt" \
                        "bwa mem -t 2 ${ref_gen} ${sample_dir}/${sample3}_tr_1P.fastq.gz ${sample_dir}/${sample4}_tr_2P.fastq.gz \
                        | samtools sort -@2 -o ${mapped_dir}/${basename2}.bam - \
                        && samtools index ${mapped_dir}/${basename2}.bam \
                        && samtools flagstat ${mapped_dir}/${basename2}.bam > ${mapped_dir}/${basename2}_map_stats.txt"

}

function peakcalling_processing {
samplelist="/sc/arion/work/arayan01/project/personalized_epigenome/results/sample_ids.txt"
mapped_dir="${scratch}/mapping"

		base_samples=($(awk -F'_R[12]_001' '{print $1}' "${samplelist}" | sort -u))
			base_sample1="${base_samples[0]}"
			base_sample2="${base_samples[1]}"

		base_sample3=("$base_sample1" "$base_sample2")

		for samples in "${base_sample3[@]}"; do
    rm -f \
        "${mapped_dir}/${samples}_SC_subset.bam" \
        "${mapped_dir}/${samples}_SC_subset_dedup.bam" \
        "${mapped_dir}/${samples}_SC_subset_dedup.bam.bai" \
        "${mapped_dir}/${samples}_markdup_metrics.txt" \
        "${mapped_dir}/${samples}_SC_subset_dedup_map_stats.txt" \
        "job2_prepeakcalling_${samples}.txt"

		bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "job2_prepeakcalling_${samples}.txt" \
			"samtools view -bh ${mapped_dir}/${samples}.bam > ${mapped_dir}/${samples}_SC_subset.bam && \
			picard MarkDuplicates I=${mapped_dir}/${samples}_SC_subset.bam O=${mapped_dir}/${samples}_SC_subset_dedup.bam M=${mapped_dir}/${samples}_markdup_metrics.txt && \
			samtools index ${mapped_dir}/${samples}_SC_subset_dedup.bam && \
			samtools flagstat ${mapped_dir}/${samples}_SC_subset_dedup.bam > ${mapped_dir}/${samples}_SC_subset_dedup_map_stats.txt"
done
}

function peakcalling {

samplelist="/sc/arion/work/arayan01/project/personalized_epigenome/results/sample_ids.txt"
mapped_dir="${scratch}/mapping"
	rm -rf "${work}/peakcalling"
	mkdir -p "${work}/peakcalling"
	peakcall="${work}/peakcalling"

		base_samples=($(awk -F'_R[12]_001' '{print $1}' "${samplelist}" | sort -u))
                        base_sample1="${base_samples[0]}"
                        base_sample2="${base_samples[1]}"

                base_sample3=("$base_sample1" "$base_sample2")

                for samples in "${base_sample3[@]}"; do

		bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "job3_${samples}.txt" \
			"macs2 callpeak -t ${mapped_dir}/${samples}_SC_subset_dedup.bam -f BAMPE -n ${peakcall}/${samples}_peak -g hs --keep-dup all"
done
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
peakcalling
