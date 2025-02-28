set -x -e

adapter_path="/sc/arion/work/arayan01/test-env/envs/atacseq/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa"
atac_scratch="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-02-19_ranjan_data/data"
atac_work="/sc/arion/work/arayan01/project/personalized_epigenome/results/2025-02-24_atac"

function atac_samplelist {
	sampledir="${atac_scratch}/2022-11-22_data_from_Ranjan_Xiang/atac"
	output="/sc/arion/work/arayan01/project/personalized_epigenome/results/2025-02-24_atac/sample_ids.txt"
		ls ${sampledir} | cut -f 1 -d "." > sample_ids.txt
}

function atac_fastqc_initial {
samplelist="/sc/arion/work/arayan01/project/personalized_epigenome/results/sample_ids.txt"
#once job is finished run multiqc in directory!
	sampledir="${atac_scratch}/2022-11-22_data_from_Ranjan_Xiang/atac"
	qc_out="${atac_work}/qc/init_fastqc"

	while read -r sample ; do
		bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "${atac_scratch}/job_atacfastqc_${sample}.txt" \
    				fastqc -t 2 ${sampledir}/${sample}.fastq.gz -o ${qc_out}

done < "${samplelist}"
}

function atac_trimming {
rm -rf "${atac_scratch}/2025-02-25_atac/trimming"
mkdir -p "${atac_scratch}/2025-02-25_atac/trimming"
samplelist="/sc/arion/work/arayan01/project/personalized_epigenome/results/sample_ids.txt"
sampledir="${atac_scratch}/2022-11-22_data_from_Ranjan_Xiang/atac"

adapter_path="/sc/arion/work/arayan01/test-env/envs/atacseq/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa"

	sample1="7_SEM-1_S53_L004"
	output1="${atac_scratch}/2025-02-25_atac/trimming/${sample1}_R1_001"
	output2="${atac_scratch}/2025-02-25_atac/trimming/${sample1}_R2_001"

	bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "${atac_scratch}/job_atactrimming_${sample1}.txt" \
		trimmomatic PE -threads 2 \
		${sampledir}/${sample1}_R1_001.fastq.gz ${sampledir}/${sample1}_R2_001.fastq.gz \
		${output1}_tr_1P.fastq.gz ${output1}_tr_1U.fastq.gz \
		${output2}_tr_2P.fastq.gz ${output2}_tr_2U.fastq.gz \
		ILLUMINACLIP:${adapter_path}:2:30:10 MINLEN:87
}

function atac_trimming2 {
samplelist="/sc/arion/work/arayan01/project/personalized_epigenome/results/sample_ids.txt"
sampledir="${atac_scratch}/2022-11-22_data_from_Ranjan_Xiang/atac"

        adapter_path="/sc/arion/work/arayan01/test-env/envs/atacseq/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa"
	sample2="8_SEM-2_S54_L004"
	output1="${atac_scratch}/2025-02-25_atac/trimming/${sample2}_R1_001"
	output2="${atac_scratch}/2025-02-25_atac/trimming/${sample2}_R2_001"

        bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "${atac_scratch}/job_atactrimming_${sample2}.txt" \
		trimmomatic PE -threads 2 \
                ${sampledir}/${sample2}_R1_001.fastq.gz ${sampledir}/${sample2}_R2_001.fastq.gz \
                ${output1}_tr_1P.fastq.gz ${output1}_tr_1U.fastq.gz \
                ${output2}_tr_2P.fastq.gz ${output2}_tr_2U.fastq.gz \
                ILLUMINACLIP:${adapter_path}:2:30:10 MINLEN:87
}

function atac_fastqc_trimming {
rm -rf "${atac_work}/qc/trimming"
mkdir -p "${atac_work}/qc/trimming"

samplelist="/sc/arion/work/arayan01/project/personalized_epigenome/results/sample_ids.txt"

	qc_output="${atac_work}/qc/trimming"
	trimmed_dir="${atac_scratch}/2025-02-25_atac/trimming"

	while read -r sample; do
		fastqc -t 2 ${trimmed_dir}/*_SEM-*_S**_L004_R*_001_tr_*P.fastq.gz -o ${qc_output}

done < "${samplelist}"
}

function atac_refgen_indexing {
	gen_dir="${atac_scratch}/2025-02-25_atac/refgenome"

		bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "${atac_scratch}/job_atacgenomeindexing.txt" \
			bwa index -p ${gen_dir}/hg38  ${gen_dir}/hg38.fasta
}

function atac_mapping1 {
ref_gen="${atac_scratch}/2025-02-25_atac/refgenome/hg38"
mapped_dir="${atac_scratch}/2025-02-25_atac/mapping"
sample_dir="${atac_scratch}/2025-02-25_atac/trimming"

sample1="7_SEM-1_S53_L004_R1_001"
sample2="7_SEM-1_S53_L004_R2_001"
basename1="7_SEM-1_S53_L004"

		bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "job_atacmapping_${basename1}.txt" \
			"bwa mem -t 2 ${ref_gen} ${sample_dir}/${sample1}_tr_1P.fastq.gz ${sample_dir}/${sample2}_tr_2P.fastq.gz \
			| samtools sort -@2 -o ${mapped_dir}/${basename1}.bam - \
			&& samtools index ${mapped_dir}/${basename1}.bam \
			&& samtools flagstat ${mapped_dir}/${basename1}.bam > ${mapped_dir}/${basename1}_map_stats.txt"

}

function atac_mapping2 {
ref_gen="${atac_scratch}/2025-02-25_atac/refgenome/hg38"
mapped_dir="${atac_scratch}/2025-02-25_atac/mapping"
sample_dir="${atac_scratch}/2025-02-25_atac/trimming"

sample3="8_SEM-2_S54_L004_R1_001"
sample4="8_SEM-2_S54_L004_R2_001"
basename2="8_SEM-2_S54_L004"

                bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "job_atacmapping_${basename2}.txt" \
                        "bwa mem -t 2 ${ref_gen} ${sample_dir}/${sample3}_tr_1P.fastq.gz ${sample_dir}/${sample4}_tr_2P.fastq.gz \
                        | samtools sort -@2 -o ${mapped_dir}/${basename2}.bam - \
                        && samtools index ${mapped_dir}/${basename2}.bam \
                        && samtools flagstat ${mapped_dir}/${basename2}.bam > ${mapped_dir}/${basename2}_map_stats.txt"

}

function atac_pre_peakcalling_processing {
samplelist="/sc/arion/work/arayan01/project/personalized_epigenome/results/sample_ids.txt"
mapped_dir="${atac_scratch}/2025-02-25_atac/mapping"

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

function atac_peakcalling {

samplelist="/sc/arion/work/arayan01/project/personalized_epigenome/results/sample_ids.txt"
mapped_dir="${atac_scratch}/2025-02-25_atac/mapping"
	rm -rf "${atac_work}/peakcalling"
	mkdir -p "${atac_work}/peakcalling"
	peakcall="${atac_work}/peakcalling"

		base_samples=($(awk -F'_R[12]_001' '{print $1}' "${samplelist}" | sort -u))
                        base_sample1="${base_samples[0]}"
                        base_sample2="${base_samples[1]}"

                base_sample3=("$base_sample1" "$base_sample2")

                for samples in "${base_sample3[@]}"; do

		bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "job3_${samples}.txt" \
			"macs2 callpeak -t ${mapped_dir}/${samples}_SC_subset_dedup.bam -f BAMPE -n ${peakcall}/${samples}_peak -g hs --keep-dup all"
done
}


#atac_samplelist
#atac_fastqc_initial
#atac_trimming
#atac_trimming2
#atac_fastqc_trimming
#atac_refgen_indexing
#atac_mapping1
#atac_mapping2
#atac_pre_peakcalling_processing
atac_peakcalling
