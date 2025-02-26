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
		${output2}_tr_2P.fastq.gz ${output2}_tr_2P.fastq.gz \
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
                ${output2}_tr_2P.fastq.gz ${output2}_tr_2P.fastq.gz \
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

function atac_ref_indexing {

	ref_dir="${atac_scratch}/refgenome"
	ref_gen="${atac_scratch}/refgenome/hg38.fa.fai"


	bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "${atac_scratch}/job_refgenome_index.txt" \
		bwa index -p hg38 ${ref_dir}/${ref_gen}

}

function atac_mapping {

samplelist="/sc/arion/work/arayan01/project/personalized_epigenome/results/sample_ids.txt"
ref_gen="${atac_scratch}/refgenome/indexed/hg"
sample_dir="${atac_scratch}/2025-02-25_atac/trimming"
mapped_dir="${atac_scratch}/2025-02-25_atac/mapping"
	while read sample ; do
		bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "job1_${sample}.txt" \
			"bwa mem -t 2 ${ref_gen} ${sample_dir}/${sample}_tr_1P.fastq.gz ${sample_dir}/${sample}_tr_2P.fastq.gz \
			| samtools sort -@2 -o ${mapped_dir}/${sample}.bam - \
			&& samtools index ${mapped_dir}/${sample}.bam \
			&& samtools flagstat ${mapped_dir}/${sample}.bam > ${mapped_dir}/${sample}_map_stats.txt"


done < "${samplelist}"
}
#atac_samplelist
#atac_fastqc_initial
#atac_trimming
#atac_trimming2
atac_fastqc_trimming
#atac_ref_indexing
