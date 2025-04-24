set -x -e

scratch="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-04-07_chipsseqhg38"
work="/sc/arion/work/arayan01/project/personalized_epigenome/results/2025-04-07_chipsseqhg38"

sampledir="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-02-19_ranjan_data/data/2022-11-22_data_from_Ranjan_Xiang/chip/RSen042022_24"
samplelist="${work}/sample_ids.txt"

function samplelist {
    rm -rf "${work}/sample_ids.txt"
    mkdir -p "${work}"
    
    output="${work}/sample_ids.txt"

    
    for filepath in "${sampledir}"/*; do
        [ -f "$filepath" ] || continue
        basename=$(basename "$filepath")
        echo -e "${basename}\t${filepath}" >> "${output}"
    done
}

function fastqc_initial {
    rm -rf "${work}/qc/init_fastqc"
    mkdir "${work}/qc/init_fastqc"
    
	qc_out="${work}/qc/init_fastqc" #once job is finished run multiqc in directory!

	while IFS=$'\t' read -r base path; do
		bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "${work}/job_fastqc.txt" \
		    fastqc --noextract -t 2 "${path}" -o "${qc_out}"
	done < "${samplelist}"
}

function trimming {
    rm -rf "${scratch}/trimming"
    mkdir -p "${scratch}/trimming"

   
    out="${scratch}/trimming"

    while IFS=$'\t' read -r base path; do
        bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "${work}/trimming_job.txt" \
         trim_galore --quality 15 --stringency 3 --report ${path} -o "${out}/${base}_trimmed"
    done < "${samplelist}"
}

function bowtiemapping { #LONG RUN TIME!, prior to running this build the ref genome off the filtered file using bowtie (LONG RUN TIME ASWELL)
    rm -rf ${scratch}/mapping
    mkdir -p ${scratch}/mapping
        
    refgen="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/hg38_rfgen/bwt2/hg38_filtered"
    out="${scratch}/mapping"
    
    while IFS=$'\t' read -r base path; do
        bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "${work}/bowtiemapping_job.txt"  -eo "${work}/bowtiemapping_job_err.txt" \
            bowtie2 -x ${refgen} \
                -U ${path} \
                -S ${out}/${base}.sam \
                --met-file ${out}/bowtie2_mapping_stats.txt
    done < "${samplelist}"
}

function mappingindex {
    rm -rf ${scratch}/mapping/indexed
    mkdir -p ${scratch}/mapping/indexed

    output="${scratch}/mapping/indexed"
    mapdir="${work}/mapping" #change to scratch once done

    while IFS=$'\t' read -r base path; do
      bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "${work}/bowtieindex_job.txt"  -eo "${work}/bowtieindex_job_err.txt" \
                "samtools sort ${mapdir}/${base}.sam -o ${output}/${base}.bam && \
                samtools index ${output}/${base}.bam"
    done < "${samplelist}"
}

function chipqc {
    rm -rf ${scratch}/results/chipqc
    mkdir -p ${scratch}/results/chipqc
    
    output="${scratch}/results/chipqc"
    data="${scratch}/mapping/indexed"

     bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "${work}/chipqc_job.txt" -eo "${work}/chipqc_eo_job.txt" \
        multiBamSummary bins -b ${data}/*.bam \
        -o ${output}/samples_correlation_matrix.npz \
        -bs 1000 \
        -n 500
}

function corplot {
    rm -rf "${work}/plot/correlationplot"
    mkdir -p "${work}/plot/correlationplot"
    
        data="${scratch}/results/chipqc/samples_correlation_matrix.npz"
        output="${work}/results/plot/correlationplot.png"

        bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "corplot_job.txt" \
            plotCorrelation -in ${data} \
                -c pearson \
                -p heatmap \
                --colorMap RdBu_r \
                --zMin 0.80 --zMax 1.0 \
                -o ${output}
 }

function ipplot {
    output1="${work}/results/plot/S7_ip_plot.png"
    output2="${work}/results/plot/S8_ip_plot.png"

    s7l1="${scratch}/mapping/indexed/SEM-CTCF-1_S7_L001_R1_001.fastq.gz.bam"
    s7l2="${scratch}/mapping/indexed/SEM-CTCF-1_S7_L002_R1_001.fastq.gz.bam"

    s8l1="${scratch}/mapping/indexed/SEM-CTCF-2_S8_L001_R1_001.fastq.gz.bam"
    s8l2="${scratch}/mapping/indexed/SEM-CTCF-2_S8_L002_R1_001.fastq.gz.bam"

    bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "S7_ipplot_job.txt" -eo "eo_S7_ipplot_job.txt" \
        iplot.py "${s7l1}" "${s7l2}"
    bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "S8_ipplot_job.txt" -eo "eo_S8_ipplot_job.txt" \
        iplot.py "${s8l1}" "${s8l2}"    
    
    
    
    bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "S7_ipplot_job.txt" -eo "eo_S7_ipplot_job.txt" \
        plotFingerprint -b "${s7l1}" "${s7l2}" \
        -l "SEM-CTCF-1_S7_L001" "SEM-CTCF-1_S7_L002" \
        -plot "${output1}" \
        -v \
        -bs 10000

    bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "S8_ipplot_job.txt" -eo "eo_S8_ipplot_job.txt" \
        plotFingerprint -b "${s8l1}" "${s8l2}" \
        -l "SEM-CTCF-2_S8_LOO1" "SEM-CTCF-2_S8_LOO2" \
        -plot "${output2}" \
        -v \
        -bs 10000
}

#fastqc_initial
#trimming
#bowtiemapping
#mappingindex
#chipqc
#corplot
ipplot
