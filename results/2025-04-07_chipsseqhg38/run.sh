set -x -e

scratch="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-04-07_chipsseqhg38"
work="/sc/arion/work/arayan01/project/personalized_epigenome/results/2025-04-07_chipsseqhg38"
job="${work}/job"

sampledir1="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-02-19_ranjan_data/2022-11-22_data_from_Ranjan_Xiang/chip/RSen042022_24"
sampledir2="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-02-19_ranjan_data/2022-11-22_data_from_Ranjan_Xiang/chip/RSen072522_10"

samplelist="${work}/sample_ids.txt"


function samplelist_generate {
    output="${work}/sample_ids.txt"
    for filepath in ${sampledir1}/*; do
        [ -f "$filepath" ] || continue
        basename=$(basename "$filepath")
        echo -e "${basename}\t${filepath}" >> "$output"
    done

    for filepath in ${sampledir2}/*; do
        [ -f "$filepath" ] || continue
        basename=$(basename "$filepath")
        echo -e "${basename}\t${filepath}" >> "$output"
    done
}

function fastqc_initial {
    rm -rf "${work}/qc/init_fastqc"
    mkdir "${work}/qc/init_fastqc"
    
	qc_out="${work}/qc/init_fastqc" #once job is finished run multiqc in directory!

	while IFS=$'\t' read -r base path; do
		bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "${job}/job_fastqc.txt" \
		    fastqc --noextract -t 2 "${path}" -o "${qc_out}"
	done < "${samplelist}"
}

function trimming {
    rm -rf "${scratch}/trimming"
    mkdir -p "${scratch}/trimming"

    out="${scratch}/trimming"

    while IFS=$'\t' read -r base path; do
        bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "${job}/trimming_job.txt" \
         trim_galore --quality 15 --stringency 3 --report ${path} -o "${out}"
    done < "${samplelist}"
}

function fastqc_posttrimming {
    rm -rf "${work}/qc/trimming_fastqc"
    mkdir "${work}/qc/trimming_fastqc"
    
    trimdir="${scratch}/trimming"

	qc_out="${work}/qc/trimming_fastqc" #once job is finished run multiqc in directory!

	while IFS=$'\t' read -r base path; do
        clean_base="${base%.fastq.gz}"
        fq_file="${trimdir}/${clean_base}_trimmed.fq.gz"

        if [[ -f "$fq_file" ]]; then
            bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "${job}/posttrimming_job_fastqc.txt" -eo "eo_posttrimming_job_fastqc.txt" \
                fastqc --noextract -t 2 "$fq_file" -o "$qc_out"   
        fi   
    done < "${samplelist}"
}
function bowtiemapping { #LONG RUN TIME!, prior to running this build the ref genome off the filtered file using bowtie (LONG RUN TIME ASWELL)
    rm -rf ${scratch}/mapping
    mkdir -p ${scratch}/mapping
        
    refgen="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/hg38_rfgen/bwt/hg38_filtered"
    out="${scratch}/mapping"
    trimdir="${scratch}/trimming"

    while IFS=$'\t' read -r base path; do
        clean_base="${base%.fastq.gz}"
        trimmed_file="${trimdir}/${clean_base}_trimmed.fq.gz"
        if [[ -f "$trimmed_file" ]]; then
            bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" \
                -o "${job}/bowtiemapping_job.txt" \
                -eo "bowtiemapping_job_err.txt" \
                bowtie2 -x "${refgen}" \
                        -U "${trimmed_file}" \
                        -S "${out}/${base}.sam" \
                        --met-file "${out}/bowtie2_mapping_stats.txt"
        fi
    done < "${samplelist}"
}

function mappingindex {
    rm -rf ${scratch}/mapping/indexed
    mkdir -p ${scratch}/mapping/indexed

    output="${scratch}/mapping/indexed"
    mapdir="${scratch}/mapping" #change to scratch once done

    while IFS=$'\t' read -r base path; do
      bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "${job}/bowtieindex_job.txt"  -eo "${job}/bowtieindex_job_err.txt" \
                "samtools sort ${mapdir}/${base}.sam -o ${output}/${base}.bam && \
                samtools index ${output}/${base}.bam"
    done < "${samplelist}"
}

function chipqc {
    rm -rf ${scratch}/results/chipqc
    mkdir -p ${scratch}/results/chipqc
    
    output="${scratch}/results/chipqc"
    data="${scratch}/mapping/indexed"

     bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "${job}/chipqc_job.txt" -eo "${job}/chipqc_eo_job.txt" \
        multiBamSummary bins -b ${data}/*.bam \
        -o ${output}/samples_correlation_matrix.npz \
        -bs 1000 \
        -n 500
}

function corplot {
    rm -rf "${work}/plot/correlationplot"
    mkdir -p "${work}/plot/correlationplot"
    
        data="${scratch}/results/chipqc/samples_correlation_matrix.npz"
        output="${work}/plot/correlationplot.png"

        bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "${job}/corplot_job.txt" -eo "corplot_eo.txt"\
            plotCorrelation -in ${data} \
                -c pearson \
                -p heatmap \
                --colorMap RdBu_r \
                --zMin 0.80 --zMax 1.0 \
                -o ${output}
 }

function callpeaks {
    rm -rf "${scratch}/peaks"
    mkdir -p "${scratch}/peaks"

    input_dir="${scratch}/mapping/indexed"
    output_dir="${scratch}/peaks"

    while IFS=$'\t' read -r base path; do
        bam="${input_dir}/${base}.bam"

        bsub -P acc_oscarlr -q premium -n 2 -W 12:00 -R "rusage[mem=8000]" \
            -o "${work}/peakcall_${base}.out" -eo "${work}/peakcall_${base}.err" \
            macs2 callpeak -t "${bam}" \
                -f BAM \
                -g hs \
                -n "${base}_peaks" \
                --outdir "${output_dir}" \
                --keep-dup all \
                --nomodel \
                --extsize 200 \
                -q 0.01
    done < "${samplelist}"
}

function bamcov_merge {
    rm -rf "${scratch}/heatplot/prep"
    rm -rf "${scratch}/heatplot/merged_bam"
    
    mkdir -p "${scratch}/heatplot/prep"
    mkdir -p "${scratch}/heatplot/merged_bam"

    bamdir="${scratch}/mapping/indexed"

    merged_names=(
        "SEM-CTCF-1_S7_merged.bam"
        "SEM-CTCF-2_S8_merged.bam"
        "RS411-CTCF_merged.bam"
        "RS411-Rad21_merged.bam"
        "SEM-Rad21_merged.bam"
    )

    for i in {0..4}; do
        start_line=$((i * 2 + 1))
        end_line=$((start_line + 1))
        pair=$(sed -n "${start_line}p;${end_line}p" "$samplelist")

        bams=()
        while IFS=$'\t' read -r base path; do
            bams+=("${bamdir}/${base}.bam")
        done <<< "$pair"

        merged="${scratch}/heatplot/merged_bam/${merged_names[$i]}"
        bsub -P acc_oscarlr -q premium -n 2 -W 8:00 -R "rusage[mem=8000]" \
            -o "${job}/pair${i}_merge_bam_job.txt" -eo "${job}/pair${i}_merge_bam_err.txt" \
            samtools merge -f "$merged" "${bams[@]}"
    done
}

function index_merged_bam {
    mergeddir="${scratch}/heatplot/merged_bam"

    for bam in "${mergeddir}"/*.bam; do
        bsub -P acc_oscarlr -q premium -n 1 -W 01:00 -R "rusage[mem=2000]" \
            -o "${job}/$(basename "$bam").index.out" \
            -eo "${job}/$(basename "$bam").index.err" \
            samtools index "$bam"
    done
}
function mergepeaks {
    peakdir="${scratch}/peaks"
    outdir="${scratch}/heatplot/merged_peaks"
    mkdir -p "${outdir}"

    declare -A peak_groups

    peak_groups["RS411-CTCF"]="RS411-CTCF-1_S1_R1_001.fastq.gz_peaks_peaks.narrowPeak RS411-CTCF-2_S2_R1_001.fastq.gz_peaks_peaks.narrowPeak"
    peak_groups["SEM-CTCF_S7"]="SEM-CTCF-1_S7_L001_R1_001.fastq.gz_peaks_peaks.narrowPeak SEM-CTCF-1_S7_L002_R1_001.fastq.gz_peaks_peaks.narrowPeak"
    peak_groups["SEM-CTCF_S8"]="SEM-CTCF-2_S8_L001_R1_001.fastq.gz_peaks_peaks.narrowPeak SEM-CTCF-2_S8_L002_R1_001.fastq.gz_peaks_peaks.narrowPeak"
    peak_groups["RS411-Rad21"]="RS411-Rad21-1_S3_R1_001.fastq.gz_peaks_peaks.narrowPeak RS411-Rad21-2_S4_R1_001.fastq.gz_peaks_peaks.narrowPeak"
    peak_groups["SEM-Rad21"]="SEM-Rad21-1_S5_R1_001.fastq.gz_peaks_peaks.narrowPeak SEM-Rad21-2_S6_R1_001.fastq.gz_peaks_peaks.narrowPeak"
    
    for group in "${!peak_groups[@]}"; do
        echo "Processing $group"
        files="${peak_groups[$group]}"
        inputs=""
        for file in $files; do
            inputs+="${peakdir}/${file} "
        done

        bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" \
            -o "${outdir}/merge_${group}.out" -eo "${outdir}/merge_${group}.err" \
            "cat ${inputs} | bedtools sort -i - | bedtools merge -i - > ${outdir}/merged_${group}_peaks.bed"
    done
}

function bamcov_run {
    rm -rf "${scratch}/heatplot/prep"
    mkdir -p "${scratch}/heatplot/prep"

    bamdir="${scratch}/heatplot/merged_bam"
    outdir="${scratch}/heatplot/prep"

    for bam in ${bamdir}/*.bam; do
        sample=$(basename "$bam" .bam)

        bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" \
            -o "${sample}_bamcov_job.txt" -eo "${sample}_bamcov_err.txt" \
            bamCoverage -b "$bam" \
            -o "${outdir}/${sample}.bw" \
            --outFileFormat bigwig \
            --binSize 25 \
            --normalizeUsing RPGC \
            --effectiveGenomeSize 2913022398 \
            --verbose
    done
}

function matrix {
    beddir="${scratch}/heatplot/merged_peaks"
    bwdir="${scratch}/heatplot/prep"
    outdir="${scratch}/heatplot/prep"

    declare -A bw_to_bed
    bw_to_bed["SEM-CTCF-1_S7_merged"]="merged_SEM-CTCF_S7_peaks.bed"
    bw_to_bed["SEM-CTCF-2_S8_merged"]="merged_SEM-CTCF_S8_peaks.bed"
    bw_to_bed["RS411-CTCF_merged"]="merged_RS411-CTCF_peaks.bed"
    bw_to_bed["RS411-Rad21_merged"]="merged_RS411-Rad21_peaks.bed"
    bw_to_bed["SEM-Rad21_merged"]="merged_SEM-Rad21_peaks.bed"

    for label in "${!bw_to_bed[@]}"; do
        bw="${bwdir}/${label}.bw"
        bed="${beddir}/${bw_to_bed[$label]}"
        outmatrix="${outdir}/${label}_matrix.gz"
        outlog=" ${label}_matrix_job.out"
        outerr="${label}_matrix_job.err"

        bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" \
            -o "$outlog" -eo "$outerr" \
            computeMatrix reference-point \
            -S "$bw" \
            -R "$bed" \
            --referencePoint center \
            --upstream 3000 \
            --downstream 3000 \
            --outFileName "$outmatrix" \
            --skipZeros \
            --verbose
    done
}

function heatplots {
    outdir="${work}/plot/heatlplot"

    declare -A matrixfiles
    matrixfiles["SEM-CTCF-1_S7"]="SEM-CTCF-1_S7_merged_matrix.gz"
    matrixfiles["SEM-CTCF-2_S8"]="SEM-CTCF-2_S8_merged_matrix.gz"
    matrixfiles["RS411-CTCF"]="RS411-CTCF_merged_matrix.gz"
    matrixfiles["RS411-Rad21"]="RS411-Rad21_merged_matrix.gz"
    matrixfiles["SEM-Rad21"]="SEM-Rad21_merged_matrix.gz"
    
    for label in "${!matrixfiles[@]}"; do
        matrixfile="${scratch}/heatplot/prep/${matrixfiles[$label]}"
        outpng="${outdir}/${label}_heatmap.png"

        bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" \
            -o "${label}_heatmap_job.txt" -eo "${label}_heatmap_err.txt" \
            plotHeatmap -m "$matrixfile" \
            --heatmapWidth 8 \
            --heatmapHeight 8 \
            --colorMap RdYlBu \
            -o "$outpng" \
            --kmeans 4
    done
}


#samplelist
#fastqc_initial
#trimming
#fastqc_posttrimming
#bowtiemapping
#mappingindex
#chipqc
#corplot
#callpeaks
#bamcov_merge
#index_merged_bam 
#mergepeaks
#bamcov_run
#matrix
heatplots


#05/13/2025
# left off at mapping job, finished atac MAPQ dist. chart but make sure the axis are correct compared to this also after this start hg19?
#https://www.acgt.me/blog/2014/12/16/understanding-mapq-scores-in-sam-files-does-37-42