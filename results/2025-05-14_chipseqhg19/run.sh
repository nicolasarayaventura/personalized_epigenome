set -x -e

scratch="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-05-14_chipseqhg19"
work="/sc/arion/work/arayan01/project/personalized_epigenome/results/2025-05-14_chipseqhg19"
job="${work}/job"

samplelist="/sc/arion/work/arayan01/project/personalized_epigenome/results/2025-04-07_chipsseqhg38/sample_ids.txt"


function bowtiemapping { #LONG RUN TIME!, prior to running this build the ref genome off the filtered file using bowtie (LONG RUN TIME ASWELL)
    rm -rf ${scratch}/mapping
    mkdir -p ${scratch}/mapping
        
    refgen="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/hg19_rfgen/bwt/hg19_filtered"
    out="${scratch}/mapping"
    trimdir="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-04-07_chipsseqhg38/trimming"

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
unction mergepeaks {
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
#bowtiemapping
#mappingindex
#callpeaks
bamcov_merge #05/20/2025 4:55pm
#index_merged_bam
#mergepeaks
#bamcov_run
#matrix
#heatlplot