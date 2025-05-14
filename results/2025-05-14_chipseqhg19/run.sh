set -x -e

scratch="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-05-14_chipseqhg19"
work="/sc/arion/work/arayan01/project/personalized_epigenome/results/2025-05-14_chipseqhg19"
job="${work}/job"

samplelist="/sc/arion/work/arayan01/project/personalized_epigenome/results/2025-04-07_chipsseqhg38/sample_ids.txt"


function bowtiemapping { #LONG RUN TIME!, prior to running this build the ref genome off the filtered file using bowtie (LONG RUN TIME ASWELL)
    rm -rf ${scratch}/mapping
    mkdir -p ${scratch}/mapping
        
    refgen="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/hg38_rfgen/bwt/hg19_filtered"
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
        output="${work}/results/plot/correlationplot.png"

        bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "${job}/corplot_job.txt" \
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
    mkdir -p "${scratch}/heatplot/prep"
    mkdir -p "${scratch}/heatplot/merged_bam"

    bamdir="${scratch}/mapping/indexed"
    merge_dir="${scratch}/heatplot/merged_bam"

    declare -A sample_bams_ctcf1
    declare -A sample_bams_ctcf2

    # Group BAMs by sample prefix (SEM-CTCF-1 or SEM-CTCF-2)
    while IFS=$'\t' read -r base path; do
        if [[ "$base" == SEM-CTCF-1* ]]; then
            sample_bams_ctcf1["SEM-CTCF-1"]="${sample_bams_ctcf1["SEM-CTCF-1"]} ${bamdir}/${base}.bam"
        elif [[ "$base" == SEM-CTCF-2* ]]; then
            sample_bams_ctcf2["SEM-CTCF-2"]="${sample_bams_ctcf2["SEM-CTCF-2"]} ${bamdir}/${base}.bam"
        fi
    done < "${samplelist}"

    # Function to merge BAM files for a sample
    merge_bam_files() {
        sample=$1
        sample_bams=$2
        merged_bam="${merge_dir}/${sample}.bam"

        bsub -P acc_oscarlr -q premium -n 2 -W 8:00 -R "rusage[mem=8000]" \
            -o "${sample}_merge_bam_job.txt" -eo "${sample}_merge_bam_err.txt" \
            samtools merge -f "$merged_bam" $sample_bams
    }

    merge_bam_files "SEM-CTCF-1" "${sample_bams_ctcf1["SEM-CTCF-1"]}"
    merge_bam_files "SEM-CTCF-2" "${sample_bams_ctcf2["SEM-CTCF-2"]}"
}

function index_merged_bams {
    bamdir="${scratch}/heatplot/merged_bam"

    for bam in ${bamdir}/*.bam; do
        bsub -P acc_oscarlr -q premium -n 1 -W 01:00 -R "rusage[mem=2000]" \
            -o "${bam%.bam}_index.out" -eo "${bam%.bam}_index.err" \
            samtools index "$bam"
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
function merge_ctcf1_peaks {
    outdir="${scratch}/heatplot/merged_peaks"
    peakdir="${scratch}/peaks"

    mkdir -p "${outdir}"

    rm -f "${outdir}/concatenated_ctcf1_peaks.narrowPeak"

    cat ${peakdir}/SEM-CTCF-1_S7_L001_R1_001.fastq.gz_peaks_peaks.narrowPeak \
        ${peakdir}/SEM-CTCF-1_S7_L002_R1_001.fastq.gz_peaks_peaks.narrowPeak > \
        "${outdir}/concatenated_ctcf1_peaks.narrowPeak"

    bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" \
        -o "${outdir}/merge_ctcf1_peaks_job.txt" \
        "bedtools sort -i ${outdir}/concatenated_ctcf1_peaks.narrowPeak | \
         bedtools merge -i - > ${outdir}/merged_ctcf1_peaks.bed"
}
function merge_ctcf2_peaks {
    peakdir="${scratch}/peaks"
    outdir="${scratch}/heatplot/merged_peaks"

    # Merge peaks for CTCF-2 samples
    cat ${peakdir}/SEM-CTCF-2_S8_L001_R1_001.fastq.gz_peaks_peaks.narrowPeak \
        ${peakdir}/SEM-CTCF-2_S8_L002_R1_001.fastq.gz_peaks_peaks.narrowPeak > \
        "${outdir}/concatenated_ctcf2_peaks.narrowPeak"

    bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" \
        -o "${outdir}/merge_ctcf2_peaks_job.txt" \
        "bedtools sort -i ${outdir}/concatenated_ctcf2_peaks.narrowPeak | \
         bedtools merge -i - > ${outdir}/merged_ctcf2_peaks.bed"
}

function matrix1 {
    beddir="${scratch}/heatplot/merged_peaks"
    bwdir="${scratch}/heatplot/prep"
    outdir="${scratch}/heatplot/prep"

    bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" \
        -o "matrix_job.out" -eo "matrix_job.err" \
        computeMatrix reference-point \
        -S "${bwdir}/SEM-CTCF-1.bw" \
        -R "${beddir}/merged_ctcf1_peaks.bed"  \
        --referencePoint center \
        --upstream 3000 \
        --downstream 3000 \
        --outFileName "${outdir}/CTCF_matrix.gz" \
        --skipZeros \
        --verbose
}

function matrix2 {
    beddir="${scratch}/heatplot/merged_peaks"
    bwdir="${scratch}/heatplot/prep"
    outdir="${scratch}/heatplot/prep"

    bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" \
        -o "matrix_job.out" -eo "matrix_job.err" \
        computeMatrix reference-point \
        -S "${bwdir}/SEM-CTCF-2.bw" \
        -R "${beddir}/merged_ctcf2_peaks.bed"  \
        --referencePoint center \
        --upstream 3000 \
        --downstream 3000 \
        --outFileName "${outdir}/CTCF2_matrix.gz" \
        --skipZeros \
        --verbose
}

function heatp1 {
    outdir="${work}/plot/heatlplot"
    matrixfile="${scratch}/heatplot/prep/CTCF_matrix.gz"

        bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "CTCF__heatmap_job.txt" -eo "CTCF__heatmap_err.txt" \
            plotHeatmap -m ${matrixfile} \
            --heatmapWidth 8 \
            --heatmapHeight 8 \
            --colorMap RdYlBu \
            -o ${outdir}/CTCF_k4_heatmap.png \
            --kmeans 4
                
}

function heatp2 {
    # rm -rf "${work}/plot/heatlplot"
    #mkdir -p "${work}/plot/heatlplot"

    outdir="${work}/plot/heatlplot"
    matrixfile="${scratch}/heatplot/prep/CTCF2_matrix.gz"

        bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "CTCF2_heatmap_job.txt" -eo "CTCF2__heatmap_err.txt" \
            plotHeatmap -m ${matrixfile} \
            --heatmapWidth 8 \
            --heatmapHeight 8 \
            --colorMap RdYlBu \
            -o ${outdir}/CTCF2_k4_heatmap.png \
            --kmeans 4
                
}

bowtiemapping
#mappingindex
#chipqc
#corplot
#callpeaks
#bamcov_merge
#index_merged_bams 
#merge_ctcf1_peaks
#merge_ctcf2_peaks
#bamcov_run
#matrix1
#matrix2 
#heatp1
#heatp2


#05/13/2025
# left off at mapping job, finished atac MAPQ dist. chart but make sure the axis are correct compared to this also after this start hg19?
#https://www.acgt.me/blog/2014/12/16/understanding-mapq-scores-in-sam-files-does-37-42