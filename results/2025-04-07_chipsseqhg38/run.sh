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

    peak_groups["RS411-CTCF"]="RS411-CTCF-1_S1_R1_001.fastq.gz_peaks_peaks.gz.bam RS411-CTCF-2_S2_R1_001.fastq.gz_peaks_peaks.gz.bam]"
    peak_groups["SEM-CTCF_S7"]="SEM-CTCF-1_S7_L001_R1_001.fastq.gz_peaks_peaks.gz.bam] SEM-CTCF-1_S7_L002_R1_001.fastq.gz_peaks_peaks.gz.bam]"
    peak_groups["SEM-CTCF_S8"]="SEM-CTCF-2_S8_L001_R1_001.fastq.gz_peaks_peaks.gz.bam] SEM-CTCF-2_S8_L002_R1_001.fastq.gz_peaks_peaks.gz.bam]"
    peak_groups["RS411-Rad21"]="RS411-Rad21-1_S3_R1_001.fastq.gz_peaks_peaks.gz.bam] RS411-Rad21-2_S4_R1_001.fastq.gz_peaks_peaks.gz.bam]"
    peak_groups["SEM-Rad21"]="SEM-Rad21-1_S5_R1_001.fastq.gz_peaks_peaks.gz.bam] SEM-Rad21-2_S6_R1_001.fastq.gz_peaks_peaks.gz.bam]"
    
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
            --effectiveGenomeSize 2423000000 \
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

#peak counts 05/27/2025

function merged_bam {
    bamdir="${scratch}/mapping/indexed/"
    outdir="${scratch}/mapping/merged_bam"
    mkdir -p "${outdir}"

    declare -A bamfile_groups
    bamfile_groups["RS411-CTCF"]="RS411-CTCF-1_S1_R1_001.fastq.gz.bam RS411-CTCF-2_S2_R1_001.fastq.gz.bam"
    bamfile_groups["SEM-CTCF_S7"]="SEM-CTCF-1_S7_L001_R1_001.fastq.gz.bam SEM-CTCF-1_S7_L002_R1_001.fastq.gz.bam"
    bamfile_groups["SEM-CTCF_S8"]="SEM-CTCF-2_S8_L001_R1_001.fastq.gz.bam SEM-CTCF-2_S8_L002_R1_001.fastq.gz.bam"
    bamfile_groups["RS411-Rad21"]="RS411-Rad21-1_S3_R1_001.fastq.gz.bam RS411-Rad21-2_S4_R1_001.fastq.gz.bam"
    bamfile_groups["SEM-Rad21"]="SEM-Rad21-1_S5_R1_001.fastq.gz.bam SEM-Rad21-2_S6_R1_001.fastq.gz.bam"

    for group in "${!bamfile_groups[@]}"; do
        files="${bamfile_groups[$group]}"
        inputs=""
        for file in $files; do
            inputs+=" ${bamdir}/${file}"
        done

        outfile="${outdir}/${group}.bam"

        bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=16000]" -o "${job}/${group}_merged_job.txt" \
            "samtools merge -f ${outfile} ${inputs}"
    done
}

function merge_callpeak {
    merged_bam_dir="${scratch}/mapping/merged_bam"
    outdir="${scratch}/mapping/merged_bam/peaks"
    mkdir -p "${outdir}"

    for bam in "${merged_bam_dir}"/*.bam; do
        sample_name=$(basename "${bam}" .bam)

        bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=16000]" \
            -o "${work}/${sample_name}_peakcalling_job.txt" \
            -e "${work}/${sample_name}_peakcalling_job.err" \
            "macs2 callpeak -t '${bam}' -n '${outdir}/${sample_name}' -g hs --keep-dup all"
    done
}

function mappingquality {
    rm -rf "${work}/mapping_quality"
    mkdir -p "${work}/mapping_quality"

    sampledir="${scratch}/mapping/merged_bam"

    for file in "${sampledir}"/*.bam; do
        sample_name=$(basename "${file}" .bam)
        bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" \
            -o "${work}/mapping_quality/${sample_name}_job.txt" \
            -eo "${work}/mapping_quality/${sample_name}_job_eo.txt" \
            python mapping_quality.py "${file}" "${work}/mapping_quality/${sample_name}_mapping_quality.png"
    done
}
function count_peaks {
    peaks_dir="${scratch}/mapping/merged_bam/peaks"
    output_file="${work}/peak_counts_hg38_all.txt"

    # Header for the combined output file
    printf "File Name | Chromosome | Start | End | Locus | Peak Count\n" > "$output_file"

    coordinate_windows=(
        "chr2 88826861 90237547"
        "chr14 105422420 107043718"
        "chr22 22005516 22922912"
    )

    # Loop over every .narrowPeak file
    for narrowpeak_file in "${peaks_dir}"/*.narrowPeak; do
        base_narrowpeak_file=$(basename "$narrowpeak_file")

        for window in "${coordinate_windows[@]}"; do
            chr=$(echo $window | cut -d ' ' -f1)
            start=$(echo $window | cut -d ' ' -f2)
            end=$(echo $window | cut -d ' ' -f3)

            case "$chr" in
                chr2) locus="IgK" ;;
                chr14) locus="IgH" ;;
                chr22) locus="IgL" ;;
                *) locus="Unknown" ;;
            esac

            # Count peaks overlapping the window (partial overlap)
            count=$(awk -v chr="$chr" -v start="$start" -v end="$end" \
                '$1 == chr && $3 > start && $2 < end' "$narrowpeak_file" | wc -l)

            printf "%s | %s | %s | %s | %s | %d\n" \
                "$base_narrowpeak_file" "$chr" "$start" "$end" "$locus" "$count" >> "$output_file"
        done
    done
}

function atac_chipmerge {
    atacdir="/sc/arion/work/arayan01/project/personalized_epigenome/results/2025-02-24_atachg38/peakcalling" 
    chipdir="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-04-07_chipsseqhg38/mapping/merged_bam/peaks"
    outdir="${scratch}/atac_chip_merge"
    mkdir -p "$outdir"
    
    declare -A mergegroups
    mergegroups["SEM-CTCF"]="SEM-CTCF_S7_peaks.narrowPeak SEM-CTCF_S8_peaks.narrowPeak"
    mergegroups["SEM-ATAC"]="7_SEM-1_S53_L004_peak_peaks.narrowPeak 8_SEM-2_S54_L004_peak_peaks.narrowPeak"
    
    for group in "${!mergegroups[@]}"; do
        read -r file1 file2 <<< "${mergegroups[$group]}"

        # Route to correct directories
        if [[ "$group" == *"CTCF"* ]]; then
            path1="${chipdir}/${file1}"
            path2="${chipdir}/${file2}"
        else
            path1="${atacdir}/${file1}"
            path2="${atacdir}/${file2}"
        fi

        if [[ -f "$path1" && -f "$path2" ]]; then
            bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=16000]" \
            -o "${job}/${group}_merge_job.txt" \
            -e "${job}/${group}_merge_job.err" \
            -- \
            bash -c "cat \"$path1\" \"$path2\" | sort -k1,1 -k2,2n | bedtools merge > \"${outdir}/${group}_merged.bed\""
        fi
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
#heatplots
#merged_bam
#merge_callpeak
#count_peaks
#mappingquality
#atac_chipmerge

