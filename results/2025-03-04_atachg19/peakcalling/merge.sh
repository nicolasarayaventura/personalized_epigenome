set -x -e

scratch="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-03-04_atachg19"
work="/sc/arion/work/arayan01/project/personalized_epigenome/results/2025-03-04_atachg19"

function merge {
    merged_output="${work}/peakcalling/merged_SEM.bam"

    bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=16000]" -o "${work}/peakcalling/merged_job.txt" \
        "samtools merge -f ${merged_output} \
        ${work}/peakcalling/7_SEM-1_S53_L004_SC_subset_dedup.bam \
        ${work}/peakcalling/8_SEM-2_S54_L004_SC_subset_dedup.bam"
}

function merge_callpeak {
    merged_output="${work}/peakcalling/merged_SEM.bam"
    rm -rf ${work}/peakcalling/merged
    mkdir ${work}/peakcalling/merged
    output_dir="${work}/peakcalling/merged"

     bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=16000]" -o "${work}/peakcalling/merged_peaks_job.txt" \
        "macs2 callpeak -t ${merged_output} -f BAMPE -n ${output_dir}/merged_SEM_peak -g hs --keep-dup all"
}

function count_peaks {
    narrowpeak_file="${work}/peakcalling/merged/merged_SEM_peak_peaks.narrowPeak"  # Use correct path
    output_file="${work}/peakcalling/merged/peak_counts_hg19.txt"

     # Get the base name of the narrowPeak file (without the path)
    base_narrowpeak_file=$(basename "$narrowpeak_file")

    # Write header to the output file
    printf "File Name | Chromosome | Start | End | Locus | Peak Count\n" > "$output_file"

    # Define coordinate windows for counting
    coordinate_windows=(
        "chr2  88826861 90237547"
        "chr14  105422420 107043718"
        "chr22  22005516 22922912"
    )

    # Loop over each coordinate window and count peaks within each
    for window in "${coordinate_windows[@]}"; do
        chr=$(echo $window | cut -d ' ' -f 1)
        start=$(echo $window | cut -d ' ' -f 2)
        end=$(echo $window | cut -d ' ' -f 3)

        # Assign locus name based on chromosome
        if [[ "$chr" == "chr2" ]]; then
            locus="IgK"
        elif [[ "$chr" == "chr14" ]]; then
            locus="IgH"
        elif [[ "$chr" == "chr22" ]]; then
            locus="IgL"
        else
            locus="Unknown"
        fi

        # Count the number of peaks in the specified window
        count=$(awk -v chr="$chr" -v start="$start" -v end="$end" \
            '$1 == chr && $2 >= start && $3 <= end' "$narrowpeak_file" | wc -l)

        # Append results to the output file with the locus name included
        printf "%s | %s | %s | %s | %s | %d peaks\n" "$base_narrowpeak_file" "$chr" "$start" "$end" "$locus" "$count" >> "$output_file"
    done
}

#merge
#merge_callpeak
count_peaks