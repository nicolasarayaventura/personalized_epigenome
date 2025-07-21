set -x -e

###
work="/sc/arion/work/arayan01/project/personalized_epigenome/results/2025-07-08_hichiphg38"
scratch="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-07-08_hichiphg38"
sampledir1="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-02-19_ranjan_data/data/2022-11-22_data_from_Ranjan_Xiang/hichip/RSen090622"
sampledir2="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-02-19_ranjan_data/data/2022-11-22_data_from_Ranjan_Xiang/hichip/RSen090722"
###
tmp="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-07-08_hichiphg38/tmp"
ref="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/hg38_rfgen/bwa/hg38_filtered"
chrmsize="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/hg38_rfgen/bwa/hg38_filtered.chrom.sizes"
###
jobs="/sc/arion/work/arayan01/project/personalized_epigenome/results/2025-07-08_hichiphg38/jobs"
errmsg="/sc/arion/work/arayan01/project/personalized_epigenome/data/2025-07-08_hichiphg38/jobs/eo"
###
cores=4
samplelist="${work}/sample_ids.txt"

function samplelists {
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

function bmem {
    rm -rf "${scratch}/bmem"

    outdir="${scratch}/bmem"  
    mkdir -p "${outdir}"       

    while IFS=$'\t' read -r basename filepath; do
        if [[ "$filepath" == *_R1*.fastq.gz ]]; then
            r1="$filepath"
            [[ -e "$r1" ]] || continue

            r2="${r1/_R1_/_R2_}"
            [[ -e "$r2" ]] || continue

            filename=$(basename "$r1")
            base=${filename%%_R1*}

            bsub -P acc_oscarlr -q premium -n ${cores} -W 4:00 -R "rusage[mem=16000]" \
                -o "${jobs}/bmem_${base}.log" -eo "${errmsg}/bmem_${base}.log" \
                bash -c "bwa mem -5SP -T0 -t ${cores} ${ref} ${r1} ${r2} > ${outdir}/${base}_aligned.sam"
        fi
    done < "$samplelist"
}

function bmem_chip_peakcall {
    rm -rf "${scratch}/chip_bam"
    mkdir -p "${scratch}/chip_bam"

    while IFS=$'\t' read -r basename filepath; do
        if [[ "$filepath" == *_R1*.fastq.gz ]]; then
            r1="$filepath"
            [[ -e "$r1" ]] || continue

            r2="${r1/_R1_/_R2_}"
            [[ -e "$r2" ]] || continue

            filename=$(basename "$r1")
            base=${filename%%_R1*}

            bsub -P acc_oscarlr -q premium -n ${cores} -W 4:00 -R "rusage[mem=16000]" \
                -o "${jobs}/chip_${base}.log" -eo "${errmsg}/chip_${base}.log" \
                bash -c "bwa mem -t ${cores} ${ref} ${r1} ${r2} | \
                         samtools view -bS -q 30 - | \
                         samtools sort -@ ${cores} -o ${scratch}/chip_bam/${base}.bam"
        fi
    done < "$samplelist"
}


function pt_parse {
    samdir="${scratch}/bmem"
    outdir="${scratch}/parsed"
    mkdir -p "$outdir"

    for file in "${samdir}"/*_aligned.sam; do
        [[ -e "$file" ]] || continue
        base=$(basename "$file" "_aligned.sam")

        bsub -P acc_oscarlr -q premium -n ${cores} -W 2:00 -R "rusage[mem=8000]" \
            -o "${jobs}/parse_${base}.log" -eo "${errmsg}/parse_${base}_eo.log" \
            bash -c "pairtools parse \
                --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 \
                --nproc-in ${cores} --nproc-out ${cores} --chroms-path ${chrmsize} \
                '${file}' > '${outdir}/${base}.parsed.pairs'"
    done
}

function pt_sort {
    pardir="${scratch}/parsed"
    outdir="${scratch}/sorted"
    mkdir -p "$outdir"

    for file in "${pardir}"/*.parsed.pairs; do
        [[ -e "$file" ]] || continue  # Skip if no parsed.pairs files found
        base=$(basename "$file" ".parsed.pairs")

        bsub -P acc_oscarlr -q premium -n ${cores} -W 2:00 -R "rusage[mem=8000]" \
            -o "${jobs}/sort_${base}.log" -eo "${errmsg}/sort_${base}_eo.log" \
            bash -c "pairtools sort --tmpdir=/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-07-08_hichiphg38/tmp --nproc ${cores} < ${file} > ${outdir}/${base}.sorted.pairs"
    done
}

function pt_dedup {
    sortdir="${scratch}/sorted"
    outdir="${scratch}/deduped"
    mkdir -p "$outdir"

    for file in "${sortdir}"/*.sorted.pairs; do
        [[ -e "$file" ]] || continue
        base=$(basename "$file" ".sorted.pairs")

        bsub -P acc_oscarlr -q premium -n ${cores} -W 2:00 -R "rusage[mem=8000]" \
            -o "${jobs}/dedup_${base}.log" -eo "${errmsg}/dedup_${base}_eo.log" \
            bash -c "pairtools dedup --nproc-in ${cores} --nproc-out ${cores} --mark-dups \
                --output-stats ${outdir}/${base}_stats.txt \
                --output ${outdir}/${base}.pairsam \
                ${file}"
    done
}

function pt_split {
    dedupdir="${scratch}/deduped"
    outdir="${scratch}/final_pairs"
    mkdir -p "$outdir"

    for file in "${dedupdir}"/*.pairsam; do
        [[ -e "$file" ]] || continue
        base=$(basename "$file" ".pairsam")

        bsub -P acc_oscarlr -q premium -n ${cores} -W 2:00 -R "rusage[mem=8000]" \
            -o "${jobs}/split_${base}.log" -eo "${errmsg}/split_${base}_eo.log" \
            bash -c "pairtools split \
                --nproc-in ${cores} --nproc-out ${cores} \
                --output-pairs ${outdir}/${base}.pairs \
                --output-sam ${outdir}/${base}.bam \
                ${file}"
    done
}

function samtools_sort_and_index {
    bamdir="${scratch}/final_pairs"   # directory with unsorted BAM files
    sorted_dir="${scratch}/sorted_bam"
    mkdir -p "$sorted_dir" "$jobs" "$errmsg"

    for file in "${bamdir}"/*.bam; do
        [[ -e "$file" ]] || continue
        base=$(basename "$file" ".bam")

        bsub -P acc_oscarlr -q premium -n ${cores} -W 3:00 -R "rusage[mem=8000]" \
            -o "${jobs}/sort_index_${base}.log" -eo "${errmsg}/sort_index_${base}_eo.log" \
            bash -c "samtools sort -@ ${cores} -T /sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-07-08_hichiphg38/tmp/${base}_temp \
                -o ${sorted_dir}/${base}_sorted.bam ${file} && samtools index ${sorted_dir}/${base}_sorted.bam"
    done
}

function run_get_qc {
    dedupdir="${scratch}/deduped"
    qc_outdir="${work}/qc_reports"
    mkdir -p "$qc_outdir" "$jobs" "$errmsg"

    merged_qc="${qc_outdir}/merged_qc_summary.txt"
    echo -n "" > "$merged_qc"  # clear previous content

    for statsfile in "${dedupdir}"/*_stats.txt; do
        [[ -e "$statsfile" ]] || continue
        base=$(basename "$statsfile" "_stats.txt")

        bsub -P acc_oscarlr -q premium -n 1 -W 1:00 -R "rusage[mem=4000]" \
            -o "${jobs}/qc_${base}.log" -eo "${errmsg}/qc_${base}_eo.log" \
            bash -c "python3 /sc/arion/scratch/arayan01/projects/hichip_tut/data/HiChiP/get_qc.py -p '${statsfile}' > '${qc_outdir}/${base}_qc.txt' && \
                     echo -e '=== ${base} ===' >> '${merged_qc}' && \
                     cat '${qc_outdir}/${base}_qc.txt' >> '${merged_qc}' && \
                     echo >> '${merged_qc}'"
    done
}



function enrichment {
    mapdir="${scratch}/mapping"

    for bam in "${mapdir}"/*_mapped.PT.bam; do
        [[ -e "$bam" ]] || continue  # Skip if no BAM files found

        base=$(basename "$bam" "_mapped.PT.bam")

        bsub -P acc_oscarlr -q premium -n 4 -W 1:00 -R "rusage[mem=4000]" \
            -o "${jobs}/enrich_${base}.log" -eo "${errmsg}/enrich_${base}.log" \
            bash /sc/arion/scratch/arayan01/projects/hichip_tut/data/HiChiP/enrichment_stats.sh \
            -g "$ref" \
            -b "$bam" \
            -p "$peakfile" \
            -t 16 \
            -x CTCF
    done
}

function global_enrichment {
    chromsizes="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/hg38_rfgen/hg38.chrom.sizes"
    tmpdir="/sc/arion/scratch/arayan01/projects/hichip_tut/data/tmp"
    stats="/sc/arion/scratch/arayan01/projects/hichip_tut/results/mapping_stats.txt"
    outbam="/sc/arion/scratch/arayan01/projects/hichip_tut/results/mapped.PT.bam"
    pairs="/sc/arion/scratch/arayan01/projects/hichip_tut/results/mapped.pairs"
    cores=4
    #_bed for is for bed files other one is for narrowpeak files
    bsub -P acc_oscarlr -q premium -n 4 -W 1:00 -R "rusage[mem=4000]" -o global_enrichment_job.txt -eo global_enrichment_job_eo.txt \
        python3 /sc/arion/scratch/arayan01/projects/hichip_tut/data/HiChiP/plot_chip_enrichment_bed.py \
        -bam /sc/arion/scratch/arayan01/projects/hichip_tut/results/mapped.PT.bam \
        -peaks /sc/arion/scratch/arayan01/projects/hichip_tut/data/samples/ENCFF017XLW.bed \
        -output /sc/arion/scratch/arayan01/projects/hichip_tut/results/enrichment.png

}
function contact_matrix_juicer {
        bsub -P acc_oscarlr -q premium -n 4 -W 1:00 -R "rusage[mem=4000]" -o juicer_job.txt -eo juicer_eo.txt \
            java -Xmx48000m  \
            -Djava.awt.headless=true \
            -jar /sc/arion/scratch/arayan01/projects/hichip_tut/data/HiChiP/juicertools.jar \
            pre --threads 16 \
            /sc/arion/scratch/arayan01/projects/hichip_tut/results/mapped.pairs \
            contact_map.hic \
            /sc/arion/scratch/arayan01/projects/hichip_tut/results/hg38.chrom.sizes

}

function genpairx {
    pfile=/sc/arion/scratch/arayan01/projects/hichip_tut/results/mapped.pairs
    bgzip ${pfile}
    bsub -P acc_oscarlr -q premium -n 4 -W 1:00 -R "rusage[mem=4000]" -o genpairx_job.txt -eo genpairx_job_eo.txt \
        pairix ${pfile}.gz 
}
function viscool1 {
    pfile=/sc/arion/scratch/arayan01/projects/hichip_tut/results/mapped.pairs.gz
    chrom=/sc/arion/scratch/arayan01/projects/hichip_tut/results/hg38.chrom.sizes
    coolfile=/sc/arion/scratch/arayan01/projects/hichip_tut/results/matrix_1kb.cool
    #bsub -P acc_oscarlr -q premium -n 4 -W 1:00 -R "rusage[mem=4000]" -o vis_job.txt -eo vis_job_eo.txt \
     #   "cooler cload pairix -p 16 ${chrom}:5000 ${pfile} matrix_1kb.cool"
#use cooler zoomify --balance -p 16 matrix_1kb.cool afterwards
    bsub -P acc_oscarlr -q premium -n 4 -W 1:00 -R "rusage[mem=4000]" -o vis_job.txt -eo vis_job_eo.txt \
        cooler zoomify --balance -p 16 ${coolfile}
}

function loopcalling {
    grep -v '#' /sc/arion/scratch/arayan01/projects/hichip_tut/results/mapped.pairs| awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"$7}' | gzip -c > hicpro_mapped.pairs.gz
    fithichip=/sc/arion/scratch/arayan01/projects/hichip_tut/data/FitHiChIP/FitHiChIP_Docker.sh
    config=/sc/arion/scratch/arayan01/projects/hichip_tut/data/config.txt
    bsub -P acc_oscarlr -q premium -n 4 -W 1:00 -R "rusage[mem=4000]" -o loop_job.txt -eo loop_job_eo.txt \
        bash ${fithichip} -C ${config}
}

#samplelists
#bmem
#bmem_chip_peakcall
#pt_parse
#pt_sort
#pt_dedup
#pt_split
#samtools_sort_and_index
run_get_qc
#enrichmentß