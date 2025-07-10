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
jobs="/sc/arion/work/arayan01/projects/personalized_epigenome/data/2025-07-08_hichiphg38/jobs"
errmsg="/sc/arion/work/arayan01/projects/personalized_epigenome/data/2025-07-08_hichiphg38/jobs/eo"
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

function mapping {
    outdir="${scratch}/mapping"
    mkdir -p "$outdir"

    while IFS=$'\t' read -r basename filepath; do

        if [[ "$filepath" == *_R1*.fastq.gz ]]; then
            r1="$filepath"
            [[ -e "$r1" ]] || continue

            r2="${r1/_R1_/_R2_}"
            [[ -e "$r2" ]] || continue

            filename=$(basename "$r1")
            base=${filename%%_R1*}

            bsub -P acc_oscarlr -q premium -n ${cores} -W 4:00 -R "rusage[mem=16000]" \
                -o "${jobs}/pairtools_${base}.log" -eo "${errmsg}/pairtools_${base}.log" \
                bash -c "
                    bwa mem -5SP -T0 -t ${cores} ${ref} \"${r1}\" \"${r2}\" | \
                    pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 \
                        --nproc-in ${cores} --nproc-out ${cores} --chroms-path ${chrmsize}| \
                    pairtools sort --tmpdir=${tmp} --nproc ${cores} | \
                    pairtools dedup --nproc-in ${cores} --nproc-out ${cores} --mark-dups --output-stats ${outdir}/${base}_stats.txt | \
                    pairtools split --nproc-in ${cores} --nproc-out ${cores} \
                        --output-pairs ${outdir}/${base}_mapped.pairs --output-sam - | \
                    samtools view -bS -@${cores} | \
                    samtools sort -@${cores} -o ${outdir}/${base}_mapped.PT.bam && \
                    samtools index ${outdir}/${base}_mapped.PT.bam
                "
        fi
    done < "$samplelist"
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
mapping
#enrichment