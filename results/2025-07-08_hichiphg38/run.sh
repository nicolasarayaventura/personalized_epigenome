set -x -e

###
work="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-07-08_hichiphg38"
scratch="/sc/arion/scratch/arayan01/projects/personalized_epigenome/results/2025-07-08_hichiphg38"
sampledir="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-02-19_ranjan_data/2022-11-22_data_from_Ranjan_Xiang/hichip"
ref="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/hg38_rfgen/hg38.fa"
###
jobs="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-07-08_hichiphg38/jobs/jobs.txt"
errmsg="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-07-08_hichiphg38/jobs/eo/errormsg.txt"


function run_bwa_mem {
    bsub -P acc_oscarlr -q premium -n ${cores} -W 2:00 -R "rusage[mem=8000]" -o bwa_job.txt -eo bwa_job_eo.txt \
        bash -c "bwa mem -5SP -T0 -t ${cores} ${ref} \
        ${sampledir}/HiChiP_CTCF_2M_R1.fastq.gz ${sampledir}/HiChiP_CTCF_2M_R2.fastq.gz > aligned.sam"
}

function run_pairtools_parse {
    bsub -P acc_oscarlr -q premium -n ${cores} -W 2:00 -R "rusage[mem=8000]" -o parse_job.txt -eo parse_job_eo.txt \
        bash -c "pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 \
        --nproc-in ${cores} --nproc-out ${cores} --chroms-path ${chromsizes} \
        < aligned.sam > parsed.pairs"
}

function run_pairtools_sort {
    bsub -P acc_oscarlr -q premium -n ${cores} -W 2:00 -R "rusage[mem=8000]" -o sort_job.txt -eo sort_job_eo.txt \
        bash -c "pairtools sort --tmpdir=${tmpdir} --nproc ${cores} < parsed.pairs > sorted.pairs"
}

function run_pairtools_dedup {
    bsub -P acc_oscarlr -q premium -n ${cores} -W 2:00 -R "rusage[mem=8000]" -o dedup_job.txt -eo dedup_job_eo.txt \
        bash -c "pairtools dedup --nproc-in ${cores} --nproc-out ${cores} --mark-dups \
        --output-stats ${stats} < sorted.pairs > deduped.pairs"
}

function run_pairtools_split {
    bsub -P acc_oscarlr -q premium -n ${cores} -W 2:00 -R "rusage[mem=8000]" -o split_job.txt -eo split_job_eo.txt \
        bash -c "pairtools split --nproc-in ${cores} --nproc-out ${cores} --output-pairs ${pairs} --output-sam - < deduped.pairs > split.sam"
}

function run_samtools_view {
    bsub -P acc_oscarlr -q premium -n ${cores} -W 2:00 -R "rusage[mem=8000]" -o samtools_view.txt -eo samtools_view_eo.txt \
        "samtools view -bS -@ ${cores} split.sam > split.bam"
}

function run_samtools_sort_and_index1 {
    bsub -P acc_oscarlr -q premium -n ${cores} -W 2:00 -R "rusage[mem=8000]" -o samtools_sort.txt -eo samtools_sort_eo.txt \
        "samtools sort -@ ${cores} -o ${outbam} split.bam"
}
function run_samtools_sort_and_index2 {
    bsub -P acc_oscarlr -q premium -n 1 -W 1:00 -R "rusage[mem=4000]" -o samtools_index.txt -eo samtools_index_eo.txt \
        "samtools index ${outbam}"
}
function get_qc {
    bsub -P acc_oscarlr -q premium -n 1 -W 1:00 -R "rusage[mem=4000]" -o get_qc_job.txt -eo get_qc_job_eo.txt \
        python /sc/arion/scratch/arayan01/projects/hichip_tut/data/HiChiP/get_qc.py -p /sc/arion/scratch/arayan01/projects/hichip_tut/results/mapping_stats.txt
}
function enrichment {
    bsub -P acc_oscarlr -q premium -n 4 -W 1:00 -R "rusage[mem=4000]" -o enrichment_job.txt -eo enrichment_job_eo.txt \
        bash /sc/arion/scratch/arayan01/projects/hichip_tut/data/HiChiP/enrichment_stats.sh \
        -g /sc/arion/scratch/arayan01/projects/personalized_epigenome/data/hg38_rfgen/hg38.fa \
        -b /sc/arion/scratch/arayan01/projects/hichip_tut/results/mapped.PT.bam \
        -p /sc/arion/scratch/arayan01/projects/hichip_tut/data/samples/ENCFF017XLW.bed \
        -t 16 \
        -x CTCF
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