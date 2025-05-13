set -x -e

scratch_base="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data"

scratch_atac_hg38="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-02-24_atachg38"
scratch_atac_hg19="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-03-04_atachg19"

scratch_chip_hg38="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-04-07_chipsseqhg38"
scratch_chip_hg19="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-04-07_chipsseqhg19"

function data_download {
    mkdir -p /sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-02-19_ranjan_data
    bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=16000]" -o "/sc/arion/work/arayan01/project/personalized_epigenome/data/datadownload.txt" \
        cp /sc/arion/projects/oscarlr/projects/IGH_epigenome/data/2025-02-11_transfer_from_UofL/data.tar.gz \
           /sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-02-19_ranjan_data
}

function uncompress {
        bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=16000]" -o "/sc/arion/work/arayan01/project/personalized_epigenome/data/datadownload.txt" \
                tar -xzvf /sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-02-19_ranjan_data/data.tar.gz

}

function ref_genomes {
    rm -rf "${scratch_base}/hg38_rfgen"
    rm -rf "${scratch_base}/hg19_rfgen"

    mkdir -p "${scratch_base}/hg38_rfgen"
    mkdir -p "${scratch_base}/hg19_rfgen"

    ref_dir1="${scratch_base}/hg38_rfgen"
    ref_dir2="${scratch_base}/hg19_rfgen"

    cat url_links.txt | while read refgen38 url1 refgen19 url2; do
        wget -P "${ref_dir1}" "${url1}"
        wget -P "${ref_dir2}" "${url2}"
    done < url_links.txt
}
function genome_unzip {
    bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=16000]" -o "/sc/arion/work/arayan01/project/personalized_epigenome/data/datadownload.txt" \
        "gunzip ${scratch_base}/hg38_rfgen/${refgen38}/*.gz" 

    bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=16000]" -o "/sc/arion/work/arayan01/project/personalized_epigenome/data/datadownload.txt" \
        "gunzip ${scratch_base}/hg19_rfgen/${refgen19}/*.gz"
}

function genomefilter {
#hg38
    hg38gen_dir="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/hg38_rfgen"
    samtools faidx "${hg38gen_dir}/hg38.fa"
    grep -vE "Un|random|alt" "${hg38gen_dir}/hg38.fa.fai" > "${hg38gen_dir}/hg38_filtered.fa.fai"
    cut -f1 "${hg38gen_dir}/hg38_filtered.fa.fai" > "${hg38gen_dir}/filtered_chromosomes.txt"

    rm -f "${hg38gen_dir}/hg38_filtered.fa"
    while read chr; do
        samtools faidx "${hg38gen_dir}/hg38.fa" "$chr" >> "${hg38gen_dir}/hg38_filtered.fa"
    done < "${hg38gen_dir}/filtered_chromosomes.txt"

#hg19
    hg19gen_dir="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/hg19_rfgen"
    samtools faidx "${hg19gen_dir}/hg19.fa"
    grep -vE "Un|random|alt|_" "${hg19gen_dir}/hg19.fa.fai" > "${hg19gen_dir}/hg19_filtered.fa.fai"
    cut -f1 "${hg19gen_dir}/hg19_filtered.fa.fai" > "${hg19gen_dir}/filtered_chromosomes.txt"

    rm -f "${hg19gen_dir}/hg19_filtered.fa"
    while read chr; do
        samtools faidx "${hg19gen_dir}/hg19.fa" "$chr" >> "${hg19gen_dir}/hg19_filtered.fa"
    done < "${hg19gen_dir}/filtered_chromosomes.txt"

#job submission
    bsub -P acc_oscarlr -q premium -W 24:00 -R "rusage[mem=16000]" \
        -o "job_hg38filterd_indexing.txt" -eo "job_hg38filterd_indexing_err.txt" \
        "bwa index -p ${hg38gen_dir}/hg38_filtered ${hg38gen_dir}/hg38_filtered.fa"

    bsub -P acc_oscarlr -q premium -W 24:00 -R "rusage[mem=16000]" \
        -o "job_hg19filterd_indexing.txt" -eo "job_hg19filterd_indexing_err.txt" \
        "bwa index -p ${hg19gen_dir}/hg19_filtered ${hg19gen_dir}/hg19_filtered.fa"
}

function bwt_genome {
    hg38_bwtdir="${scratch_base}/hg38_rfgen/bwt"
    hg38_genome="${scratch_base}/hg38_rfgen/hg38_filtered.fa"
    hg38_out="${bwtdir}/hg38_filtered"


    hg19_bwtdir="${scratch_base}/hg19_rfgen/bwt"
    hg19_genome="${scratch_base}/hg19_rfgen/hg19_filtered.fa"
    hg19_out="${bwtdir}/hg19_filtered"
    
    bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=16000]" -o "job_hg38bwt.txt" \
        bowtie2-build ${hg38_genome} ${hg38_out}

    bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=16000]" -o "job_hg19bwt.txt" \
        bowtie2-build ${hg19_genome} ${hg19_out}

}

# Main script execution

#data_download
#uncompress
#ref_genomes
#genome_unzip
genomefilter
#bwt_genome