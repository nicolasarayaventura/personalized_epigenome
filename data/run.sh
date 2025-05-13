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

function ref_genome_hg38 {
    rm -rf "${scratch_base}/hg38_rfgen"
            mkdir -p "${scratch_base}/hg38_rfgen"
            ref_dir="${scratch_base}/hg38_rfgen"
    cat url_links.txt | while read ref_gen38 url1 refgen19 url2 ; do
                wget -P "${ref_dir}/${refgen38}" "${url1}"
            done < url_links.txt
}
function hg38_unzip {
    bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=16000]" -o "/sc/arion/work/arayan01/project/personalized_epigenome/data/datadownload.txt" \
        gunzip "${scratch_base}/hg38_rfgen/${refgen38}/*.gz"
}
function filterhg38 {
    gen_dir="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/hg38_rfgen"
    samtools faidx "${gen_dir}/hg38.fa"
    grep -vE "Un|random|alt" "${gen_dir}/hg38.fa.fai" > "${gen_dir}/hg38_filtered.fa.fai"
    cut -f1 "${gen_dir}/hg38_filtered.fa.fai" > "${gen_dir}/filtered_chromosomes.txt"

    rm -f "${gen_dir}/hg38_filtered.fa"
    while read chr; do
        samtools faidx "${gen_dir}/hg38.fa" "$chr" >> "${gen_dir}/hg38_filtered.fa"
    done < "${gen_dir}/filtered_chromosomes.txt"

    bsub -P acc_oscarlr -q premium -W 24:00 -R "rusage[mem=16000]" \
    -o "${work}/job_hg38filterd_indexing.txt" -eo "${work}/job_hg38filterd_indexing_err.txt" \
    "bwa index -p ${gen_dir}/hg38_filtered ${gen_dir}/hg38_filtered.fa"
}

function bwt_hg38 {
    bwtdir="${scratch_base}/bwt"
    genome="${scratch_base}/hg38_rfgen/hg38_filtered.fa"
    out="${bwtdir}/hg38_filtered"

    bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=16000]" -o "${work}/job_hg38bwt.txt" \
        bowtie2-build ${genome} ${out}

}
function ref_genome_hg19 {
    rm -rf "${scratch_base}/hg19_rfgen"
            mkdir -p "${scratch_base}/hg19_rfgen"
            ref_dir="${scratch_base}/hg19_rfgen"

    cat url_links.txt | while read ref_gen38 url1 refgen19 url2 ; do
                wget -P "${ref_dir}" "${url2}"
    done < url_links.txt
}

#data_download
#uncompress
#ref_genome_hg38
#hg38_unzip
filterhg38
#bwt_hg38