set -x -e

scratch="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-04-07_chipsseqhg38"
work="/sc/arion/work/arayan01/project/personalized_epigenome/results/2025-04-07_chipsseqhg38"

sampledir="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-02-19_ranjan_data/data/2022-11-22_data_from_Ranjan_Xiang/chip/RSen042022_24"
samplelist="${work}/sample_ids.txt"

function samplelist {
    rm -rf "${work}/sample_ids.txt"
    mkdir -p "${work}"
    
    output="${work}/sample_ids.txt"

    
    for filepath in "${sampledir}"/*; do
        [ -f "$filepath" ] || continue
        basename=$(basename "$filepath")
        echo -e "${basename}\t${filepath}" >> "${output}"
    done
}

function fastqc_initial {
    rm -rf "${work}/qc/init_fastqc"
    mkdir "${work}/qc/init_fastqc"
    
	qc_out="${work}/qc/init_fastqc" #once job is finished run multiqc in directory!

	while IFS=$'\t' read -r base path; do
		bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "${work}/job_fastqc.txt" \
		    fastqc --noextract -t 2 "${path}" -o "${qc_out}"
	done < "${samplelist}"
}

function trimming {
    rm -rf "${scratch}/trimming"
    mkdir -p "${scratch}/trimming"

   
    out="${scratch}/trimming"

    while IFS=$'\t' read -r base path; do
        bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "${work}/trimming_job.txt" \
         trim_galore --quality 15 --stringency 3 --report ${path} -o "${out}/${base}_trimmed"
    done < "${samplelist}"
}

#samplelist
#fastqc_initial
trimming