#!/bin/bash
set -e -x

function samplelists {
    output="${work}/sample_ids.txt"
    for sample in SEM RS411
    do
        data_path=${sampledir1}
        if [ ${sample} == "RS411" ]
        then
            data_path=${sampledir2}
        fi
        for tf in CTCF-1 CTCF-2 RAD21-1 RAD21-2 K27ac-1 K27ac-2
        do
            if [ ${sample} == "SEM" ]
            then
                if [ ${tf} == "K27ac-1" ]
                then 
                      data_path=${sampledir2}
                fi
                if [ ${tf} == "K27ac-1=2" ]
                then 
                      data_path=${sampledir2}
                fi
            r1=`ls ${data_path}/${sample}-${tf}_*R1*.fastq.gz`
            r2=`ls ${data_path}/${sample}-${tf}_*R2*.fastq.gz`
            echo "${sample}\t${tf}\t${r1}\t${r2}"
        done
    done > "$output"
}


function bmem {
    rm -rf "${scratch}/bmem"

    outdir="${scratch}/bmem"  
    mkdir -p "${outdir}"       

    cat "${work}/sample_ids.txt" | while read sample tf r1 r2
    do
          bsub -P acc_oscarlr -q premium -n ${cores} -W 4:00 -R "rusage[mem=16000]" \
                -o "${jobs}/bmem_${base}.log" -eo "${errmsg}/bmem_${base}.log" \
                bash -c "bwa mem -5SP -T0 -t ${cores} ${ref} ${r1} ${r2} > ${outdir}/${sample}_${tf}_aligned.sam"
    done < "$samplelist"
}

