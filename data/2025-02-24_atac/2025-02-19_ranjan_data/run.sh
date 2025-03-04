set -x -e

scratch_hg38="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-02-24_atachg38"
scratch_hg19="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-03-04_atachg19"


function data_download {
	mkdir -p /sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-02-19_ranjan_data
	cp /sc/arion/projects/oscarlr/projects/IGH_epigenome/data/2025-02-11_transfer_from_UofL/data.tar.gz /sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-02-19_ranjan_data

}

function ref_genome_hg38 {
rm -rf "${scratch_hg38}/refgenome"
	mkdir -p ${scratch_hg38}/refgenome
	ref_dir="${scratch_hg38}/refgenome"

   cat url_links.txt | while read ref_gen38 url1 refgen19 url2 ; do
            wget -P "${refgen38}" "${url1}"
	done < url_links.txt
}

function ref_genome_hg19 {
rm -rf "${scratch_hg19}/refgenome"
	mkdir -p ${scratch_hg19}/refgenome
	ref_dir="${scratch_hg19}/refgenome"

   cat url_links.txt | while read ref_gen38 url1 refgen19 url2 ; do
            wget -P "${ref_dir}" "${url2}"
done < url_links.txt
}

#data_download
#ref_genome_hg38
ref_genome_hg19