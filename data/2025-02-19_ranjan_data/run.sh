set -x -e

function data_download {

mkdir -p /sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-02-19_ranjan_data
cp /sc/arion/projects/oscarlr/projects/IGH_epigenome/data/2025-02-11_transfer_from_UofL/data.tar.gz /sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-02-19_ranjan_data

}

function ref_genome {
	mkdir -p /sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-02-19_ranjan_data/data/2025-02-25_atac/refgenome
	ref_dir="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/2025-02-19_ranjan_data/data/2025-02-25_atac/refgenome"

	cat url_links.txt | while read ref_gen url ; do
		wget -P "${ref_dir}" "${url}"
done < url_links.txt
}

#data_download
ref_genome
