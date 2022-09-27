if [[ $# -lt 1 ]]; then
  echo "Please provide file with mgnify accessions"
  exit 1
fi

uhgg_ftp="http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/"
file=$1
while read p; do
  accession=$p
  genome_src="${uhgg_ftp}uhgg_catalogue/${accession::-2}/$accession/genome/"
  gene_dir="${GROUP_HOME}/uhgg/genes/"
  genome_dir="${GROUP_HOME}/uhgg/reference_genomes/"

  wget -P $gene_dir "${genome_src}${accession}.gff"

  wget -P $genome_dir "${genome_src}${accession}.fna"
done <$file
