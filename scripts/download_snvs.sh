if [[ $# -lt 1 ]]; then
  echo "Please provide file with mgnify accessions"
  exit 1
fi

uhgg_ftp="http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/"
file=$1
while read p; do
  accession=$p
  genome_src="${uhgg_ftp}uhgg_catalogue/${accession::-2}/$accession/genome/"
  snv_src="${uhgg_ftp}snv_catalogue/${accession}.tsv.tar.lz4"
  save_path="${GROUP_SCRATCH}/uhgg/snvs/"
  file="${save_path}/${accession}.tsv.tar.lz4"
  if [[ -f "$save_path" ]]; then
    echo "$file exists."
    continue
  fi
  #echo $snv_src
  wget -P $save_path $snv_src
done <$file
