dir="${GROUP_SCRATCH}/uhgg/snvs/"
file=$1
while read p; do
  accession=$p
  filename="${dir}${accession}.tsv.tar.lz4"
  echo $filename;
  lz4 -dc < $filename | tar -C $dir -xvf -
done <$file
