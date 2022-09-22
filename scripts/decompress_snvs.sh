dir="${GROUP_SCRATCH}/uhgg/snvs/"
for file in ${dir}*.tsv.tar.lz4;
do
  echo $file;
  lz4 -dc < $file | tar -C $dir -xvf -
done;
