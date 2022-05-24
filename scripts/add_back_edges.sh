#!/bin/bash

zip_directory=$1
output_file=$2
executable=$3

for file in ${zip_directory}/*;
do
  echo "start processing file ${file}"
  tmp_file="${zip_directory}/tmp"
  gzip -dkc $file > ${tmp_file}
  $executable --infile $tmp_file --outfile $output_file --max_weight 254 --num_vid_bytes 4
  rm ${tmp_file}
done

