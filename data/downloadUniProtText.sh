#!/usr/bin/bash

dir=$1
for file in $dir/*.txt.gz ; do
    filename=$(basename "$file")
    # echo $filename;
    acc="${filename##*/}"
    acc="${acc%.txt.gz}"
    # echo $acc;
    # python3 fetchVariants.py $file
    url="https://rest.uniprot.org/uniprotkb/${acc}.txt"
    wget "$url" -O $dir/$acc.txt
    echo $url
    # os.system('gzip ' + fasta_path + acc+'.txt ')
done