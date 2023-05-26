#!/usr/bin/bash

dir=$1
for file in $dir/*.txt.gz ; do
    # echo $file;
    python3 fetchVariants.py $file
done