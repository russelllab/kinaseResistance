#!/bin/bash

## bash script to run calculate_log_scores.py

for file in ../kin_seq/*
do
	if [ "${file: -6}" == ".fa.gz" ]
	then
		# echo "$file"
		IFS='/'
		read -a fileName <<< "$file"
		IFS='.'
		read -a name <<< "${fileName[2]}"
		# echo "${name[0]}"
		python3 make_aln.py "${name[0]}" > ../logs/"${name[0]}".log
		# exit 1
	fi
done
