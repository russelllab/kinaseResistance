#!/bin/bash

## bash script to run calculate_log_scores.py

for file in ../AK_hmm_models/*
do
	IFS='/'
	read -a fileName <<< "$file"
	IFS='.'
	read -a name <<< "${fileName[2]}"
	# echo "${name[0]}"
	python3 calculate_log_scores.py "$file" ../hmm_build/"$name".hmm.gz > ../log_odds_scores/"$name".scores.txt
	# exit 1
done
