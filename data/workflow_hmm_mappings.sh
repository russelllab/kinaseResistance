#!/bin/bash
# Bash script to create a workflow to make HMM

# Step 1: Extract the human kinases from the fasta file
# and split them if they are more than one domains
# Save the output in the file humanKinasesHitsSplit.fasta
python find_hits.py humanKinases.fasta

# Step 2: Run the hmmalign to align the human sequences
# with the Pkinase HMM profile
# Save the output in the file humanKinasesHitsSplit.aln
hmmalign -o humanKinasesHitsSplit.sto ../pfam/Pkinase.hmm humanKinasesHitsSplit.fasta

# Step 3: Convert the stockholm file to CLUSTAL format
# Save the output in the file humanKinasesHitsSplit.aln
sed '/^#/s/.*/CLUSTAL/' humanKinasesHitsSplit.sto > humanKinasesHitsSplit.aln

# Step 4: Move the .aln file to the alignments folder
mv humanKinasesHitsSplit.aln ../alignments/
cd ../alignments/

# Step 5: Run the trim_alignment_split_script.py to trim the alignment
# Save the output in the file humanKinasesHitsSplitTrimmed.aln/fasta
# Also save the web vsersion humanKinasesHitsSplitTrimmedWeb.aln
# The script automatically adds the keyword "trimmed" to the file name
# and saves in both fasta and aln format
# the script also saves the jalview annoptation file
# jalview_annotations.txt
python trim_alignment_split_outside.py humanKinasesHitsSplit.aln ../data/humanKinases.fasta 32177 32940 30

# Step 6: Build a hidden Markov model (HMM) from the alignment
# Use the fasta file as input because the aln file fails (!!!!)
# Save the output in the file humanKinasesHitsSplitTrimmed.hmm
# in the pfam folder
rm -rf ../pfam/humanKinasesHitsSplitTrimmed.hmm*
hmmbuild --symfrac 0.0 ../pfam/humanKinasesHitsSplitTrimmed.hmm humanKinasesHitsSplitTrimmed.fasta
# hmmbuild ../pfam/humanKinasesHitsSplitTrimmed.hmm humanKinasesHitsSplitTrimmed.fasta
hmmpress ../pfam/humanKinasesHitsSplitTrimmed.hmm

# Step 7: Do hmmsearch with the HMM profile against
# the full kinase sequences. Save the output in the file
# humanKinasesHitsSplitHmmsearchTrimmed.txt.gz (gzip compressed)
cd ../data/
hmmsearch -o humanKinasesHitsSplitHmmsearchTrimmed.txt ../pfam/humanKinasesHitsSplitTrimmed.hmm humanKinases.fasta
gzip -f humanKinasesHitsSplitHmmsearchTrimmed.txt

# Step 8: Map the UniProt AA to the new domain numbering in the
# domain name 'humanKinasesHitsSplitTrimmed' (given as input
# to the script) using the file humanKinasesHitsSplitHmmsearchTrimmed.txt.gz
# Save the output in the file humanKinasesHitsSplitHmmsearchTrimmedMapped.txt.gz
# The scipt automatically adds the keyword "Mappings" to the input file name
# python map2domain.py humanKinasesHitsSplitHmmsearchTrimmed.txt.gz humanKinasesHitsSplitTrimmed 30 793
# python map2domain.py humanKinasesHitsSplitHmmsearchTrimmed.txt.gz ../alignments/humanKinasesHitsSplitTrimmed.fasta humanKinasesHitsSplitTrimmed 30 794
python map2aln.py ../alignments/humanKinasesHitsSplitTrimmed.fasta humanKinasesHitsSplitTrimmed

# Step 9: Map the PTM sites to the new domain numbering in the
# domain name 'humanKinasesHitsSplitTrimmed' (given as input
# to the script) using the file humanKinasesHitsSplitHmmsearchTrimmed.txt.gz
# Save the output in the file humanKinasesHitsSplitHmmsearchTrimmedPTM.tsv
# python3 make_table_psp_kinase_trimmed.py humanKinasesHitsSplitTrimmed humanKinasesHitsSplitHmmsearchTrimmedMappings.tsv.gz humanKinases.fasta humanKinasesHitsSplitTrimmedPTM.tsv
python3 make_table_psp_kinase_trimmed.py humanKinasesHitsSplitTrimmed humanKinasesHitsSplitTrimmedMappings.tsv.gz humanKinases.fasta humanKinasesHitsSplitTrimmedPTM.tsv