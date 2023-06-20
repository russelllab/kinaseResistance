#!/bin/bash
# Bash script to create a workflow to make HMM

# Step 1: Extract the human kinases from the fasta file
# and split them if they are more than one domains
# Save the output in the file humanKinasesHitsSplit.fasta
python find_hits.py humanKinases.fasta
cat humanKinasesPkinaseHitsSplit.fasta \
    humanKinasesPkinaseHitsSplitInactive.fasta \
    humanKinasesPK_Tyr_Ser-ThrHitsSplit.fasta \
    humanKinasesPK_Tyr_Ser-ThrHitsSplitInactive.fasta > \
    humanKinasesPkinasePK_Tyr_Ser-ThrAll.fasta
cat humanKinasesPkinaseHitsSplitInactive.fasta \
    humanKinasesPK_Tyr_Ser-ThrHitsSplitInactive.fasta > \
    humanKinasesPkinasePK_Tyr_Ser-ThrAllInactive.fasta


# Step 2: Run the hmmalign to align the human sequences
# with the Pkinase & PK-Tyr HMM profile
# Save the output in the file humanKinasesHitsSplit.aln
# hmmalign -o humanKinasesHitsSplit.sto ../pfam/Pkinase.hmm humanKinasesHitsSplit.fasta
hmmalign -o humanKinasesPkinaseHitsSplit.sto ../pfam/Pkinase.hmm humanKinasesPkinaseHitsSplit.fasta
hmmalign -o humanKinasesPK_Tyr_Ser-ThrHitsSplit.sto ../pfam/Pkinase.hmm humanKinasesPK_Tyr_Ser-ThrHitsSplit.fasta

# python prepareAlignments.py -i humanKinasesPkinaseHitsSplit.sto > ali1.fasta
# python prepareAlignments.py -i humanKinasesPK_Tyr_Ser-ThrHitsSplit.sto > ali2.fasta

# Step 3: Convert the stockholm file to CLUSTAL format
# Save the output in the file humanKinasesHitsSplit.aln
# sed '/^#/s/.*/CLUSTAL/' humanKinasesHitsSplit.sto > humanKinasesHitsSplit.aln
sed '/^# STOCKHOLM/s/.*/CLUSTAL/' humanKinasesPkinaseHitsSplit.sto | grep -v '#' > humanKinasesPkinaseHitsSplit.aln
sed '/^# STOCKHOLM/s/.*/CLUSTAL/' humanKinasesPK_Tyr_Ser-ThrHitsSplit.sto | grep -v '#' > humanKinasesPK_Tyr_Ser-ThrHitsSplit.aln

# Step 4: Move the .aln file to the alignments folder
# mv humanKinasesHitsSplit.aln ../alignments/
mv humanKinasesPkinaseHitsSplit.aln ../alignments/
mv humanKinasesPK_Tyr_Ser-ThrHitsSplit.aln ../alignments/
cd ../alignments/
# exit 0

# Step 5: Run the trim_alignment_split_script.py to trim the alignment
# Save the output in the file humanKinasesHitsSplitTrimmed.aln/fasta
# Also save the web vsersion humanKinasesHitsSplitTrimmedWeb.aln
# The script automatically adds the keyword "trimmed" to the file name
# and saves in both fasta and aln format
# the script also saves the jalview annoptation file
# jalview_annotations.txt
# python trim_alignment_split_outside.py humanKinasesHitsSplit.aln ../data/humanKinases.fasta 32177 32940 30
python trim_alignment_split_outside.py humanKinasesPkinaseHitsSplit.aln ../data/humanKinases.fasta 32178 32816 30
python trim_alignment_split_outside.py humanKinasesPK_Tyr_Ser-ThrHitsSplit.aln ../data/humanKinases.fasta 1950 2546 30

# Step 6: Merge the two alignments
cat humanKinasesPkinaseHitsSplitTrimmed.fasta humanKinasesPK_Tyr_Ser-ThrHitsSplitTrimmed.fasta > ali.fasta
ruby ../../../makemergetable.rb humanKinasesPkinaseHitsSplitTrimmed.fasta humanKinasesPK_Tyr_Ser-ThrHitsSplitTrimmed.fasta > subMSAtable
mafft --merge subMSAtable ali.fasta > humanKinasesHitsSplitTrimmed.fasta

# Step 7: Build a hidden Markov model (HMM) from the alignment
# Use the fasta file as input because the aln file fails (!!!!)
# Save the output in the file humanKinasesHitsSplitTrimmed.hmm
# in the pfam folder
rm -rf ../pfam/humanKinasesHitsSplitTrimmed.hmm*
hmmbuild --symfrac 0.0 ../pfam/humanKinasesHitsSplitTrimmed.hmm humanKinasesHitsSplitTrimmed.fasta
# hmmbuild ../pfam/humanKinasesHitsSplitTrimmed.hmm humanKinasesHitsSplitTrimmed.fasta
hmmpress ../pfam/humanKinasesHitsSplitTrimmed.hmm

# Step 8: Do hmmsearch with the HMM profile against
# the full kinase sequences. Save the output in the file
# humanKinasesHitsSplitHmmsearchTrimmed.txt.gz (gzip compressed)
cd ../data/
hmmsearch -o humanKinasesHitsSplitHmmsearchTrimmed.txt ../pfam/humanKinasesHitsSplitTrimmed.hmm humanKinases.fasta
gzip -f humanKinasesHitsSplitHmmsearchTrimmed.txt

hmmsearch -o humanKinasesPkinasePK_Tyr_Ser-ThrAllInactiveHmmsearch.txt ../pfam/humanKinasesHitsSplitTrimmed.hmm humanKinasesPkinasePK_Tyr_Ser-ThrAllInactive.fasta
gzip -f humanKinasesPkinasePK_Tyr_Ser-ThrAllInactiveHmmsearch.txt

# Step 9: Map the UniProt AA to the new domain numbering in the
# domain name 'humanKinasesHitsSplitTrimmed' (given as input
# to the script) using the file humanKinasesHitsSplitHmmsearchTrimmed.txt.gz
# Save the output in the file humanKinasesHitsSplitHmmsearchTrimmedMapped.txt.gz
# The scipt automatically adds the keyword "Mappings" to the input file name
# python map2domain.py humanKinasesHitsSplitHmmsearchTrimmed.txt.gz humanKinasesHitsSplitTrimmed 30 793
# python map2domain.py humanKinasesHitsSplitHmmsearchTrimmed.txt.gz ../alignments/humanKinasesHitsSplitTrimmed.fasta humanKinasesHitsSplitTrimmed 30 794
python map2aln.py ../alignments/humanKinasesHitsSplitTrimmed.fasta humanKinasesHitsSplitTrimmed

# Step 10: Map the PTM sites to the new domain numbering in the
# domain name 'humanKinasesHitsSplitTrimmed' (given as input
# to the script) using the file humanKinasesHitsSplitHmmsearchTrimmed.txt.gz
# Save the output in the file humanKinasesHitsSplitHmmsearchTrimmedPTM.tsv
# python3 make_table_psp_kinase_trimmed.py humanKinasesHitsSplitTrimmed humanKinasesHitsSplitHmmsearchTrimmedMappings.tsv.gz humanKinases.fasta humanKinasesHitsSplitTrimmedPTM.tsv
python3 make_table_psp_kinase_trimmed.py humanKinasesHitsSplitTrimmed humanKinasesHitsSplitTrimmedMappings.tsv.gz humanKinases.fasta humanKinasesHitsSplitTrimmedPTM.tsv

# Step 11: Align all kinases (including inactive and pseudokinases) with the new HMM profile
# Save the output in the file humanKinasesPkinasePK_Tyr_Ser-ThrAll.sto
# Convert the stockholm file to CLUSTAL format
hmmalign -o humanKinasesPkinasePK_Tyr_Ser-ThrAll.sto ../pfam/humanKinasesHitsSplitTrimmed.hmm humanKinasesPkinasePK_Tyr_Ser-ThrAll.fasta
sed '/^# STOCKHOLM/s/.*/CLUSTAL/' humanKinasesPkinasePK_Tyr_Ser-ThrAll.sto | grep -v '#' > humanKinasesPkinasePK_Tyr_Ser-ThrAll.aln
mv humanKinasesPkinasePK_Tyr_Ser-ThrAll.aln ../alignments/
cd ../alignments/
python trim_alignment_split_outside.py \
        humanKinasesPkinasePK_Tyr_Ser-ThrAll.aln \
        ../data/humanKinases.fasta \
        32155 \
        33025 \
        30