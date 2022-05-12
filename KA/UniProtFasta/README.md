## Mappings of UniProt and COSMIC sequences

### script extractUniprotFasta.py
1. Fetch all UniProt sequences (canonical) using REST API by UniProt (see **_UniProt.fasta** files) with keyword as the gene symbol provided by COSMIC

### script doClustalo.py
2. Prepare input for clustalo by putting COSMIC and UniProt sequences in one file (see **_clustalo_input.fasta** files)
3. Run clustalo to align all the input and store the output (see **_clustalo.aln** files)
4. Map COSMIC sequences to UniProt sequences (see **_mappings.tsv** files) (remember, numbering starts from 0, and not 1)
