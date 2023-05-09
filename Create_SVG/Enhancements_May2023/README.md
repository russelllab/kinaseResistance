To use the script we execute the following command:
`python3 create_svg_20230509_kinases_GS.py P46734 84 9 5 humanKinasesTrimmed.clustal Mutational_Infofile_Kinases_V2.txt GenerelleKonservierung_May-09-2023.txt none SeqIdentity_Matrix_May-09-2023.txt 2` 
This command has several fields after calling the script:

| Field        | Example           | Description  |
| ------------- |:-------------:| -----:|
|0     | P46734 | The uniprot ID of the protein we are interested in |
|1     | 84 | The position to be highlighted |
|2     | 9 | The Windowsize, we show +/- the windowsize of residues around the highlighted position. **WARNING** Does not work if the starting position of the sequences is != 1.|
|3      |5 | Shows the top 10 sequences, based on information content|
|4     | humanKinasesTrimmed.clustal | The alignment file |
|5     | Mutational_Infofile_Kinases_V2.txt | The file containing positional information |
|6      |GenerelleKonservierung_May-09-2023.txt|A dictionary containing the 3 most frequent characters for each alignment position, and uniprotIDs that have functional info at that position|
|7     | none | A file containing structural/domain features, numbering based on **protein of interest** Hasn't been used in a while, might break |
|8     |SeqIdentity_Matrix_May-09-2023.txt|A ditionary containing a sequence identity matrix based on the provided alignment|
|9     |2|Sorting value. **1** Sorts by available functional information based on windowsize and shown sequences, **2** sorts by sequence identity instead.




## Prep Work

The dictionaries **GenerelleKonservierung_May-09-2023.txt** and **SeqIdentity_Matrix_May-09-2023.txt** can be created using the script 'Precalculate.py'.
 
Use the script like this:

'python3 Precalculate.py humanKinasesTrimmed.clustal Mutational_Infofile_Kinases_V2.txt 2 P15056 600 20'.
The last 4 fields ("2 P15056 600 20") are only for debugging purposes. Simply provide the alignment and the information and the aforementioned dictionaries will be written to disk.
