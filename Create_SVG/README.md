___
# Introduction
## This script
The script **Enhancements_May2023/27June/create_svg_20230706_kinases_GS.py** comes with several python3 functions that can annotate or highlight positions in a Clustal-formatted alignment.

___
## General use case
In general, this script needs two things: 
* An alignment in clustal format 
* a python dictionary formatted as {Protein:{Feature:\[Residues]}}
* an optional protein feature dictionary, formatted as (and based on protein of interest) {Feature:[Startposition, Endposition]}

to work. The dictionary can be created elsewhere and could contain different features than the one I included here, so it is **versatile**.

___
## Required libraries/software

Python 3.6+

import svgwrite

import ast

import sys


___
## Features
- tooltip functionality to the svgs with hyperlinks to sources (Pubmed, COSMIC, Phosphositeplus etc.)
- option to use a parameter to show only the Top X sequences (+ sequence of interest)
- Upgraded the code to work under Python3.6
- Using transparent rectangles to highlight a sequence conservation (= identity) over >= 70 %, based on the sequence of interest. The colors for this are taken from CLUSTAL/Jalview
- Showing basic heatmapping above the alignment, indicating how many highlights per position & per category we have
- start and end positions for each displayed sequence
- Command line functionality

___
## Using the script
# Prep Work
The dictionaries **GenerelleKonservierung_MONTH-DAY-YEAR.txt** and **SeqIdentity_Matrix_MONTH-DAY-YEAR.txt** can be created using the script **Enhancements_May2023/27June/Precalculate_20230627.py**.

Use the script like this:

`python3 Precalculate_20230627.py Clustal_Alignment Mutation_Infofile 2 P15056 600 20`.

The last 4 fields ("2 P15056 600 20") are only for debugging purposes. Simply provide the alignment and the information and the aforementioned dictionaries will be written to disk.

# Main Script
To use the script we can now execute the following command:

`python3 create_svg_20230706_kinases_GS.py PIM1=P11309=1-313 97 10 20 humanKinasesPkinasePK_Tyr_Ser-ThrAll_no_gapsTrimmedWeb.aln sample_dic_mutation_info.txt GenerelleKonservierung_Jun-28-2023.txt ss.tsv SeqIdentity_Matrix_Jun-28-2023.txt 1`   

This command has several fields after calling the script:

| Field        | Example           | Description  |
| ------------- |:-------------:| -----:|
| 0     | PIM1=P11309=1-313 | The sequence identifier. In this case "=" will be replaced with "\|" later |
| 1     | 97 | The position to be highlighted |
| 2     | 10 | The Windowsize, we show +/- the windowsize of residues around the highlighted position|
| 3     |20 | Shows the top 20 sequences, based on information content|
| 4     | humanKinasesPkinasePK_Tyr_Ser-ThrAll_no_gapsTrimmedWeb.a | The alignment file |
| 5     | sample_dic_mutation_info.txt | The file containing positional information |
| 6     | GenerelleKonservierung_Jun-28-2023.txt | A file containing conservation values for each alignment position. We use this to remove columns with mostly gaps |
| 7     |ss.tsv|A file containing structural features to show above the alignment|
| 8     |SeqIdentity_Matrix_Jun-28-2023.txt|A file similar to field number 6|
| 9     |1|Determines the way of filtering that is used|

The filtering values are:

| Value        | Filter Mechanism | 
| ------------- |:-------------:| 
| 1     | Number of functional infos per kinase| 
| 2     | Overall sequence similarity | 
| 3     | Activating | 
| 4     | Deactivating| 
| 5     | Resistance | 
| 6     | Phosphorylation | 
| 7     | Acetylation | 
| 8     | Ubiquitination | 
| 9     | Sumoylation | 
| 10     | O-GlcNAc | 
| 11     | Methylation | 
| 12     | Activating | 

Information of the given type (1,3-12) is weighted stronger the closer it is mapped to the position of interest for the queried protein.

**Note**: The script allows for a little hack here. If you want a (large) .svg containing the whole alignment just give a big number in field 2, for example 20000. The script will then produce a complete alignment view.

**Note**: The script can also be imported from within python via 'import create_svg_20230706_kinases_GS.py as create_svg'

___
## The most recent type of results
MAP2K3, position 84, with tooltips and showing MAP2K3 + top30 sequences sorted by known "Activating" variants. Open in new tab to use tooltips.
<img src="https://github.com/russelllab/kinaseResistance/blob/9ef9ee144c2e443a89dd56c4ab4f5e8033c5bbec/Create_SVG/Enhancements_May2023/27June/MAP2K3%7CP46734%7C28-347_Position84_Windowsize30_Topguns30_Sorting_3.svg?sanitize=true">

___
## FAQs
# Is there a more general version of this script available?
Yes, there will be a more general version of this alignment script available together with automated data collection for any protein in uniprot. For more information contact torsten.schmenger(at)bioquant.uni-heidelberg.de.

___
## Contact
Alignment: torsten.schmenger(at)bioquant.uni-heidelberg.de

Website: gurdeep.singh(at)bioquant.uni-heidelberg.de
