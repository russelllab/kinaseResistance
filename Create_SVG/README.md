# Introduction (small one)

## This script
The file **create_svg_20230210.py** comes with 3 python functions and a dictionary with annotations (Line 229 ff. taken from CONNECTORα).
You will also have to use the second file, which contains an alignment, as input for the script.

`# change the following line (223) to the location of the downloaded alignmentfile.`

`alignmentfile = "/INSERT_PATH_HERE/clustalo-I20230119-084041-0905-58774427-p1m.clustal_num"`	

## General use case
In general, this script needs two things: 
* An alignment in clustal format 
* a python dictionary formatted as {Protein:{Feature:\[Residues]}} 

to work. The dictionary can be created elsewhere and could contain different features than the one I included here, so it is **versatile**.

## Example
The same functions are part of my CONNECTOR workflow, presented in the lab meeting on Thursday 09. Feb. 2023, where the features are read through the Uniprot API and the alignment is read in from a file & parsed via Biopython. Both are then handed over to the functions included in this script.


## Required libraries/software

Python 2.7. (sorry folks)

import svgwrite

from Bio import AlignIO


## Features (latest Version)

- Conservation: Gives a black rectangle as an indicator of sequence identity (top).

- Residue Numbering: Gives the residue number (every 10 residues), based on sequence of interest, and highlights the input residue in red.

- Feature Annotation: Gives a colored background based on type of annotation (taken from uniprot) to the respective residue.

- Sorted Sequences: Sequences with fewer uniprot annotations are sorted to the bottom of the alignment.

- GAPs removed: Gaps are printed with white color (i.e. invisible on a white background). Additionally, columns with more than 90 % GAPs are removed from the alignment. Sequences affected by this (i.e. the up to 10 % of sequences that did not have a gap at that position) **are kept and not removed**. 



<img src="https://github.com/russelllab/kinaseResistance/blob/bb35c9fb85daf276d0cffc44edae6f7622b676ca/Create_SVG/Version2/sequence.svg?sanitize=true">

