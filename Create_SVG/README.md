___
# Introduction
## This script
The file **create_svg_20230228.py** comes with 3 python functions and two dictionaries with annotations (Line 229 ff. taken from CONNECTORα and Uniprot/PFAM/Interpro).
You will also have to use the second file, which contains an alignment, as input for the script.

`# change the following line (247) to the location of the downloaded alignmentfile.`

`alignmentfile = "/INSERT_PATH_HERE/clustalo-I20230119-084041-0905-58774427-p1m.clustal_num"`	
___
## General use case
In general, this script needs two things: 
* An alignment in clustal format 
* a python dictionary formatted as {Protein:{Feature:\[Residues]}} 
* an optional protein feature dictionary, formatted as (and based on protein of interest) {Feature:[Startposition, Endposition]}

to work. The dictionary can be created elsewhere and could contain different features than the one I included here, so it is **versatile**.

___
## Example
The same functions are part of my CONNECTOR workflow, presented in the lab meeting on Thursday 09. Feb. 2023.
I read the features for the protein of interest and homologs (from the alignment) directly through from uniprot (for example: https://rest.uniprot.org/uniprotkb/P04637.txt).
The alignment is read in from a file & parsed via Biopython. Both are then handed over to the functions included in this script.

___
## Required libraries/software

Python 2.7. (sorry folks)

import svgwrite

from Bio import AlignIO

**Note**: Libraries are also available for >= Python 3.6, print statements and dictionary handling could be converted 

___
## Features (Version 3)

- Conservation: Gives a black rectangle as an indicator of sequence identity (top).

- Residue Numbering: Gives the residue number (every 10 residues), based on sequence of interest, and highlights the input residue in red.

- Feature Annotation: Gives a colored background based on type of annotation (taken from uniprot) to the respective residue.

- Sorted Sequences: Sequences with fewer uniprot annotations are sorted to the bottom of the alignment.

- GAPs removed: Gaps are printed with white color (i.e. invisible on a white background). Additionally, columns with more than 90 % GAPs are removed from the alignment. Sequences affected by this (i.e. the up to 10 % of sequences that did not have a gap at that position) **are kept and not removed**. 

- **NEW** Highlighting protein features, *here* for example p-loop, Switch I and the Effector region of RHOA. We currently support the displaying of up to 9 features (dependent on the given colors in *featurecolors* on line 2518 of this example script).


<img src="https://github.com/russelllab/kinaseResistance/blob/68b6218879d1c1e53a2bc3c0b8605b125be59fb2/Create_SVG/Version3/sequence_20230228.svg?sanitize=true">

