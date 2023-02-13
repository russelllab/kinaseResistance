# Introduction (small one)

The file **create_svg_20230208.py** comes with 3 python functions and a dictionary with annotations.
You will also have to use the second file, which contains an alignment, as input for the script.

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

