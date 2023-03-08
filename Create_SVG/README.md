___
# Introduction
## This script
The file **create_svg_20230306_kinases.py** comes with several python2 functions that can annotate or highlight positions in a Clustal-formatted alignment.

___
## General use case
In general, this script needs two things: 
* An alignment in clustal format 
* a python dictionary formatted as {Protein:{Feature:\[Residues]}}, **also see Version3/AlignmentAnnotationdDictionary_Example.txt**
* an optional protein feature dictionary, formatted as (and based on protein of interest) {Feature:[Startposition, Endposition]}, **also see Version3/RegionAnnotationsDictionary_Example.txt**

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
import ast
imposrt sys
from Bio import AlignIO

**Note**: Libraries are also available for >= Python 3.6, print statements and dictionary handling could be converted 

In case we have to deal with .pkl files I did a dirty solution to reformat them into a script that Python2.7 can use. Ironically this script is obviously in Python3. **See reformatter.py**

___
## Features (Version K)
- **NEW** Command line functionality. 
To use the script we can now execute the following command:
`python create_svg_20230306_kinases.py P46734 210 30 MAP2K3.aln Mutational_Infofile.txt Features_Infofile.txt` 
This command has several field after calling the script:

| Field        | Example           | Description  |
| ------------- |:-------------:| -----:|
| 0     | P46734 | The uniprot ID of the protein we are interested in |
| 1     | 210 | The position to be highlighted |
| 2     | 30 | The Windowsize, we show +/- the windowsize of residues around the highlighted position|
| 3     | MAP2K3.aln | The alignment file |
| 4     | Mutational_Infofile.txt | The file containing positional information |
| 5     | Features_Infofile.txt | A file containing structural/domain features, numbering based on **protein of interest** |

**Note**: The script allows for a little hack here. If you want a (large) .svg containing the whole alignment just give a big number in field 2, for example 20000. The script will then produce a complete alignment view.

- **NEW** Named Output files. The resultfile will already be named depending on the input settings, so one can easily try different settings. The name follows this format: 
`poi+"_Position"+str(startposition)+"_Windowsize"+str(windowsize)+".svg"`

- Conservation: Gives a black rectangle as an indicator of sequence identity (top) for the POI residue at that position.

- Residue Numbering: Gives the residue number (every 10 residues), based on sequence of interest, and highlights the input residue in red.

- Feature Annotation: Gives a colored background based on type of annotation (taken from uniprot) to the respective residue.

- Sorted Sequences: Sequences with fewer uniprot annotations are sorted to the bottom of the alignment.

- GAPs removed: Gaps are printed with white color (i.e. invisible on a white background). Additionally, columns with more than 90 % GAPs are removed from the alignment. Sequences affected by this (i.e. the up to 10 % of sequences that did not have a gap at that position) **are kept and not removed**. 

- Highlighting protein features, *here* for example p-loop, Switch I and the Effector region of RHOA. We currently support the displaying of up to 9 features (dependent on the given colors in *featurecolors* on line 2518 of this example script).


MAP2K3, Position 84, Windowsize 30

<img src="https://github.com/russelllab/kinaseResistance/blob/61b2365956c9f8157cf562a5827d359d837e5f74/Create_SVG/Version_K(inases)/MAP2K3_Position84_Windowsize30.svg?sanitize=true">

MAP2K3, Position 110, Windowsize 30

<img src="https://github.com/russelllab/kinaseResistance/blob/dd3f1a3b9d234ea833e8641500c198141347fa88/Create_SVG/Version_K(inases)/MAP2K3_Position110_Windowsize30.svg?sanitize=true">

MAP2K3, Position 210, Windowsize 30

<img src="https://github.com/russelllab/kinaseResistance/blob/61b2365956c9f8157cf562a5827d359d837e5f74/Create_SVG/Version_K(inases)/MAP2K3_Position210_Windowsize30.svg?sanitize=true">

