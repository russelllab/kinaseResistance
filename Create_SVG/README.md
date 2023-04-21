___
# Introduction
## This script
The file **create_svg_20230314_kinases.py** comes with several python2 functions that can annotate or highlight positions in a Clustal-formatted alignment.

The file **create_svg_20230421_kinases.py** has similar (see below) functionalities and was upgraded to work under Python3.6.
Further work will only be conducted on the Python3+ script.

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

Python 2.7. or Python 3.6

import svgwrite

import ast

imposrt sys

from Bio import AlignIO (only when using Python 2.7)


In case we have to deal with .pkl files I did a quick solution to reformat them into a script that Python2.7 can use. Ironically this script is obviously in Python3. **See reformatter.py**

___
## Features (Vlatest)
- **New** Added tooltip functionality to the svgs.
- **New** Added a parameter to show only the Top X sequences (+ sequence of interest).
- **New** Upgraded the code to work under Python3.6
- **New** Added transparent rectangles to highlight a sequence conservation (= identity) over >= 70 %, based on the sequence of interest. The colors for this are taken from CLUSTAL/Jalview.
- **New** Change highlighting to circles. Circle radius can later be adjusted based on evidence.
- **New** Added basic heatmapping above the alignment, showing how many highlights per position & per category we have.
- **New** Added start and end positions for each displayed sequence.
- Command line functionality. 

To use the script we can now execute the following command:
`python3 create_svg_20230421_kinases.py O96017 372 30 10 humanKinasesTrimmed.clustal Mutational_Infofile_Kinases.txt` 
This command has several fields after calling the script:

| Field        | Example           | Description  |
| ------------- |:-------------:| -----:|
| 0     | P46734 | The uniprot ID of the protein we are interested in |
| 1     | 210 | The position to be highlighted |
| 2     | 30 | The Windowsize, we show +/- the windowsize of residues around the highlighted position. **WARNING** Does not work if the starting position of the sequences is != 1.|
|3      |10 | Shows the top 10 sequences, based on information content|
| 4     | MAP2K3.aln | The alignment file |
| 5     | Mutational_Infofile.txt | The file containing positional information |
| 6     | Features_Infofile.txt | A file containing structural/domain features, numbering based on **protein of interest** |

**Note**: The script allows for a little hack here. If you want a (large) .svg containing the whole alignment just give a big number in field 2, for example 20000. The script will then produce a complete alignment view. **New** Giving "none" instead of a position to be highlighted (field 1) works the same + it removes the position specific rectangle.

- Named Output files. The resultfile will already be named depending on the input settings, so one can easily try different settings. The name follows this format: 
`poi+"_Position"+str(startposition)+"_Windowsize"+str(windowsize)+".svg"`

- Conservation: Gives a black rectangle as an indicator of sequence identity (top) for the POI residue at that position.

- Residue Numbering: Gives the residue number (every 10 residues), based on sequence of interest, and highlights the input residue in red.

- Feature Annotation: Gives a colored background based on type of annotation (taken from uniprot) to the respective residue.

- Sorted Sequences: Sequences with fewer uniprot annotations are sorted to the bottom of the alignment.

- GAPs removed: Gaps are printed with white color (i.e. invisible on a white background). Additionally, columns with more than 90 % GAPs are removed from the alignment. Sequences affected by this (i.e. the up to 10 % of sequences that did not have a gap at that position) **are kept and not removed**. 

- Highlighting protein features, *here* for example p-loop, Switch I and the Effector region of RHOA. We currently support the displaying of up to 9 features (dependent on the given colors in *featurecolors* on line 2518 of this example script).

___
## What is next?
- Make circle size adjustable by evidence.
- add mouseover tooltips

___
## The most recent type of results
MAP2K3, A84 with the most recent testdata (cleaned and trimmed kinase alignment made by Gurdeep). Kinases with no highlights were removed from the SVG visualization. Columns with > 90 % gaps were removed, unless there is a position with an annotation, then the column is kept despite having many gaps (see for example BRAF V600).
<img src="https://github.com/russelllab/kinaseResistance/blob/646b21fcfc6729faf1219a352a9ac8e0679d4a1a/Create_SVG/Version_K(inases)/MP2K3_HUMAN_Position84_Windowsize30000.svg?sanitize=true">


PLK3, a random position close to the activation loop, windowsize 30
<img src="https://github.com/russelllab/kinaseResistance/blob/5e63c6a6d701c3e0062cb9bb3e4cb74b80e7bdb1/Create_SVG/Version_K(inases)/PLK3_Position220_Windowsize30.svg?sanitize=true">


___
## Legacy Results
MAP2K3, Position 84, Windowsize 30

<img src="https://github.com/russelllab/kinaseResistance/blob/ac8fb82c5fbf26a116d23f3b84c61e7c543108b2/Create_SVG/Version_K(inases)/MAP2K3_Position84_Windowsize30.svg?sanitize=true">

MAP2K3, Position 84, Windowsize 30, highlights done via CONNECTOR, region annotations taken from Interpro

<img src="https://github.com/russelllab/kinaseResistance/blob/1f81be16c8fd62121950b02d21f4da526a8962cc/Create_SVG/Version_K(inases)/AnnotatedAlignment.svg?sanitize=true">


MAP2K3, Position 210, Windowsize 30

<img src="https://github.com/russelllab/kinaseResistance/blob/61b2365956c9f8157cf562a5827d359d837e5f74/Create_SVG/Version_K(inases)/MAP2K3_Position210_Windowsize30.svg?sanitize=true">

MAP2K3, Position 210, Windowsize 30, highlights done via CONNECTOR, region annotations taken from Interpro

<img src="https://github.com/russelllab/kinaseResistance/blob/59022d9441ee7207c8f14cf474dafa269d612416/Create_SVG/Version_K(inases)/AnnotatedAlignment_G210C.svg?sanitize=true">
