## Changelog 05. June 2023
Updated the feature annotation to work with the contents of the file ss.tsv. My script reads the lines without "#" and replaces any number of spaces with a single space, then reads the feature and the start & end positions. This file can therefore be updated later and should not break the script as long as the format remains the same.

Do either

`create_svg.main("1", identitydictionary, overallconser,"humanKinasesHitsSplitTrimmedWeb.aln", "BRAF|P15056|430-762", 600,30, 10, possis, feature_dict)`

or

`python3 create_svg_20230605_kinases_GS.py MAP2K3=P46734=37-347 84 30000 30000 humanKinasesHitsSplitTrimmedWeb.aln sample_dic_mutation_info.txt GenerelleKonservierung_May-24-2023.txt ss.tsv SeqIdentity_Matrix_May-24-2023.txt 1`

## Changelog 24. May 2023
- Cosmetic fixes. Removed yellow highlighting of sequence of interest and shortened horizontal rectangle to avoid inconsistencies with sequence identifiers of significantly different number of characters.
- Adjusted the font size of labels on the left (imaginary y-axis) axis to the same size.
- Modified the mouse-over information. In line 1 in only says Genename|UniprotID|Mutation, the sequence range was removed to reduce clutter of the infobox.


## Changelog 16. May 2023
- Adjusted code to accomodate for the new alignment, using different identifiers in case of multiple kinase domains per protein. Changed the position of the bootstrap arrow box and fixed the colors.
Everything can be found in [Create_SVG/Enhancements_May2023/16thMay](https://github.com/russelllab/kinaseResistance/tree/main/Create_SVG/Enhancements_May2023/16thMay)

## Changelog 12. May 2023
- Added [Bootstrap Box arrow up-right](https://icons.getbootstrap.com/icons/box-arrow-up-right/) icon as .svg to the .svg.
You will find it in [Create_SVG/Enhancements_May2023/May11th](https://github.com/russelllab/kinaseResistance/tree/main/Create_SVG/Enhancements_May2023/May11th) with the name **create_svg_20230512_kinases_GS.py**

## Changelog 11. May 2023
- See [Create_SVG/Enhancements_May2023/May11th](https://github.com/russelllab/kinaseResistance/tree/main/Create_SVG/Enhancements_May2023/May11th)
- Added Clickable hyperlinks
- Cosmetic changes 
  - changed position of the mouseover box to make hyperlinks more easily clickable
  - changed rule how circle radii are determined. **Now**: Circle radius is determined for categories with a decreasing number of entries, zum Beispiel categories with the most entries will have the biggest circle radius


## Changelog 10. May 2023
- See [Create_SVG/Enhancements_May2023/May10th](https://github.com/russelllab/kinaseResistance/tree/main/Create_SVG/Enhancements_May2023/May10th)
- Changed some of the positionings (heatmap, color legend, conservation barplots) to be dynamically placed based on the length of provided functional annotations
- Added new sorting rules

| Sorting Value        | Sorting by ...  | 
| ------------- |:-------------:| 
| 1 | No. of functional infos |
| 2 | Sequence Identity |
| 3| Activating |
| 4 |Deactivating  |
| 5 | Resistance |
| 6 |Phosphorylation  |
| 7 | Acetylation |
| 8 | Ubiquitination |
| 9 | Sumoylation |
| 10 | O-GlcNAc |
| 11 |  Methylation|
     
- Please use the dictionaries **GenerelleKonservierung_May-10-2023.txt** and **SeqIdentity_Matrix_May-10-2023.txt**
- **No changes to how the script or main is called!**



## Changelog 09. May 2023
- Conservation barplots now go from bottom to top
  - Conservation shown is based on the first non-gap character frequency. This fixes the issue that BRAF V600E was previously shown as highly conserved, although that came from this position being ~ 99 % GAPs.
  - Conservation is now calculated based on the complete alignment, not only on the shown sequences
- Heatmap is now scaled based on the complete alignment, not only on shown sequences. Additionally the colors are scaled based on the highest absolute value for a given category, and are not based on the visible window
- Split the mouseover information into two lines, first line showing Mechismo-Style format, second line showing information
  - This enables the user in principle to copy & paste the information, as the user mouse can now access the mouseover box
- Adjusted the Activating/Deactivating/Resistance labels and boxes
- Introduced a second sorting rule. Default is **1**, that is sorting by available functional information. Rule **2** sorts by sequence identity.

## Using the updated script

To use the script we execute the following command:

`python3 create_svg_20230509_kinases_GS.py P46734 84 9 5 humanKinasesTrimmed.clustal Mutational_Infofile_Kinases_V2.txt GenerelleKonservierung_May-09-2023.txt none SeqIdentity_Matrix_May-09-2023.txt 2` 

This command has several fields after calling the script, explained below:

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

`python3 Precalculate.py humanKinasesTrimmed.clustal Mutational_Infofile_Kinases_V2.txt 2 P15056 600 20`.

The last 4 fields ("2 P15056 600 20") are only for debugging purposes. Simply provide the alignment and the information and the aforementioned dictionaries will be written to disk.
