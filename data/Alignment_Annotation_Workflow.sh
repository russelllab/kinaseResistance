#!/bin/bash
# Bash script to prepare the files required to create annotated alignments and running some tests





########################################################### 	###	#####	#####	#####		  ##	###########################################################
###########################################################	#	  #	#	#   #		 # #	###########################################################
###########################################################	 #	  #	####	#####		#  #	###########################################################
###########################################################	  #       #	#	#		   #	###########################################################
###########################################################	###       #	#####	#		   #	###########################################################

# Collect the alignment file and the dictionary holding the annotating information and store everything together with the script "Precalculate_20230524.py" in the same directory
# alignment file: humanKinasesHitsSplitTrimmedWeb.aln
# annotations: sample_dic_mutation_info.txt

python3 Precalculate_20230524.py humanKinasesHitsSplitTrimmedWeb.aln sample_dic_mutation_info.txt 1 P00533 790 20

# This will create two dictionaries "GenerelleKonservierung" and "SeqIdentityMatrix" with the date of execution
# Resultfile 1: GenerelleKonservierung_May-24-2023.txt
# Resultfile 2: SeqIdentity_Matrix_May-24-2023.txt


########################################################### 	###	#####	#####	#####		#####	###########################################################
###########################################################	#	  #	#	#   #		    #	###########################################################
###########################################################	 #	  #	####	#####		 ####	###########################################################
###########################################################	  #       #	#	#		#	###########################################################
###########################################################	###       #	#####	#		#####	###########################################################

# Run some tests!
# Open the file "createsvg_examples_20230524.py" and edit the following lines to suit the latest version of the above mentioned files

# Example:
#script = "create_svg_20230524_kinases_GS.py"
#mutinfos =  "sample_dic_mutation_info.txt"
#konserve = "GenerelleKonservierung_May-24-2023.txt"
#seqident =  "SeqIdentity_Matrix_May-24-2023.txt"
#import create_svg_20230524_kinases_GS as create_svg

# Make sure that the function "def test_createSVG():" is okay
# a) Check if the correct alignment file is provided
# b) Check that the sequence identifier are valid, i.e. check if "BRAF|P15056|430-762" is in the alignmentfile (instead of say BRAF|P15056|430-743)

python3 createsvg_examples_20230524.py

########################################################### 	###	#####	#####	#####		#####	###########################################################
###########################################################	#	  #	#	#   #		    #	###########################################################
###########################################################	 #	  #	####	#####		 ####	###########################################################
###########################################################	  #       #	#	#		    #	###########################################################
###########################################################	###       #	#####	#		#####	###########################################################

# Annotate Alignments!
# Please see https://github.com/russelllab/kinaseResistance/blob/main/Create_SVG/Enhancements_May2023/README.md
# Example command
python3 create_svg_20230524_kinases_GS.py CHEK2=O96017=193-536 373 15 10 humanKinasesHitsSplitTrimmedWeb.aln sample_dic_mutation_info.txt GenerelleKonservierung_May-24-2023.txt none SeqIdentity_Matrix_May-24-2023.txt 1











