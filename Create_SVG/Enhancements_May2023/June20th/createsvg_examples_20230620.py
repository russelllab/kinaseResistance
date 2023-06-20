#!/usr/bin/env/Python 3.6.8
import sys
import os
import logging
import ast
import re
whitespace_killer=re.compile(r"\s+")

script = "create_svg_202306020_kinases_GS.py"
mutinfos =  "sample_dic_mutation_info.txt"
konserve = "GenerelleKonservierung_Jun-20-2023.txt"
seqident =  "SeqIdentity_Matrix_Jun-20-2023.txt"
featurefile = "ss.tsv"
import create_svg_20230620_kinases_GS as create_svg


with open(mutinfos) as f:
	data_align = f.read()
possis = ast.literal_eval(data_align)

with open(konserve) as g:
	overconserv = g.read()
overallconser = ast.literal_eval(overconserv)

with open(seqident) as h:
	data_ident = h.read()
identitydictionary = ast.literal_eval(data_ident)

feature_dict = {}
with open(featurefile) as ff:
	for line in ff:
		if "#" not in line:
			newline = whitespace_killer.sub(" ",line)
			feature = newline.split(" ")[0]
			start = int(newline.split(" ")[2].split("-")[0])
			end = int(newline.split(" ")[2].split("-")[1].replace("\n",""))
			if feature not in feature_dict:
				feature_dict[feature]=[]
			for i in range(start, end+1):
				if i not in feature_dict[feature]:
					feature_dict[feature].append(i)

###create_svg.main(sortingvalue, identitydictionary, overallconservation,\alignmentFile, kinase, mutation_position,\int(ws), int(topN), dic_mutations_info)
def test_createSVG():
	print("1")
	create_svg.main("1", identitydictionary, overallconser,"humanKinasesPkinasePK_Tyr_Ser-ThrAllTrimmedWeb.aln", "BRAF|P15056|405-748", 600,30, 10, possis, feature_dict)
	print("2")
	create_svg.main("2", identitydictionary, overallconser,"humanKinasesPkinasePK_Tyr_Ser-ThrAllTrimmedWeb.aln", "BRAF|P15056|405-748", 600,30, 10, possis, feature_dict)
	print("3")
	create_svg.main("1", identitydictionary, overallconser,"humanKinasesPkinasePK_Tyr_Ser-ThrAllTrimmedWeb.aln", "EGFR|P00533|671-1000", 790,30, 10, possis, feature_dict)
	print("4")
	create_svg.main("2", identitydictionary, overallconser,"humanKinasesPkinasePK_Tyr_Ser-ThrAllTrimmedWeb.aln", "EGFR|P00533|671-1000", 790,20, 5, possis, feature_dict)
	print("5")
	create_svg.main("3", identitydictionary, overallconser,"humanKinasesPkinasePK_Tyr_Ser-ThrAllTrimmedWeb.aln", "EGFR|P00533|671-1000", 790,10, 2, possis, feature_dict)
	print("6")
	create_svg.main("2", identitydictionary, overallconser,"humanKinasesPkinasePK_Tyr_Ser-ThrAllTrimmedWeb.aln", "MAP2K3|P46734|29-347", 84,30, 10, possis, feature_dict)
	print("7")
	create_svg.main("3", identitydictionary, overallconser,"humanKinasesPkinasePK_Tyr_Ser-ThrAllTrimmedWeb.aln", "MAP2K3|P46734|29-347", 84,30, 10, possis, feature_dict)
	print("8")

if __name__ == '__main__':
    test_createSVG()

