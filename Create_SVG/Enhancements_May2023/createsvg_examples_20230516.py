#!/usr/bin/env/Python 3.6.8
import sys
import os
import logging
import ast

script = "create_svg_20230516_kinases_GS.py"
mutinfos =  "MutationInfo_Final.txt"
konserve = "GenerelleKonservierung_May-16-2023.txt"
seqident =  "SeqIdentity_Matrix_May-16-2023.txt"
import create_svg_20230516_kinases_GS as create_svg


with open(mutinfos) as f:
	data_align = f.read()
possis = ast.literal_eval(data_align)

with open(konserve) as g:
	overconserv = g.read()
overallconser = ast.literal_eval(overconserv)

with open(seqident) as h:
	data_ident = h.read()
identitydictionary = ast.literal_eval(data_ident)

###create_svg.main(sortingvalue, identitydictionary, overallconservation,\alignmentFile, kinase, mutation_position,\int(ws), int(topN), dic_mutations_info)
def test_createSVG():
	print("1")
	create_svg.main("1", identitydictionary, overallconser,"humanKinasesHitsSplitTrimmedWeb.aln", "BRAF|P15056|430-743", 600,30, 10, possis)
	print("2")
	create_svg.main("2", identitydictionary, overallconser,"humanKinasesHitsSplitTrimmedWeb.aln", "BRAF|P15056|430-743", 600,30, 10, possis)
	print("3")
	create_svg.main("1", identitydictionary, overallconser,"humanKinasesHitsSplitTrimmedWeb.aln", "EGFR|P00533|683-995", 790,30, 10, possis)
	print("4")
	create_svg.main("2", identitydictionary, overallconser,"humanKinasesHitsSplitTrimmedWeb.aln", "EGFR|P00533|683-995", 790,20, 5, possis)
	print("5")
	create_svg.main("3", identitydictionary, overallconser,"humanKinasesHitsSplitTrimmedWeb.aln", "EGFR|P00533|683-995", 790,10, 2, possis)
	print("6")
	create_svg.main("2", identitydictionary, overallconser,"humanKinasesHitsSplitTrimmedWeb.aln", "MAP2K3|P46734|37-347", 84,30, 10, possis)
	print("7")
	create_svg.main("3", identitydictionary, overallconser,"humanKinasesHitsSplitTrimmedWeb.aln", "MAP2K3|P46734|37-347", 84,30, 10, possis)
	print("8")

if __name__ == '__main__':
    test_createSVG()





