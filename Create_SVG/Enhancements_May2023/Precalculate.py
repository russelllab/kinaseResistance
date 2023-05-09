#!/usr/bin/env/Python 3.6.8
import sys
sys.path.append('/home/bq_tschmenger/svgwrite/')	### https://pypi.org/project/svgwrite/
import svgwrite	### preferably install it using pip
#from Bio import AlignIO
import ast
import logging
import re
nums = re.compile(r"[+-]?\d+(?:\.\d+)?")
whitespace_killer=re.compile(r"\s+")
import math
import heapq
from datetime import date
today = date.today()
import logging ### use this in a try: exception: setup, write: logging.exception("message") under except:
### THIS VERSION WORKS FOR PYTHON3


feature_dict = {}
gapletters = [".","-"]
translator = {}
def CUSTOM_ALIGN(targetfile):
	alignments_to_keep = {}
	beginnerdict = {}
	with open(targetfile,"r") as alignfile:
		for line in alignfile:
			newline = whitespace_killer.sub("=",line).replace("\n","")
			idcontainer 	= newline.split("=")[0]
			seq		= newline.split("=")[1].replace("\n","")
			checkvar = str(idcontainer).count("|")		### sp|Q8NEV4|MYO3A_HUMAN|1
			if checkvar == 0:
					try:
						uniprotler = idcontainer.replace(" ","")
						seq_name = uniprotler
						if uniprotler not in translator:
							translator[uniprotler]=seq_name
						final_name = uniprotler
					except:
						final_name = idcontainer
					if final_name not in alignments_to_keep:
						alignments_to_keep[final_name]=str(seq)

			elif checkvar == 1:
					try:
						uniprotler = idcontainer.split("|")[0]
						seq_name = idcontainer.split("|")[1]
						if uniprotler not in translator:
							translator[uniprotler]=seq_name
						final_name = uniprotler
					except:
						final_name = idcontainer
					if final_name not in alignments_to_keep:
						alignments_to_keep[final_name]=str(seq)
			elif checkvar == 2:
					try:
						uniprotler = idcontainer.split("|")[1]
						seq_name = idcontainer.split("|")[2]
						if uniprotler not in translator:
							translator[uniprotler]=seq_name
						final_name = uniprotler
					except:
						final_name = idcontainer
					if final_name not in alignments_to_keep:
						alignments_to_keep[final_name]=str(seq)
			elif checkvar == 3:
					try:
						uniprotler = idcontainer.split("|")[1]
						seq_name = idcontainer.split("|")[2]
						if uniprotler not in translator:
							translator[uniprotler]=seq_name
						final_name = uniprotler
	
					except:
						final_name = idcontainer
					beginpos = idcontainer.split("|")[3]
					if uniprotler not in beginnerdict:
						beginnerdict[uniprotler]=[beginpos]
					if final_name not in alignments_to_keep:
						alignments_to_keep[final_name]=str(seq)						
			else:
					if final_name not in alignments_to_keep:
						alignments_to_keep[final_name]=str(seq)
	return alignments_to_keep, beginnerdict
# ---------------------------------------------------------------------------------------------------------------------------------------------
def most_frequent(listus):
    from collections import Counter
    c = Counter(listus)
    return c.most_common(3)
# ---------------------------------------------------------------------------------------------------------------------------------------------
def overallconserve(seqdict, beginnerdict, relevantpositions):
    overalldictionary = {}
    sequencelength = 1433  ### the sequence length is the same for every sequence in this alignment anyway
    for candidate in seqdict:
        nongapcounter = int(beginnerdict[candidate][0])-1###will correspond to the starting protein position of the respective sequence
        alignposition = 1
        for i in range(0,sequencelength):
            keepITin = "false"
	### i corresponds to the alignment position!
            identitycontainer = []
            buchstabe = seqdict[candidate][i]
            if buchstabe not in gapletters:
                nongapcounter+=1             
            if alignposition not in overalldictionary:
                overalldictionary[alignposition]=[[buchstabe],[]]
            else:
                overalldictionary[alignposition][0].append(buchstabe)
            try:
                for categ in relevantpositions[candidate]:
                    #if candidate == "P15056":
                    #    print(i,"\t",candidate,"\t",nongapcounter,"\t",buchstabe,"\t",relevantpositions[candidate][categ])
                    if str(nongapcounter) in relevantpositions[candidate][categ]:
                        ###if candidate == "P15056":
                            ###print("yes","\t",nongapcounter,"\t",i,"\t",alignposition)
                        if candidate not in overalldictionary[alignposition][1]:
                                            overalldictionary[alignposition][1].append(candidate)
            except:
                pass
            alignposition+=1
    ###print(overalldictionary)
    for a in overalldictionary:
        oldval = overalldictionary[a][0]
        mostrecurrent_raw = most_frequent(oldval) ### this [('.', 518), ('g', 2)] OR [('E', 137), ('Q', 108), ('-', 90)]
        ###print(mostrecurrent_raw)           ### DEBUGGING
        try:
            eins = str(mostrecurrent_raw[0]).split(",")[1].replace(")","").replace("]","") ### this is the raw count
            eins_char = str(mostrecurrent_raw[0]).split(",")[0].replace(")","").replace("]","").replace("(","").replace("'","")
            eins_identitypercent = float(eins)/float(len(oldval))
        except:
            logging.exception("message")
            eins = "null"
        try:
            zwei = str(mostrecurrent_raw[1]).split(",")[1].replace(")","").replace("]","")
            zwei_char = str(mostrecurrent_raw[1]).split(",")[0].replace(")","").replace("]","").replace("(","").replace("'","")
            zwei_identitypercent = float(zwei)/float(len(oldval))
        except:
            zwei = "null"
        try:
            third = str(mostrecurrent_raw[2]).split(",")[1].replace(")","").replace("]","")
            third_char = str(mostrecurrent_raw[2]).split(",")[0].replace(")","").replace("]","").replace("(","").replace("'","")
            third_identitypercent = float(third)/float(len(oldval))
        except:
            third = "null"

        if eins != "null":
                overalldictionary[a] = {}
                overalldictionary[a]["First"]=[eins_identitypercent, eins_char]
                overalldictionary[a]["Allowance"]=[]
        else:
                overalldictionary[a] = {}
                overalldictionary[a]["First"]=["null", "null"]
        if zwei != "null":
                overalldictionary[a]["Second"]=[zwei_identitypercent, zwei_char]
        else:
                overalldictionary[a]["Second"]=["null", "null"]
        if third != "null":
                overalldictionary[a]["Third"]=[third_identitypercent, third_char]
        else:
                overalldictionary[a]["Third"]=["null", "null"]

    return overalldictionary
# ---------------------------------------------------------------------------------------------------------------------------------------------
def seqidentitybestimmung(seqdict, beginnerdict):
	diebestimmung = {}
	for protein in seqdict:
		tempdict = {}
		sequenzler = seqdict[protein]
		tempcharacter = 0
		for i, letter in enumerate(sequenzler,start = 1):
			if letter not in gapletters:
				tempcharacter += 1
			if i not in tempdict:
				tempdict[i]=letter
		for another in seqdict:
			tempcounter = 0
			anotherseq = seqdict[another]
			for i, anotherletter in enumerate(anotherseq,start = 1):
				if anotherletter in tempdict[i]:
					if anotherletter not in gapletters:
						if tempdict[i] not in gapletters:
							tempcounter+=1
			gleichheit = float(tempcounter)/float(tempcharacter)
			if protein not in diebestimmung:
				diebestimmung[protein]={}
				diebestimmung[protein][gleichheit]=[another]
			elif gleichheit not in diebestimmung[protein]:
				diebestimmung[protein][gleichheit]=[another]
			else:
				diebestimmung[protein][gleichheit].append(another)
	return diebestimmung
	#print(diebestimmung)
# ---------------------------------------------------------------------------------------------------------------------------------------------
def SHOWORDER(tophits, identdict, sortregel, seqs, doi, starti, endi, goi, indix_posis):
	# dictionary of interest, residue of interest, windowsize, gene of interest
	showtime = {}
	if sortregel == "1":
		for k in doi:	### uniprot ID = k
			featurecount = []
			sequenzler = seqs[k]
			try:
				residue = int(indix_posis[k][0])-1	
			except:
				residue = 0
			for i, letter in enumerate(sequenzler,start = 1):
				if letter not in gapletters:
					residue += 1
					if i >= starti:
						if i <= endi:
							for v in doi[k]: ### categories, i.e. VARIANT = v
								for vv in doi[k][v]:	### residue number = vv
									if int(vv) == residue:
										if int(vv) not in featurecount:
											featurecount.append(vv)
			if k != goi:							
				if k not in showtime:
					showtime[k]=len(featurecount)
		ranking = sorted(showtime, key=lambda x: (-showtime[x], x))
		ranking.insert(0,goi)
		rankingcomplete = ranking
		for entry in seqs:
			if entry not in rankingcomplete:
				rankingcomplete.append(entry)
	elif sortregel == "2":
		ranking = []
		maxidentities = heapq.nlargest(len(identdict[goi]),identdict[goi]) ### this will retrieve the highest identities
		for val in maxidentities:
			for candidate in identdict[goi][val]:
				if candidate != goi:
					if len(ranking)<= tophits:
						ranking.append(candidate)	
		ranking.insert(0,goi)
		rankingcomplete = ranking
		for entry in seqs:
			if entry not in rankingcomplete:
				rankingcomplete.append(entry)
	else:
		pass
	return ranking, rankingcomplete
# ---------------------------------------------------------------------------------------------------------------------------------------------
alignmentfile = sys.argv[1]
sequences, trackstart 	= CUSTOM_ALIGN(alignmentfile)
with open(sys.argv[2]) as f:
	data_align = f.read()
positions = ast.literal_eval(data_align)
sortingvalue = sys.argv[3]
poi = sys.argv[4]			### "P46734"
startposition = int(sys.argv[5])	### 84
windowsize = int(sys.argv[6])		### 10
###########################################################################################################
###print(sequences)  ### 'Q14289': 'sdiyaeipdet ...	### DEBUGGING
###print(trackstart) ### Q9H4A3': ['197'], ...		### DEBUGGING

SeqIdentity_Dictionary = seqidentitybestimmung(sequences, trackstart)
###print(SeqIdentity_Dictionary) ### {UniprotID: { 0.20634920634920634: ['Q00534', 'P31749', ...	### DEBUGGING
seqiddict_name = "SeqIdentity_Matrix_"+today.strftime("%b-%d-%Y")+".txt"
with open(seqiddict_name, 'w') as f:
    f.write(str(SeqIdentity_Dictionary))


############## PREPARING THE DICTIONARY THAT HAS THE 3 MOST RECURRENT CHARACTERS, THEIR SEQUENCE IDENTITY VALUE AND WHICH PROTEINS HAVE FUNCTIONAL INFORMATION AT THIS POSITION
konservierung = overallconserve(sequences, trackstart, positions)
#print(konservierung)
for seque in sequences:
    identification = seque
    sequenz        = sequences[identification]
    nongapcounter = int(trackstart[seque][0])-1	###will correspond to the starting protein position of the respective sequence
    alignposition = 1
    for i in range(0,len(sequenz)):
        buchstabe = sequenz[i]
        if buchstabe not in gapletters:
            nongapcounter+=1
            try:
                for categ in positions[identification]:
                 ###if identification == "P15056":		### DEBUGGING
                    ###print(i,"\t",identification,"\t",nongapcounter,"\t",buchstabe,"\t",positions[identification][categ],"\t",alignposition)		### DEBUGGING
                 if str(nongapcounter) in positions[identification][categ]:                                              
                    ###print(i,"\t",identification,"\t",nongapcounter,"\t",buchstabe,"\t",positions[identification][categ],"\t",alignposition,"\t","worked")		### DEBUGGING
                    if identification not in konservierung[alignposition]["Allowance"]:
                        konservierung[alignposition]["Allowance"].append(identification)
            except:
                ###logging.exception("message")	### DEBUGGING
                pass
        alignposition += 1
conserv_name = "GenerelleKonservierung_"+today.strftime("%b-%d-%Y")+".txt"
with open(conserv_name, 'w') as f:
        f.write(str(konservierung))
###print(konservierung)		### DEBUGGING



sequence_of_interest = sequences[poi]
non_minus_count = 0
distance_end = len(sequence_of_interest)+100		### to make sure it gets weeded out below, if none of the if statements directly below trigger
distance_start = 0					### to make sure it gets weeded out below, if none of the if statements directly below trigger
try:
    sequencestart = int(wheretobegin[poi][0])
    non_minus_count = sequencestart
except:
    pass
alignmentstart = 0
lettercounter = 0
for i, letter in enumerate(sequence_of_interest,start = 1):
    alignmentstart += 1
    if letter not in gapletters:
        lettercounter += 1
        non_minus_count += 1
        if non_minus_count == startposition:
            startpos = i	### this is the alignment position that corresponds to the residue of interest. alignment position includes "-"
        if non_minus_count == startposition+windowsize+1:
            distance_end = i
        if non_minus_count == startposition-windowsize+1:
            distance_start = i
maxcharactercnt = lettercounter		### should capture the true length of the sequence of interest
### make sure the windowsize does not conflict with positions close to the start or end of the sequence
if distance_start <= 0:
    distance_start = 1
if distance_end > len(sequence_of_interest):
    distance_end = len(sequence_of_interest)
roworder, roworder_complete = SHOWORDER(10, SeqIdentity_Dictionary, sortingvalue, sequences, positions, distance_start, distance_end, poi, trackstart)
roworder = roworder[0:10+1]
###print(roworder) 			### DEBUGGING
###print(roworder_complete) 		### DEBUGGING

###print(konservierung)			### DEBUGGING
###for k in konservierung:		### DEBUGGING
    ###print(k,"\t", konservierung[k])	### DEBUGGING
###print("\n")				### DEBUGGING
forbidden = []
for alignposition in konservierung: ### overalldictionary[a]["First"]=[eins_identitypercent, eins_char]
    firstchar = konservierung[alignposition]["First"][1]
    firstval  = float(konservierung[alignposition]["First"][0])
    listofitems = konservierung[alignposition]["Allowance"]
    if firstval >= 0.90:
        if str(firstchar) in gapletters:
            if bool(set(roworder).intersection(listofitems)) == False:
                ###print(alignposition,"\t", roworder,"\t", listofitems, "\t", "False")		### DEBUGGING
                forbidden.append(alignposition)
            else:
                ###print(alignposition,"\t", roworder,"\t", listofitems, "\t", "True")		### DEBUGGING
                pass
###print(forbidden)			### DEBUGGING






