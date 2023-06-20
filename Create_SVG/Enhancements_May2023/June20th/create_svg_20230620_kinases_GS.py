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
### THIS VERSION WORKS FOR PYTHON3

feature_dict = {}
gapletters = [".","-"]
translator = {}
def CUSTOM_ALIGN(targetfile):	### RPS6KA1|Q15418|33-320 qpskdegvlk
	alignments_to_keep = {}
	beginnerdict = {}
	with open(targetfile,"r") as alignfile:
		for line in alignfile:
			if line.startswith("#") or line.startswith("\n"): continue
			idcontainer 	= line.split(" ")[0]
			seq		= line.split(" ")[1].replace("\n","")
			checkvar 	= str(idcontainer).count("|")		### RPS6KA1|Q15418|33-320 qpskdegvlk
			if idcontainer not in alignments_to_keep:
						alignments_to_keep[idcontainer]=str(seq)
			beginpos 	= idcontainer.split("|")[2].split("-")[0]
			endpos		= idcontainer.split("|")[2].split("-")[1]
			if idcontainer not in beginnerdict:
				beginnerdict[idcontainer]=[beginpos,endpos]
	###print(alignments_to_keep)
	###print(beginnerdict)
	return alignments_to_keep, beginnerdict
# ---------------------------------------------------------------------------------------------------------------------------------------------
def featuresort(tophits, identdict, sortregel, seqs, doi, starti, endi, goi, indix_posis, feature):
	showtime = {}
	for k in doi:	### uniprot ID = k
		featurecount = []
		generalfeatures = []
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
						for klass in doi[k]:
							if feature in klass:
								for vv in doi[k][feature]:	### {'P00519': {'Activating': [['337', '<a href="https://pubmed.ncbi
									if str(residue) in vv:
										if int(residue) not in featurecount:
											featurecount.append(residue)
							else:
								for vv in doi[k][klass]:
									if str(residue) in vv:
										if int(residue) not in generalfeatures:
											generalfeatures.append(residue)	### we do this to count features other than the requested one
		if k != goi:							
			if k not in showtime:
				if len(featurecount) == 0:
					showtime[k]=len(generalfeatures)/1000
				else:
					showtime[k]=len(featurecount)
	ranking = sorted(showtime, key=lambda x: (-showtime[x], x))
	ranking.insert(0,goi)
	rankingcomplete = ranking
	for entry in seqs:
		if entry not in rankingcomplete:
			rankingcomplete.append(entry)
	return ranking, rankingcomplete
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
									if str(residue) in vv:
										if int(residue) not in featurecount:
											featurecount.append(residue)
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
	elif sortregel == "3":
		ranking,rankingcomplete= featuresort(tophits, identdict, sortregel, seqs, doi, starti, endi, goi, indix_posis, "Activating")
	elif sortregel == "4":
		ranking,rankingcomplete= featuresort(tophits, identdict, sortregel, seqs, doi, starti, endi, goi, indix_posis, "Deactivating")
	elif sortregel == "5":
		ranking,rankingcomplete= featuresort(tophits, identdict, sortregel, seqs, doi, starti, endi, goi, indix_posis, "Resistance")
	elif sortregel == "6":
		ranking,rankingcomplete= featuresort(tophits, identdict, sortregel, seqs, doi, starti, endi, goi, indix_posis, "Phosphorylation")
	elif sortregel == "7":
		ranking,rankingcomplete= featuresort(tophits, identdict, sortregel, seqs, doi, starti, endi, goi, indix_posis, "Acetylation")
	elif sortregel == "8":
		ranking,rankingcomplete= featuresort(tophits, identdict, sortregel, seqs, doi, starti, endi, goi, indix_posis, "Ubiquitination")
	elif sortregel == "9":
		ranking,rankingcomplete= featuresort(tophits, identdict, sortregel, seqs, doi, starti, endi, goi, indix_posis, "Sumoylation")
	elif sortregel == "10":
		ranking,rankingcomplete= featuresort(tophits, identdict, sortregel, seqs, doi, starti, endi, goi, indix_posis, "O-GlcNAc")
	elif sortregel == "11":
		ranking,rankingcomplete= featuresort(tophits, identdict, sortregel, seqs, doi, starti, endi, goi, indix_posis, "Methylation")
	else:
		ranking,rankingcomplete= featuresort(tophits, identdict, sortregel, seqs, doi, starti, endi, goi, indix_posis, "Activating")
	return ranking, rankingcomplete
# ---------------------------------------------------------------------------------------------------------------------------------------------
# https://www.jalview.org/help/html/colourSchemes/clustal.html
Clustalcolors = {"A":"hydrophobic",
		"I":"hydrophobic",
		"L":"hydrophobic",
		"M":"hydrophobic",
		"F":"hydrophobic",
		"W":"hydrophobic",
		"V":"hydrophobic",
		"C":"hydrophobic",
		"K":"positive",
		"R":"positive",
		"E":"negative",
		"D":"negative",
		"N":"polar",
		"Q":"polar",
		"S":"polar",
		"T":"polar",
		"G":"glycine",
		"P":"proline",
		"H":"aromatic",
		"Y":"aromatic"}
clustaltypes = {"hydrophobic":"blue",			
		"positive":"red",
		"negative":"magenta",
		"polar":"green",
		"glycine":"black",
		"proline":"orange",
		"aromatic":"cyan"}
# ---------------------------------------------------------------------------------------------------------------------------------------------
def create_svg(sortrule, seqgleichheit, konservierung, sequences, positions, colordict, coloringcategories, featurecolors, startposition, windowsize, poi, forbidden, proteinfeatures, wheretobegin, topguns, path=''):
    heatmapper = {}
    startposition_checker = startposition	
    #### do this when havng constructed the dictionary with interesting positions
    #### here it is supplied as is, but needs to be further modified
    Heatmapstart = 100-(len(coloringcategories)+1)*10
    Konservierungsypsilon = Heatmapstart - 20
    Categoryypsilon = Heatmapstart - 120
    featureypsilon = Heatmapstart - 45
    for item in positions:
        for categ in colordict:
            if categ not in positions[item]:
                positions[item][categ]=[]
    ### adjust the proteinfeatures dictionary
    for feat in proteinfeatures:
        mini = min(proteinfeatures[feat])
        maxi = max(proteinfeatures[feat])
        for i in range(mini,maxi+1):
            if i not in proteinfeatures[feat]:
                proteinfeatures[feat].append(i)
    if startposition == "none":
        startposition = 1
    try:
    	filename = translator[poi]+"_Position"+str(startposition)+"_Windowsize"+str(windowsize)+"_Topguns"+str(topguns)+"_Sorting_"+str(sortrule)+".svg"
    except:
        filename = poi+"_Position"+str(startposition)+"_Windowsize"+str(windowsize)+"_Topguns"+str(topguns)+"_Sorting_"+str(sortrule)+".svg"
    dwg = svgwrite.Drawing(path+filename, profile='full')
    x = 50
    y = 120
    for interesting in sequences:
        if poi in interesting:
            sequence_of_interest = sequences[interesting]
            break
    non_minus_count = 0
    distance_end = len(sequence_of_interest)+100	### to make sure it gets weeded out below, if none of the if statements directly below trigger
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
    roworder, roworder_complete = SHOWORDER(topguns, seqgleichheit, sortrule, sequences, positions, distance_start, distance_end, poi, wheretobegin)
    if topguns != 30000:
        roworder = roworder[0:topguns+1]
    forbidden = []
    for alignposition in konservierung: ### overalldictionary[a]["First"]=[eins_identitypercent, eins_char]
       firstchar = konservierung[alignposition]["First"][1]
       firstval  = float(konservierung[alignposition]["First"][0])
       listofitems = konservierung[alignposition]["Allowance"]
       if firstval >= 0.89:
          if str(firstchar) in gapletters:
              if bool(set(roworder).intersection(listofitems)) == False:
                  forbidden.append(int(alignposition))
    maximumdistance = distance_end - distance_start
    viewboxcounter = 1
    all_x_vals = []
    highlightingID = 0
    highlightsaver = {}
    proteincounter = 0
    for uniprot in roworder_complete:
        textcounter = 0
        proteincounter += 1
        seq 	= sequences[uniprot]
        namus 	= uniprot
        startingpoint = startposition - windowsize	### this is required for the correct labeling according to the sequence of interest
        
        try:
            drawname = translator[namus]
        except:
            drawname = namus
        if poi in namus:	##### this if/else conditional can probably be put in yet another function to reduce the amount of code being used here
            old_x = x
            old_y = y
            x = 50
            y = 100

            if len(drawname) < 8:
                #dwg.add(dwg.rect((x-70, y), (85, 14), fill="yellow"))
                dwg.add(dwg.text(drawname, insert = (x-40,y+7), text_anchor='end', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill='black'))
                dwg.add(dwg.text(drawname, insert = (x-40,Konservierungsypsilon+5), text_anchor='end', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill='black'))

            else:
                #dwg.add(dwg.rect((x-90, y), (85, 14), fill="yellow"))
                dwg.add(dwg.text(drawname, insert = (x-40,y+7), text_anchor='end', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill='black'))
                dwg.add(dwg.text(drawname, insert = (x-40,Konservierungsypsilon+5), text_anchor='end', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill='black'))

        else:
             if proteincounter <= topguns+1: 
                 if len(drawname) < 8:
                     dwg.add(dwg.text(drawname, insert = (x-40,y+7), text_anchor='end', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill='black'))
                 else:
                     dwg.add(dwg.text(drawname, insert = (x-40,y+7), text_anchor='end', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill='black'))
	#charactercount = 0
        totalcount = 0 ### is equivalent with the alignment position
        if startingpoint <= 0:
            startnumberlabel = 1
        elif startingpoint >= maxcharactercnt:
            startnumberlabel = maxcharactercnt
        else:
            startnumberlabel = startingpoint 
        try:
            charactercount = int(wheretobegin[uniprot][0])-1
            startingpoint = charactercount+1
            startnumberlabel = charactercount+1		
        except:
            charactercount = 0

        tempfeat = {}
        featcount = 0
        firstdone = "false"
        lastdone = "false"
        forbidden_start = "false"
        forbidden_end = "false"
        gapcounter = 0
        for i, letter in enumerate(seq, start=1):
            totalcount += 1		#### gives the alignment position, including gaps
            letter = seq[i-1]
            if x not in all_x_vals:
                all_x_vals.append(x)
            if letter not in gapletters:
                charactercount += 1
                if totalcount <= distance_end:	### distance_end refers to the last alignment position that will be considered, which is +windowsize non-gap residues from the input position
                                                ### our total count should always be smaller than the final end position, both refer to the alignment position
                    endcounter = charactercount ### this will be used to draw the last residue number at the end of the row
                    testlenge = int(distance_end)-int(totalcount)
                    if testlenge <= maximumdistance:	### checks that we still operate around the position of interest +/- residues only
                                                        ### for example, if the starting alignment position is 400, and the ending alignment position is 1000, then distance_end(=1000) - totalcount should be smaller than 600 (=maximumdistance)
                        if totalcount >= distance_start: ### the total count must be greater or equal to the distance_start, referring to the alignment position
                            if firstdone == "false":
                                forbidden_start = "true"
                                startcounter = charactercount
                                if proteincounter <= topguns+1:
                                    dwg.add(dwg.text(startcounter, insert=(35, y+8), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill="black"))					
                                firstdone = "true"
                            if totalcount not in forbidden:
                                if poi in namus:
                                    ### {1: {'First': [0.20961538461538462, '.'],
                                    konserv_char = konservierung[totalcount]["First"][1]
                                    konserv_val = konservierung[totalcount]["First"][0]
                                    if konserv_char in gapletters:
                                        konserv_char = konservierung[totalcount]["Second"][1]
                                        konserv_val = konservierung[totalcount]["Second"][0]
                                    if konserv_char in gapletters:
                                        konserv_char = konservierung[totalcount]["Third"][1]
                                        konserv_val = konservierung[totalcount]["Third"][0]
                                    if float(konserv_val)>= 0.7:
                                        dwg.add(dwg.rect((x,y),(10,len(roworder)*20), fill= clustaltypes[Clustalcolors[letter.upper()]], opacity=0.2))

                                    viewboxcounter += 1
                                    if int(charactercount) >= 100:
                                            drawfontsize = "6px"
                                    else:
                                            drawfontsize = "8px"
                                    if int(charactercount) == int(startposition):
                                        position_interest_x = x
                                        position_interest_y = y
                                        dwg.add(dwg.rect((x, Konservierungsypsilon), (10, 14), fill="black"))
                                        toproof = 1-float(konserv_val)
                                        dwg.add(dwg.rect((x, Konservierungsypsilon), (10, 14*toproof), fill="white"))
                                        if startposition_checker != "none":
                                            dwg.add(dwg.text(str(charactercount), insert=(x+5, Konservierungsypsilon-3), text_anchor='middle', dominant_baseline='central', font_size=drawfontsize, font_family='Arial', font_weight='bold', fill='red'))			
                                    else:
                                        if int(charactercount)>= startingpoint:
                                            dwg.add(dwg.rect((x, Konservierungsypsilon), (10, 14), fill="black"))
                                            toproof = 1-float(konserv_val)
                                            dwg.add(dwg.rect((x, Konservierungsypsilon), (10, 14*toproof), fill="white"))
                                            if int(charactercount)%10 == False:
                                                dwg.add(dwg.text(str(charactercount), insert=(x+5, Konservierungsypsilon-3), text_anchor='middle', dominant_baseline='central', font_size=drawfontsize, font_family='Arial', font_weight='bold', fill='black'))
                                    for feat in proteinfeatures:
                                        if int(totalcount) in proteinfeatures[feat]:
                                            if feat not in tempfeat:
                                                tempfeat[feat]=[featurecolors[featcount],featcount]
                                                featcount+=1
                                            elevator = tempfeat[feat][1]
                                            elevator_floor = 0
                                            if elevator >= 3:
                                                if elevator_floor <= 3:
                                                    elevator = elevator_floor
                                                    elevator_floor += 1
                                                else:
                                                    elevator_floor = 0
                                                    elevator = elevator_floor
                                            y_level = featureypsilon + (elevator*5)
                                            if textcounter > 4:
                                                textcounter = 0
                                            textlevel = featureypsilon + (textcounter*6.5)
                                            dwg.add(dwg.rect((x, y_level), (10, 2), fill=tempfeat[feat][0]))
                                            if "done" not in tempfeat[feat]:                                               
                                                dwg.add(dwg.text(str(feat), insert=(x+10+len(str(feat)), textlevel-30), text_anchor='middle', dominant_baseline='central', font_size='8px', font_family='Arial', font_weight='bold', fill=tempfeat[feat][0]))
                                                tempfeat[feat].append("done")	
                                                textcounter+=1		
                                    startnumberlabel+=1
                                try:
                                    drawn = 0
                                    radi = 7
                                    hightlightstring = ""
                                    bigtosmall = sorted(positions[namus], key=lambda kk: len(positions[namus][kk]),reverse=True)
                                    for colorcateg in bigtosmall:
                                        for lister in positions[namus][colorcateg]:
                                            if str(charactercount) in lister:
                                                hyperinfo = lister[2]
                                                mutantinfo = lister[1]
                                                if drawn != 1:
                                                    highlightingID += 1
                                                    drawstring = drawname.split("|")[0]+"|"+drawname.split("|")[1]
                                                    hightlightstring = drawstring.replace(">","")+"|"+letter+str(charactercount) + "++" + mutantinfo+"{"+colorcateg+"{"+hyperinfo
                                                else:
                                                    hightlightstring = hightlightstring + "}" + mutantinfo+"{"+colorcateg+"{"+hyperinfo  
                                                radius = radi
                                                if proteincounter <= topguns+1:
                                                    dwg.add(dwg.circle((x+5, y+7.5), (radius), fill=colordict[colorcateg]))
                                                    dwg.add(dwg.text(letter, insert=(x+5, y+8), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill="black", id=str(highlightingID)))
                                                    drawn = 1
                                                    radi -= 1
                                                if x not in heatmapper:
                                                    heatmapper[x]={}
                                                    heatmapper[x][colorcateg]=1
                                                elif colorcateg not in heatmapper[x]:
                                                    heatmapper[x][colorcateg]=1
                                                else:
                                                    heatmapper[x][colorcateg]+=1
                                        #radi -= 0.5
                                    if drawn == 1:
                                        if str(highlightingID) not in highlightsaver:
                                             highlightsaver[str(highlightingID)]=[x+5,y+7.5,hightlightstring]
                                    if drawn == 0:
                                        if proteincounter <= topguns+1:
                                             dwg.add(dwg.text(letter, insert=(x+5, y+8), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill="black"))
                                except:
                                    logging.exception("message")
                                    if proteincounter <= topguns+1: 
                                        dwg.add(dwg.text(letter, insert=(x+5, y+8), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill="black"))	
                                x += 10
                            else:
                                gapcounter += 1
            else:	### will draw just a "-" for a gap in the alignment
                if totalcount >= distance_start:
                    if totalcount <= distance_end:
                        if totalcount not in forbidden:
                                if poi in namus: 
                                    for feat in proteinfeatures:
                                        if int(totalcount) in proteinfeatures[feat]:
                                            if feat not in tempfeat:
                                                tempfeat[feat]=[featurecolors[featcount],featcount]
                                                featcount+=1
                                            elevator = tempfeat[feat][1]
                                            elevator_floor = 0
                                            if elevator >= 3:
                                                if elevator_floor <= 3:
                                                    elevator = elevator_floor
                                                    elevator_floor += 1
                                                else:
                                                    elevator_floor = 0
                                                    elevator = elevator_floor
                                            y_level = featureypsilon + (elevator*5)
                                            if textcounter > 4:
                                                textcounter = 0
                                            textlevel = featureypsilon + (textcounter*6.5)
                                            dwg.add(dwg.rect((x, y_level), (10, 2), fill=tempfeat[feat][0]))
                                            if "done" not in tempfeat[feat]:                                               
                                                dwg.add(dwg.text(str(feat), insert=(x+10+len(str(feat)), textlevel-30), text_anchor='middle', dominant_baseline='central', font_size='8px', font_family='Arial', font_weight='bold', fill=tempfeat[feat][0]))
                                                tempfeat[feat].append("done")	
                                                textcounter+=1	       
                                x += 10
        if proteincounter <= topguns+1:
            last_display_x = x
            last_display_y = x
        viewboxcounter = last_display_x
        lastx = last_display_x
        lasty = last_display_y
        finalresidue = startcounter+gapcounter+(2*windowsize)
        if proteincounter <= topguns+1:
            dwg.add(dwg.text(endcounter, insert=(lastx+20, y+8), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill="black"))	
							
        if poi in namus:
            dwg.add(dwg.rect((45,y),(lastx-40,14), fill="none",stroke="black",stroke_width=1))	
            x = 50
            y = old_y
        else:
            x = 50
            y += 20
    viewboxwidth = (viewboxcounter+250)
    viewboxheight = len(roworder)*20+100+(100-Categoryypsilon)
    dwg.viewbox(-120, Categoryypsilon-20,viewboxwidth,viewboxheight)

    if startposition_checker != "none":
        dwg.add(dwg.rect((position_interest_x, position_interest_y), (10, len(roworder)*20),fill="none",stroke="black",stroke_width=1)) ### Note: This will fail if the requested position is outside the available alignment, this should not be an issue when given full alignments


    x = 50
    y = 0
    maxfinderdict = {}
    for protein in positions:
        for cat in positions[protein]:
            lenval = len(positions[protein][cat])
            if cat not in maxfinderdict:
                maxfinderdict[cat]=[lenval]
            else:
                maxfinderdict[cat].append(lenval)

    maxfinder = {}
    for xval in heatmapper:
        for category in colordict:
            if category not in heatmapper[xval]:
                heatmapper[xval][category]=0
        for categ in heatmapper[xval]:
            if categ not in maxfinder:
                maxfinder[categ]=[int(heatmapper[xval][categ])]
            else:
                maxfinder[categ].append(int(heatmapper[xval][categ]))
    for allxval in all_x_vals:
       if allxval not in heatmapper:
           heatmapper[allxval]={}
           for category in colordict:
               if category not in heatmapper[allxval]:		
                   heatmapper[allxval][category]=0.0

    mapx = 40
    mapy = Heatmapstart
    for category in colordict:
        try:
            heatmap_maximum = max(maxfinderdict[category]) ### use this if we want to scale the colors based on the complete available information
                                                           ### use maxfinder[category] instead if we want to scale the colors based on the chosen windowsize (still considering everything, even from sequences
							   ### that are not shown
        except:
            heatmap_maximum = 1
        dwg.add(dwg.text(category, insert=(45, mapy+5), text_anchor='end', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill='black'))
        for xval in heatmapper:
            try:
                opac = float(heatmapper[xval][category])/float(heatmap_maximum)
            except: ### this might happen if we do not have anything for one of the 3 categories A, D, R. In this care heatmap_maximum might be 0
                opac = 0.0
            if float(opac) == 0.0:
                dwg.add(dwg.rect((xval, mapy), (10, 10), fill="white", opacity = 0.15 ))
            else:
                dwg.add(dwg.rect((xval, mapy), (10, 10), fill=colordict[category], opacity = opac ))
            if mapy == 60:
                pass		#### I need to get all necessary xvals first, dammit
        dwg.add(dwg.rect((50, mapy), (lastx-mapx-10, 10),fill="none",stroke="black",stroke_width=0.5))	### <<<<
        mapy += 10
    for i in range(50,lastx-10,10):
        dwg.add(dwg.rect((i, Heatmapstart), (10, len(colordict)*10),fill="none",stroke="black",stroke_width=0.5))

    x = 50
    y = 0
    categocounter = 0
    for category in colordict:
        categocounter += 1
        if categocounter > 5:
            categocounter = 1
            Categoryypsilon += 10
            x = 50
        dwg.add(dwg.rect((x-30, Categoryypsilon), (80, 10), fill=colordict[category]))
        dwg.add(dwg.text(category, insert=(x+10, Categoryypsilon+5), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill='black'))
        x += 80


    dwg.save()

    styletext = """<style>
   <![CDATA[
    text.moo {
         font-family: "arial";
         fill: black;
         font-size: 50%;
    }
    text.hyper {
         font-family: "arial";
         
         font-size: 0.45em;
    }
    rect.hiss {
         fill:white;
    }
   ]]>
   .bootstrap {display: none;}
   svg text.moo {display: none;}
   svg text.hyper {display: none;}
   svg rect.hiss {display: none;}
   svg g:hover text {display: block;}
   svg g:hover rect {display: block;}
   svg g:hover .bootstrap {display: block;}
 </style>"""

    imagefile = open(path+filename,"r", encoding='utf-8')
    data = imagefile.read()
    newdata = data.replace("</svg>", styletext+"</svg>")
    imagefile.close()
    writeFile = open(path+filename, "w", encoding='utf-8')
    writeFile.write(newdata)
    writeFile.close()

    circletext = ""
    for hlid in highlightsaver:
        txt = highlightsaver[hlid][2]
        uppertext = txt.split("++")[0]
        lowertext = txt.split("++")[1]
        cx = highlightsaver[hlid][0]
        cy = highlightsaver[hlid][1]
        delty = len(lowertext.split("}"))*10
        tspany = cy+15
        whiteboxheight = len(lowertext.split("}"))*20+30
        tspanner = ""
        ### mutation type source
        for showfeature in lowertext.split("}"):
            showmutation = showfeature.split("{")[0]
            showfeaturetext = showfeature.split("{")[1]

            tspanfeaturetext = showmutation+"\t"+showfeaturetext
            tspanfeaturetextcolor = colordict[showfeaturetext]
            linkinformation = showfeature.split("{")[2]
            bootstrapper = """<svg xmlns="http://www.w3.org/2000/svg" class="bootstrap" width="8" height="8" x='"""+str(cx)+"""' y='"""+str(tspany-35-delty)+"""' fill='blue'  viewBox="0 0 16 16">
  <path fill-rule="evenodd" d="M8.636 3.5a.5.5 0 0 0-.5-.5H1.5A1.5 1.5 0 0 0 0 4.5v10A1.5 1.5 0 0 0 1.5 16h10a1.5 1.5 0 0 0 1.5-1.5V7.864a.5.5 0 0 0-1 0V14.5a.5.5 0 0 1-.5.5h-10a.5.5 0 0 1-.5-.5v-10a.5.5 0 0 1 .5-.5h6.636a.5.5 0 0 0 .5-.5z"/>
  <path fill-rule="evenodd" d="M16 .5a.5.5 0 0 0-.5-.5h-5a.5.5 0 0 0 0 1h3.793L6.146 9.146a.5.5 0 1 0 .708.708L15 1.707V5.5a.5.5 0 0 0 1 0v-5z"/></svg>"""


            tspanner = tspanner + """<a href=\""""+linkinformation+"""\" target="_blank">"""+bootstrapper+"""<text class="hyper" fill='blue' x='"""+str(cx+12)+"""' y='"""+str(tspany-28-delty)+"""'><tspan class="text">"""+str(tspanfeaturetext)+"""</tspan></text></a>"""
            tspany += 15

        circletext = circletext+"""<g xmlns="http://www.w3.org/2000/svg">
          <circle xmlns="http://www.w3.org/2000/svg" cx='"""+str(cx)+"""' cy='"""+str(cy)+"""' r="7" style="fill:transparent;stroke:transparent;stroke-width:0.5;fill-opacity:0.25;stroke-opacity:0.25"/>      
          <rect class="hiss" x='"""+str(cx-5)+"""' y='"""+str(cy-40-delty)+"""' height='"""+str(whiteboxheight)+"""' width='"""+str(len(uppertext)+80)+"""'></rect>
          <text class="moo" x='"""+str(cx)+"""' y='"""+str(cy-28-delty)+"""'><tspan class="text">"""+uppertext+"""</tspan></text>"""+tspanner+"""</g>"""

    imagefile = open(path+filename,"r", encoding='utf-8')
    imagefile.seek(0)
    data= imagefile.read()
    imagefile.close()
    newdata = data.replace("</svg>", circletext+"</svg>")
    writeFile = open(path+filename, "w", encoding='utf-8')
    writeFile.write(newdata)
    writeFile.close()
    return filename

# ---------------------------------------------------------------------------------------------------------------------------------------------
def main(sortingvalue, identitydictionary,overallconservation, alignmentfile, protein_of_interest, position_of_interest, window, topguns, positions, feature_dict, path=''):
	sequences, trackstart 	= CUSTOM_ALIGN(alignmentfile)
	TheForbiddenPositions = []
	
	### to define the annotation colors we want to use

	positioncolors = ["#009e73","#d55e00","#0072b2","lightgreen","salmon","yellow","blueviolet","deeppink","olive","dodgerblue","palegreen"]
	generalcategories = ["Activating","Deactivating","Resistance","Phosphorylation","Acetylation","Ubiquitination","Sumoylation","O-GlcNAc","Methylation"]
	colors = {}
	coloringcategories = []
	counter = 0
	for k in positions:
		for v in positions[k]:
			if v not in coloringcategories: 
				coloringcategories.append(v)
	for item in generalcategories:
		if item in coloringcategories:
			colors[item]=positioncolors[counter]
			counter+=1
			
	for seqident in sequences:
		if seqident not in positions:
			positions[seqident]={}
			for colcateg in colors:
				positions[seqident][colcateg]=[]

	featurecolors = ["firebrick","tomato","orange","olive","black","teal","dodgerblue","blueviolet","deeppink",
			"firebrick","tomato","orange","olive","black","teal","dodgerblue","blueviolet","deeppink",
			"firebrick","tomato","orange","olive","black","teal","dodgerblue","blueviolet","deeppink",
			"firebrick","tomato","orange","olive","black","teal","dodgerblue","blueviolet","deeppink",
			"firebrick","tomato","orange","olive","black","teal","dodgerblue","blueviolet","deeppink",
			"firebrick","tomato","orange","olive","black","teal","dodgerblue","blueviolet","deeppink",
			"firebrick","tomato","orange","olive","black","teal","dodgerblue","blueviolet","deeppink"]
	###

	filename = create_svg(sortingvalue, identitydictionary, overallconservation, sequences, positions, colors, coloringcategories, featurecolors, position_of_interest, window, protein_of_interest, TheForbiddenPositions, feature_dict, trackstart, topguns, path)	
	return filename

if __name__ == "__main__":
	alignmentfile = sys.argv[5]	#### change this to the location of the alignmentfile.
	###
	protein_of_interest = sys.argv[1].replace("=","|")
	try:
		position_of_interest = int(sys.argv[2])
	except:
		position_of_interest = str(sys.argv[2])
	window = int(sys.argv[3])
	topguns = int(sys.argv[4])

	if position_of_interest == "none":
		window = 30000

	with open(sys.argv[6]) as f:
		data_align = f.read()
	positions = ast.literal_eval(data_align)
	with open(sys.argv[7]) as g:
		overconserv = g.read()
	overallconservation = ast.literal_eval(overconserv)
	with open(sys.argv[9]) as h:
		data_ident = h.read()
	identitydictionary = ast.literal_eval(data_ident)
	sortingvalue = sys.argv[10]
	feature_dict = {}
	try:
		if sys.argv[8] != "none":
			with open(sys.argv[8]) as ff:
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
	except:
		#print(logging.exception("message"))
		pass
	main(sortingvalue, identitydictionary, overallconservation, alignmentfile, protein_of_interest, position_of_interest, window, topguns, positions, feature_dict)
