#!/usr/bin/env/Python 3.6.8
## https://github.com/googlesamples/assistant-sdk-python/issues/236
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
### THIS VERSION WORKS FOR PYTHON3

feature_dict = {}
gapletters = [".","-"]
translator = {}
def CUSTOM_ALIGN(targetfile):
	alignments_to_keep = {}
	beginnerdict = {}
	with open(targetfile,"r") as alignfile:
		for line in alignfile:
			#print(line.replace("\n",""))
			newline = whitespace_killer.sub("=",line).replace("\n","")
			idcontainer 	= newline.split("=")[0]
			#print(idcontainer," except")
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
#	print(alignments_to_keep)
#	print(beginnerdict)
	return alignments_to_keep, beginnerdict
# ---------------------------------------------------------------------------------------------------------------------------------------------
def conservation_checker(identifier, seqdict, beginnerdict, relevantpositions, rankedlist):
	for protein in rankedlist:
		#print(protein,"\t",identifier)
		sequence = seqdict[protein]
		if identifier in protein:
			truesequenzler = sequence
	sequencelength = len(truesequenzler) ### including gaps, meaning this is the alignment length
	positionalcounter = int(beginnerdict[identifier][0])
	ident_startposition = positionalcounter
	conservational_dictionary = {}
	theforbiddenalignpos = []
	for i in range(0,sequencelength):
		keepITin = "false"
		### i corresponds to the alignment position!
		identitycontainer = []
		for ident in rankedlist:
			if identifier in ident:
				orires = seqdict[ident][i]
				identitycontainer.append(seqdict[ident][i].upper())
			else:
				identitycontainer.append(seqdict[ident][i].upper())
		identitypercentage = float(identitycontainer.count(orires.upper()))/float(len(identitycontainer))	### so far this also includes "-" as the original truesequence residue, be cautious
		oritype = "none"

		if orires not in gapletters:
#			print positionalcounter,"\t", identitycontainer, "\t", orires,"\t",identitypercentage,"\t",typuspercentage ,"\t", aminotypes,"\t",oritype
			if int(positionalcounter) not in conservational_dictionary:
				conservational_dictionary[int(positionalcounter)] = [float(identitypercentage), orires]
			positionalcounter+=1
		elif orires in gapletters:
			if identitypercentage >= 0.90:	### more than 90% GAPs at this alignment position
				for uniprotID in rankedlist:
					relevantstartingposition = int(beginnerdict[uniprotID][0])
					relevantsequence = seqdict[uniprotID]
					residue_counter = 0
					for a in range(0,i):	### going to the correct alignment position
						residueletter = relevantsequence[a]
						if residueletter not in gapletters:
							residue_counter += 1
					if residueletter not in gapletters:					
						for category in relevantpositions[uniprotID]:
							if str(residue_counter+relevantstartingposition) in relevantpositions[uniprotID][category]:
								keepITin = "true"
						
				if keepITin == "false":
					theforbiddenalignpos.append(i+1)	
		else:
			pass
	return conservational_dictionary, theforbiddenalignpos	
# ---------------------------------------------------------------------------------------------------------------------------------------------
def SHOWORDER(seqs, doi, starti, endi, goi, indix_posis):
	# dictionary of interest, residue of interest, windowsize, gene of interest
	showtime = {}
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
	#print showtime
	ranking = sorted(showtime, key=lambda x: (-showtime[x], x))
	ranking.insert(0,goi)
	#print ranking
	#for ranked in ranking:
	#	print ranked, "\t", showtime[ranked]
	return ranking
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
def create_svg(sequences_dict, positions, colordict, coloringcategories, featurecolors, startposition, windowsize, poi, forbidden, proteinfeatures, wheretobegin, topguns, path=''):
    heatmapper = {}
    startposition_checker = startposition	
    #### do this when havng constructed the dictionary with interesting positions
    #### here it is supplied as is, but needs to be further modified
    for item in positions:
        for categ in colordict:
            if categ not in positions[item]:
                positions[item][categ]=[]
    ### adjust the proteinfeatures dictionary
    for feat in proteinfeatures:
        oldstuff = proteinfeatures[feat]
        proteinfeatures[feat] = map(int, oldstuff)
        maxi = max(proteinfeatures[feat])
        mini = min(proteinfeatures[feat])
        for i in range(mini,maxi+1):
            proteinfeatures[feat].append(i)
    if startposition == "none":
        startposition = 1
    try:
    	filename = translator[poi]+"_Position"+str(startposition)+"_Windowsize"+str(windowsize)+".svg"
    except:
        filename = poi+"_Position"+str(startposition)+"_Windowsize"+str(windowsize)+".svg"
    dwg = svgwrite.Drawing(filename, profile='full')
    x = 50
    y = 120
    sequence_of_interest = sequences_dict[poi]
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

    #print(non_minus_count)
    for i, letter in enumerate(sequence_of_interest,start = 1):
        alignmentstart += 1
        #print(i)
        #print(letter)
        if letter not in gapletters:
            lettercounter += 1
            non_minus_count += 1
            if non_minus_count == startposition:
                startpos = i	### this is the alignment position that corresponds to the residue of interest. alignment position includes "-"
            if non_minus_count == startposition+windowsize:
                distance_end = i
            if non_minus_count == startposition-windowsize:
                distance_start = i
    maxcharactercnt = lettercounter		### should capture the true length of the sequence of interest
    ### make sure the windowsize does not conflict with positions close to the start or end of the sequence

    #print(distance_start)
    #print(distance_end)
    #print(len(sequence_of_interest))

    if distance_start <= 0:
        distance_start = 1
    if distance_end > len(sequence_of_interest):
        distance_end = len(sequence_of_interest)

    #print(distance_start)
    #print(distance_end)
    #print(len(sequence_of_interest))

    roworder = SHOWORDER(sequences_dict, positions, distance_start, distance_end, poi, wheretobegin)
    if topguns != 30000:
        roworder = roworder[0:topguns+1]
    Konserve, forbidden 	= conservation_checker(poi ,sequences_dict, wheretobegin, positions, roworder)
    maximumdistance = distance_end - distance_start
    viewboxcounter = 1
    all_x_vals = []

    #print(distance_start)#560
    #print(distance_end)#1152
    #print(maximumdistance)#592
    #print(startposition)#372
    #print(windowsize)#30
    #print(startpos)#912
    #print(maxcharactercnt)#517
    highlightingID = 0
    highlightsaver = {}
    for uniprot in roworder:
        seq 	= sequences_dict[uniprot]
        namus 	= uniprot
        startingpoint = startposition - windowsize	### this is required for the correct labeling according to the sequence of interest
        
        try:
            drawname = translator[namus]
        except:
            drawname = namus

	#print drawname
        if poi in namus:	##### this if/else conditional can probably be put in yet another function to reduce the amount of code being used here
            old_x = x
            old_y = y
            x = 50
            y = 100

            dwg.add(dwg.rect((x-90, y), (80, 14), fill="yellow"))
            if len(drawname) < 8:
                dwg.add(dwg.text(drawname, insert = (x-60,y+7), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill='black'))
            else:
                dwg.add(dwg.text(drawname, insert = (x-60,y+7), text_anchor='middle', dominant_baseline='central', font_size='7px', font_family='Arial', font_weight='bold', fill='black'))
        else:
             if len(drawname) < 8:
                 dwg.add(dwg.text(drawname, insert = (x-60,y+7), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill='black'))
             else:
                 dwg.add(dwg.text(drawname, insert = (x-60,y+7), text_anchor='middle', dominant_baseline='central', font_size='7px', font_family='Arial', font_weight='bold', fill='black'))
	#charactercount = 0
        totalcount = 0
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

        #print(uniprot)
        #print(charactercount)
        #print(startingpoint)

        #print(distance_start)
        #print(distance_end)
        #print(maximumdistance)
        #print(startnumberlabel) 
        #print(startposition)
    
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
                                dwg.add(dwg.text(startcounter, insert=(35, y+8), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill="black"))					
                                firstdone = "true"
                            if totalcount not in forbidden:
                                #print(uniprot,"\t",charactercount,"\t",startnumberlabel)
                                if poi in namus:
                                    if float(Konserve[charactercount][0])>= 0.7:
                                        dwg.add(dwg.rect((x,y),(10,len(roworder)*20), fill= clustaltypes[Clustalcolors[letter.upper()]], opacity=0.2))

                                    viewboxcounter += 1
                                    if int(charactercount) == int(startposition):
                                        position_interest_x = x
                                        position_interest_y = y
                                        dwg.add(dwg.rect((x, 40), (10, 14*float(Konserve[charactercount][0])), fill="black"))
                                        if startposition_checker != "none":
                                            dwg.add(dwg.text(str(charactercount), insert=(x+5, 37), text_anchor='middle', dominant_baseline='central', font_size='8px', font_family='Arial', font_weight='bold', fill='red'))			
                                    else:
                                        if int(charactercount)>= startingpoint:
                                            dwg.add(dwg.rect((x, 40), (10, 14*float(Konserve[charactercount][0])), fill="black"))
                                            if int(charactercount)%10 == False:
                                                dwg.add(dwg.text(str(charactercount), insert=(x+5, 37), text_anchor='middle', dominant_baseline='central', font_size='8px', font_family='Arial', font_weight='bold', fill='black'))
                                    for feat in proteinfeatures:
                                        if charactercount in proteinfeatures[feat]:
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
                                            y_level = 22 + (elevator*3)
                                            dwg.add(dwg.rect((x, y_level), (10, 2), fill=tempfeat[feat][0]))
                                            if "done" not in tempfeat[feat]:
                                                dwg.add(dwg.text(str(feat), insert=(x+15, 20), text_anchor='middle', dominant_baseline='central', font_size='8px', font_family='Arial', font_weight='bold', fill=tempfeat[feat][0]))
                                                tempfeat[feat].append("done")			
                                    startnumberlabel+=1
                                try:
                                    drawn = 0
                                    radi = 7
                                    hightlightstring = ""
                                    for colorcateg in coloringcategories:
                                        if str(charactercount) in positions[namus][colorcateg]:
                                            if drawn != 1:
                                                highlightingID += 1
                                                hightlightstring = drawname+"/"+str(charactercount) + "/" + colorcateg
                                            else:
                                                hightlightstring = hightlightstring + "/" + colorcateg  
                                            try:
                                                radius = math.log10(int(positions[namus][colorcateg][str(charactercount)]))
						### assuming a new nested dictionary like this: {unip:{colorcateg/A/R/D:{position:evidence}}}
                                            except:
                                                radius = radi
                                            dwg.add(dwg.circle((x+5, y+7.5), (radius), fill=colordict[colorcateg]))
                                            dwg.add(dwg.text(letter, insert=(x+5, y+8), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill="black", id=str(highlightingID)))
                                            drawn = 1
                                            if x not in heatmapper:
                                                heatmapper[x]={}
                                                heatmapper[x][colorcateg]=1
                                            elif colorcateg not in heatmapper[x]:
                                                heatmapper[x][colorcateg]=1
                                            else:
                                                heatmapper[x][colorcateg]+=1
                                        radi -= 1
                                    if drawn == 1:
                                        if str(highlightingID) not in highlightsaver:
                                            highlightsaver[str(highlightingID)]=[x+5,y+7.5,hightlightstring]
                                    if drawn == 0:
                                        dwg.add(dwg.text(letter, insert=(x+5, y+8), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill="black"))
                                except:
                                    logging.exception("message") 
                                    dwg.add(dwg.text(letter, insert=(x+5, y+8), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill="black"))	
                                x += 10
                            else:
                                gapcounter += 1
            else:	### will draw just a "-" for a gap in the alignment
                if totalcount >= distance_start:
                    if totalcount <= distance_end:
                        if totalcount not in forbidden:
                            x += 10
        viewboxcounter = x
        lastx = x
        lasty = y
        finalresidue = startcounter+gapcounter+(2*windowsize)
        dwg.add(dwg.text(endcounter, insert=(lastx+20, y+8), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill="black"))	
							
        if poi in namus:
            dwg.add(dwg.rect((-40,y),(x+50,14), fill="none",stroke="black",stroke_width=1))	
            x = 50
            y = old_y
        else:
            x = 50
            y += 20
    viewboxwidth = (viewboxcounter+200)
    viewboxheight = len(roworder)*20+100
    dwg.viewbox(-70, 0,viewboxwidth,viewboxheight)

    if startposition_checker != "none":
    	dwg.add(dwg.rect((position_interest_x, position_interest_y), (10, len(roworder)*20),fill="none",stroke="black",stroke_width=1))

    x = 50
    y = 0

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
    mapy = 60
    for category in colordict:
        heatmap_maximum = max(maxfinder[category])
        dwg.add(dwg.text(category, insert=(40, mapy+5), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill='black'))
        for xval in heatmapper:
	    #print xval, "\t", mapy
            opac = float(heatmapper[xval][category])/float(heatmap_maximum)
            if float(opac) == 0.0:
                dwg.add(dwg.rect((xval, mapy), (10, 10), fill="lightblue", opacity = 0.15 ))
            else:
                dwg.add(dwg.rect((xval, mapy), (10, 10), fill=colordict[category], opacity = opac ))
            if mapy == 60:
                pass		#### I need to get all necessary xvals first, dammit
        dwg.add(dwg.rect((50, mapy), (lastx-mapx-10, 10),fill="none",stroke="black",stroke_width=0.5))	### <<<<
        mapy += 10
    for i in range(50,lastx-10,10):
        dwg.add(dwg.rect((i, 60), (10, 30),fill="none",stroke="black",stroke_width=0.5))
    x = 50
    y = 0
    for category in colordict:
        dwg.add(dwg.rect((x-30, y), (60, 10), fill=colordict[category]))
        dwg.add(dwg.text(category, insert=(x, y+5), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill='black'))
        x += 60


    dwg.save()
    #print(highlightsaver) 


    styletext = """<style>
   <![CDATA[
    text.moo {
         font-family: "arial";
         fill: black;
         font-size: 100%;
    }
    rect.hiss {
         fill:white;
    }
   ]]>
   svg text.moo {display: none;}
   svg rect.hiss {display: none;}
   svg g:hover text {display: block;}
   svg g:hover rect {display: block;}
 </style>"""

    imagefile = open(filename,"r")
    data = imagefile.read()
    newdata = data.replace("</svg>", styletext+"</svg>")
    imagefile.close()
    writeFile = open(filename, "w")
    writeFile.write(newdata)
    writeFile.close()

    circletext = ""
    for hlid in highlightsaver:
        cx = highlightsaver[hlid][0]
        cy = highlightsaver[hlid][1]
        txt = highlightsaver[hlid][2]
        circletext = circletext+"""<g xmlns="http://www.w3.org/2000/svg">
          <circle xmlns="http://www.w3.org/2000/svg" cx='"""+str(cx)+"""' cy='"""+str(cy)+"""' r="7" style="fill:transparent;stroke:transparent;stroke-width:0.5;fill-opacity:0.25;stroke-opacity:0.25"/>      
          <rect class="hiss" x='"""+str(cx)+"""' y='"""+str(cy-32)+"""' height='20' width='"""+str(len(txt)*10)+"""'></rect>
          <text class="moo" x='"""+str(cx)+"""' y='"""+str(cy-20)+"""'>"""+txt+"""</text>
          </g>"""

    #print(circletext) 
    imagefile = open(filename,"r")
    imagefile.seek(0)
    data= imagefile.read()
    #print(data)
    imagefile.close()
    newdata = data.replace("</svg>", circletext+"</svg>")
    writeFile = open(filename, "w")
    writeFile.write(newdata)
    writeFile.close()

# ---------------------------------------------------------------------------------------------------------------------------------------------
def main(alignmentfile, protein_of_interest, position_of_interest, window, topguns, positions, path=''):
	# print (positions['Q9NYV4'])
	# sys.exit()
	sequences, trackstart 	= CUSTOM_ALIGN(alignmentfile)
	Konserve = {}
	TheForbiddenPositions = []
	#Konserve, TheForbiddenPositions 	= conservation_checker(protein_of_interest,sequences, trackstart, positions)

	### note: This is the real CONNECTORalpha dictionary I could fetch from uniprot on 31.01.2023

	feature_dict = {}
	try:
		with open(sys.argv[7]) as ff:
			data_feat = ff.read()
		feature_dict = ast.literal_eval(data_feat)
		for k in feature_dict:
			start 	= feature_dict[k][0]
			end	= feature_dict[k][1]
			for i in range(start, end+1):
				if i not in feature_dict[k]:
					feature_dict[k].append(i)
	except:
		pass
	### to define the annotation colors we want to use

	positioncolors = ["lightgreen","salmon","yellow","purple","lightblue"]
	colors = {}
	coloringcategories = []
	for k in positions:
		counter = 0
		for v in positions[k]:
			if v in colors: continue
			if v not in coloringcategories:
				coloringcategories.append(v)
			colors[v]=positioncolors[counter]
			counter+=1
	for seqident in sequences:
		if seqident not in positions:
			positions[seqident]={}
			for colcateg in colors:
				positions[seqident][colcateg]=[]

	featurecolors = ["firebrick","tomato","orange","olive","palegreen","teal","dodgerblue","blueviolet","deeppink"]
	###

	filename = create_svg(sequences, positions, colors, coloringcategories, featurecolors, position_of_interest, window, protein_of_interest, TheForbiddenPositions, feature_dict, trackstart, topguns, path)	### def create_svg(sequences_dict, positions, colors, startposition, windowsize, poi):
	return filename

if __name__ == "__main__":
	alignmentfile = sys.argv[5]	#### change this to the location of the alignmentfile.
	###
	protein_of_interest = sys.argv[1].replace('=', '|')
	# print (protein_of_interest)
	# sys.exit()
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
	main(alignmentfile, protein_of_interest, position_of_interest, window, topguns, positions)
