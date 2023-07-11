import psycopg2
import logging



### fetching uniprotIDs and genenames from the database
### storing everything in 2x lists and 1x dict
conn = psycopg2.connect(database="kinase_project2", 
	host="pevolution2.bioquant.uni-heidelberg.de",
	user="gurdeep",
	password="hellokitty")

cursor = conn.cursor()
commander_farsight = "select gene, acc from kinases;"
cursor.execute(commander_farsight)
results = cursor.fetchall()

genes = []
accessions = []
genemapper = {}
drugcounter = 0

for result in results:
	genes.append(result[0])
	accessions.append(result[1])
	if result[1] not in genemapper:
		genemapper[result[1]]=result[0]

### 520 kinases



### crossing our kinase set with what is present in the drug targets file
### looking for cases where the appoval status says "Y"[es] and then check if the target uniprotIDs match to our set of kinases
### likewise if the status neither says "Y"[es] not "N"[o] then we assume the drug is in development for the respective targets



# counting approved drugs & their targets
approved_drugs = []
approved_targets = []
# counting in-development drugs & their targets
development_drugs = []
development_targets = [] 	

with open("pkidb_2023-06-30.tsv","r") as drugfile:	### NOTE: the file has a line-break error that one has to fix manually after downloading, no idea what they are doing...
	for line in drugfile:
		if "BrandName" not in line:
			drugcounter += 1
			try:	
				drug		= line.split("\t")[1].replace("\n","")
				raw_targets 	= line.split("\t")[22].replace("\n","")
				status 		= line.split("\t")[33].replace("\n","")
				if "Y" in status.replace(" ",""):
					for item in raw_targets.split("|"):
						if item.replace(" ","") in genemapper:
							genus = genemapper[item.replace(" ","")]
							if genus not in approved_targets:	### to prevent counting targets several times
								approved_targets.append(genus)	
					if drug not in approved_drugs:
						approved_drugs.append(drug)
				elif "N" in status.replace(" ",""):
					pass
				else:
					for item in raw_targets.split("|"):
						if item.replace(" ","") in genemapper:
							genus = genemapper[item.replace(" ","")]
							if genus not in development_targets:	### to prevent counting targets several times
								development_targets.append(genus)	
					if drug not in development_drugs:
						development_drugs.append(drug)
			except:
				print(line)
				logging.exception("message")
				pass


### removing drugtargets from the in-development container if there is already an approved drug available
for drugtarget in approved_targets:
	try:
		development_targets.remove(drugtarget)
	except:
		pass

print("Number of approved drug targets", len(approved_targets))
print("Number of approved drugs", len(approved_drugs))
print("Number of drug targets in development", len(development_targets))
print("Number of drugs in development", len(development_drugs))
print("The number of drugs in the file is", drugcounter)

#Number of approved drug targets 294
#Number of approved drugs 66
#Number of drug targets in development 7
#Number of drugs in development 10
#The number of drugs in the file is 369





