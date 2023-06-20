import os, sys, gzip

'''
Takes the human kinases fasta file, and finds the best domain hit
for each kinase. The best domain hit is the one with the highest
score.
'''

class Kinase:
    def __init__(self, name):
        self.name = name
        self.seq = ''
        self.domains = []
        self.inactive = False

def runHmmsearch(acc, path2fasta):
    pdomains = ['Pkinase', 'PK_Tyr_Ser-Thr']
    row = []
    for pdomain in pdomains:
        os.system('hmmsearch --tblout hh '+ ' ../pfam/' + pdomain +'.hmm '\
                  + path2fasta + acc + '.fasta.gz')
        for line in open('hh', 'r'):
            if line.startswith('#'): continue
            if line.split()[0] == acc: continue
            print (line.split())
            # print (line.split('\t')[4], acc)
            if float(line.split()[4]) > 1e-5: continue
            row.append(line.split()[3])
    if len(row) == 0:
        print (acc, 'nothing works')
        sys.exit()
    else:
        return row

def getKinaseDomain(acc):
    domains = []
    
    if acc in ['O60885', 'Q9UIG0', 'P21675', 'Q13263',\
               'O15164', 'Q9NRL2', 'P53004', 'Q9Y5P4',\
               'Q5VZY9', 'Q9UPN9', 'P11274', 'Q8NI60',\
               'Q12979', 'Q12979', 'Q15059', 'Q58F21',\
               'P25440', 'Q58F21', 'Q8IZX4']:
        domains = []
    else:
        path2fasta = '../KA/UniProtFasta2/'
        for line in gzip.open(path2fasta + acc +'.txt.gz', 'rt'):
            if line.startswith('DR') == False: continue
            if line.split(';')[0].split()[1] != 'Pfam': continue
            if line.split(';')[1].split()[0] in list_kinase_domains:
                domains.append(line.split(';')[1].split()[0].split('.')[0])
        if domains == []:
            domains = runHmmsearch(acc, path2fasta)
    return domains

inactive_domains = ['O43187', 'O15197', 'Q13308', 'Q58A45',\
                    'Q5JZY3', 'Q6VAB6', 'Q7RTN6', 'Q7Z7A4',\
                    'Q86YV5', 'Q8IVT5', 'Q8IWB6', 'Q8IZE3',\
                    'Q8NB16', 'Q8NCB2', 'Q8NE28', 'Q8TEA7',\
                    'Q92519', 'Q96C45', 'Q96KG9', 'Q96RU7',\
                    'Q96RU8', 'Q9C0K7', 'Q9H792', 'Q9NSY0',\
                    'Q9UHY1', 'Q9Y616']
pseudokinases = ['Q58A45', 'Q96S38_1', 'P29597_1', 'P23458_1', 'P52333_1', 'O60674_1']
list_kinase_domains = []
for line in open('kinase_domains.txt', 'r'):
    list_kinase_domains.append(line.rstrip())
kinases = {}
flag = 0
for line in open(sys.argv[1], 'r'):
    if line[0] == '>':
        kinase = line.split('>')[1].split()[0].strip()
        kinases[kinase] = Kinase(kinase)
        acc = kinase.split('|')[1]
        kinases[kinase].domains = getKinaseDomain(acc)
        flag = 1
        ## Mark the inactive kinases
        if line.split('>')[1].split()[1].strip() == 'Inactive':
            kinases[kinase].inactive = True
            continue
        ## Mark the one that have catalytically inactive domain/pseduokinases
        ## It's hard to know which one is the catalytically inactive domain
        ## but try to make a list as you go along (see above)
        ## Pseudokinases with multiple domains are treated later
        ## after running the hmmsearch
        if line.split('|')[1] in inactive_domains or line.split('|')[1] in pseudokinases:
            kinases[kinase].inactive = True
            continue
    else:
        kinases[kinase].seq += line.strip()

# Create a fasta file per domain
domain2name = {'PF00069': 'Pkinase', 'PF07714': 'PK_Tyr_Ser-Thr'}
for domain in list_kinase_domains:
    if domain not in ['PF00069', 'PF07714']: continue
    domainName = domain2name[domain]
    fasta = ''
    os.system('hmmsearch --domE 1.0 --noali\
              --domtblout domains/humanKinases'+domainName+'Hits.hmmsearch \
              ../pfam/'+domainName+'.hmm '+\
              sys.argv[1])

class KinaseDomain:
    def __init__(self, name):
        self.name = name
        self.domainEvalues = {}

kinaseDomains = {}
pseudokinasesHits = []
for domain in list_kinase_domains:
    if domain not in ['PF00069', 'PF07714']: continue
    domainName = domain2name[domain]
    hits = []
    for line in open('domains/humanKinases'+domainName+'Hits.hmmsearch', 'r'):
        if line.startswith('#'):
            continue
        total_hits = int(line.split()[10])
        eValue = float(line.split()[6])
        name = line.split()[0]
        if total_hits == 1:
            hit = name + '|' + 'start' + '-' + 'end'
            if name not in kinaseDomains: kinaseDomains[name] = KinaseDomain(name)
            if domain not in kinaseDomains[name].domainEvalues: kinaseDomains[name].domainEvalues[domain]=[]
            kinaseDomains[name].domainEvalues[domain] += [eValue]
            hits.append(hit)
            continue
        current_hit = int(line.split()[9])
        
        current_start = str(line.split()[17])
        current_end = str(line.split()[18])
        if current_hit == 1:
            hit = name + '|' + 'start' + '-' + current_end
        elif current_hit == total_hits:
            hit = name + '|' + current_start + '-' + 'end'
        else:
            hit = name + '|' + current_start + '-' + current_end
        if name not in kinaseDomains: kinaseDomains[name] = KinaseDomain(name)
        if domain not in kinaseDomains[name].domainEvalues: kinaseDomains[name].domainEvalues[domain]=[]
        kinaseDomains[name].domainEvalues[domain] += [eValue]
        hits.append(hit)
        ## Save names of pseudokinases which multiple domains
        ## The ones with single domain were already marked as inactive
        ## before running the hmmsearch
        acc = line.split()[0].split('|')[1]
        ignore = False
        for entry in pseudokinases:
            if '_' not in entry: continue
            pseudokinase = entry.split('_')[0]
            dom = entry.split('_')[1]
            if acc == pseudokinase and int(dom) == int(current_hit):
                # print (acc)
                pseudokinasesHits.append(hit)
                # ignore = True
                break
        # if ignore: continue
        ##

# print (hits)

# for kinase in kinases:
#     found = 0
#     for hit in hits:
#         name = hit.split('|')[0]+'|'+hit.split('|')[1]+'|'+hit.split('|')[2].lstrip().rstrip()
#         if kinase == name:
#             found = 1
#             break
#     if found == 0:
#         print (kinase+'+')

# sys.exit()

def getLikelyDomain(name):
    maxEvalue = None
    likelyDomain = None
    for domain in kinaseDomains[name].domainEvalues:
        if domain not in ['PF00069', 'PF07714']: continue
        for eValue in kinaseDomains[name].domainEvalues[domain]:
            if maxEvalue == None:
                maxEvalue = eValue
                likelyDomain = domain
            elif eValue < maxEvalue:
                maxEvalue = eValue
                likelyDomain = domain
    # print (maxEvalue, likelyDomain, hit, kinaseDomains[hit].domainEvalues)
    return likelyDomain

for domain in list_kinase_domains:
    if domain not in ['PF00069', 'PF07714']: continue
    domainName = domain2name[domain]
    l = ''; m = ''
    for hit in hits:
        kinase = '|'.join(hit.split('|')[:-1])
        # fetch the best domain
        likelyDomain = getLikelyDomain(kinase)
        # ignore if the best domain is not the one we are looking for
        if likelyDomain != domain: continue
        start = hit.split('|')[-1].split('-')[0]
        end = hit.split('|')[-1].split('-')[1]
        if start == 'start' and end == 'end':
            seq = kinases[kinase].seq
        elif start == 'start' and end != 'end':
            end = int(end)
            seq = kinases[kinase].seq[:end]
        elif start != 'start' and end == 'end':
            start = int(start)-1
            seq = kinases[kinase].seq[start:]
        elif start != 'start' and end != 'end':
            start = int(start)-1
            end = int(end)
            seq = kinases[kinase].seq[start:end]
        # print ('>'+hit+'\n'+seq)
        ## Eliminate the ones that are not present in the kinases dictionary
        ## These would be the inactive kinases
        if kinases[kinase].inactive is False and hit not in pseudokinasesHits:
            l += '>'+hit+'\n'+seq+'\n'
            continue
        m += '>'+hit+'\n'+seq+'\n'

    open('humanKinases'+domainName+'HitsSplit.fasta', 'w').write(l)
    open('humanKinases'+domainName+'HitsSplitInactive.fasta', 'w').write(m)