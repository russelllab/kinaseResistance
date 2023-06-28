#!/usr/bin/env python
# coding: utf-8

import os, sys, gzip
from tqdm import tqdm
import psycopg2

'''
List of functions that fetch data from
the database
'''

PTM_TYPES = ['ac', 'gl', 'm1', 'm2', 'm3', 'me', 'p', 'sm', 'ub']
AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',\
      'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

def connection(db_name='kinase_project'):
    '''Function to connect to postgresql database'''
    mydb = psycopg2.connect(
                            database = db_name,
                            user = "gurdeep",
                            password = "hellokitty",
                            host = "localhost",
                            port = "5432")
    return mydb

def checkGZfile(filename):
    '''Function to check if the file is gzipped'''
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    else:
        return open(filename, 'r')

def mutTypes(mycursor, acc, mutation):
    '''Function to fetch mutation types from the DB'''
    mycursor.execute("select mut_type from mutations where \
                     acc=%s and mutation=%s", (acc, mutation,))
    hits = mycursor.fetchall()
    mut_types = [hit[0] for hit in hits]
    return mut_types

def getAlnPos(mycursor, hmmpos):
    '''Function to fetch aln pos from the DB for a given hmmpos'''
    mycursor.execute("select alnpos from hmm where \
                     pfampos=%s", (hmmpos,))
    hits = mycursor.fetchone()
    if hits is None:
        return '-'
    else:
        return hits[0]

def getRegion(mycursor, acc, mutation):
    '''Function to fetch region from the DB'''
    mycursor.execute("select pfampos from positions where \
                     acc=%s and uniprotpos=%s", (acc, mutation[1:-1],))
    hits = mycursor.fetchone()
    if hits is None:
        return '-'
    else:
        alnpos = hits[0]
        # print (alnpos)
        if alnpos == '-': return '-'
        else: alnpos = int(alnpos)
        region = []
        with gzip.open('../alignments/humanKinasesHitsSplitTrimmed_ss.tsv.gz', 'rt', encoding='utf-8') as regionFile:
            for line in regionFile:
                if line.startswith('#'): continue
                region_name = line.split()[0]
                start, end = int(line.split()[1].split('-')[0]), int(line.split()[1].split('-')[1])
                if alnpos in range(start, end+1):
                    region.append(region_name)
        
        # if no specific region found 
        # return a hyphen
        if len(region) == 0:
            return '-'
        else:
            return ';'.join(region)

def getSSdic():
    dic_region = {}
    with gzip.open('../alignments/humanKinasesHitsSplitTrimmed_ss.tsv.gz', 'rt', encoding='utf-8') as regionFile:
        for line in regionFile:
            if line.startswith('#'): continue
            region_name = line.split()[0]
            start, end = int(line.split()[1].split('-')[0]), int(line.split()[1].split('-')[1])
            for position in range(start, end+1):
                if position not in dic_region: dic_region[position] = []
                dic_region[position].append(region_name)
    return dic_region

def getRegion2(mycursor, acc, mutation, dic_region):
    '''Function to fetch region from the DB'''
    mycursor.execute("select pfampos from positions where \
                     acc=%s and uniprotpos=%s", (acc, mutation[1:-1],))
    hits = mycursor.fetchone()
    if hits is None:
        return '-'
    else:
        alnpos = hits[0]
        # print (alnpos)
        if alnpos == '-': return '-'
        else: alnpos = int(alnpos)
        region = []
        if alnpos in dic_region:
            region = dic_region[alnpos]
        
        # if no specific region found 
        # return a hyphen
        if len(region) == 0:
            return '-'
        else:
            return ';'.join(region)

def getAccGene(mycursor, name):
    '''Function to fetch acc and gene from the DB'''
    check_with = ['acc', 'gene', 'uniprot_id']
    for check in check_with:
        mycursor.execute("select acc, gene, uniprot_id, protein_name, fasta\
                         from kinases where "+check+" = %s", (name,))
        hits = mycursor.fetchone()
        if hits is not None: break
    if hits is None:
        print (f'Neither acc nor gene with name {name} found')
        return None, None, None, None, None
    acc, gene, uniprot_id, protein_name, protein_length = hits[0], hits[1], hits[2], hits[3], len(hits[4])
    return acc, gene, uniprot_id, protein_name, protein_length

def validateMutation(mycursor, acc, wtPos, wtAA, mutAA):
    mycursor.execute("select fasta from kinases where \
                        acc=%s", (acc,))
    hits = mycursor.fetchone()
    if hits is None:
        print (f'No sequence found for {acc}')
        return False
    fasta = hits[0]
    try:
        if wtPos > len(fasta):
            print (f'wtPos {wtPos} is greater than length of fasta {len(fasta)} of {acc}')
            return False
    except:
        print (f'wtPos {wtPos} is not an integer {acc}')
        return False
    if wtAA != fasta[wtPos-1]:
        print (f'wtAA {wtAA} does not match with fasta {fasta[wtPos-1]} at position {wtPos} of {acc}')
        print (fasta[wtPos-3:wtPos+3])
        return False
    return True

def checkInputPositionAA(acc, mutation, mycursor):
	'''
	Checks if the input position and amino acid are valid
	'''
	mycursor.execute("select uniprotaa from positions \
		  				where acc=%s and uniprotpos=%s", (acc, mutation[1:-1],))
	hits = mycursor.fetchone()
	# print (hits)
	if hits is None:
		return 1, 'Position not found'
	if hits[0] != mutation[0]:
		return 2, hits[0]
	return 0, 'OK'

def getAdjacentSites(mycursor, acc, position, WS):
    mycursor.execute("select fasta from kinases where acc=%s", (acc,))
    hits = mycursor.fetchone()
    if hits is None:
        return None
    else:
        seq = hits[0]
        startSite = position-WS-1
        if startSite < 0: startSite = 0
        endSite = position+WS+1
        if endSite > len(seq): endSite = len(seq)
        return seq[startSite:position-1]+seq[position-1].lower()+seq[position+1:endSite]

def retrieve_entries(mycursor, acc, mutation):
	'''For a give acc/gene/uniprot ID, retireve known information'''
	## fetch ptm_types
	mycursor.execute(\
					'select ptmtype from ptms \
					where acc=%s and uniprotpos=%s', \
					(acc, mutation[1:-1],)\
					)
	hits = mycursor.fetchone()
	if hits is not None: ptmType = hits[0]
	else: ptmType = 'None'
	
	## fetch mut_types
	mycursor.execute(\
					'select mut_type from mutations \
					where acc=%s and mutation=%s', \
					(acc, mutation,)\
					)
	hits = mycursor.fetchone()
	if hits is not None: mutType = hits[0]
	else: mutType = 'None'

	print (ptmType, mutType)
	return (ptmType, mutType)

def fetchFasta(kinases, Kinase, mycursor):
    '''
    for line in open('../data/humanKinases.fasta', 'r'):
        #print (line)
        if line[0] == '>':
            acc = line.split('|')[1].replace('\n', '')
            gene = line.split('GN=')[1].split()[0]
            kinases[acc] = Kinase(acc, gene)
            # flag = 0
            # if acc not in exceptions:
            #     kinases[acc] = kinase(acc, gene)
            #     flag = 1
        else:
            # if flag == 1:
            kinases[acc].fasta += line.replace('\n', '')
    '''
    mycursor.execute("select acc, gene, fasta from kinases")
    for acc, gene, fasta in mycursor.fetchall():
        # print (acc, gene, fasta)
        kinases[acc] = Kinase(acc, gene)
        kinases[acc].fasta = fasta

def fetchGroup(kinases, Kinase):
    for line in open('../data/kinases.tsv', 'r'):
        acc = line.split('\t')[7].split('>')[1].split('<')[0]
        if acc in kinases:
            kinases[acc].group = line.split('\t')[4]

def fetchPkinaseHMM(mycursor):
    dic_ss = {'G': 1, 'H': 1, 'B': 2, 'C': 3, 'E': 4, 'S': 5, 'T': 6, '-':7}
    hmm = {} # hmmPosition > AA > bit-score
    '''
    # for line in open('../pfam/Pkinase.hmm'):
    for line in open('../pfam/humanKinasesHitsSplitTrimmed.hmm'):
        if len(line.split()) > 2:
            if line.split()[-2] == '-' and line.split()[-3] == '-':
                #print (line.split())
                position = int(line.split()[0])
                ss = dic_ss[line.split()[-1].replace('\n', '')]
                hmm[position] = {'ss': ss}
                for value, aa in zip(line.split()[1:-5], AA):
                    hmm[position][aa] = float(value)
            elif line.split()[0] == 'HMM':
                AA = line.replace('\n', '').split()[1:]
    '''
    mycursor.execute("select * from hmm")
    for row in mycursor.fetchall():
        pfampos = row[0]
        if pfampos == '-': continue
        pfamaa = row[1]
        for aa, bitscore in zip(AA, row[2:-1]):
            if pfampos not in hmm: hmm[pfampos] = {'ss': 'ss'}
            hmm[pfampos][aa] = float(bitscore)
    return hmm

def fetchHmmsearch(kinases, Kinase):
    '''
    Function to do an hmmsearch of all kinases against Pkinase.hmm
    and store the mappings. Note that some kinases may have more than
    one Pkinase domain.
    '''
    # os.system('hmmsearch -o out.txt ../pfam/Pkinase.hmm ../data/humanKinases.fasta')
    os.system('hmmsearch -o out.txt ../pfam/humanKinasesHitsSplitTrimmed.hmm\
                                    ../data/humanKinases.fasta')
    flag = 0
    with open('out.txt', 'r') as hmmsearchFile:
        for line in hmmsearchFile:
            if line[:2] == '>>':
                acc = line.split('|')[1]
                flag = 1
                #print (acc)
            if flag == 1 and line.split()!= [] and acc in kinases:
                if '== domain' in line:
                    domainNum = int(line.split('domain')[1].split()[0])
                    kinases[acc].domains[domainNum] = {}
                    kinases[acc].seq2pfam[domainNum] = {}
                elif line.split()[0] == 'humanKinasesHitsSplitTrimmed':
                    hmmStart = int(line.split()[1])
                    hmmSeq = line.split()[2]
                    hmmEnd = int(line.split()[3])
                elif acc in line.split()[0]:
                    kinaseStart = line.split()[1]
                    if kinaseStart == '-': continue
                    kinaseStart = int(kinaseStart)
                    kinaseSeq = line.split()[2]
                    kinaseEnd = int(line.split()[3])
                    for hmmChar, kinaseChar in zip(hmmSeq, kinaseSeq):
                        if hmmChar not in ['.', '-'] and kinaseChar not in ['.', '-']:
                            #kinases[acc].domains[domainNum][kinaseStart] = hmmStart
                            kinases[acc].domains[domainNum][hmmStart] = kinaseStart
                            kinases[acc].seq2pfam[domainNum][kinaseStart] = hmmStart
                            hmmStart += 1
                            kinaseStart += 1
                        elif hmmChar in ['.', '-']:
                            kinaseStart += 1
                        elif kinaseChar in ['.', '-']:
                            hmmStart += 1
            #print (kinases[acc].domains)
            #sys.exit()
    # print (kinases['Q96NX5'].domains)
    # print (kinases['Q96NX5'].domains[1][1384])

def createDicForDSSP(dic, position, mutation, value):
    if position not in dic: dic[position] = {}
    dic[position][mutation] = float(value)

def dsspScores(kinases, Kinase):
    remove = []
    dir = '/net/home.isilon/ds-russell/mechismoX/analysis/features/data/VLatest/'
    for num, acc in enumerate(kinases):
        for line in gzip.open(dir + acc[:4] + '/AF-' + acc + '-F1-model_v1.dssp-scores.gz', 'rt'):
            if line.split()[0] == 'length':
                continue
            #print (acc, line.split())
            position = int(line.split()[0])
            mutation = line.split()[2] + line.split()[0] + line.split()[10]
            ## Dihedral angles
            torsional = line.split()[18]
            #print (kinases[acc].dihedral)
            createDicForDSSP(kinases[acc].dihedral, position, mutation, torsional)
            ## Secondary structures
            secondary = line.split()[22]
            createDicForDSSP(kinases[acc].sec, position, mutation, secondary)
            ## Accessibility
            accessibility = line.split()[26]
            createDicForDSSP(kinases[acc].access, position, mutation, accessibility)
            ## Buried
            buried = line.split()[30]
            createDicForDSSP(kinases[acc].burr, position, mutation, buried)
            #break
        #break
    print (list(set(remove)))

def iupredScores(kinases, Kinase):
    remove = []
    dir = '/net/home.isilon/ds-russell/mechismoX/analysis/features/data/VLatest/'
    for num, acc in enumerate(kinases):
        if ((num+1)%50 == 0):
            print (num+1)
        for line in gzip.open(dir + acc[:4] + '/AF-' + acc + '-F1-model_v1.iupred.gz', 'rt'):
            #print (acc, line.split())
            position = int(line.split()[0])
            mutation = line.split()[2] + line.split()[0] + line.split()[3]
            ## IUPred
            iupred = float(line.split()[9])
            createDicForDSSP(kinases[acc].iupred, position, mutation, iupred)


def mechismoScores(kinases, Kinase):
    remove = []
    dir = '/net/home.isilon/ds-russell/mechismoX/analysis/features/data/VLatest/'
    for num, acc in enumerate(kinases):
        if ((num+1)%50 == 0):
            print (num+1)
        for line in gzip.open(dir + acc[:4] + '/AF-' + acc + '-F1-model_v1.mech_intra.gz', 'rt'):
            if line.split()[0] == 'MECH':
                #print (acc, line.split())
                #sys.exit()
                position = int(line.split()[1])
                mutation = line.split()[2] + line.split()[1] + line.split()[3]
                ## Mechismo score
                mechismo = float(line.split()[6])
                createDicForDSSP(kinases[acc].mechismo, position, mutation, mechismo)

def homologyScores(kinases, Kinase):
    remove = []
    path = '/net/home.isilon/ds-russell/mechismoX/analysis/alignments/data/HUMAN/orthologs_only/'
    for acc in tqdm(kinases):
        for dic, fileEnd in zip([
                            kinases[acc].allHomologs,
                            kinases[acc].orthologs,
                            kinases[acc].exclParalogs,
                            kinases[acc].specParalogs,
                            kinases[acc].bpso,
                            kinases[acc].bpsh],
                            [
                            '_all_homs.scores.txt.gz',
                             '_orth.scores.txt.gz',
                             '_excl_para.scores.txt.gz',
                             '_spec_para.scores.txt.gz',
                             '_bpso.scores.txt.gz',
                             '_bpsh.scores.txt.gz'
                             ]):
            if os.path.isfile(path + acc[:4] + '/' + acc + fileEnd) is False:
                print (path + acc[:4] + '/' + acc + fileEnd, 'does not exist')
                continue
            for line in gzip.open(path + acc[:4] + '/' + acc + fileEnd, 'rt'):
                #print (acc, line.split())
                #sys.exit()
                value = line.split()[0].split('/')[1]
                position = int(value[1:-1])
                residue = value[-1]
                #print (mutation, position)
                ## Mechismo score
                score = float(line.split()[4])
                createDicForDSSP(dic, position, residue, score)
                #if acc == 'Q9NYV4' and position == 877:
                #    print (dic[position])

def getHomologyScores(mycursor, acc, wtAA, position, mutAA):
    row = []
    for position in range(position-2, position+3):
        for homology in ['all_homs','excl_para','spec_para','orth','bpso','bpsh']:
            mycursor.execute("SELECT * FROM "+homology+" \
                            WHERE acc=%s and position=%s", (acc, str(position),))
            hit = mycursor.fetchone()
            if hit is None:
                row = [] # empty row if no hits found
                break
            # print (hit)
            acc = hit[0]
            AA = 'ACDEFGHIKLMNPQRSTVWY'
            homology_score_wt, homology_score_mut = None, None
            for logodd, aa in zip(hit[3:-1], AA):
                if aa == mutAA:
                    # row.append(float(logodd))
                    homology_score_mut = float(logodd)
                if aa == wtAA:
                    homology_score_wt = float(logodd)
            if homology_score_wt is None or homology_score_mut is None:
                return None
            row.append(homology_score_mut - homology_score_wt)
    return row

def getTaxonScores(mycursor, acc, wtAA, position, mutAA):
    '''
    Returns a list of taxon scores for a given position and mutation
    '''
    row = []
    for position in range(position-2, position+3):
        for taxon in [
                    'eukaryotes',
                    'mammals',
                    'metazoa',
                    'vertebrates'
                    ]:
            mycursor.execute("SELECT * FROM "+taxon+" \
                            WHERE acc=%s and position=%s", (acc, str(position),))
            hit = mycursor.fetchone()
            if hit is None:
                # row = [] # empty row if no hits found
                return None
            # print (hit)
            acc = hit[0]
            AA = 'ACDEFGHIKLMNPQRSTVWY'
            taxon_score_wt, taxon_score_mut = None, None
            for logodd, aa in zip(hit[3:-1], AA):
                if aa == mutAA:
                    # row.append(float(logodd))
                    taxon_score_mut = float(logodd)
                if aa == wtAA:
                    taxon_score_wt = float(logodd)
            if taxon_score_wt is None or taxon_score_mut is None:
                return None
            row.append(taxon_score_mut - taxon_score_wt)
    return row

def getIUPredScore(mycursor, acc, wtAA, position, mutAA):
    iupred_score_wt, iupred_score_mut = None, None
    mycursor.execute("SELECT * FROM iupred \
                        WHERE acc=%s and position=%s", (acc, str(position),))
    hit = mycursor.fetchone()
    if hit is not None:
        acc = hit[0]
        AA = 'ACDEFGHIKLMNPQRSTVWY'
        for logodd, aa in zip(hit[3:-1], AA):
            if aa == mutAA:
                # row.append(float(logodd))
                iupred_score_mut = float(logodd)
            if aa == wtAA:
                iupred_score_wt = float(logodd)
    # return row
    if iupred_score_wt is None or iupred_score_mut is None:
        return None
    return iupred_score_mut - iupred_score_wt

def getMechIntraScores(mycursor, acc, wtAA, position, mutAA):
    row = []
    mycursor.execute("SELECT * FROM mech_intra \
                        WHERE acc=%s and position=%s", (acc, str(position),))
    hit = mycursor.fetchone()
    if hit is not None:
        acc = hit[0]
        ncontacts = hit[3]
        nresidues = hit[4]
        row = [int(ncontacts), int(nresidues)]
        AA = 'ACDEFGHIKLMNPQRSTVWY'
        mech_score_wt, mech_score_mut = None, None
        for logodd, aa in zip(hit[5:-1], AA):
            if aa == mutAA:
                # row.append(float(logodd))
                mech_score_mut = float(logodd)
            if aa == wtAA:
                mech_score_wt = float(logodd)
        row.append(mech_score_mut - mech_score_wt)
    if len(row) == 0: row = None
    return row

def getDSSPScores(mycursor, acc, wtAA, position, mutAA):
    row = []
    dsspTypes = ['phi_psi', 'sec', 'burr', 'acc']
    for dsspType in dsspTypes:
        mycursor.execute("SELECT * FROM "+dsspType+" \
                            WHERE acc=%s and position=%s", (acc, str(position),))
        hit = mycursor.fetchone()
        if hit is not None:
            acc = hit[0]
            AA = 'ACDEFGHIKLMNPQRSTVWY'
            dssp_score_wt, dssp_score_mut = None, None
            for logodd, aa in zip(hit[3:-1], AA):
                if aa == mutAA:
                    # row.append(float(logodd))
                    dssp_score_mut = float(logodd)
                if aa == wtAA:
                    dssp_score_wt = float(logodd)
            row.append(dssp_score_mut - dssp_score_wt)
        else:
            # You can't take any entry if at least one of the dsspTypes is missing
            row = []
            break
    if len(row) == 0: row = None
    return row

def getHmmPkinaseScore(mycursor, acc, wtAA, position, mutAA):
    # print (f'HMMscore in {acc} for {wtAA}{position}{mutAA}')
    mycursor.execute("SELECT pfampos, alnpos FROM positions \
                     WHERE acc = %s and uniprotpos = %s", (acc, str(position)))
    hits = mycursor.fetchone()
    # print (hits)
    # print (acc, wtAA, position, mutAA, hits)
    if hits == None:
        print (acc, wtAA, position, mutAA, hits)
        return '-', None, None, None, None
    pfampos = hits[0]
    alnpos = hits[1]
    # print (f'pfampos of {acc}/{wtAA}{position}{mutAA} is {pfampos}')
    mycursor.execute("SELECT * FROM hmm \
                        WHERE pfampos = %s", (str(pfampos),))
    row = mycursor.fetchone()
    
    AA = 'ACDEFGHIKLMNPQRSTVWY'
    for aa, bitscore in zip(AA, row[4:]):
        if aa == mutAA:
            mut_bitscore = bitscore
        elif aa == wtAA:
            wt_bitscore = bitscore
    return alnpos, pfampos, wt_bitscore, mut_bitscore, 0

def getATPbindingScores(mycursor, acc, position):
    count_atp = []
    mycursor.execute("SELECT pfampos FROM positions \
                     WHERE acc = %s and uniprotpos = %s", (acc, str(position)))
    hits = mycursor.fetchone()
    # print (acc, wtAA, position, mutAA, hits)
    if hits == None:
        print (acc, position, hits)
        return 0
    pfampos = hits[0]
    # print (f'pfampos of {acc}/{wtAA}{position}{mutAA} is {pfampos}')
    mycursor.execute("SELECT ligand, acc FROM ligands \
                        WHERE pfampos = %s", (str(pfampos),))
    hits = mycursor.fetchall()
    if hits is None:
        return count_atp
    for hit in hits:
        ligand, kinase = hit
        if 'ATP' in ligand: count_atp.append(kinase)
    count_atp = list(set(count_atp))
    count_atp = 1 if len(count_atp) > 0 else 0
    # return len(count_atp)
    return count_atp

def getPTMscore(mycursor, acc, mutation_position, ws=0):
    if ws > 0: ws -= 1
    ws = int(ws/2)
    # print (acc, kinases[acc].gene, mutation_position)
    
    row = []
    for position in range(mutation_position-ws, mutation_position+ws+1):
        mycursor.execute("SELECT ptmtype FROM ptms \
                        WHERE acc = %s and uniprotpos = %s", (acc, str(position)))
        hit = mycursor.fetchone()
        if hit is not None:
            ptm_type_at_position = hit[0]
        else:
            ptm_type_at_position = 'None'
        ## prepare vector for known information
        for ptm_type in PTM_TYPES:
            if ptm_type_at_position == 'None':
                row.append('0')
            elif ptm_type == ptm_type_at_position:
                row.append('1')
            else:
                row.append('0')
            # if ptm_type not in kinases[acc].ptm:
            #     row.append('0')
            # elif position in kinases[acc].ptm[ptm_type]:
            #     row.append('1')
            # else:
            #     row.append('0')
        
        ## prepare vector for inference
        mycursor.execute("SELECT pfampos FROM ptms \
                        WHERE acc = %s and uniprotpos = %s", (acc, str(position)))
        hit = mycursor.fetchone()
        if hit is not None:
            hmmPos = hit[0]
        else:
            hmmPos = 'None'
        if hmmPos not in ['-', 'None']:
            mycursor.execute("SELECT ptmtype FROM ptms \
                            WHERE pfampos = %s", (str(hmmPos),))
            hits = mycursor.fetchall()
            ptms_type_at_hmmPos = []
            for hit in hits:
                ptms_type_at_hmmPos.append(hit[0])
        # hmm_position = kinases[acc].returnhmmPos(position)
        if hmmPos in ['-', 'None']:
            for ptm_type in PTM_TYPES:
                row.append('0')
        else:
            for ptm_type in PTM_TYPES:
                count_ptm_type = ptms_type_at_hmmPos.count(ptm_type)
                row.append( '0' if count_ptm_type==0 else str(count_ptm_type) )
                # row.append( '0' if count_ptm_type<3 else '1' )
    
    # if row.count('1')>=5:
    #     print (row)
    #     sys.exit()
    return row

def getCountAAchange(mycursor, acc, mutation_position, kinases, ws=0):
    if ws > 0: ws -= 1
    ws = int(ws/2)
    mut_types = ['A', 'D', 'R']
    dic_mut_types = {'A': ['activating', 'increase'], 'D': ['loss', 'decrease'], 'R': ['resistance']}
    # mut_types = ['activating', 'increase', 'loss', 'decrease', 'resistance']
    adr_row = []
    for position in range(mutation_position-ws, mutation_position+ws+1):
        '''mycursor.execute("SELECT mut_type FROM mutations \
                            WHERE acc = %s and wtpos = %s", (acc, str(position)))
        hits = mycursor.fetchall()
        if hits is None:
            adr_row.extend([0, 0, 0])
            continue
        mut_types_at_position = []
        for entry in hits:
            mut_types_at_position.append(entry[0])
        # print (acc, position, mut_types_at_position)
        for mut_type in mut_types:
            categories = dic_mut_types[mut_type]
            # flag = 0
            count = 0
            for category in categories:
                if category in mut_types_at_position:
                    # adr_row.append(1)
                    # flag = 1
                    # break
                    count += 1
            # if flag == 0:
            #     adr_row.append(0)
            adr_row.append(count)'''

        mycursor.execute("SELECT pfampos FROM positions \
                            WHERE acc = %s and uniprotpos = %s", (acc, str(position)))
        hmmPos = mycursor.fetchone()
        if hmmPos is None:
            adr_row.extend([0, 0, 0])
            continue
        else:
            hmmPos = hmmPos[0]
        mycursor.execute("SELECT wtaa, mutaa, mut_type FROM mutations \
                        WHERE pfampos = %s", (str(hmmPos),))
        hits = mycursor.fetchall()
        # print (hits)
        mut_types_at_hmmpos = []
        wtaas_at_hmmpos = []
        mutaas_at_hmmpos = []
        for entry in hits:
            wtaas_at_hmmpos.append(entry[0])
            mutaas_at_hmmpos.append(entry[1])
            mut_types_at_hmmpos.append(entry[2])
        for mut_type in mut_types:
            dic_count = {}
            categories = dic_mut_types[mut_type]
            for wtaa, mutaa, mut_type_at_hmmpos in zip(wtaas_at_hmmpos, mutaas_at_hmmpos, mut_types_at_hmmpos):
                # flag = 0
                count = 0
                for category in categories:
                    if category == mut_type_at_hmmpos:
                        # print (wtaa, mutaa, mut_type_at_hmmpos, category)
                        if wtaa not in dic_count: dic_count[wtaa] = 0
                        if mutaa not in dic_count: dic_count[mutaa] = 0
                        dic_count[wtaa] -= 1
                        dic_count[mutaa] += 1
            
            for aa in AA:
                if aa not in dic_count:
                    adr_row.append(0)
                else:
                    # print (aa, dic_count[aa], mut_type, dic_count)
                    adr_row.append(dic_count[aa])
            # adr_row.append(count)
            # if mut_type in mut_types_at_hmmpos: adr_row.append(1)
            # else: adr_row.append(0)
    # print (acc, kinases[acc].gene, mutation_position, adr_row)
    # sys.exit()

    return adr_row

def getADRvector(mycursor, acc, mutation_position, kinases, ws=0):
    if ws > 0: ws -= 1
    ws = int(ws/2)
    mut_types = ['A', 'D', 'R']
    dic_mut_types = {'A': ['activating', 'increase'], 'D': ['loss', 'decrease'], 'R': ['resistance']}
    # mut_types = ['activating', 'increase', 'loss', 'decrease', 'resistance']
    adr_row = []
    for position in range(mutation_position-ws, mutation_position+ws+1):
        mycursor.execute("SELECT mut_type FROM mutations \
                            WHERE acc = %s and wtpos = %s", (acc, str(position)))
        hits = mycursor.fetchall()
        if hits is None:
            adr_row.extend([0, 0, 0])
            continue
        mut_types_at_position = []
        for entry in hits:
            mut_types_at_position.append(entry[0])
        # print (acc, position, mut_types_at_position)
        for mut_type in mut_types:
            categories = dic_mut_types[mut_type]
            # flag = 0
            count = 0
            for category in categories:
                count += mut_types_at_position.count(category)
                # if category in mut_types_at_position:
                #     # adr_row.append(1)
                #     # flag = 1
                #     # break
                #     count += 1
            # if flag == 0:
            #     adr_row.append(0)
            adr_row.append(count)

        mycursor.execute("SELECT pfampos FROM positions \
                            WHERE acc = %s and uniprotpos = %s", (acc, str(position)))
        hmmPos = mycursor.fetchone()
        if hmmPos is None:
            adr_row.extend([0, 0, 0])
            continue
        else:
            hmmPos = hmmPos[0]
        mycursor.execute("SELECT mut_type FROM mutations \
                        WHERE pfampos = %s", (str(hmmPos),))
        hits = mycursor.fetchall()
        mut_types_at_hmmpos = []
        for entry in hits:
            mut_types_at_hmmpos.append(entry[0])
        for mut_type in mut_types:
            categories = dic_mut_types[mut_type]
            # flag = 0
            count = 0
            for category in categories:
                count += mut_types_at_hmmpos.count(category)
                # if category in mut_types_at_hmmpos:
                #     print (acc, position, hmmPos, mut_type, category)
                #     # adr_row.append(1)
                #     # flag = 1
                #     # break
                #     count += 1
            # if flag == 0:
            #     adr_row.append(0)
            adr_row.append(count)
            # if mut_type in mut_types_at_hmmpos: adr_row.append(1)
            # else: adr_row.append(0)
    # print (acc, kinases[acc].gene, mutation_position, adr_row)
    # sys.exit()

    return adr_row

def getAAvector(wtAA, mutAA):
    row = []
    for amino_acid in [wtAA, mutAA]:
        for aa in 'ACDEFGHIKLMNPQRSTVWY':
        # for aa in 'DEKQRSTY':
            if aa == amino_acid:
                row.append('1')
            else:
                row.append('0')
    return row
