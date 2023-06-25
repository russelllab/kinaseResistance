#!/usr/bin/env python3

'''
Script to make a database for
kinase mutations project
'''

# import mysql.connector
import os, gzip, sys
from tqdm import tqdm
import psycopg2
# from backports import zoneinfo
import pandas as pd
sys.path.insert(1, '../ML/')
import fetchData

def connection():
    '''Function to connect to postgresql database'''
    '''mydb = mysql.connector.connect(
    host="localhost",
    user="kinase_user",
    password=""
    )'''
    mydb = psycopg2.connect(
                            database = "kinase_project2",
                            user = "gurdeep",
                            password = "hellokitty",
                            host = "localhost",
                            port = "5432")
    return mydb

def create_alignment_table(mycursor)->None:
    aln2pfam = {}
    for line in open('../pfam/humanKinasesHitsSplitTrimmed.hmm', 'r'):
        if len(line.split()) <= 2:
            continue
        if line.split()[-2] == '-' and line.split()[-3] == '-':
            #print (line.split())
            pfam_position = int(line.split()[0])
            aln_position = int(line.split()[-5])
            aln2pfam[aln_position] = pfam_position
    
    dic_aln = {}
    for line in open('../alignments/humanKinasesHitsSplitTrimmed.fasta', 'r'):
        if line[0] == '>':
            name = line.split('>')[1].rstrip().lstrip()
            if name not in dic_aln: dic_aln[name] = ''
            continue
        dic_aln[name] += line.replace('\n', '')
    
    dic = {}
    for name in dic_aln:
        sequence = dic_aln[name]
        # sequence = line.replace('\n', '')
        for num, residue in enumerate(sequence, start=1):
            if num not in dic: dic[num] = []
            dic[num].append(residue)
    data = []
    for num in range(1, len(dic)+1):
        row = []
        row.append(num)
        row.append(aln2pfam[num] if num in aln2pfam else '-')
        row.append(''.join(dic[num]))
        data.append(row)
    df = pd.DataFrame(data, columns=['alnPos', 'pfampos', 'alnAA'])
    # df.to_sql('alignment', con=mycursor, if_exists='replace', index=False)
    mycursor.execute("DROP TABLE IF EXISTS alignment CASCADE")
    mycursor.execute("CREATE TABLE alignment (\
                     alnPos VARCHAR(5) PRIMARY KEY, \
                     pfamPos VARCHAR(10) REFERENCES hmm DEFERRABLE, \
                     alnAA TEXT\
                        )")
    tmp_df = "./tmp_dataframe.csv"
    df.to_csv(tmp_df, index=False, header=False)
    f = open(tmp_df, 'r')
    mycursor.copy_from(f, 'alignment', sep=',')
    

def create_hmm_table(mycursor)->None:
    '''Function to create HMM table'''
    AA = 'ACDEFGHIKLMNPQRSTVWY'
    # mycursor.execute("SET FOREIGN_KEY_CHECKS=0")
    mycursor.execute("DROP TABLE IF EXISTS hmm CASCADE")
    hmm_score_names = ['pfam'+aa + 'VARCHAR(1)' for aa in AA]
    mycursor.execute("CREATE TABLE hmm (\
                     pfamPos VARCHAR(5) PRIMARY KEY, pfamAA VARCHAR(5),\
                     pfamSS VARCHAR(5), alnPos VARCHAR(5),\
                     pfamA FLOAT, pfamC FLOAT, pfamD FLOAT, pfamE FLOAT,\
                     pfamF FLOAT, pfamG FLOAT, pfamH FLOAT, pfamI FLOAT,\
                     pfamK FLOAT, pfamL FLOAT, pfamM FLOAT, pfamN FLOAT,\
                     pfamP FLOAT, pfamQ FLOAT, pfamR FLOAT, pfamS FLOAT,\
                     pfamT FLOAT, pfamV FLOAT, pfamW FLOAT, pfamY FLOAT\
                     )")
    dic_ss = {'G': 1, 'H': 1, 'B': 2, 'C': 3, 'E': 4, 'S': 5, 'T': 6, '-':7}
    hmm = {} # hmmPosition > AA > bit-score
    # for line in open('../pfam/Pkinase.hmm'):
    for line in open('../pfam/humanKinasesHitsSplitTrimmed.hmm', 'r'):
        if len(line.split()) > 2:
            if line.split()[-2] == '-' and line.split()[-3] == '-':
                #print (line.split())
                position = int(line.split()[0])
                aln_position = int(line.split()[-5])
                pfam_aa = line.split()[-4]
                pfam_ss = dic_ss[line.split()[-1].replace('\n', '')]
                hmm[position] = {'ss': pfam_ss}
                for value, aa in zip(line.split()[1:-5], AA):
                    hmm[position][aa] = float(value)
                mycursor.execute("INSERT INTO hmm (pfamPos, pfamAA, \
                                pfamSS, alnPos,\
                                pfamA, pfamC, pfamD, pfamE, pfamF, pfamG, pfamH, \
                                pfamI, pfamK, pfamL, pfamM, pfamN, pfamP, pfamQ, \
                                pfamR, pfamS, pfamT, pfamV, pfamW, pfamY) \
                                VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, \
                                        %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, \
                                        %s, %s)", \
                                (position, pfam_aa, \
                                pfam_ss, aln_position, \
                                hmm[position]['A'], hmm[position]['C'], \
                                hmm[position]['D'], hmm[position]['E'], \
                                hmm[position]['F'], hmm[position]['G'], \
                                hmm[position]['H'], hmm[position]['I'], \
                                hmm[position]['K'], hmm[position]['L'], \
                                hmm[position]['M'], hmm[position]['N'], \
                                hmm[position]['P'], hmm[position]['Q'], \
                                hmm[position]['R'], hmm[position]['S'], \
                                hmm[position]['T'], hmm[position]['V'], \
                                hmm[position]['W'], hmm[position]['Y']))
                                 
            elif line.split()[0] == 'HMM':
                AA = line.replace('\n', '').split()[1:]
    
    # Inserting the '-' position
    mycursor.execute("INSERT INTO hmm (pfamPos, pfamAA, \
                                pfamSS, alnPos,\
                                pfamA, pfamC, pfamD, pfamE, pfamF, pfamG, pfamH, \
                                pfamI, pfamK, pfamL, pfamM, pfamN, pfamP, pfamQ, \
                                pfamR, pfamS, pfamT, pfamV, pfamW, pfamY) \
                                VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, \
                                        %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, \
                                        %s, %s)", \
                                ('-', '-', \
                                '-', '-', \
                                0, 0, \
                                0, 0, \
                                0, 0, \
                                0, 0, \
                                0, 0, \
                                0, 0, \
                                0, 0, \
                                0, 0, \
                                0, 0, \
                                0, 0))
    # return hmm

def create_ptm_table(mycursor)->None:
    '''Function to create the PTM table'''
    mycursor.execute("DROP TABLE IF EXISTS ptms CASCADE")
    mycursor.execute("CREATE TABLE ptms (\
                     uniprotPos INT, uniprotAA VARCHAR(1),\
                     pfamPos VARCHAR(10) REFERENCES hmm DEFERRABLE, \
                     pfamAA VARCHAR(1),\
                     acc VARCHAR(10) REFERENCES kinases(acc) DEFERRABLE, \
                     gene VARCHAR(25), ptmType VARCHAR(10), name VARCHAR(50), \
                     CONSTRAINT name FOREIGN KEY (name) REFERENCES positions(name) \
                     )")
    # for line in tqdm(open('../data/Kinase_psites_hits_split_trimmed.tsv', 'r')):
    for line in tqdm(open('../data/humanKinasesHitsSplitTrimmedPTM.tsv', 'r')):
        if line.startswith('#'): continue
        line = line.rstrip().split('\t')
        acc = line[0]
        ## Neglect the isoforms
        if '-' in acc: continue
        gene = line[1]
        ptmType = line[3].split('-')[1]
        residue = line[3].split('-')[0]
        uniprotPos = int(residue[1:])
        uniprotAA = residue[0]
        pfamAA = (line[5])[0]
        pfamPos = line[4]
        name = acc + '/' + uniprotAA + str(uniprotPos)

        mycursor.execute("select * from positions where name=%s", (name,))
        if mycursor.fetchone() is None:
            print (f'{name} not found in the positions table but in the PTMs file, \
                hence excluded')
            continue

        mycursor.execute("INSERT INTO ptms (\
                            uniprotPos, uniprotAA, pfamPos, pfamAA, acc, gene, ptmType, name) \
                            VALUES (%s, %s, %s, %s, %s, %s, %s, %s)", \
                            (uniprotPos, uniprotAA, pfamPos, pfamAA, acc, gene, ptmType, name))

def fetch_mappings_dic():
    '''Function to create the positions table'''

    '''
    mycursor.execute("DROP TABLE IF EXISTS positions")
    mycursor.execute("CREATE TABLE positions (id SERIAL PRIMARY KEY, \
                     uniprotPos INT, uniprotAA VARCHAR(1),\
                     pfamPos INT, pfamAA VARCHAR(1),\
                     acc VARCHAR(20), uniprot_id VARCHAR(25) \
                     )")
                    # UNIQUE(mutation, wtAA, wtPos, mutAA, mut_type, acc, gene, info, source)\
    '''
    mappings = {}
    # for line in tqdm(gzip.open('../data/humanKinasesHitsHmmsearchMappings.tsv.gz', 'rt')):
    # for line in tqdm(gzip.open('../data/humanKinasesHitsSplitHmmsearchTrimmedMappings.tsv.gz', 'rt')):
    for line in tqdm(gzip.open('../data/humanKinasesHitsSplitTrimmedMappings.tsv.gz', 'rt')):
        if line[0] == '#': continue
        acc = line.split('\t')[0].split('|')[1]
        uniprot_id = line.split('\t')[0].split('|')[2]
        if acc not in mappings: mappings[acc] = {'uniprot_id': uniprot_id, 'positions': {}}
        uniprotAA = line.split('\t')[1]
        uniprotPos = int(line.split('\t')[2])
        pfamAA = line.split('\t')[3]
        pfamPos = line.split('\t')[4].rstrip()
        alnPos = line.split('\t')[5].rstrip()
        mappings[acc]['positions'][uniprotPos] = {'uniprotAA': uniprotAA,
                                                'pfamAA': pfamAA,
                                                'pfamPos': pfamPos,
                                                'alnPos': alnPos}
        
        # print (mutation, wtAA, wtPos, mutAA, mut_type, acc, gene, info, source)
        '''
        mycursor.execute("INSERT INTO positions (uniprotP->dictos, uniprotAA, pfamPos, pfamAA, \
                         acc, uniprot_id) \
                         VALUES (%s, %s, %s, %s, %s, %s)", \
                            (uniprotPos, uniprotAA, pfamPos, pfamAA, acc, uniprot_id)
                            )
        '''
    # Extract information of inactive kinases from the
    # hmmsearch output file
    domain = None
    flag = False
    for line in gzip.open('../data/humanKinasesPkinasePK_Tyr_Ser-ThrAllInactiveHmmsearch.txt.gz', 'rt'):
        if line.rstrip().lstrip() == '': continue
        if line.split()[0] == 'Query:': domain = line.split()[1]
        if domain is None: continue
        if line.startswith('>>'):
            name = line.split()[1]
            acc = name.split('|')[1]
            if acc not in mappings: mappings[acc] = {'uniprot_id': uniprot_id, 'positions': {}}
            uniprot_start = 0 if name.split('|')[3].split('-')[0] == 'start' else int(name.split('|')[3].split('-')[0])-1
            uniprot_id = name.split('|')[2]
            flag = True
            continue
        if flag == False: continue
        if line.split()[0] == domain:
            dom_start = int(line.split()[1])
            dom_end = int(line.split()[3])
            dom_seq = line.split()[2]
        elif line.split()[0] == name:
            if line.split()[1] == '-': continue
            name_start = int(line.split()[1]) + uniprot_start
            name_end = int(line.split()[3]) + uniprot_start
            name_seq = line.split()[2]
            name_counter = 0; dom_counter = 0
            for i in range(len(name_seq)):
                if name_seq[i] not in ['-', '.'] and dom_seq[i] not in ['-', '.']:
                    uniprotPos = name_start + name_counter
                    pfamPos = dom_start + dom_counter
                    uniprotAA = name_seq[i]
                    pfamAA = dom_seq[i]
                    # You need to check if the position is already present
                    # in the dictionary. In that case just skip it                    
                    if uniprotPos in mappings[acc]['positions']: continue
                    # If the uniprotPos is not present in the dictionary
                    # then add it to the dictionary (this means that it is
                    # inactive/pseudokinase)
                    mappings[acc]['positions'][uniprotPos] = {'uniprotAA': uniprotAA,
                                                            'pfamAA': pfamAA,
                                                            'pfamPos': pfamPos,
                                                            'alnPos': '-'}
                    dom_counter += 1; name_counter += 1
                elif name_seq[i] in ['-', '.']: dom_counter += 1
                elif dom_seq[i] in ['-', '.']: name_counter += 1

    return mappings

def find_pfampos(mycursor, acc, uniprotPos)->None:
    mycursor.execute("select pfampos from positions \
                     where acc=%s and uniprotpos=%s", (acc, str(uniprotPos)))
    pfamPos = mycursor.fetchone()
    if pfamPos is not None: pfamPos = pfamPos[0]
    else: pfamPos = '-'
    return pfamPos

def create_mutations_table(mycursor)->None:
    '''Function to create the mutations table'''
    mycursor.execute("DROP TABLE IF EXISTS mutations CASCADE")
    mycursor.execute("CREATE TABLE mutations (id SERIAL PRIMARY KEY, \
                     mutation VARCHAR(10), wtAA VARCHAR(1), wtPos INT, mutAA VARCHAR(1), \
                     pfamPos VARCHAR(10) REFERENCES hmm DEFERRABLE, \
                     mut_type VARCHAR(10), \
                     acc VARCHAR(10) REFERENCES kinases(acc) DEFERRABLE, \
                     gene VARCHAR(10), \
                     info TEXT, pubmed TEXT, source VARCHAR(200) \
                     )")
                    # UNIQUE(mutation, wtAA, wtPos, mutAA, mut_type, acc, gene, info, source)\
    # for line in open('../AK_mut_w_sc_feb2023/act_deact_v2.tsv', 'r'):
    # for line in open('../data/mutations/act_deact_mutlist_may2023.tsv', 'r'):
    for line in gzip.open('../data/mutations/ad_mutations.tsv.gz', 'rt'):
        if line.split()[0] == 'UniProtAcc': continue
        # print (line)
        # gene = line.split('\t')[0]
        acc = line.split('\t')[0]
        gene = line.split('\t')[1]
        mutation = line.split('\t')[2]
        wtAA = mutation[0]
        mutAA = mutation[-1]
        if len(wtAA) > 1 or len(mutAA) > 1: continue
        wtPos = int(mutation[1:-1])
        if fetchData.validateMutation(mycursor, acc, wtPos, wtAA, mutAA) is False:
            # sys.exit()
            continue
        pfamPos = find_pfampos(mycursor, acc, wtPos)
        acc, gene, uniprot_id, protein_name, protein_length = fetchData.getAccGene(mycursor, acc)
        # print (gene, acc, mutation, pfamPos)
        # if pfamPos is None: continue
        mut_type = line.split('\t')[6].lstrip().rstrip()
        if mut_type not in ['increase', 'activating', 'decrease', 'loss']:
            print (mutation, gene, 'is unknown', mut_type)
            continue
        # mut_type = 'A' if mut_type in ['increase', 'activation', 'activating'] else 'D'
        # print (acc, kinases[acc].gene, wtAA, position, mutAA)
        # print (mutation, mut_type, gene)
        # mutation = wtAA + wtPos + mutAA
        info = line.split('\t')[4]
        # info = 'info'
        # source = line.split('\t')[-1]
        source = line.split('\t')[7].rstrip()
        pubmed = line.split('\t')[5]
        # print (mutation, wtAA, wtPos, mutAA, mut_type, acc, gene, info, source)
        mycursor.execute("INSERT INTO mutations (mutation, wtAA, wtPos, mutAA, \
                         pfamPos, mut_type, acc, gene, info, pubmed, source) \
                         VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)", \
                        (mutation, wtAA, wtPos, mutAA, pfamPos, mut_type,\
                         acc, gene, info, pubmed, source))

    '''Fetch resistant mutation data'''
    # for line in open('../AK_mut_w_sc_feb2023/res_mut_v3_only_subs_KD_neighb.tsv', 'r'):
    for line in open('../data/mutations/res_muts_apr_2023_from_cosmic.tsv', 'r'):
        if line.split()[0] == 'Mutation_uniprot': continue
        acc = line.split('/')[0]
        mutation = line.split('/')[1].rstrip()
        acc, gene, uniprot_id, protein_name, protein_length = fetchData.getAccGene(mycursor, acc)
        wtAA = mutation[0]
        mutAA = mutation[-1]
        if mutAA == 'X': continue
        if len(wtAA) > 1 or len(mutAA) > 1: continue
        wtPos = mutation[1:-1]
        wtPos = int(wtPos)
        if fetchData.validateMutation(mycursor, acc, wtPos, wtAA, mutAA) is False:
            # sys.exit()
            continue
        pfamPos = find_pfampos(mycursor, acc, wtPos)
        # mutation = wtAA + wtPos + mutAA
        mut_type = 'resistance'
        source = 'COSMIC'
        info = '-'
        pubmed = '-'
        # if wtPos not in seq2pfam[acc]:
        #     print (f'{uniprot_position} seems to be outside the domain in {acc} and reported {mut_type}')
        #     continue
        mycursor.execute("INSERT INTO mutations (mutation, wtAA, wtPos, mutAA, \
                         pfamPos, mut_type, acc, gene, info, pubmed, source) \
                         VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)", \
                        (mutation, wtAA, wtPos, mutAA, pfamPos, mut_type,\
                         acc, gene, info, pubmed, source))
    
    '''Fetch neutral mutation data'''
    # for line in open('../AK_mut_w_sc_feb2023/nat_mut_tidy_v2_march2023.tsv', 'r'):
    #for line in open('../data/mutations/nat_mut_tidy_all_samples_tidy_2.tsv', 'r'):
    # for line in open('../data/mutations/20230525_AllNeutrals_1pSNP_5pHom.txt', 'r'):
    for line in gzip.open('../data/mutations/neutral_nodis_AFgeq0.001_homgeq2.txt.gz', 'rt'):
        # if line.split('\t')[0] == 'Uniprot': continue
        acc = line.split()[14].lstrip().rstrip()
        acc, gene, uniprot_id, protein_name, protein_length = fetchData.getAccGene(mycursor, acc)
        mutation = line.split()[20].lstrip().rstrip()
        wtAA = mutation[0]
        mutAA = mutation[-1]
        if mutAA == 'X': continue
        if len(wtAA) > 1 or len(mutAA) > 1: continue
        wtPos = mutation[1:-1]
        if wtPos.isdigit() == False: continue
        wtPos = int(wtPos)
        if fetchData.validateMutation(mycursor, acc, wtPos, wtAA, mutAA) is False:
            # sys.exit()
            continue
        # print (acc, wtPos)
        pfamPos = find_pfampos(mycursor, acc, wtPos)
        mut_type = 'neutral'
        source = 'gnomAD'
        # info = 'AC/AN:'+str(line.split('\t')[7]) + '; Hom/AC:'+str(line.split('\t')[8].rstrip())
        info = '-'
        pubmed = '-'
        # if acc not in seq2pfam:
        #     continue
        # if uniprot_position not in seq2pfam[acc]:
        #     print (f'{uniprot_position} seems to be outside the domain and reported {mut_type}')
        #     print (seq2pfam[acc])
        #     continue
        mycursor.execute("INSERT INTO mutations (mutation, wtAA, wtPos, mutAA, \
                         pfamPos, mut_type, acc, gene, info, pubmed, source) \
                         VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)", \
                        (mutation, wtAA, wtPos, mutAA, pfamPos, mut_type, acc,\
                         gene, info, pubmed, source))

def create_ligands_table(mycursor)->None:
    '''Function to create the ligands table'''
    mycursor.execute("DROP TABLE IF EXISTS ligands CASCADE")
    mycursor.execute("CREATE TABLE ligands (id SERIAL PRIMARY KEY, \
                     ligand VARCHAR(100), ligand_id VARCHAR(100), uniprotpos INT, \
                     pfamPos VARCHAR(10) REFERENCES hmm DEFERRABLE, \
                     acc VARCHAR(10) REFERENCES kinases(acc) DEFERRABLE, \
                     gene VARCHAR(10), \
                     PDB TEXT, pubmed TEXT \
                     )")
                    # UNIQUE(mutation, wtAA, wtPos, mutAA, mut_type, acc, gene, info, source)\
    for line in gzip.open('../data/ATP_binding_sites/ATP_binding_sites_UniProt.tsv.gz', 'rt'):
        if line.split()[0] == 'UniProtAcc': continue
        acc = line.split('\t')[0]
        gene = line.split('\t')[1]
        uniprotpos = line.split('\t')[2]
        pfamPos = find_pfampos(mycursor, acc, uniprotpos)
        acc, gene, uniprot_id, protein_name, protein_length = fetchData.getAccGene(mycursor, acc)
        ligand = line.split('\t')[3]
        ligand_id = line.split('\t')[4]
        if ligand_id == '': ligand_id = '-'
        pdbs = line.split('\t')[6]
        if pdbs == '': pdbs = '-'
        pubmeds = line.split('\t')[7]
        if pubmeds == '': pubmeds = '-'
        # print (ligand, ligand_id, uniprotpos, pfamPos, acc, gene, pdbs, pubmeds)
        mycursor.execute("INSERT INTO ligands (ligand, ligand_id, uniprotpos, \
                         pfamPos, acc, gene, PDB, pubmed) \
                         VALUES (%s, %s, %s, %s, %s, %s, %s, %s)", \
                        (ligand, ligand_id, uniprotpos, pfamPos, acc, gene, pdbs, pubmeds))

def createDicForDSSP(dic, position, mutation, value):
    if position not in dic: dic[position] = {}
    dic[position][mutation] = float(value)

def create_homology_table(mycursor) -> None:
    '''Function to create the homology  table'''
    # path = '/net/home.isilon/ds-russell/mechismoX/analysis/alignments/data/HUMAN/orthologs_only/'
    path = '../data/homFiles/'
    mycursor.execute("DROP TABLE IF EXISTS homology CASCADE")
    mycursor.execute('select acc from kinases')
    accs = mycursor.fetchall()
    for fileEnd in tqdm([
                    '_all_homs.scores.txt.gz',
                    '_orth.scores.txt.gz',
                    '_excl_para.scores.txt.gz',
                    '_spec_para.scores.txt.gz',
                    '_bpso.scores.txt.gz',
                    '_bpsh.scores.txt.gz'
                    ]):
        homology = fileEnd.split('.scores')[0]
        homology = homology[1:]
        mycursor.execute("DROP TABLE IF EXISTS "+homology+" CASCADE")
        AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        mycursor.execute("CREATE TABLE "+homology+" (\
                    acc VARCHAR(10) REFERENCES kinases(acc) DEFERRABLE, \
                    wtaa VARCHAR(5), position INT, \
                    A_score FLOAT, C_score FLOAT, D_score FLOAT, E_score FLOAT, \
                    F_score FLOAT, G_score FLOAT, H_score FLOAT, I_score FLOAT, \
                    K_score FLOAT, L_score FLOAT, M_score FLOAT, N_score FLOAT, \
                    P_score FLOAT, Q_score FLOAT, R_score FLOAT, S_score FLOAT, \
                    T_score FLOAT, V_score FLOAT, W_score FLOAT, Y_score FLOAT, \
                    info TEXT) \
                    ")
        data = []
        num = 0
        for acc_tuple in tqdm(accs):
            dic = {}
            num += 1
            # if num == 10: break
            acc = acc_tuple[0]
            if os.path.isfile(path + '/' + acc + fileEnd) is False:
                print (path + '/' + acc + fileEnd, 'does not exist')
                continue
            for line in gzip.open(path + '/' + acc + fileEnd, 'rt'):
                #print (acc, line.split())
                #sys.exit()
                mutation = line.split()[0].split('/')[1]
                position = int(mutation[1:-1])
                wtaa = mutation[0]
                mutaa = mutation[-1]
                wtscore = float(line.split()[2])
                mutscore = float(line.split()[3])
                diffscore = float(line.split()[4])
                # info = line.split()[5].rstrip()
                info = '-'
                if position not in dic: dic[position] = {'wtaa': wtaa, 'info': info}
                dic[position][mutaa+'_score'] = mutscore
            for position in range(1, len(dic)+1):
                row = []
                try:
                    row = [acc, dic[position]['wtaa'], position]
                except:
                    print (acc, position, fileEnd, 'does not exist')
                    sys.exit()
                for aa in AA:
                    if aa+'_score' in dic[position]:
                        row.append(dic[position][aa+'_score'])
                    else:
                        row.append(None)
                row.append(dic[position]['info'] if 'info' in dic[position] else None)
                data.append(row)
        if len(data) == 0: continue
        df = pd.DataFrame(data)
        # print (df)
        tmp_df = "./tmp_dataframe.csv"
        df.to_csv(tmp_df, index=False, header=False)
        f = open(tmp_df, 'r')
        mycursor.copy_from(f, homology, sep=',')

def create_taxon_tables(mycursor) -> None:
    '''Function to create the homology  table'''
    # path = '/net/home.isilon/ds-russell/mechismoX/analysis/alignments/data/HUMAN/orthologs_only/'
    path = '../data/taxFiles/'
    mycursor.execute("DROP TABLE IF EXISTS homology CASCADE")
    mycursor.execute('select acc from kinases')
    accs = mycursor.fetchall()
    for taxon in tqdm([
                    'arthropods',
                    'eukaryotes',
                    'mammals',
                    'metazoa',
                    'vertebrates'
                        ]):
    # for fileEnd in tqdm([
    #                 '_all_homs.scores.txt.gz',
    #                 '_orth.scores.txt.gz',
    #                 '_excl_para.scores.txt.gz',
    #                 '_spec_para.scores.txt.gz',
    #                 '_bpso.scores.txt.gz',
    #                 '_bpsh.scores.txt.gz'
    #                 ]):
        mycursor.execute("DROP TABLE IF EXISTS "+taxon+" CASCADE")
        AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        mycursor.execute("CREATE TABLE "+taxon+" (\
                    acc VARCHAR(10) REFERENCES kinases(acc) DEFERRABLE, \
                    wtaa VARCHAR(5), position INT, \
                    A_score FLOAT, C_score FLOAT, D_score FLOAT, E_score FLOAT, \
                    F_score FLOAT, G_score FLOAT, H_score FLOAT, I_score FLOAT, \
                    K_score FLOAT, L_score FLOAT, M_score FLOAT, N_score FLOAT, \
                    P_score FLOAT, Q_score FLOAT, R_score FLOAT, S_score FLOAT, \
                    T_score FLOAT, V_score FLOAT, W_score FLOAT, Y_score FLOAT, \
                    info TEXT) \
                    ")
        data = []
        num = 0
        for acc_tuple in tqdm(accs):
            dic = {}
            num += 1
            # if num == 10: break
            acc = acc_tuple[0]
            if os.path.isfile(path + '/' + taxon + '/' + acc + '.txt.gz') is False:
                print (path + '/' + taxon + '/' + acc + '.txt.gz', 'does not exist')
                continue
            for line in gzip.open(path + '/' + taxon + '/' + acc + '.txt.gz', 'rt'):
                #print (acc, line.split())
                #sys.exit()
                mutation = line.split()[0].split('/')[1]
                position = int(mutation[1:-1])
                wtaa = mutation[0]
                mutaa = mutation[-1]
                wtscore = float(line.split()[2])
                mutscore = float(line.split()[3])
                diffscore = float(line.split()[4])
                # info = line.split()[5].rstrip()
                info = '-'
                if position not in dic: dic[position] = {'wtaa': wtaa, 'info': info}
                dic[position][mutaa+'_score'] = mutscore
            for position in range(1, len(dic)+1):
                row = []
                try:
                    row = [acc, dic[position]['wtaa'], position]
                except:
                    print (acc, position, taxon, acc, 'does not exist')
                    sys.exit()
                for aa in AA:
                    if aa+'_score' in dic[position]:
                        row.append(dic[position][aa+'_score'])
                    else:
                        row.append(None)
                row.append(dic[position]['info'] if 'info' in dic[position] else None)
                data.append(row)
        if len(data) == 0: continue
        df = pd.DataFrame(data)
        # print (df)
        tmp_df = "./tmp_dataframe.csv"
        df.to_csv(tmp_df, index=False, header=False)
        f = open(tmp_df, 'r')
        mycursor.copy_from(f, taxon, sep=',')

def create_iupred_table(mycursor) -> None:
    '''Function to create the iupred table'''
    # path = '/net/home.isilon/ds-russell/mechismoX/analysis/alignments/data/HUMAN/orthologs_only/'
    path = '../data/strucFiles/'
    mycursor.execute("DROP TABLE IF EXISTS structure CASCADE")
    mycursor.execute('select acc from kinases')
    accs = mycursor.fetchall()
    fileEnd = 'iupred.gz'
    strucType = fileEnd.split('.gz')[0]
    mycursor.execute("DROP TABLE IF EXISTS "+strucType+" CASCADE")
    AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    mycursor.execute("CREATE TABLE "+strucType+" (\
                acc VARCHAR(10) REFERENCES kinases(acc) DEFERRABLE, \
                wtaa VARCHAR(5), position INT, \
                A_score FLOAT, C_score FLOAT, D_score FLOAT, E_score FLOAT, \
                F_score FLOAT, G_score FLOAT, H_score FLOAT, I_score FLOAT, \
                K_score FLOAT, L_score FLOAT, M_score FLOAT, N_score FLOAT, \
                P_score FLOAT, Q_score FLOAT, R_score FLOAT, S_score FLOAT, \
                T_score FLOAT, V_score FLOAT, W_score FLOAT, Y_score FLOAT, \
                info TEXT) \
                ")
    data = []
    num = 0
    for acc_tuple in tqdm(accs):
        dic = {}
        num += 1
        # if num == 10: break
        acc = acc_tuple[0]
        filePath = path + '/AF-' + acc + '-F1-model_v1.' + fileEnd
        if os.path.isfile(filePath) is False:
            print (filePath, 'does not exist')
            continue
        for line in gzip.open(filePath, 'rt'):
            #print (acc, line.split())
            #sys.exit()
            # mutation = line.split()[0].split('/')[1]
            position = int(line.split()[0])
            wtaa = line.split()[2]
            mutaa = line.split()[3]
            wtscore = float(line.split()[5])
            mutscore = float(line.split()[7])
            diffscore = float(line.split()[9].rstrip())
            # info = line.split()[5].rstrip()
            info = '-'
            if position not in dic: dic[position] = {'wtaa': wtaa, 'info': info}
            dic[position][mutaa+'_score'] = mutscore
        for position in range(1, len(dic)+1):
            row = []
            row = [acc, dic[position]['wtaa'], position]
            for aa in AA:
                if aa+'_score' in dic[position]:
                    row.append(dic[position][aa+'_score'])
                else:
                    row.append(None)
            row.append(dic[position]['info'] if 'info' in dic[position] else None)
            data.append(row)
    # if len(data) == 0: continue
    df = pd.DataFrame(data)
    # print (df)
    tmp_df = "./tmp_dataframe.csv"
    df.to_csv(tmp_df, index=False, header=False)
    f = open(tmp_df, 'r')
    mycursor.copy_from(f, strucType, sep=',')

def create_dssp_tables(mycursor) -> None:
    '''Function to create dssp tables'''
    # path = '/net/home.isilon/ds-russell/mechismoX/analysis/alignments/data/HUMAN/orthologs_only/'
    path = '../data/strucFiles/'
    mycursor.execute("DROP TABLE IF EXISTS structure CASCADE")
    mycursor.execute('select acc from kinases')
    accs = mycursor.fetchall()
    fileEnd = 'dssp-scores.gz'
    AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    data_phi_psi = []
    data_sec = []
    data_acc = []
    data_burr = []
    num = 0
    for acc_tuple in tqdm(accs):
        dic_phi_psi = {}
        dic_sec = {}
        dic_acc = {}
        dic_burr = {}
        num += 1
        # if num == 10: break
        acc = acc_tuple[0]
        filePath = path + '/AF-' + acc + '-F1-model_v1.' + fileEnd
        if os.path.isfile(filePath) is False:
            print (filePath, 'does not exist')
            continue
        for line in gzip.open(filePath, 'rt'):
            if line.split()[0] == 'length': continue
            #sys.exit()
            # mutation = line.split()[0].split('/')[1]
            position = int(line.split()[0])
            wtaa = line.split()[2]
            mutaa = line.split('PHI/PSI')[0].split()[-1]
            content_after_phi_psi = line.split('PHI/PSI')[1]
            # print (acc, line.split())
            for scores, dic in [[(4,5),dic_phi_psi],
                                [(8,9), dic_sec],
                                [(12,13), dic_burr],
                                [(16,17), dic_acc]
                                ]:
                wtscore, mutscore = float(content_after_phi_psi.split()[scores[0]]), float(content_after_phi_psi.split()[scores[1]])
                if position not in dic: dic[position] = {'wtaa': wtaa}
                dic[position][mutaa+'_score'] = mutscore
        for data, dic in [
                            [data_phi_psi, dic_phi_psi],
                            [data_sec, dic_sec],
                            [data_burr, dic_burr],
                            [data_acc, dic_acc]
                        ]:
            for position in range(1, len(dic)+1):
                row = []
                row = [acc, dic[position]['wtaa'], position]
                for aa in AA:
                    if aa+'_score' in dic[position]:
                        row.append(dic[position][aa+'_score'])
                    else:
                        row.append(None)
                row.append(dic[position]['info'] if 'info' in dic[position] else None)
                data.append(row)
    for data, dsspType in [
                    [data_phi_psi, 'phi_psi'],
                    [data_sec, 'sec'],
                    [data_burr, 'burr'],
                    [data_acc, 'acc']
                    ]:
        mycursor.execute("DROP TABLE IF EXISTS "+dsspType+" CASCADE")
        mycursor.execute("CREATE TABLE "+dsspType+" (\
                    acc VARCHAR(10) REFERENCES kinases(acc) DEFERRABLE, \
                    wtaa VARCHAR(5), position INT, \
                    A_score FLOAT, C_score FLOAT, D_score FLOAT, E_score FLOAT, \
                    F_score FLOAT, G_score FLOAT, H_score FLOAT, I_score FLOAT, \
                    K_score FLOAT, L_score FLOAT, M_score FLOAT, N_score FLOAT, \
                    P_score FLOAT, Q_score FLOAT, R_score FLOAT, S_score FLOAT, \
                    T_score FLOAT, V_score FLOAT, W_score FLOAT, Y_score FLOAT, \
                    info TEXT) \
                    ")
        # if len(data) == 0: continue
        df = pd.DataFrame(data)
        # print (df)
        tmp_df = "./tmp_dataframe.csv"
        df.to_csv(tmp_df, index=False, header=False)
        f = open(tmp_df, 'r')
        mycursor.copy_from(f, dsspType, sep=',')

def create_mechismo_table(mycursor) -> None:
    '''Function to create the mechismo table'''
    # path = '/net/home.isilon/ds-russell/mechismoX/analysis/alignments/data/HUMAN/orthologs_only/'
    path = '../data/strucFiles/'
    mycursor.execute("DROP TABLE IF EXISTS mechismo CASCADE")
    mycursor.execute('select acc from kinases')
    accs = mycursor.fetchall()
    fileEnd = 'mech_intra.gz'
    strucType = fileEnd.split('.gz')[0]
    mycursor.execute("DROP TABLE IF EXISTS "+strucType+" CASCADE")
    AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    mycursor.execute("CREATE TABLE "+strucType+" (\
                acc VARCHAR(10) REFERENCES kinases(acc) DEFERRABLE, \
                wtaa VARCHAR(5), position INT, \
                ncontacts INT, nresidues INT, \
                A_score FLOAT, C_score FLOAT, D_score FLOAT, E_score FLOAT, \
                F_score FLOAT, G_score FLOAT, H_score FLOAT, I_score FLOAT, \
                K_score FLOAT, L_score FLOAT, M_score FLOAT, N_score FLOAT, \
                P_score FLOAT, Q_score FLOAT, R_score FLOAT, S_score FLOAT, \
                T_score FLOAT, V_score FLOAT, W_score FLOAT, Y_score FLOAT, \
                info TEXT) \
                ")
    data = []
    num = 0
    for acc_tuple in tqdm(accs):
        dic = {}
        num += 1
        # if num == 10: break
        acc = acc_tuple[0]
        filePath = path + '/AF-' + acc + '-F1-model_v1.' + fileEnd
        if os.path.isfile(filePath) is False:
            print (filePath, 'does not exist')
            continue
        for line in gzip.open(filePath, 'rt'):
            if line.split()[0] == 'CO':
                position = int(line.split()[1])
                wtaa = line.split()[2]
                ncontacts = int(line.split()[3])
                nresidues = int(line.split()[4])
                if position not in dic: dic[position] = {'wtaa': wtaa}
                dic[position]['ncontacts'] = ncontacts
                dic[position]['nresidues'] = nresidues
                # dic[position][mutaa+'_score'] = mutscore
            elif line.split()[0] == 'MECH':
                position = int(line.split()[1])
                wtaa = line.split()[2]
                mutaa = line.split()[3]
                wtscore = float(line.split()[4])
                mutscore = float(line.split()[5])
                if len(line.split()) == 8:
                    info = line.split()[7].rstrip()
                else:
                    info = '-'
                if position not in dic: dic[position] = {'wtaa': wtaa, 'info': info}
                if 'info' not in dic[position]: dic[position]['info'] = info
                dic[position][mutaa+'_score'] = mutscore
                # dic[position][mutaa+'_score'] = mutscore
        for position in range(1, len(dic)+1):
            row = []
            row = [
                    acc,
                    dic[position]['wtaa'],
                    position,
                    dic[position]['ncontacts'],
                    dic[position]['nresidues']
                ]
            for aa in AA:
                if aa+'_score' in dic[position]:
                    row.append(dic[position][aa+'_score'])
                else:
                    row.append(None)
            row.append(dic[position]['info'] if 'info' in dic[position] else None)
            data.append(row)
    # if len(data) == 0: continue
    df = pd.DataFrame(data)
    # print (df)
    tmp_df = "./tmp_dataframe.csv"
    df.to_csv(tmp_df, index=False, header=False)
    f = open(tmp_df, 'r')
    mycursor.copy_from(f, strucType, sep=',')

'''
def fetch_map_fasta2aln():
    # Function to map from fasta to alignment
    dic_aln = {}
    for line in open('../alignments/humanKinasesHitsSplitTrimmed.fasta', 'r'):
        if line[0] == '>':
            name = line.split('>')[1].rstrip().lstrip()
            if name not in dic_aln: dic_aln[name] = ''
            continue
        dic_aln[name] += line.replace('\n', '')

    map_fasta2aln = {}
    # for line in open('../alignments/humanKinasesHitsSplitTrimmed.fasta', 'r'):
    for name in dic_aln:
        sequence = dic_aln[name]
        # if line.startswith('>'):
        # acc = line.split('|')[1]
        acc = name.split('|')[1]
        # print (acc)
        # start = line.split('|')[3].split('-')[0]
        start = name.split('|')[3].split('-')[0]
        if start == 'start': start = 0
        else: start = int(start)
        # start += int(line.split('|')[4].split()[0])
        counter = int(name.split('|')[4].split()[0])
        # if counter != 1:
        start += counter
        # continue
        if acc not in map_fasta2aln: map_fasta2aln[acc] = {}
        current_position = start
        for alnpos, aa in enumerate(sequence.rstrip(), start=1):
            if aa not in ['-', '.']:
                map_fasta2aln[acc][current_position] = alnpos
                current_position += 1
    print (map_fasta2aln['P08581'][1100])
    # sys.exit()
    return map_fasta2aln
'''

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

def fetch_pfam_domains():
    path2fasta = '../KA/UniProtFasta2/'
    dic_kinase = {}
    for line in open('../data/humanKinases.fasta', 'r'):
        if line[0] != '>': continue
        acc = line.split('|')[1]
        if '-' in acc: continue
        if acc not in dic_kinase: dic_kinase[acc] = []
        for line in gzip.open(path2fasta + acc +'.txt.gz', 'rt'):
            if line.startswith('DR') == False: continue
            if line.split(';')[0].split()[1] != 'Pfam': continue
            if line.split(';')[1].split()[0] in ['PF00069', 'PF07714', 'PF00454',\
                                                'PF02518', 'PF02816', 'PF06743',\
                                                'PF01163', 'PF10494', 'PF03109']:
                dic_kinase[acc].append(line.split(';')[1].split()[0])

    for acc in dic_kinase:
        if len(dic_kinase[acc]) == 0:
            if acc in ['O60885', 'Q9UIG0', 'P21675', 'Q13263',\
                    'O15164', 'Q9NRL2', 'P53004', 'Q9Y5P4',\
                        'Q5VZY9', 'Q9UPN9', 'P11274', 'Q8NI60',\
                        'Q12979', 'Q12979', 'Q15059', 'Q58F21',\
                        'P25440', 'Q58F21', 'Q8IZX4']:
                dic_kinase[acc] = '-'
            else:
                dic_kinase[acc] = runHmmsearch(acc, path2fasta)
    return dic_kinase

def create_kinases_table(mycursor)->None:
    '''Function to create the kinases table'''
    mycursor.execute("DROP TABLE IF EXISTS kinases CASCADE")
    mycursor.execute("CREATE TABLE kinases (\
                     acc VARCHAR(10) PRIMARY KEY, \
                     gene VARCHAR(10), uniprot_id VARCHAR(25), \
                     pfam_domain VARCHAR(50), protein_name VARCHAR(100), \
                     fasta TEXT) \
                     ")
                    #  UNIQUE(acc, gene, uniprot_id, fasta))\
    mycursor.execute("DROP TABLE IF EXISTS positions CASCADE")
    mycursor.execute("CREATE TABLE positions (\
                     uniprotPos INT, uniprotAA VARCHAR(1), alnPos VARCHAR(50),\
                     pfamPos VARCHAR(10) REFERENCES hmm DEFERRABLE, \
                     pfamAA VARCHAR(1), \
                     acc VARCHAR(10) REFERENCES kinases DEFERRABLE, \
                     uniprot_id VARCHAR(25), \
                     name VARCHAR(50) PRIMARY KEY \
                     )")
    # CONSTRAINT pfamPos FOREIGN KEY (pfamPos) REFERENCES hmm(pfamPos) DEFERRABLE\
                     
    kinases = {}
    for line in open('../data/humanKinases.fasta', 'r'):
        if line[0] == '>':
            # print (line)
            acc = line.split('|')[1]
            uniprot_id = line.split('|')[2]
            gene = line.split('GN=')[1].split()[0]
            protein_name = line.split('_HUMAN')[1].split('OS=')[0].rstrip().lstrip()
            name = acc + '|' + uniprot_id + '|' + gene + '|' + protein_name
            kinases[name] = ''
        else:
            kinases[name] += line.rstrip()

    pfam_domains = fetch_pfam_domains()
    mappings = fetch_mappings_dic()
    # map_fasta2aln = fetch_map_fasta2aln()
    num = 0
    data = []
    for kinase in tqdm(kinases):
        num += 1
        # if num == 3:
        #    break
        acc = kinase.split('|')[0]
        uniprot_id = kinase.split('|')[1].split()[0]
        gene = kinase.split('|')[2]
        protein_name = kinase.split('|')[3]
        fasta = str(kinases[kinase])
        pfam_domain = ';'.join(pfam_domains[acc])
        # print (acc, pfam_domain)
        # print (acc, gene, uniprot_id, fasta)
        mycursor.execute("INSERT INTO kinases (acc, gene, uniprot_id, pfam_domain, protein_name, fasta) \
                            VALUES (%s, %s, %s, %s, %s, %s)", \
                            (acc, gene, uniprot_id, pfam_domain, protein_name, fasta))
        # if acc not in mappings: continue
        # if acc not in map_fasta2aln: continue
        for uniprotPos, uniprotAA in enumerate(fasta, start=1):
            pfamPos, pfamAA, alnPos = '-', '-', '-'
            if acc in mappings:
                if uniprotPos in mappings[acc]['positions']:
                    if mappings[acc]['positions'][uniprotPos]['uniprotAA'] != uniprotAA:
                        print ('ERROR', acc, gene, uniprotPos, uniprotAA, mappings[acc]['positions'][uniprotPos]['uniprotAA'])
                        sys.exit()
                    pfamPos = mappings[acc]['positions'][uniprotPos]['pfamPos']
                    pfamAA = mappings[acc]['positions'][uniprotPos]['pfamAA']
                    alnPos = mappings[acc]['positions'][uniprotPos]['alnPos']
                    # print (uniprotPos, uniprotAA, pfamPos, pfamAA, acc, uniprot_id)
            
            # alnPos = '-'
            # if acc in map_fasta2aln:
            #     if uniprotPos in map_fasta2aln[acc]:
            #         alnPos = map_fasta2aln[acc][uniprotPos]
                
            name = acc+'/'+uniprotAA+str(uniprotPos)
            # if acc == 'P08581' and int(uniprotPos) == 1100:
            #     print (name, uniprotPos, uniprotAA, alnPos, pfamPos, pfamAA, acc, uniprot_id)
            #     sys.exit()
            row = []
            row.append(uniprotPos)
            row.append(uniprotAA)
            row.append(alnPos)
            row.append(pfamPos)
            row.append(pfamAA)
            row.append(acc)
            row.append(uniprot_id)
            row.append(name)
            data.append(row)
            '''
            mycursor.execute("INSERT INTO positions (uniprotPos, uniprotAA, alnPos, pfamPos, pfamAA, \
                                acc, uniprot_id, name) \
                                VALUES (%s, %s, %s, %s, %s, %s, %s, %s)", \
                                (uniprotPos, uniprotAA, alnPos, pfamPos, pfamAA, acc, uniprot_id, name)
                                )
            '''
    df = pd.DataFrame(data, columns=['uniprotPos', 'uniprotAA', 'alnPos', 'pfamPos', 'pfamAA', 'acc', 'uniprot_id', 'name'])
    # print (df)
    tmp_df = "./tmp_dataframe.csv"
    df.to_csv(tmp_df, index=False, header=False)
    f = open(tmp_df, 'r')
    mycursor.copy_from(f, 'positions', sep=',')

if __name__ == '__main__':
    mydb = connection()
    mydb.autocommit = True
    mycursor = mydb.cursor()
    # Execute SQL query to check if database exists
    # mycursor.execute("SHOW DATABASES")
    # sql = '''CREATE DATABASE kinase_project''';
    # mycursor.execute(sql)

    '''
    # Loop through results and check if database exists
    db_name = 'kinase_project'
    db_exists = False
    for db in mycursor:
        if db[0] == db_name:
            db_exists = True
            break
    if db_exists == False: mycursor.execute("CREATE DATABASE "+db_name)
    mycursor.execute("use "+db_name)
    '''
    # Create tables
    # print ('Creating HMM table')
    # create_hmm_table(mycursor)
    # print ('Creating kinases table')
    # create_kinases_table(mycursor)
    # print ('Creating mutation table')
    # create_mutations_table(mycursor)
    # print ('Creating ligands table')
    # create_ligands_table(mycursor)
    # print ('Creating homology table')
    # create_homology_table(mycursor)
    print ('Creating taxon tables')
    create_taxon_tables(mycursor)
    # print ('Creating IUPRED table')
    # create_iupred_table(mycursor)
    # print ('Creating Mechismo table')
    # create_mechismo_table(mycursor)
    # print ('Creating DSSP tables')
    # create_dssp_tables(mycursor)
    # print ('Creating PTM table')
    # create_ptm_table(mycursor)
    # print ('Creating alignment table')
    # create_alignment_table(mycursor)
    mydb.commit()

    # Use mysqldump to create backup file
    # backup_file = "kinaseDB.sql"
    # os.system(f"mysqldump -u {mydb.user} {mydb.database} > {backup_file}")

    # Close MySQL connection
    mydb.close()
    # main()
