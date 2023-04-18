#!/usr/bin/env python3

'''
Script to make a database for
kinase mutations project
'''

import mysql.connector
import os, gzip, sys
from tqdm import tqdm
import psycopg2

def connection():
    '''Function to connect to postgresql database'''
    '''mydb = mysql.connector.connect(
    host="localhost",
    user="kinase_user",
    password=""
    )'''
    mydb = psycopg2.connect(
                            database = "kinase_project",
                            user = "gurdeep",
                            password = "hellokitty",
                            host = "localhost",
                            port = "5432")
    return mydb
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
                     pfamPos VARCHAR(10), pfamAA VARCHAR(1),\
                     acc VARCHAR(20), gene VARCHAR(25), \
                     ptmType VARCHAR(10),\
                     name VARCHAR(50), \
                     CONSTRAINT name FOREIGN KEY (name) REFERENCES positions(name) \
                     )")
    for line in open('../data/Kinase_psites_hits_split_trimmed.tsv', 'r'):
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
    for line in tqdm(gzip.open('../data/humanKinasesHitsHmmsearchMappings.tsv.gz', 'rt')):
        if line[0] == '#': continue
        acc = line.split('\t')[0].split('|')[1]
        uniprot_id = line.split('\t')[0].split('|')[2]
        if acc not in mappings: mappings[acc] = {'uniprot_id': uniprot_id, 'positions': {}}
        uniprotAA = line.split('\t')[1]
        uniprotPos = int(line.split('\t')[2])
        pfamAA = line.split('\t')[3]
        pfamPos = line.split('\t')[4].rstrip()
        mappings[acc]['positions'][uniprotPos] = {'uniprotAA': uniprotAA,
                                                'pfamAA': pfamAA,
                                                'pfamPos': pfamPos}
        
        # print (mutation, wtAA, wtPos, mutAA, mut_type, acc, gene, info, source)
        '''
        mycursor.execute("INSERT INTO positions (uniprotP->dictos, uniprotAA, pfamPos, pfamAA, \
                         acc, uniprot_id) \
                         VALUES (%s, %s, %s, %s, %s, %s)", \
                            (uniprotPos, uniprotAA, pfamPos, pfamAA, acc, uniprot_id)
                            )
        '''
    return mappings

def create_mutations_table(mycursor)->None:
    '''Function to create the mutations table'''
    mycursor.execute("DROP TABLE IF EXISTS mutations CASCADE")
    mycursor.execute("CREATE TABLE mutations (id SERIAL PRIMARY KEY, \
                     mutation VARCHAR(10), wtAA VARCHAR(1), wtPos INT, \
                     mutAA VARCHAR(1), mut_type VARCHAR(10), \
                     acc VARCHAR(10), gene VARCHAR(10), \
                     info TEXT, source VARCHAR(200) \
                     )")
                    # UNIQUE(mutation, wtAA, wtPos, mutAA, mut_type, acc, gene, info, source)\
    for line in open('../AK_mut_w_sc_feb2023/act_deact_v2.tsv', 'r'):
        if line.split()[0] == 'uniprot_name': continue
        gene = line.split('\t')[0]
        acc = line.split('\t')[1]
        wtAA = line.split('\t')[2].replace(',', '')
        mutAA = line.split('\t')[4].replace(',', '')
        if len(wtAA) > 1 or len(mutAA) > 1: continue
        wtPos = str(line.split('\t')[3])
        mut_type = line.split('\t')[5]
        # print (acc, kinases[acc].gene, wtAA, position, mutAA)
        mutation = wtAA + wtPos + mutAA
        info = line.split('\t')[-2]
        # info = 'info'
        source = line.split('\t')[-1]
        # print (mutation, wtAA, wtPos, mutAA, mut_type, acc, gene, info, source)
        mycursor.execute("INSERT INTO mutations (mutation, wtAA, wtPos, mutAA, mut_type, \
                         acc, gene, info, source) \
                         VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)", \
                        (mutation, wtAA, wtPos, mutAA, mut_type, acc, gene, info, source))

def createDicForDSSP(dic, position, mutation, value):
    if position not in dic: dic[position] = {}
    dic[position][mutation] = float(value)

def create_homology_table(mycursor) -> None:
    '''Function to create the homology  table'''
    mycursor.execute("DROP TABLE IF EXISTS homology CASCADE")
    mycursor.execute("CREATE TABLE homology (id SERIAL PRIMARY KEY, \
                     acc VARCHAR(10), mutation VARCHAR(10), \
                     wtaa VARCHAR(5), position INT, mutaa VARCHAR(5), \
                     wtscore FLOAT, mutscore FLOAT, diffscore FLOAT, \
                     info TEXT) \
                     ")
    path = '/net/home.isilon/ds-russell/mechismoX/analysis/alignments/data/HUMAN/orthologs_only/'
    mycursor.execute('select acc from kinases')
    print (mycursor.fetchall())
    for row in tqdm(mycursor.fetchall()):
        acc = row[0]
        for fileEnd in [
                        '_all_homs.scores.txt.gz',
                        '_orth.scores.txt.gz',
                        '_excl_para.scores.txt.gz',
                        '_spec_para.scores.txt.gz',
                        '_bpso.scores.txt.gz',
                        '_bpsh.scores.txt.gz'
                        ]:
            if os.path.isfile(path + acc[:4] + '/' + acc + fileEnd) is False:
                print (path + acc[:4] + '/' + acc + fileEnd, 'does not exist')
                continue
            for line in gzip.open(path + acc[:4] + '/' + acc + fileEnd, 'rt'):
                #print (acc, line.split())
                #sys.exit()
                mutation = line.split()[0].split('/')[1]
                position = int(mutation[1:-1])
                wtaa = mutation[0]
                mutaa = mutation[1]
                wtscore = float(line.split()[2])
                mutscore = float(line.split()[3])
                diffscore = float(line.split()[4])
                info = line.split()[5].rstrip()
                mycursor.execute('INSERT INTO homology (acc, mutation, \
                                 wtaa, position, mutaa, \
                                 wtscore, mutscore, diffscore, \
                                 info) \
                                 VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)', \
                                (acc, mutation, wtaa, position, mutaa, \
                                wtscore, mutscore, diffscore, \
                                info))

def create_kinases_table(mycursor)->None:
    '''Function to create the kinases table'''
    mycursor.execute("DROP TABLE IF EXISTS kinases CASCADE")
    mycursor.execute("CREATE TABLE kinases (id SERIAL PRIMARY KEY, \
                     acc VARCHAR(10), gene VARCHAR(10), uniprot_id VARCHAR(25), \
                     fasta TEXT) \
                     ")
                    #  UNIQUE(acc, gene, uniprot_id, fasta))\
    mycursor.execute("DROP TABLE IF EXISTS positions CASCADE")
    mycursor.execute("CREATE TABLE positions (\
                     uniprotPos INT, uniprotAA VARCHAR(1),\
                     pfamPos VARCHAR(10) REFERENCES hmm DEFERRABLE, \
                     pfamAA VARCHAR(1), acc VARCHAR(20), \
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
            name = acc + '|' + uniprot_id + '|' + gene
            kinases[name] = ''
        else:
            kinases[name] += line.rstrip()

    mappings = fetch_mappings_dic()
    num = 0
    for kinase in tqdm(kinases):
        num += 1
        if num == 10:
            break
        acc = kinase.split('|')[0]
        uniprot_id = kinase.split('|')[1].split()[0]
        gene = kinase.split('|')[2]
        fasta = str(kinases[kinase])
        # print (acc, gene, uniprot_id, fasta)
        mycursor.execute("INSERT INTO kinases (acc, gene, uniprot_id, fasta) \
                            VALUES (%s, %s, %s, %s)", \
                            (acc, gene, uniprot_id, fasta))
        if acc not in mappings: continue
        for uniprotPos, uniprotAA in enumerate(fasta, start=1):
            if uniprotPos in mappings[acc]['positions']:
                if mappings[acc]['positions'][uniprotPos]['uniprotAA'] != uniprotAA:
                    print ('ERROR', uniprotPos, uniprotAA, mappings[acc]['positions'][uniprotPos]['uniprotAA'])
                    sys.exit()
                pfamPos = mappings[acc]['positions'][uniprotPos]['pfamPos']
                pfamAA = mappings[acc]['positions'][uniprotPos]['pfamAA']
                # print (uniprotPos, uniprotAA, pfamPos, pfamAA, acc, uniprot_id)
            else:
                pfamPos, pfamAA = '-', '-'
            name = acc+'/'+uniprotAA+str(uniprotPos)
            mycursor.execute("INSERT INTO positions (uniprotPos, uniprotAA, pfamPos, pfamAA, \
                                acc, uniprot_id, name) \
                                VALUES (%s, %s, %s, %s, %s, %s, %s)", \
                                (uniprotPos, uniprotAA, pfamPos, pfamAA, acc, uniprot_id, name)
                                )

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
    create_hmm_table(mycursor)
    create_mutations_table(mycursor)
    create_kinases_table(mycursor)
    create_homology_table(mycursor)
    sys.exit()
    create_ptm_table(mycursor)
    mydb.commit()

    # Use mysqldump to create backup file
    # backup_file = "kinaseDB.sql"
    # os.system(f"mysqldump -u {mydb.user} {mydb.database} > {backup_file}")

    # Close MySQL connection
    mydb.close()
    # main()