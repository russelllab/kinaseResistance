from Bio.PDB import *
from Bio.SeqUtils import seq1
import os, sys, gzip
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

parser = MMCIFParser()
pdb = '1gag.cif'
structure = parser.get_structure(pdb, ""+pdb)
dic_canonical2pdb = {}
fasta_count = 0
for model in structure:
    for chain in model:
        if (chain.id == 'A'):
            fasta = '>'+pdb+'_'+chain.id+'\n'
            map = []
            for residue in chain:
                for atom in residue:
                    if (atom.id == 'CA'):
                        if residue.id[0] == ' ':
                            fasta_count += 1
                            #print (pdb, chain.id, residue.get_resname(), residue.id)
                            # AA = protein_letters_3to1[residue.get_resname()]
                            AA = seq1(residue.get_resname())
                            fasta += AA
                            dic_canonical2pdb[fasta_count] = residue.id[1]
                            print (residue.id[1], AA, fasta_count)
            print (fasta)
            open(pdb+'_'+chain.id+'.fasta', 'w').write(fasta)
            break
# sys.exit()
# run hmmsearch
os.system('hmmsearch -o hmmsearch_INSR.txt ../pfam/Pkinase.hmm '+pdb+'_A.fasta')

# read hmmsearch output
go = 0
dic_canonical2domain = {}
dic_domain2canonical = {}
for line in open('hmmsearch_INSR.txt', 'r'):
    if len(line.split()) <= 1: continue
    if line[0] == '>':
        go = 1
    elif go == 1 and line.split()[0] == 'Pkinase':
        qs = int(line.split()[1])
        qe = line.split()[3]
        query = line.split()[2]
    elif go == 1 and line.split()[0] == pdb+'_A':
        ss = int(line.split()[1])
        se = line.split()[3]
        sbjct = line.split()[2]

        q = qs
        s = ss
        for i in range(0, len(sbjct)):
            if sbjct[i] not in ['-', '.'] and query[i] not in ['-', '.']:
                dic_canonical2domain[s] = q
                dic_domain2canonical[q] = s
                s += 1
                q += 1
            if sbjct[i] in ['-', '.'] and query[i] not in ['-', '.']:
                q += 1
            if sbjct[i] not in ['-', '.'] and query[i] in ['-', '.']:
                s += 1
# print (dic_domain2canonical[442])
# print (dic_canonical2pdb[dic_domain2canonical[442]])

def make_colors(num_colors):
    # Exclude colors: red, green, light green, coral, blue
    excluded_colors = ['red', 'green', 'lightgreen', 'coral', 'blue']

    # Generate a colormap with gradually changing colors
    colormap = plt.cm.get_cmap('viridis')

    # Create an array to store the colors
    colors = []

    # Generate unique colors
    for i in range(num_colors):
        while True:
            # Generate a random color from the colormap
            color = colormap(np.random.uniform())
            r , g, b, _ = color

            # Convert the color to its hexadecimal representation
            hex_color = mcolors.rgb2hex(color)

            # Check if the color is not in the excluded colors
            if hex_color not in excluded_colors:
                # colors.append(hex_color)
                colors.append([r, g, b])
                break

    # Print the array of colors
    print(colors)
    colors = ['cyan', 'brightorange', 'brown', 'deepolive', 'deepsalmon',
                'deepteal', 'gray', 'lightblue', 'dirtyviolet', 'purple',
                'aquamarine', 'skyblue', 'yellow', 'sulfur', 'wheat',
                'smudge', 'sand', 'violetpurple', 'tv_blue', 'marine',
                'greencyan', 'palecyan', 'magenta', 'salmon', 'orange',
                'greencyan', 'palecyan', 'magenta', 'salmon', 'orange'
            ]
    colors = []
    for line in open('ss_colors.tsv', 'r'):
        if line.startswith('#'): continue
        colors.append(line.split()[1].rstrip())
    return colors

dic_ss = {}
for line in open('ss.tsv', 'r'):
    if line.startswith('#'): continue
    if line.split()[0] in ['DFG-motif', 'HrD-motif', 'APE-motif',
                            'Catalytic-Lys', 'Gly-rich-loop']:
        continue
    name = line.split()[0]
    start, end = line.split()[1].rstrip().split('-')
    print (line)
    print (name, start, end)
    start = int(start)
    end = int(end)
    dic_ss[name] = [start, end]
colors = make_colors(len(dic_ss.keys()))
for color, name in zip(colors, dic_ss):
    dic_ss[name].append(color)
# sys.exit()

pdb_colors = {}
for line in open('INSR_colors.txt', 'r'):
    if line.startswith('#'): continue
    start, end = line.split()[0].rstrip().split('-')
    start = int(start)
    end = int(end)
    color = line.split()[1].rstrip()
    for i in range(start, end+1):
        pdb_colors[i] = color

chain = 'A'
cmd.load(pdb)
cmd.hide('everything')
cmd.show('cartoon', 'chain '+str(chain))
cmd.color('hydrogen', 'chain '+chain)
# cmd.color('wheat', 'chain B')
# cmd.show('surface',  'chain B')
cmd.set('cartoon_fancy_helices', 1)
cmd.set('cartoon_fancy_sheets', 1)
cmd.set('antialias', 2)
# cmd.select('chain', 'chain '+chain)
cmd.set('cartoon_oval_width', 0.5)
cmd.set('cartoon_oval_length', 0.9)
cmd.set('cartoon_oval_quality', 100000)
cmd.bg_color('white')
for pdb_position in pdb_colors:
    cmd.color(pdb_colors[pdb_position], chain+'/'+str(pdb_position)+'/ca')
    # cmd.show('spheres', chain+'/'+str(pdb_position)+'/ca')
for pdb_position in range(978, 994):
    cmd.select('sele', 'chain '+chain+' and res '+str(pdb_position))
    cmd.hide('everything', 'sele')
    cmd.deselect()
    # cmd.show('spheres', chain+'/'+str(pdb_position)+'/ca')
'''
for domain_position in dic_domain2canonical:
    domain_color = 'cyan'
    for row in dic_ss:
        start, end, color = dic_ss[row]
        if start <= domain_position <= end:
            domain_color = color
            break
    # print (domain_position, domain_color)
    canonical_position = dic_domain2canonical[domain_position]
    if canonical_position not in dic_canonical2pdb: continue
    pdb_position = dic_canonical2pdb[canonical_position]
    if int(domain_position) == 141:
        print (canonical_position, pdb_position, domain_color)
    cmd.color(domain_color, chain+'/'+str(pdb_position)+'/ca')
    '''
    # if dic[position] in diseased:
    #     cmd.color('green', chain.id+'/'+str(dic[position])+'/ca')
    #     cmd.show('spheres', chain.id+'/'+str(dic[position])+'/ca')
    # elif dic[position] in neutral:
    #     cmd.color('red', chain.id+'/'+str(dic[position])+'/ca')
    #     cmd.show('spheres', chain.id+'/'+str(dic[position])+'/ca')
    #else:
    #    cmd.color('gold', chain.id+'/'+str(dic[pos])+'/ca')
