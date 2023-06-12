import os, sys, gzip

inactive_domains = ['O43187', 'O15197', 'Q13308', 'Q58A45',\
                    'Q5JZY3', 'Q6VAB6', 'Q7RTN6', 'Q7Z7A4',\
                    'Q86YV5', 'Q8IVT5', 'Q8IWB6', 'Q8IZE3',\
                    'Q8NB16', 'Q8NCB2', 'Q8NE28', 'Q8TEA7',\
                    'Q92519', 'Q96C45', 'Q96KG9', 'Q96RU7',\
                    'Q96RU8', 'Q9C0K7', 'Q9H792', 'Q9NSY0',\
                    'Q9UHY1', 'Q9Y616']
pseudokinases = ['Q58A45', 'Q96S38_1', 'P29597_1', 'P23458_1', 'P52333_1', 'O60674_1']
kinases = {}
flag = 0
for line in open(sys.argv[1], 'r'):
    if line[0] == '>':
        ## Get rid of the inactive kinases
        if line.split('>')[1].split()[1].strip() == 'Inactive':
            flag = 0
            continue
        ## Get rid of the one that have catalytically inactive domain/pseduokinases
        ## It's hard to know which one is the catalytically inactive domain
        ## but try to make a list as you go along (see above)
        ## Pseudokinases with multiple domains are removed later
        ## after running the hmmsearch
        if line.split('|')[1] in inactive_domains or line.split('|')[1] in pseudokinases:
            flag = 0
            continue
        kinase = line.split('>')[1].split()[0].strip()
        kinases[kinase] = ''
        flag = 1
    else:
        if flag == 1:
            kinases[kinase] += line.strip()

os.system('hmmsearch --domE 1.0 --noali --domtblout humanKinasesHits.hmmsearch ../pfam/Pkinase.hmm '+sys.argv[1])

hits = []
for line in open('humanKinasesHits.hmmsearch', 'r'):
    if line.startswith('#'):
        continue
    total_hits = int(line.split()[10])
    if total_hits == 1:
        hit = line.split()[0] + '|' + 'start' + '-' + 'end'
        hits.append(hit)
        continue
    current_hit = int(line.split()[9])
    ## Get rid of the pseudokinases which multiple domains
    ## The ones with single domain were already removed
    ## before running the hmmsearch
    acc = line.split()[0].split('|')[1]
    ignore = False
    for entry in pseudokinases:
        if '_' not in entry: continue
        pseudokinase = entry.split('_')[0]
        dom = entry.split('_')[1]
        if acc == pseudokinase and int(dom) == int(current_hit):
            # print (acc)
            ignore = True
            break
    if ignore: continue
    ##
    current_start = str(line.split()[17])
    current_end = str(line.split()[18])
    if current_hit == 1:
        hit = line.split()[0] + '|' + 'start' + '-' + current_end
    elif current_hit == total_hits:
        hit = line.split()[0] + '|' + current_start + '-' + 'end'
    else:
        hit = line.split()[0] + '|' + current_start + '-' + current_end
    hits.append(hit)

l = ''
for hit in hits:
    kinase = '|'.join(hit.split('|')[:-1])
    ## Eliminate the ones that are not present in the kinases dictionary
    ## These would be the inactive kinases
    if kinase not in kinases:
        continue
    start = hit.split('|')[-1].split('-')[0]
    end = hit.split('|')[-1].split('-')[1]
    if start == 'start' and end == 'end':
        seq = kinases[kinase]
    elif start == 'start' and end != 'end':
        end = int(end)
        seq = kinases[kinase][:end]
    elif start != 'start' and end == 'end':
        start = int(start)-1
        seq = kinases[kinase][start:]
    elif start != 'start' and end != 'end':
        start = int(start)-1
        end = int(end)
        seq = kinases[kinase][start:end]
    # print ('>'+hit+'\n'+seq)
    l += '>'+hit+'\n'+seq+'\n'

open('humanKinasesHitsSplit.fasta', 'w').write(l)