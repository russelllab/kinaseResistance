import os, sys, gzip

kinases = {}
flag = 0
for line in open(sys.argv[1], 'r'):
    if line[0] == '>':
        ## Get rid of the inactive kinases
        if line.split('>')[1].split()[1].strip() == 'Inactive':
            flag = 0
            continue
        kinase = line.split('>')[1].split()[0].strip()
        kinases[kinase] = ''
        flag = 1
    else:
        if flag == 1:
            kinases[kinase] += line.strip()

os.system('hmmsearch --noali --domtblout humanKinasesHits.hmmsearch ../pfam/Pkinase.hmm '+sys.argv[1])

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