import gzip, os


kinases = []
for line in open('humanKinases.fasta', 'r'):
    if line.startswith('>'):
        acc = line.split('|')[1]
        kinases.append(acc)

path2fasta = '../KA/UniProtFasta2/'
text = ''
for acc in kinases:
    fastaFile = path2fasta + acc + '.fasta.gz'
    if os.path.exists(fastaFile) is False:
        os.system('wget https://www.uniprot.org/uniprot/' + acc + '.fasta -P ' + path2fasta)
        os.system('gzip ' + path2fasta + acc + '.fasta')
    for line in gzip.open(fastaFile, 'rt'):
        text += line

open('humanKinases2.fasta', 'w').write(text)
