'''
Script to make alignment files by splitting
given aln types in 25, 50, 75 and 100
'''
import os, gzip, sys
import argparse

# acc = 'P15056'
# aln_type = '_orth.aln.gz'
path_to_files = '../kin_seq/'
path_to_output = '../GS_output/'
pfam_dom = 'Pkinase'
pfam_path = '../pfam/'

class Kinase:
    def __init__(self, name) -> None:
        self.name = name
        self.seq = ''

def extract_fasta(acc, aln_type):
    '''
    Extract FASTA from *aln_type
    '''
    kinase_dic = {}
    for line in gzip.open(path_to_files+acc+aln_type, 'rt'):
        if len(line.split()) == 0:
                continue
        if line.split()[0] in ['CLUSTAL', '//']:
                continue
        name = line.split(' ')[0].replace('\n', '')
        if name not in kinase_dic: kinase_dic[name] = Kinase(name)
        kinase_dic[name].seq += line.split()[1].replace('\n', '').replace('-', '')
    return kinase_dic

def main(acc):
    '''
    Extract order from the given alignment
    and realign using hmmalign against Pkinase
    '''
    # for aln_type in ['_bpsh.aln.gz']:
    for aln_type in ['_all_homs.aln.gz', '_bpsh.aln.gz', '_bpso.aln.gz', '_excl_para.aln.gz', '_orth.aln.gz' , '_spec_para.aln.gz']:
        kinase_dic = extract_fasta(acc, aln_type)
        order = []
        for line in gzip.open(path_to_files+acc+aln_type, 'rt'):
            if line.split() == []:
                continue
            if line.split()[0] in ['CLUSTAL', '//']:
                continue
            name = line.split(' ')[0]
            if name not in order: order.append(name)

        dic_arr = {25: '', 50: '', 75: '', 100: ''}
        count = 0
        for name in order:
            count += 1
            # print (name)
            for num in dic_arr:
                if float(count)/len(order)*100 <= num:
                    dic_arr[num] += '>' + name + '\n' + kinase_dic[name].seq + '\n'

        for num in dic_arr:
            fasta_file_name = acc+aln_type.split('.')[0]+'_'+str(num)+'_seq.fa'
            if dic_arr[num] == '':
                print (fasta_file_name, 'cannot be created coz it will be empty')
                continue
            open(path_to_output+fasta_file_name, 'w').write(dic_arr[num])
            ## hmmalign
            aln_file_name = acc+aln_type.split('.')[0]+'_'+str(num)+'_seq_hmmalign.aln'
            os.system(# 'zcat ' + path_to_output + fasta_file_name + '|'
                'hmmalign -o '+
                path_to_output + aln_file_name + ' ' +
                pfam_path+pfam_dom+'.hmm' + ' ' +
                path_to_output + fasta_file_name
                )
            # sys.exit()
            '''
            ## mafft
            aln_file_name = acc+aln_type.split('.')[0]+'_'+str(num)+'_seq_mafft.aln.gz'
            os.system(# 'zcat ' + path_to_output + fasta_file_name + '|'
                'mafft --thread -1 --clustalout '+ path_to_output + fasta_file_name + '>' + 
                path_to_output + aln_file_name)
            '''
            ## gzip the fasta
            # os.system('gzip '+path_to_output + fasta_file_name)
            # break

if __name__ == '__main__':
    '''
    Execute this when file runs as script
    and not when run as module
    '''
    parser = argparse.ArgumentParser(description='Generate log-odds score from *.aln.gz and corresponding *.hmm.gz files',
                                        epilog='gurdeep.singh[at]bioquant[.]uni[-]heidelberg[.]de',
                                        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('acc', help='UniProt accession')
    args = parser.parse_args()
    acc = args.acc
    # kinase_dic = extract_fasta(acc)
    main(acc)
