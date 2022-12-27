'''
Script to make alignment files by splitting
given aln types in 25, 50, 75 and 100
'''
import os, gzip, sys
import argparse, threading

# acc = 'P15056'
# aln_type = '_orth.aln.gz'
path_to_files = '../alignments/'
path_to_output = '../GS_output/'
path_to_blastdb = '../alignments_blastdb/'
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
    order = []
    for line in gzip.open(path_to_files+acc[:4]+'/'+acc+aln_type, 'rt'):
        if len(line.split()) == 0:
                continue
        if line.split()[0] in ['CLUSTAL', '//']:
                continue
        name = line.split(' ')[0].replace('\n', '')
        if name not in order: order.append(name)
        if name not in kinase_dic: kinase_dic[name] = Kinase(name)
        kinase_dic[name].seq += line.split()[1].replace('\n', '').replace('-', '')
    return kinase_dic

def do_blast(acc, aln_type, kinase_dic):
    fasta_db = ''
    for name in kinase_dic:
        if name.split('_')[1] == acc:
            continue
        fasta_db += '>' + name + '\n' + kinase_dic[name].seq + '\n'
    fasta_file_name = acc+aln_type.split('.')[0]+'_seq.fa'
    open(path_to_blastdb+fasta_file_name, 'w').write(fasta_db)
    if os.path.isfile(path_to_blastdb+acc+'.fa') is False:
        open(path_to_blastdb+acc+'.fa', 'w').write('>HUMAN_'+acc+'\n'+kinase_dic['HUMAN_'+acc].seq)
    os.system('makeblastdb -dbtype prot -in '+
                path_to_blastdb+fasta_file_name +
                ' -input_type fasta ' +
                '-out '+path_to_blastdb+fasta_file_name.split('.fa')[0]
                )
    os.system('blastp -query '+path_to_blastdb+acc+'.fa' +
                ' -db '+ path_to_blastdb + fasta_file_name.split('.fa')[0] +
                ' -outfmt 6' +
                ' -out ' + path_to_blastdb + fasta_file_name.split('_seq')[0] + '_blastp.txt'
            )
    os.system('gzip '+path_to_blastdb + fasta_file_name.split('_seq')[0] + '_blastp.txt')
    order = []
    for line in gzip.open(path_to_blastdb + fasta_file_name.split('_seq')[0] + '_blastp.txt.gz', 'rt'):
        name = line.split('\t')[1]
        if name not in order: order.append(name)
    ## Delete all just created files
    # os.system('rm -rf '+path_to_blastdb + fasta_file_name.split('_seq')[0] + '_blastp.txt.gz')
    os.system('rm -rf '+path_to_blastdb+fasta_file_name.split('.fa')[0] + '*')
    # os.system('rm -rf '+path_to_blastdb+acc+'.fa')
    return order

def main(acc, aln_type):
    '''
    Extract order from the given alignment
    and realign using hmmalign against Pkinase
    '''
    # for aln_type in ['_bpsh.aln.gz']:
    # for aln_type in ['_all_homs.aln.gz', '_bpsh.aln.gz', '_bpso.aln.gz', '_excl_para.aln.gz', '_orth.aln.gz' , '_spec_para.aln.gz']:
    if True:
        # print (aln_type)
        kinase_dic = extract_fasta(acc, aln_type)
        order = do_blast(acc, aln_type, kinase_dic)
        # order = []
        # for line in gzip.open(path_to_files+acc[:4]+'/'+acc+aln_type, 'rt'):
        #     if line.split() == []:
        #         continue
        #     if line.split()[0] in ['CLUSTAL', '//']:
        #         continue
        #     name = line.split(' ')[0]
        #     if name not in order: order.append(name)

        ## Split and create files based on identity
        dic_arr = {25: '', 50: '', 75: '', 100: ''}
        count = 0
        for name in order:
            count += 1
            # print (name)
            for num in dic_arr:
                if float(count)/len(order)*100 <= num:
                    dic_arr[num] += '>' + name + '\n' + kinase_dic[name].seq + '\n'

        ## Create FASTA and run hmmalign
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
            os.system('gzip '+path_to_output + fasta_file_name)
            os.system('gzip '+path_to_output + aln_file_name)
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

class myThread (threading.Thread):
   def __init__(self, threadID, acc, aln_type):
      threading.Thread.__init__(self)
      self.threadID = threadID
      self.acc = acc
      self.aln_type = aln_type
   def run(self):
      print ("Starting ", self.acc, self.aln_type)
      main(self.acc, self.aln_type)
      print ("Exiting ", self.acc, self.aln_type)

if __name__ == '__main__':
    '''
    Execute this when file runs as script
    and not when run as module
    '''
    parser = argparse.ArgumentParser(description='Generate log-odds score from *.aln.gz and corresponding *.hmm.gz files',
                                        epilog='gurdeep.singh[at]bioquant[.]uni[-]heidelberg[.]de',
                                        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('threads', help='# threads')
    args = parser.parse_args()
    num_threads = int(args.threads)
    dic_threads = {}
    # kinase_dic = extract_fasta(acc)
    # for files in os.listdir('../kin_seq/'):
    for line in open('humanKinases.fasta', 'r'):
        # if files.endswith('.fa.gz') is False:
        #     continue
        if line[0] != '>':
            continue
        acc = line.split('|')[1]
        if os.path.isdir(path_to_files+acc[:4]) is False:
            continue
        for aln_type in [
                        #'_all_homs.aln.gz',
                        '_bpsh.aln.gz',
                        '_bpso.aln.gz',
                        # '_excl_para.aln.gz',
                        '_orth.aln.gz' ,
                        '_spec_para.aln.gz']:
            # main(acc, aln_type)
            dic_threads[acc+aln_type] = myThread(0, acc, aln_type)
            dic_threads[acc+aln_type].start()
            print (threading.active_count())
            while threading.active_count() >= num_threads:
                continue
        break
    print ('done')

            # dic_threads = myThread(2, "Thread-2", 2)

