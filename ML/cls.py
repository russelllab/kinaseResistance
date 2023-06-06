#!/usr/bin/env python
# coding: utf-8

# Class for kinase information
class Kinase:
    def __init__(self, acc, gene):
        self.acc = acc
        self.gene = gene
        self.fasta = ''
        self.group = ''
        self.mutations = {}
        self.hmm = {}
        self.ptm = {}
        self.oneHotEncoding = {}
        self.domains = {}
        self.seq2pfam = {}
        self.hmmsearch = []
        self.access = {}
        self.dihedral = {}
        self.sec = {}
        self.burr = {}
        self.iupred = {}
        self.mechismo = {}
        self.allHomologs = {}
        self.exclParalogs = {}
        self.specParalogs = {}
        self.orthologs = {}
        self.bpso = {}
        self.bpsh = {}
    
    def returnhmmPos(self, seqPos):
        domainNum = 1
        while domainNum > 0:
            if domainNum not in self.domains:
                break
            for hmmPos in self.domains[domainNum]:
                if seqPos == self.domains[domainNum][hmmPos]:
                    return hmmPos
            domainNum += 1
        return None

class Mutation:
    def __init__(self, mutation, mut_type, acc, dataset):
        self.position = int(mutation[1:-1])
        self.dataset = dataset
        self.wtAA = mutation[0]
        self.mutAA = mutation[-1]
        self.positionHmm = None
        self.mut_types = [mut_type]
    
    def findChangeInCharge(self):
        '''
        Function to calculate the change in charge of the mutation
        and return the array of charges for the wildtype and mutant
        '''
        positiveAA = ['R', 'K', 'H']
        negativeAA = ['D', 'E']
        neutralAA = ['S', 'T', 'N', 'Q', 'C', 'G', 'P', 'A', 'V', 'I', 'L', 'M', 'F', 'W', 'Y']
        charges_AA = {}
        for aa in positiveAA:
            charges_AA[aa] = 1.0
        for aa in negativeAA:
            charges_AA[aa] = -1.0
        for aa in neutralAA:
            charges_AA[aa] = 0.0
        charges =[]
        charges.append(charges_AA[self.wtAA])
        charges.append(charges_AA[self.mutAA])
        charges.append(charges_AA[self.mutAA] - charges_AA[self.wtAA])
        return charges
        
    def checkPhosphomimic(self):
        if self.wtAA in ['S', 'T', 'Y'] and self.mutAA in ['D', 'E']:
            return 1
        elif self.wtAA in ['S', 'T', 'Y']:
            return -1
        else:
            return 0
    
    def checkAcetylmimic(self):
        if self.wtAA in ['K'] and self.mutAA in ['Q']:
            return 1
        elif self.wtAA in ['K']:
            return -1
        else:
            return 0