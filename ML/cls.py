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
            if domainNum in self.domains is False:
                break
            for hmmPos in self.domains[domainNum]:
                if seqPos == self.domains[domainNum][hmmPos]:
                    return hmmPos
            domainNum += 1

class Mutation:
    def __init__(self, mutation, mut_type, acc):
        self.position = int(mutation[1:-1])
        self.wtAA = mutation[0]
        self.mutAA = mutation[-1]
        self.positionHmm = None
        self.mut_types = [mut_type]