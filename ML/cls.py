#!/usr/bin/env python
# coding: utf-8

# Class for kinase information
class Kinase:
    def __init__(self, acc, gene):
        self.acc = acc
        self.gene = gene
        self.fasta = ''
        self.group = ''
        self.hmm = {}
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