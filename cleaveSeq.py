#!/usr/bin/env python3
##read in reference protein and variant protein, then cleave the protein into peptide, next remove peptides in variant which are already in reference.
from pyteomics import parser
import sys, re

file_fasta_ref = sys.argv[1]
file_fasta_vat = sys.argv[2]
file_out = sys.argv[3]

#read in reference protein
dict_protein_ref = {}
with open(file_fasta_ref, 'r') as fin:
    for l in fin:
        l = l.rstrip()
        if l.startswith('>'):
            name = re.search('(ENST\d+)\.\d+', name).group(1)
            dict_protein_ref[name] = ''
        else:
            dict_protein_ref[name] += l

#cleave reference protein to peptide
dict_pep_ref = {}
for key in dict_protein_ref:
    peps = parser.cleave(dict_protein_ref[key], parser.expasy_rules['trypsin'], 2)
    i = 0
    for pep in peps:
        if len(pep) >= 6 and len(pep) <= 40:
            dict_pep_ref[key + "|pepID:" + str(i)] = pep
            i += 1

#read in variant protein
dict_protein_vat = {}
with open(file_fasta_vat, 'r') as fin:
    for l in fin:
        l = l.rstrip()
        if l.startswith('>'):
            name = re.search('(ENST\d+)\.\d+', l).group(1)
            dict_protein_vat[name] = ''
        else:
            dict_protein_vat[name] += l

#cleave variant protein to peptide
dict_pep_vat = {}
for key in dict_protein_vat:
    peps = parser.cleave(dict_protein_vat[key], parser.expasy_rules['trypsin'], 2)
    i = 0
    for pep in peps:
        if len(pep) >= 6 and len(pep) <= 40:
            dict_pep_vat[key + '|pepID:' + str(i)] = pep
            i += 1

#remove peptide in variant which are already in the reference
dict_pep_vat_cp = dict_pep_vat.copy()
for key in dict_pep_vat:
    if dict_pep_vat[key] in dict_pep_ref.values():
        dict_pep_vat_cp.pop(key, None)

with open(file_out, 'w') as fout:
    for key in dict_pep_vat:
        fout.write(">" + key + "\n" + dict_pep_vat[key] + "\n")