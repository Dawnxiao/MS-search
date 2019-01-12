#!/usr/bin/env python3
##change pt sequence according to vep results of dbSNP or cosmic
import sys
import re

file_fasta = sys.argv[1]
file_vep = sys.argv[2]
file_out = sys.argv[3]

#a function to substitute charcter in a string
def sub(string, p, c):
    new = [s for s in string]
    new[p] = c
    return ''.join(new)

#read in reference pt
dict_pt = {}
with open(file_fasta, "r") as fin:
    for l in fin:
        l = l.rstrip()
        if l.startswith('>'):
            name = l.split()[4].split(':')[1]
            name = re.search('ENST\d+', name).group()
            dict_pt[name] = ''
        else:
            dict_pt[name] += l

#read in vep results file and write down the postion of amino_acid and which amino_acid to change into
list_tp = []
with open(file_vep, "r") as fin:
    for l in fin:
        if not l.startswith("#"):
            ls = l.rstrip().split()
            ptID = ls[4]
            list_tp.append(ptID)
            ptPos = ls[9]
            amino_acid = ls[10].split('/')[0]
            ref_aa = ls[10].split('/')[1]
            if ptID in dict_pt:
                if re.search('-', ptPos):
                    start, end = int(re.search('(\d+)-(\d+)', ptPos).group(1)), int(re.search('(\d+)-(\d+)', ptPos).group(2))
                    for i in range(end - start + 1):
                        if len(dict_pt[ptID]) >= end:
                            if ref_aa[i] == dict_pt[ptID][start - 1 + i]
                                dict_pt[ptID] = sub(dict_pt[ptID], start - 1 + i, amino_acid[i])
                            else:
                                print('amino_acid in vep file is not consistant with amino_acid in reference protein for %s' % ptID)
                        else:
                            print(l)
                else:
                    if len(dict_pt[ptID]) >= int(ptPos):
                        if ref_aa == dict_pt[ptID][int(ptPos)]
                            dict_pt[ptID] = sub(dict_pt[ptID], int(ptPos) - 1, amino_acid)
                        else:
                            print('amino_acid in vep file is not consistant with amino_acid in reference protein for %s' % ptID)
                    else:
                        print(l)
            else:
                print(ptID + " has no conresponding transcript in ensembl")

#save the changed pt sequence to a file
with open(file_out, 'w') as fout:
    for key in set(list_tp):
        fout.write('>' + key + '\n' + dict_pt[key] + '\n')