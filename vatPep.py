#!/usr/bin/env python3
import sys, re

file_fasta = sys.argv[1] #reference protein fasta
file_dbSNP = sys.argv[2] #dbSNP file
file_out = sys.argv[3]

#a function to substitute charcter in a string
def sub(string, p, c):
    new = [s for s in string]
    new[p] = c
    return ''.join(new)

#find out all of the cleave site
def fdAll(string, pattern = 'KR'):
    return [i for i in range(len(string)) if string[i] in pattern]

#return all the conbinations of SNP in a peptide
def cbOfSNP(list_rs, l, res, tmp = []):
    if len(tmp) == l:
        res.append(tmp)
        return
    else:
        tmp1 = tmp.copy()
        tmp1.append([list_rs[0][0], list_rs[0][1], list_rs[0][3], 0])
        cbOfSNP(list_rs[1:], l, res, tmp = tmp1)
        tmp2 = tmp.copy()
        tmp2.append([list_rs[0][0], list_rs[0][1], list_rs[0][2], 1])
        cbOfSNP(list_rs[1:], l, res, tmp = tmp2)

#find out all the snp in a peptide
def fdRsInAReg(list_clv_pos, i, list_rs, ptLen, ptID):
    list_rs_reg = []
    list_rs.sort(key = lambda ele: ele[1], reverse = False)
    if i == 0:
        for rs in list_rs:
            if rs[1] >= 1 and rs[1] <= (list_clv_pos[1] + 1):
                list_rs_reg.append(rs)
            elif rs[1] > (list_clv_pos[1] + 1):
                return list_rs_reg
        return list_rs_reg
    elif i == (len(list_clv_pos) - 1):
        for rs in list_rs:
            if rs[1] > list_clv_pos[i] + 1 and rs[1] <= ptLen:
                list_rs_reg.append(rs)
            elif rs[1] > ptLen:
                print("%s is out of the range of protein %s" % (rs[0], ptID))
                sys.exit()
        return list_rs_reg
    else:
        for rs in list_rs:
            if rs[1] > list_clv_pos[i] + 1 and rs[1] <= (list_clv_pos[i + 1] + 1):
                list_rs_reg.append(rs)
            elif rs[1] > (list_clv_pos[i + 1] + 1):
                return list_rs_reg
        return list_rs_reg

#output all the result
def cbOtpt(fl, cbs, seq_pt, clv_reg, list_clv_pos, list_rs_pt, clv_reg_num, ptID, rsName = []):
    rsName_tmp = rsName.copy()
    if len(cbs) == 0:
        return
    for cb in cbs:
        seq_tmp = seq_pt
        cb.sort(key = lambda ele: ele[1], reverse = False)
        i = 0
        rsName = rsName_tmp.copy()
        start, end = clv_reg #start with 0
        while True:
            old = seq_tmp[cb[i][1] - 1]
            seq_tmp = sub(seq_tmp, cb[i][1] - 1, cb[i][2])
            if cb[i][2] in 'KR' and old not in 'KR':
                rsName.append(cb[i][0])
                print(str(cb[i]) + str(i) + '1')
                if (cb[i][1] - start) >= 6 and (cb[i][1] - start) <= 40:
                    fl.write(">" + ptID + '|' + '|'.join(rsName) + '\n' + seq_tmp[start : cb[i][1]] + '\n')
                    print(">" + ptID + '|' + '|'.join(rsName) + '\n' + seq_tmp[start : cb[i][1]] + '\n')
                if i >= (len(cb) - 1):
                    if (end - cb[i][1] + 1) >= 6 and (end - cb[i][1] + 1) <= 40:
                        fl.write(">" + ptID + '|' + '|'.join(rsName.append(cb[i][0])) + '\n' + seq_tmp[cb[i][1] : end + 1] + '\n')
                        print(">" + ptID + '|' + '|'.join(rsName.append(cb[i][0])) + '\n' + seq_tmp[cb[i][1] : end + 1] + '\n')
                    break
                start = cb[i][1]
                rsName = [cb[i][0]]
            elif cb[i][2] not in 'KR' and old in 'KR':
                print(str(cb[i]) + str(i) + '2')
                rsName.append(cb[i][0])
                next_rs = fdRsInAReg(list_clv_pos, clv_reg_num + 1, list_rs_pt, len(seq_tmp), ptID)
                next_cbs = []
                cbOfSNP(next_rs, len(next_rs), next_cbs)
                cbOtpt(fl, next_cbs, seq_tmp, (start, list_clv_pos[clv_reg_num +2]), list_clv_pos, list_rs_pt, clv_reg_num + 1, ptID, rsName = rsName)
                break
            elif cb[i][2] not in 'KR' and old not in 'KR':
                print(str(cb[i]) + str(i) + '3')
                if cb[i][3] == 1:
                    rsName.append(cb[i][0])
                if i >= len(cb) - 1:
                    if (end - start + 1) >= 6 and (end - start + 1) <= 40:
                        fl.write(">" + ptID + '|' + '|'.join(rsName) + '\n' + seq_tmp[start : end + 1] + '\n')
                        print((">" + ptID + '|' + '|'.join(rsName) + '\n' + seq_tmp[start : end + 1] + '\n'))
                    break
            else:
                print(str(cb[i]) + str(i) + '4')
                if cb[i][3] == 1:
                    rsName.append(cb[i][0])
                if (cb[i][1] - start) >= 6 and (cb[i][1] - start) <= 40 and len(rsName) != 0:
                    fl.write(">" + ptID + '|' + '|'.join(rsName) + '\n' + seq_tmp[start : cb[i][1]] + '\n')
                    print(">" + ptID + '|' + '|'.join(rsName) + '\n' + seq_tmp[start : cb[i][1]] + '\n')
                break
            i += 1
        
        

#read in reference pt
dict_pt = {}
with open(file_fasta, "r") as fin:
    for l in fin:
        l = l.rstrip()
        if l.startswith('>'):
            name = re.search('ENST\d+', l).group()
            dict_pt[name] = ''
        else:
            dict_pt[name] += l

#read in dbSNP file
dict_SNP = {}
with open(file_dbSNP, 'r') as fin:
    for l in fin:
        l = l.rstrip()
        
        
#cleave protein into peptide and change reference amino acid into alt amino acid
fout = open(file_out, 'w')
for key in dict_SNP:
    list_clv_pos = fdAll(dict_pt[key])#find out all the cleave site
	list_clv_pos.insert(0, 0)
    list_rs_pt = dict_SNP[key]
    i = 0 #cleave region number
    for i in range(len(list_clv_pos)):
        list_rs_reg = fdRsInAReg(list_clv_pos, i, list_rs_pt, len(dict_pt[key]), key)
        cbs = []
        cbOfSNP(list_rs_reg, len(list_rs_reg), cbs)
        if i == 0:
            cbOtpt(fout, cbs, dict_pt[key], (0, list_clv_pos[i + 1]), list_clv_pos, list_rs_pt, i, key)
        elif i == len(list_clv_pos) - 1:
            cbOtpt(fout, cbs, dict_pt[key], (list_clv_pos[i] + 1, len(dict_pt[key]) - 1), list_clv_pos, list_rs_pt, i, key)
        else:
            cbOtpt(fout, cbs, dict_pt[key], (list_clv_pos[i] + 1, list_clv_pos[i + 1]), list_clv_pos, list_rs_pt, i, key)

