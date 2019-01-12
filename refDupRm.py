#!/usr/bin/env python3
import sys

def fasRead(fname):
	with open(fname, 'r') as fin:
		dict_seq = {}
		for l in fin:
			l = l.rstrip()
			if l.startswith('>'):
				ls = l.split()
				key = ls[0][1:]
				dict_seq[key] = ''
			else:
				dict_seq[key] += l
		return dict_seq

i = 0
for f in sys.argv[: -1]:
	if i == 0:
		dict_total = fasRead(f)
	else:
		dict_tmp = fasRead(f)
		for key1 in dict_total:
			for key2 in dict_tmp:
				if dict_total[key1] == dict_tmp[key2]:
					dict_total.pop(key1)
					dict_total[key1 + '|' + key2] = dict_tmp[key2]
				else:
					if key1 != key2:
						dict_total[key2] = dict_tmp[key1]
					else:
						print("The two key are the same!!!")


with open(sys.argv[-1], 'w') as fout:
	for key in dict_total:
		fout.write(">%s\n%s\n" % (key, dict_total[key]))