import numpy as np
import csv
import argparse
import itertools

parser=argparse.ArgumentParser()
parser.add_argument("truth", help="set of true GTs from unphased VCF",type=str)
parser.add_argument("test_VCF", help="first VCF to compare against the true set", type=str)
args=parser.parse_args()

truf={}

def change(val):
	if '.' not in val[0]:
		new=sum([int(i) for i in val])
	else:
		new='|'.join(val)
	return new

with open(args.truth) as t:
	for line in csv.reader(t,delimiter="\t"):
		if "#" in line[0]:
			continue
		else:
			nline=[change([j for j in i.split(':')[0].split('/')]) for i in line[9:]]
			truf[line[1]]=nline

wng=0
corr=0

with open(args.test_VCF) as t:
	for line in csv.reader(t,delimiter='\t'):
		if "#" in line[0]:
			continue
		else:
			nline=[change([j for j in i.split(':')[0].split('|')]) for i in line[9:]]
			for j in range(len(nline)):
				if truf[line[1]][j]=='.|.':
					continue
				elif nline[j]==truf[line[1]][j]:
					corr+=1
				else:
					wng+=1

print "Total Compared:", wng+corr
print "Total Correct:", corr
print "Percent Correct:", float(corr)/(wng+corr)