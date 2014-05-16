import numpy
import argparse
import random
import csv

parser = argparse.ArgumentParser()
parser.add_argument("VCF", help="Give the knockout tool a VCF",type=str)
parser.add_argument("rate", help="Rate of Knocked out Genotypes between 0-1",type=float)
args = parser.parse_args()
with open(args.VCF) as t:
	for line in csv.reader(t,delimiter="\t"):
		if "#" in line:
			continue
		else:
			vars=line[9:]
			for i in range(len(vars)):
				if random.random()<=args.rate:
					gt=vars[i].split(':')
					gt[0]="./."
					gt=':'.join(gt)
					vars[i]=gt
			line=line[:9]+vars
			#print '\t'.join(line)
