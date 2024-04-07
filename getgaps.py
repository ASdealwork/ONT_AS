#!/usr/bin/env python3

# Import necessary packages
import argparse
import re
from Bio import SeqIO

# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument("fasta")
args = parser.parse_args()

with open(args.fasta) as handle:
	i = 0
	for record in SeqIO.parse(handle, "fasta"):
		for match in re.finditer('N+', str(record.seq)):
			i = i+1
			print (record.id, ".", "gap", match.start() + 1, match.end(), ".", ".", ".", "Name=gap" + str(i) + ";size=" + str(match.end()-match.start()), sep='\t')

#use the following at CMD: FILENAME.py FILENAME.fasta >> FILENAME.gff3  here
