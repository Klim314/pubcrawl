#!/usr/bin/env python3

from Bio import Entrez
from Bio import Medline
import sys

target = sys.argv[1]

def pubmedSearch(term1, term2):
	query = ' '.join([term1, "AND", term2])
	Entrez.email = "kmklim@gis.a-star.edu.sg"
	handle = Entrez.esearch(db = "pubmed", term = query , usehistory = "y")	
	record = Entrez.read(handle)

	count = int(record["Count"])
	batchSize = 10
	outName = '_'.join((term1 + '_' + term2).split(' ')) + ".compiled"
	output = open(outName, 'w')
	for start in range(0,count, batchSize):
		end = min(count, start+batchSize)
		print("Going to download record %i to %i" % (start+1, end))
		fetch_handle = Entrez.efetch(db = "pubmed",
										rettype = "medline", retmode= "text",
										retstart = start, retmax = batchSize,
										webenv = record["WebEnv"], query_key = record["QueryKey"] 
										)
		data = fetch_handle.read()
		fetch_handle.close()
		output.write(data)
	output.close()

with open(target) as f:
	lst = [i.split('\t') for i in f.read().split('\n')]

for i in lst:
	pubmedSearch(i[0], i[1])
