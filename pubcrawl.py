#!/usr/bin/env python3
"""
*pubcrawl.py
	Workhorse script of the pubcrawl package. Takes in a TSV of word 2-pairs. Searches Pubmed for "AND"-joined results and downloads abstracts and Titles.
	Abstracts and titles are saved with a prepended 5 character identifier tag
	Downloaded papers are stored in files named after the terms in the following format:
	SPECIES_1#SPECIES_2.compiled
"""


from Bio import Entrez, Medline
from time import strftime, sleep
import os
from urllib.error import HTTPError
import multiprocessing
import argparse
import re

####################
# FIXED VARIABLES  #
####################
terms = ["PMID", "TI  ", "AB  "]
Entrez.email = "kmklim@gis.a-star.edu.sg"
searcher = re.compile("(^(?:" + "|".join(terms) +  ").*\n(?:[ \t]+.*\n)*)", flags =re.M)


#/////////////#
#    CHECKS   #
#/////////////#

if not Entrez.email:
	print("NO EMAIL DETECTED")
	raise

class medExtract():
	def __init__(self, rawPaper):
		self.data = {i: i for i in terms}
		for chunk in rawPaper:
			for term in terms:
				if chunk.startswith(term):
					self.data[term] = chunk

	def export(self):
		holder = [self.data[i] for i in terms]
		return "".join([self.checkEnd(i) for i in holder])

	def checkEnd(self, line):
		if line[-1] == '\n':
			return line
		return line + '\n'


"""
grabTerm
	Takes in string in the MEDLINE format containing all acquired papers and a list of terms to extract.

	returns a list of medExtract objects

	Note:
		The first term determines breakpoints for data separation.
"""
def grabTerm(medline, terms):
	result = []
	holder =  [i for i in medline.split('\n\n')]
	holder = [searcher.findall(i) for i in holder]
	holder = [[" ".join(chunk.split("\n      ")) for chunk in paper]for paper in holder]
	for paper in holder:
		result.append(medExtract(paper))
	return result
	
"""
pubmedSearch
	Main search function. Searches the terms and feeds the result into grabTerm

	Outputs all title-abstract pairs into files.
"""
def pubmedSearch(term1, term2, outDir, retryCount = 0):
	query = ' '.join([term1, "AND", term2])

	try:
		handle = Entrez.esearch(db = "pubmed", term = query , usehistory = "y")	
		record = Entrez.read(handle)
	except:
		if retryCount <3:
			return pubmedSearch(term1, term2, retryCount +1)
		else:
			print("ERROR#" + query + " FAILED TO DOWNLOAD")
	

	count = int(record["Count"])
	batchSize = 10
	outName = outDir+'_'.join((term1 + '#' + term2).split(' ')) + ".compiled" 
	holder = []

	#download the papers
	for start in range(0,count, batchSize):
		end = min(count, start+batchSize)
		print("Going to download record %i to %i" % (start+1, end))
		print(term1, term2)
		try:
			fetch_handle = Entrez.efetch(db = "pubmed",
											rettype = "medline", retmode= "text",
											retstart = start, retmax = batchSize,
											webenv = record["WebEnv"], query_key = record["QueryKey"] 
											)
		except:
			if retryCount < 3:
				return pubmedSearch(term1, term2, retryCount + 1)
			else:
				print("#ERROR#"+ query+ "FAILED TO DOWNLOAD")

		#Stepwise Raw data preprocessing
		data = fetch_handle.read()
		data = grabTerm(data, terms)
		fetch_handle.close()
		holder.extend(data)
		sleep(0.5)
	with open(outName, 'w') as f:
		for i in holder:
			f.write(i.export())

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument( "target", help ="Target file. File must be a line-separated list of tab separated term pairs. eg: Escherichia coli  Pseudomonas aeruginosa")
	parser.add_argument( "-c", "--cores", help ="number of cores", default = 4, type = int)
	parser.add_argument("-o", "--outdir", help ="Choose output directory. Default = output/pubcrawl/", default = "output/pubcrawl/")
	args = parser.parse_args()

	args.outdir = args.outdir + strftime("%Y-%m-%d-%H_%M") + "/"
	if not os.path.exists(args.outdir):
		os.makedirs(args.outdir)

	with open(args.target) as f:
		pairIn = [i.strip().split('\t') + [args.outdir] for i in f if i != '\n']

	pool = multiprocessing.Pool(args.cores)
	mappedRuns = pool.starmap(pubmedSearch, pairIn)	
		