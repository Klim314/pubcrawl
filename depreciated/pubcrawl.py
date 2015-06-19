#!/usr/bin/env python3
"""
pubcrawl.py
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

####################
# FIXED VARIABLES  #
####################
terms = ["PMID", "TI", "AB"]
Entrez.email = "kmklim@gis.a-star.edu.sg"


#/////////////#
#  ARGUMENTS  #
#/////////////#

parser = argparse.ArgumentParser()
parser.add_argument( "target", help ="Target file. File must be a line-separated list of tab separated term pairs. eg: Escherichia coli  Pseudomonas aeruginosa")
parser.add_argument( "-c", "--cores", help ="number of cores", default = 4, type = int)
args = parser.parse_args()
target = args.target
cores = args.cores

prefix = "output/" + strftime("%Y-%m-%d-%H_%M") + "/"
if not os.path.exists(prefix):
	os.makedirs(prefix)

#/////////////#
#    CHECKS   #
#/////////////#

if not Entrez.email:
	print("NO EMAIL DETECTED")
	raise

"""
grabTerm
	Takes in string in the MEDLINE format containing all acquired papers and a list of terms to extract.

	Returns a list of papers 

	Note:
		The first term determines breakpoints for data separation. 
		ALL ENTRIES MUST HAVE THE FIRST TERM
"""
def grabTerm(medline, terms):
	medline = medline.split('\n')
	#Pad out the terms with spaces so as to render them identical to the dictionary
	terms =[i + (4-len(i)) * ' ' for i in terms]

	#output variables
	holder = []
	res = []

	curIndex = 0
	flag = -1
	#store the indexes of the term entries
	termDict = dict([(j, i) for i, j in enumerate(terms)])
	#print(terms)
	while curIndex != len(medline):
		whole = medline[curIndex]
		cur = whole[:4]
		#very start or start of a new entry
		if cur == terms[0]:
			if holder != []:
				#check if all data fields required are present
				temp = [i[:4] for i in holder]
				for i in terms:
					if i not in temp:
						holder.insert(termDict[i], i)
				res.append(holder)
				holder = []
			holder.append(whole)
			#flag to add all tabbed lines after this    
			flag = termDict[cur]
		
		elif cur in terms:
			holder.append(whole)
			flag = termDict[cur]

		curIndex +=1 

		#grab all tabbed lines belonging to the flagged line
		if flag != -1 and medline[curIndex] != '':
			while curIndex < medLen and medline[curIndex][0] == ' ':
				if flag > len(holder) or curIndex > len(holder):
					print("HOLDER: ", holder )
					print("FLAG: ", flag)
					pritn("MEDLINE", medline)
				try:
					holder[flag] += ' ' + medline[curIndex].strip()
				except:
					print("HOLDER: ", holder )
					print("FLAG: ", flag)
					raise
				curIndex+=1
			flag = -1

	temp = [i[:4] for i in holder]
	for i in terms:
		if i not in temp:
			holder.insert(termDict[i], i)
	res.append(holder)
	return res


"""
pubmedSearch
	Main search function. Searches the terms and feeds the result into grabTerm

	Outputs all title-abstract pairs into files.
"""
def pubmedSearch(term1, term2, retryCount = 0):
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
	outName = prefix + '_'.join((term1 + '#' + term2).split(' ')) + ".compiled" 
	output = ''

	#download the papers
	for start in range(0,count, batchSize):
		end = min(count, start+batchSize)
		# print("Going to download record %i to %i" % (start+1, end))
		# print(term1, term2)
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
		data = sum(data, [])
	
		data = "\n".join(data)

		
		output+= data + '\n'
		sleep(0.5)
	with open(outName, 'w') as f:
		f.write(output)

if __name__ == "__main__":

	with open(target) as f:
		lst = [i.strip().split('\t') for i in f if i != '\n']

	pool = multiprocessing.Pool(cores)
	print(lst)
	mappedRuns = pool.starmap(pubmedSearch, lst)