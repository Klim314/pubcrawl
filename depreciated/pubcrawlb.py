#!/usr/bin/env python3
# Takes in a tsv consisting of multiple lines of TSV'd terms. Searches an "AND" joined result and collects abstracts and Titles
# ./pubcrawl.py <FILENAME>
# current bugs:
# grabterm may bug up if repeated entries present
# AB  -
# AB  -
# TODO: Modify for dynamic inputs in CLI


from Bio import Entrez, Medline

from time import strftime, sleep
import os
import sys
from urllib.error import HTTPError

####################
#USER SET VARIABLES#
####################

terms = ["TI", "AB"]



target = sys.argv[1]
prefix = "output/" + strftime("%Y-%m-%d-%H_%M") + "/"

#grabs specified terms from from each article
def grabTerm(medline, terms):
    #list of grabbable terms
    #First (Start) term determines breakpoints for categorization
    #ALL ENTRIES MUST HAVE A UNIQUE START TERM
    terms =[i + (4-len(i)) * ' ' for i in terms]
    #output variables
    holder = []
    res = []

    curIndex = 0
    flag = -1
    termDict = dict([(j, i) for i, j in enumerate(terms)])
    #print(terms)
    while curIndex != len(medline):
        whole = medline[curIndex]
        cur = whole[:4]
        #very start or start of a new entry
        if cur == terms[0]:
            if holder != []:
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
            while medline[curIndex][0] == ' ':
                #print(medline[curIndex])
                holder[flag] += ' ' + medline[curIndex].strip()
                curIndex+=1
            flag = -1

    res.append(holder)
    return res


def preProc(data):
	return data.split('\n')

def pubmedSearch(term1, term2, retryCount = 0):
	query = ' '.join([term1, "AND", term2])
	Entrez.email = "kmklim@gis.a-star.edu.sg"
	handle = Entrez.esearch(db = "pubmed", term = query , usehistory = "y")	
	record = Entrez.read(handle)
	count = int(record["Count"])
	batchSize = 10
	outName = prefix + '_'.join((term1 + '_' + term2).split(' ')) + ".compiled" 
	if not os.path.exists(prefix):
		os.makedirs(prefix)
	output = open(outName, 'w')
	for start in range(0,count, batchSize):
		end = min(count, start+batchSize)
		print("Going to download record %i to %i" % (start+1, end))
		try:
			fetch_handle = Entrez.efetch(db = "pubmed",
											rettype = "medline", retmode= "text",
											retstart = start, retmax = batchSize,
											webenv = record["WebEnv"], query_key = record["QueryKey"] 
											)
		#may create false entries
		# except HTTPError:
		# 	print("ERROR")
		# 	print("retrying, try: ", count)
		# 	print("searching : ", term1, " | ", term2)
		# 	if retryCount <3:
		# 		return pubmedSearch(term1, term2, retryCount+1)
		# 	else:
		# 		raise

		#Stepwise Raw data preprocessing
		data = fetch_handle.read()
		data = preProc(data)
		data = grabTerm(data, terms)
		
		fetch_handle.close()
		data = sum(data, [])
	
		data = "\n".join(data)

		
		output.write(data + '\n')
		sleep(1)
	output.close()

with open(target) as f:
	lst = [i.split('\t') for i in f.read().split('\n')]

for i in lst:
	pubmedSearch(i[0], i[1])
	


