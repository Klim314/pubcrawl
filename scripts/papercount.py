#!/usr/bin/env python
import os

targetFolder = "files"
outFile = targetFolder + ".output"

def lc(file):
	with open(file) as f:
		return(sum(1 for line in f))

def getNames(file):	
	file = os.path.basename(file)
	name = os.path.splitext(file)[0]
	name = name.split('_')
	return (" ".join(name[0:2]).strip(), " ".join(name[2:]).strip())

def writeTSV(mat, fileName):
	length = len(mat)
	with open (fileName, 'w') as f:
		for lst in mat:
			f.write("\t".join(lst) + '\n')

def matInt2Str(mat):
	for i in range(len(mat)):
		for j in range(len(mat[i])):
			if i != 0 and j != 0:
				mat[i][j] = str(mat[i][j])



files = os.listdir(targetFolder)
holder =[]
species = set()
for i in files:
	fileName = targetFolder + '/' + i
	lines = lc(fileName)
	spNames = getNames(fileName)
	species.add(spNames[0])
	species.add(spNames[1])
	holder.append((spNames, lines))

species = sorted(list(species))
#account for first row
maDict = {j:i+1 for i,j in enumerate(species)}
print(maDict)

#create empty matrix
matrix = [["sp/sp"] + [i for i in species]]



for i in species:
	matrix.append([i]+ [0 for i in species])

for i in holder:
	first, second = maDict[i[0][0]], maDict[i[0][1]]
	matrix[first][second] = i[1]
	matrix[second][first] = i[1]

matInt2Str(matrix)
writeTSV(matrix, outFile)