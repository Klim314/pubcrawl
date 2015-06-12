from time import strftime


def readTsv(fileName):
	holder = []
	with open(fileName, 'r') as f:
		for i in f:
			holder.append(i.strip().split('\t'))
	return holder

def writeTSV(mat, fileName):
	length = len(mat)
	with open (fileName, 'w') as f:
		for lst in mat:
			f.write("\t".join(lst) + '\n')

a = readTsv("target_data.tsv")
b= [[ i[0], 'AND', i[1]] for i in a]

writeTSV(b, "and_target.in")