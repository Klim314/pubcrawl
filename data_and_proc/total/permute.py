#!usr/bin/env python3
import itertools
a = ''
with open ('totalSpeciesList.tsv') as f:
	a = f.read()

res  = itertools.combinations(a.split('\t'), 2)

res = [i for i in res]
print(len(res))