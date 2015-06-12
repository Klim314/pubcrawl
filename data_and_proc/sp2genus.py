#!/usr/bin/env python3

import sys

target = sys.argv[1]
output = sys.argv[2]

holder =set()
with open(target) as f:
	for i in f:
		splat = i.strip().split("\t")
		holder.add((splat[0].split(' ')[0], splat[1].split(' ')[0]))

with open(output, 'w') as f:
	for i in holder:
		f.write("	".join(i) + '\n')


