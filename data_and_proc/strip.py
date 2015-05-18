#!/usr/bin/env python3
import sys

target = sys.argv[1]


with open(target) as f:
	holder = ['\t'.join(i.split('\t')[:2]) for i in f.read().split('\n')]

holder = [i.split('\t') for i in holder]
for i in holder:
	if i == ['']:
		continue
	i[0] = ' '.join(i[0].split(" ")[:2])
	i[1] = ' '.join(i[1].split(" ")[:2])

holder = ["	".join(i) for i in holder]
print(holder)

with open("out" + target, 'w') as f:
	for i in holder:
		f.write(i +'\n')