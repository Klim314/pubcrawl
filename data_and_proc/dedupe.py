#!/usr/bin/env python3
import sys

target = sys.argv[1]

holder = set()

with open(target) as f:
	for i in f:
		temp = i.strip().split('\t')[::-1]
		temp1 = "\t".join(temp) 
		if temp1 in holder or temp[0] == temp[-1]:
			# a = "\t".join(i.split('\t')[::-1])
			# print(a in holder)
			# print(holder)
			print("FAIL")

			continue
		else:
			holder.add(i.strip())

holder = sorted(list(holder))

with open("out" + target, 'w') as f:
	for i in holder:
		f.write(i +'\n')