import sys
import os

operations = [0.1, 0.5, 1]
model = sys.argv[1]
results = []

for op in operations:
	file = "results_%s_%s.stats" % (op, model)
	with open(file) as f:
		lines = f.readlines()
		for line in lines:
			line = line.strip()
			if line != "":
				data = line.split("\t")
				operations = int(data[0]) * op
				data = [int(data[0]), int(operations)] + data[1:len(data)]
				results.append(data)

results.sort(key = (lambda x: x[1]))
results.sort(key = (lambda x: x[0]))

for line in results:
	print("\t".join([str(x) for x in line]))