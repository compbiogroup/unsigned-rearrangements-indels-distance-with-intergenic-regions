import sys
import os

operations = [0.1, 0.3, 0.5, 1]
qtd = 1000

def write(file, line):
	with open(file, "a") as file_object:
	    file_object.write(line + "\n")

for n in range(50, 501, 50):
	for op in operations:
		for model in ["r","t","rt"]:
			file = "output/%s_%s_%s.out" % (model, n, op)
			with open(file) as f:
				lines = f.readlines()
				results = lines[len(lines)-1]
				stats_file = "results_%s_%s.stats" % (op, model)
				write(stats_file, str(n) + "\t" + results)
				# print(stats_file)