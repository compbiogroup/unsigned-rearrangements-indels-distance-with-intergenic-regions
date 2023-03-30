import sys
import os

operations = [0.3, 1]
qtd = 1000

for n in range(100, 101, 50):
	for op in operations:
		for model in ["rt"]:
			file = "Instances/input/%s_%s_%s.in" % (model, n, op)
			output = "output/%s_%s_%s.out" % (model, n, op)
			command = "python3 main.py %s %s > %s" % (file, model, output)
			print(command)
			os.system(command)