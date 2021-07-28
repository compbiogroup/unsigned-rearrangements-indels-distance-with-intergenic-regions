import sys
import os

operations = [0.1, 0.3, 0.5, 1]
qtd = 1000

for n in range(50, 501, 50):
	for op in operations:
		for model in ["r","t","rt"]:
			file = "Instances/input/%s_%s_%s.in" % (model, n, op)
			command = "python3 main.py %s %s" % (file, model)
			print(command)
			os.system(command)