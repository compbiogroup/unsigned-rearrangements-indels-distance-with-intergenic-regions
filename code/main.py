import sys
import math
import getopt
import itertools
import statistics
import utils as utils
import operations as op
import alg_reversal_indel as alg_ri
import alg_transposition_indel as alg_ti
import alg_reversal_transposition_indel as alg_rti

from genome import Genome

utils.debug = False

fileinput = sys.argv[1]
problem = sys.argv[2]


print("Lower Bound", "Approximation", "Operations", "Indels", "Reversals", "Transpositions", sep = "\t")

data_lb = []
data_approx = []
data_operations = []
with open(fileinput) as file:
	for line in file:
		pi, bpi, biota = line.split()
		pi = [abs(int(x)) for x in pi.split(",")]
		for i in range(len(pi)):
			if pi[i] == 0:
				pi[i] = utils.alpha
		bpi = [abs(int(x)) for x in bpi.split(",")]
		biota = [abs(int(x)) for x in biota.split(",")]
		
		genome = Genome(pi, bpi, len(biota) - 1, biota)

		if problem == "r":
			genome.updateRBreakpoints()
			score = genome.score()
			lower_bound = math.ceil(score/2)

			alg_ri.sortByRearrangements(genome)
			alg_distance = genome.operationsApplied
		elif problem == "t":
			genome.updateTBreakpoints()
			score = genome.score()
			lower_bound = math.ceil(score/3)

			alg_ti.sortByRearrangements(genome)
			alg_distance = genome.operationsApplied
		elif problem == "rt":
			genome.updateRBreakpoints()
			score = genome.score()
			lower_bound = math.ceil(score/3)

			alg_rti.sortByRearrangements(genome)
			alg_distance = genome.operationsApplied

		approximation = alg_distance / lower_bound

		data_lb.append(lower_bound)
		data_approx.append(approximation)
		data_operations.append(alg_distance)
		nof_indels = genome.descOperationsApplied.count("d") + genome.descOperationsApplied.count("i")
		nof_reversals = genome.descOperationsApplied.count("r")
		nof_transpositions = genome.descOperationsApplied.count("t")
		
		print(lower_bound, approximation, alg_distance, nof_indels, nof_reversals, nof_transpositions, sep="\t")

print("Distance", max(data_operations), min(data_operations), statistics.mean(data_operations), sep="\t")
print("Approximation", max(data_approx), min(data_approx), statistics.mean(data_approx), sep="\t")
print(max(data_operations), min(data_operations), statistics.mean(data_operations), max(data_approx), min(data_approx), statistics.mean(data_approx), sep = "\t")


