#cython: language_level=3

import math
import utils as utils
import operations as op
from genome import Genome

def findStrip(genome, i):
	j = i
	while genome.hasTBreakpoint(j) == utils.BLOCK:
		j += 1
	while genome.hasTBreakpoint(i-1) == utils.BLOCK:
		i -= 1
	return (i,j)

def findStripNonIntergenic(genome, i):
	j = i
	while genome.hasTBreakpoint(j) in [utils.BLOCK, utils.UNDERCHARGED, utils.OVERCHARGED]:
		j += 1
	while genome.hasTBreakpoint(i-1) in [utils.BLOCK, utils.UNDERCHARGED, utils.OVERCHARGED]:
		i -= 1
	return (i,j)

def removeOvercharged(genome, i):
	#we need to remove nucleotides from region (i + 1)
	a = genome.pi[i]
	b = genome.pi[i+1]
	x = genome.biota[max(a,b)]
	y = genome.bpi[i+1]
	op.deletion(genome, i+1, i+1, x, y)

def removeUndercharged(genome, i):
	#we need to insert nucleotides in region (i + 1)
	x = 0
	sigma = []
	a = genome.pi[i]
	b = genome.pi[i+1]
	bsigma = [genome.biota[max(a,b)] - genome.bpi[i+1]]
	op.insertion(genome, i, x, sigma, bsigma)


def sortByRearrangements(genome):
	genome.operationsApplied = 0 #restart number of rearrangements used
	genome.updateTBreakpoints()
	score = genome.score()
	nof_op = 0
	while not genome.isIdentity():
		if len(genome.alphabetPhi) > 0:
			newelement = min(genome.alphabetPhi)
			i = genome.inverse[(newelement - 1)][0]
			i, j = findStripNonIntergenic(genome, i)
			sigma = [newelement]
			#TODO optimize inserted intergenic regions
			bsigma = [0,0]
			if genome.pi[j] == newelement - 1:
				x = 0
				bsigma[0] = genome.biota[newelement]
				op.insertion(genome, j, x, sigma, bsigma)
			else:
				x = genome.bpi[i]
				bsigma[1] = genome.biota[newelement]
				op.insertion(genome, i-1, x, sigma, bsigma)
		elif len(genome.alphabetPsi) > 0:
			i = min(genome.inverse[utils.alpha])
			i, j = findStrip(genome, i)
			#TODO optimize removed intergenic regions
			x = genome.bpi[i]
			y = 0
			op.deletion(genome, i, j+1, x, y)
		elif len(genome.overcharged) > 0:
			i = genome.overcharged[0]
			removeOvercharged(genome, i)
		elif len(genome.undercharged) > 0:
			i = genome.undercharged[0]
			removeUndercharged(genome, i)
		else:
			#the list of intergenic breakpoints is sorted
			i_minus, j = genome.intergenicBreakpoints[0:2] #first two intergenic breakpoints
			i = i_minus + 1
			m = genome.inverse[genome.pi[j] + 1][0]
			n = genome.inverse[genome.pi[j+1] - 1][0]

			if genome.bpi[j+1] + genome.bpi[m] >= genome.biota[genome.pi[m]]:
				#exchange segments (pi_i ... pi_j) and (pi_{j+1} ... pi_{m-1})		
				if genome.bpi[m] < genome.biota[genome.pi[m]]:
					x = genome.bpi[i]
					y = genome.biota[genome.pi[m]] - genome.bpi[m]
					z = 0
					op.transposition(genome, i, j+1, m, x, y, z)
				else:
					x = genome.bpi[i]
					y = 0
					z = genome.bpi[m] - genome.biota[genome.pi[m]]
					op.transposition(genome, i, j+1, m, x, y, z)
			else:
				nucleotides = genome.biota[genome.pi[m]] - (genome.bpi[j+1] + genome.bpi[m]) 
				nucleotides += max(0, genome.biota[genome.pi[j+1]] - genome.bpi[n+1])
				op.insertion(genome, j, 0, [], [nucleotides])

				if utils.debug:
					genome.printGenome()

				#exchange segments (pi_i ... pi_j) and (pi_{j+1} ... pi_{m-1})
				if genome.bpi[m] < genome.biota[genome.pi[m]]:
					x = genome.bpi[i]
					y = genome.biota[genome.pi[m]] - genome.bpi[m]
					z = 0
					op.transposition(genome, i, j+1, m, x, y, z)
				else:
					x = genome.bpi[i]
					y = 0
					z = genome.bpi[m] - genome.biota[genome.pi[m]]
					op.transposition(genome, i, j+1, m, x, y, z)

				if utils.debug:
					genome.printGenome()

				genome.updateInverse()
				#the initial element pi_{j+1} is at position i in the new string pi'
				n = genome.inverse[genome.pi[i] - 1][0]
				_, k = findStrip(genome, i)
				#if (n+1) == i the second transposition is not necessary but a deletion may be necessary
				if (n+1) != i:
					if n+1 < i:
						#exchange segments (pi'_{n+1} ... pi'_{i-1}) and (pi'_{i} ... pi'_k)
						if genome.bpi[n+1] < genome.biota[genome.pi[i]]:
							x = genome.bpi[n+1]
							y = genome.bpi[i] - (genome.biota[genome.pi[i]] - genome.bpi[n+1])
							z = 0
							op.transposition(genome, n+1, i, k+1, x, y, z)
						else:
							x = genome.biota[genome.pi[i]]
							y = genome.bpi[i]
							z = 0
							op.transposition(genome, n+1, i, k+1, x, y, z)
					else:
						#exchange segments (pi'_i ... pi'_{k}) (pi'_{k+1} ... pi_n)
						if genome.bpi[i] > genome.biota[genome.pi[i]]:
							x = genome.bpi[i] - genome.biota[genome.pi[i]]
							y = 0
							z = 0
							op.transposition(genome, i, k+1, n+1, x, y, z)
						else:
							x = 0
							y = 0 
							z = genome.biota[genome.pi[i]] - genome.bpi[i]
							op.transposition(genome, i, k+1, n+1, x, y, z)

				elif genome.bpi[i] != genome.biota[genome.pi[i]]:
					removeOvercharged(genome, i-1)


		# update genome information
		genome.updateTBreakpoints()
		score_new = genome.score()

		if utils.debug:
			genome.printGenome()
			genome.printIota()
			print("=="*100)
		if (score-score_new)/(genome.operationsApplied-nof_op) < 2/3:
			print((score-score_new)/(genome.operationsApplied-nof_op))
			print("score", score)
			print("score_new", score_new)
			print("nof_op", nof_op)
			print("nof_op_new", genome.operationsApplied)
			raise "ERROR Number of Breakpoints"
		else:
			score = score_new
			nof_op = genome.operationsApplied
	return genome.operationsApplied