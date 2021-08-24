#cython: language_level=3

import math
import utils as utils
import operations as op
from genome import Genome

def findStrip(genome, i):
	j = i
	while genome.hasRBreakpoint(j) == utils.BLOCK:
		j += 1
	while genome.hasRBreakpoint(i-1) == utils.BLOCK:
		i -= 1
	return (i,j)

def findStripNonIntergenic(genome, i):
	j = i
	while genome.hasRBreakpoint(j) in [utils.BLOCK, utils.UNDERCHARGED, utils.OVERCHARGED]:
		j += 1
	while genome.hasRBreakpoint(i-1) in [utils.BLOCK, utils.UNDERCHARGED, utils.OVERCHARGED]:
		i -= 1
	return (i,j)

def findDecreasingStrip(genome):
	i = 0
	for j in genome.intergenicBreakpoints:
		if j == i and (i > 0) and (j < len(genome.pi)-1):
			return (i,j)
		elif genome.pi[i] > genome.pi[j]:
			return (i,j)
		else:
			i = j + 1
	return (-1,-1) # there is no decreasing strips

def findDecreasingStripMinimum(genome):
	minimum = len(genome.bpi)
	i_minimum = -1
	j_minimum = -1
	i = 0
	for j in genome.intergenicBreakpoints:
		if j == i and (i > 0) and (j < len(genome.pi)-1):
			if genome.pi[j] < minimum:
				minimum = genome.pi[j]
				i_minimum = i
				j_minimum = j
		elif genome.pi[i] > genome.pi[j]:
			if genome.pi[j] < minimum:
				minimum = genome.pi[j]
				i_minimum = i
				j_minimum = j
		i = j + 1
	return (i_minimum, j_minimum)

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
	genome.updateRBreakpoints()
	score = genome.score()
	while not genome.isIdentity():
		flag_inverse_increasing_strip = False
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
			i, j = findDecreasingStripMinimum(genome)
			if i != -1:
				#there is at least one decreasing strip
				i_prime = genome.inverse[genome.pi[j] - 1][0]
				_, j_prime = findStrip(genome, i_prime)
				if i_prime < i:
					x = genome.bpi[j_prime + 1]
					y = 0
					#inserting enough nucleotides to remove intergenic breakpoints
					#this also handles the case where no decreasing strips are left after the reversal
					if genome.bpi[j_prime + 1] + genome.bpi[j + 1] < genome.biota[genome.pi[j]]:
						bsigma = (genome.biota[genome.pi[j]] + genome.biota[genome.pi[j+1]]) - (genome.bpi[j_prime + 1] + genome.bpi[j + 1])
						bsigma = [bsigma]
						op.insertion(genome, j_prime, 0, [], bsigma)
					if genome.bpi[j_prime + 1] >= genome.biota[genome.pi[j]]:
						x = genome.biota[genome.pi[j]]
						y = 0
					else:
						x = genome.bpi[j_prime + 1]
						y = genome.biota[genome.pi[j]] - genome.bpi[j_prime + 1]
					op.reversal(genome, j_prime + 1, j, x, y)
				else:
					x = genome.bpi[j+1]
					y = 0
					#inserting enough nucleotides to remove intergenic breakpoints
					if genome.bpi[j + 1] + genome.bpi[j_prime + 1] < genome.biota[genome.pi[j]]:
						aux = max(genome.pi[j+1], genome.pi[j_prime+1])
						bsigma = (genome.biota[genome.pi[j]] + genome.biota[aux]) - (genome.bpi[j + 1] + genome.bpi[j_prime + 1])
					if genome.bpi[j+1] >= genome.biota[genome.pi[j]]:
						x = genome.biota[genome.pi[j]]
						y = 0
					else:
						x = genome.bpi[j+1]
						y = genome.biota[genome.pi[j]] - genome.bpi[j+1]
					op.reversal(genome, j+1, j_prime, x, y)
			else:
				#there is no decreasing strips
				flag_inverse_increasing_strip = True
				first_breakpoint = genome.intergenicBreakpoints[0]
				#this is an increasing strip that is going to be transformed into a decreasing strip
				i_prime, j_prime = findStrip(genome, first_breakpoint + 1)
				x = genome.bpi[i_prime]
				y = 0
				op.reversal(genome, i_prime, j_prime, x, y)
		#update genome information
		genome.updateRBreakpoints()
		# genome.printGenome()
		score_new = genome.score()
		if score <= score_new and flag_inverse_increasing_strip == False:
			raise "ERROR Number of Breakpoints"
		else:
			score = score_new
	return genome.operationsApplied