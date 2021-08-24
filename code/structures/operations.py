#cython: language_level=3

import utils as utils
from genome import Genome

def deletion(genome, i, j, x, y):
	if utils.debug:
		print("deletion", i, j, x, y)
	# x_prime = genome.bpi[i] - x
	y_prime = genome.bpi[j] - y
	genome.pi = genome.pi[0:i] + genome.pi[j:len(genome.pi)]
	genome.bpi = genome.bpi[0:i] + [x + y_prime] + genome.bpi[j+1:len(genome.bpi)]
	genome.operationsApplied += 1

def insertion(genome, i, x, sigma, bsigma):
	if utils.debug:
		print("insertion", i, x, sigma, bsigma)
	x_prime = genome.bpi[i+1] - x
	genome.pi = genome.pi[0:i+1] + sigma + genome.pi[i+1:len(genome.pi)]
	bsigma[len(bsigma) - 1] += x_prime
	genome.bpi = genome.bpi[0:i+1] + [x + bsigma[0]] + bsigma[1:len(bsigma)] + genome.bpi[i+2:len(genome.bpi)]
	genome.operationsApplied += 1

def reversal(genome, i, j, x, y):
	if utils.debug:
		print("reversal", i, j, x, y)
	x_prime = genome.bpi[i] - x
	y_prime = genome.bpi[j+1] - y
	segment = genome.pi[i:j+1]
	segment.reverse()
	genome.pi = genome.pi[0:i] + segment + genome.pi[j+1:len(genome.pi)]
	bsegment = genome.bpi[i+1:j+1]
	bsegment.reverse()
	genome.bpi = genome.bpi[0:i] + [x + y] + bsegment + [x_prime + y_prime] + genome.bpi[j+2:len(genome.bpi)]
	genome.operationsApplied += 1

def transposition(genome, i, j, k, x, y, z):
	if utils.debug:
		print("transposition", i, j, k, x, y, z)
	x_prime = genome.bpi[i] - x
	y_prime = genome.bpi[j] - y
	z_prime = genome.bpi[k] - z
	genome.pi = genome.pi[0:i] + genome.pi[j:k] + genome.pi[i:j] + genome.pi[k:len(genome.pi)]
	genome.bpi = genome.bpi[0:i] + [x + y_prime] + genome.bpi[j+1:k] + [z + x_prime] + genome.bpi[i+1:j] + [y + z_prime] + genome.bpi[k+1:len(genome.bpi)]
	genome.operationsApplied += 1