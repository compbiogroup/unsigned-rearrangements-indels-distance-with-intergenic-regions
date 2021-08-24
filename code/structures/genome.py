#cython: language_level=3

import utils as utils

class Genome:
    def __init__(self, pi, bpi, iota, biota):
        self.pi = []
        self.pi.append(0)
        self.pi = self.pi + pi
        self.pi.append(iota + 1)
        self.bpi = [0] + bpi
        
        self.iota = iota
        self.biota = [0] + biota
        self.updateAlphabets()
        self.updateInverse()
        self.operationsApplied = 0

    def sizeNonExtended(self):
        return len(self.pi) - 2

    def updateInverse(self):
        self.inverse = {}
        for i in range(0, len(self.pi)):
            if self.pi[i] in self.inverse:
                self.inverse[self.pi[i]].append(i)
            else:
                self.inverse[self.pi[i]] = [i]

    def updateAlphabets(self):
        self.alphabetA = set()
        for c in self.pi[1:len(self.pi)-1]:
            if isinstance(c, int):
                self.alphabetA.add(abs(c))
            else:
                self.alphabetA.add(c)
        self.alphabetB = set(list(range(1, self.iota+1)))
        self.common = self.alphabetA.intersection(self.alphabetB)
        self.alphabetPhi = self.alphabetB.difference(self.alphabetA)
        self.alphabetPsi = self.alphabetA.difference(self.alphabetB)

    def copy(self):
        return Genome(self.pi[1:len(self.pi)-1], self.bpi[1:len(self.bpi)], self.iota, self.biota[1:len(self.biota)])

    def isIdentity(self):
        #extended permutation
        if len(self.pi) != (self.iota + 2):
            return False
        for i in range(0, self.iota + 2):
            if self.pi[i] != i:
                return False
        for i in range(0, len(self.biota)):
            if self.biota[i] != self.bpi[i]:
                return False
        return True

    def isValid(self):
        g = sorted([abs(x) for x in self.pi if isinstance(x, int)])
        # print(g)
        for x in range(len(g)-1):
            if int(g[x+1]) - int(g[x]) > 2:
                print("Missing value from range:", g[x]+1, g[x+1]-1)
                return False 
        return True

    def printGenome(self):
        print("0", end = ' - ')
        for i in range(1, len(self.pi)):
            y = "(" + str(self.bpi[i]) + ")"
            x = self.pi[i]
            print(y, x, end = ' - ', sep = " - ")
        print()

    def printIota(self):
        print("0", end = ' - ')
        for i in range(1, len(self.biota)):
            y = "(" + str(self.biota[i]) + ")"
            print(y, i, end = " - ", sep = " - ")
        print()

    def compare(self, g2):
        if (len(g2.pi) != len(self.pi)):
            return False
        
        for i in range(0, len(self.pi)):
            if (self.pi[i] != g2.pi[i]):
                return False

        for i in range(0, len(self.bpi)):
            if (self.bpi[i] != g2.bpi[i]):
                return False
        
        return True

    def hasRBreakpoint(self, i):
        if i < 0 or i == len(self.pi) - 1:
            return utils.END
        if self.pi[i] == utils.alpha and self.pi[i+1] == utils.alpha:
            return utils.BLOCK
        elif self.pi[i] == utils.alpha or self.pi[i+1] == utils.alpha:
            return utils.BREAKPOINT

        if abs(self.pi[i] - self.pi[i+1]) != 1:
            return utils.BREAKPOINT
        if self.bpi[i+1] != self.biota[max(self.pi[i], self.pi[i+1])]:
            if self.bpi[i+1] > self.biota[max(self.pi[i], self.pi[i+1])]:
                return utils.OVERCHARGED
            else:
                return utils.UNDERCHARGED
        return utils.BLOCK

    def updateRBreakpoints(self):
        self.breakpoints = []
        self.overcharged = []
        self.undercharged = []
        self.intergenicBreakpoints = []
        for i in range(len(self.pi) - 1):
            adj = self.hasRBreakpoint(i)
            if adj == utils.BREAKPOINT:
                self.breakpoints.append(i)
                self.intergenicBreakpoints.append(i)
            elif adj == utils.OVERCHARGED:
                self.overcharged.append(i)
                self.intergenicBreakpoints.append(i)
            elif adj == utils.UNDERCHARGED:
                self.undercharged.append(i)
                self.intergenicBreakpoints.append(i)

        self.updateAlphabets()
        self.updateInverse()

    def hasTBreakpoint(self, i):
        if i < 0 or i == len(self.pi) - 1:
            return utils.END
        if self.pi[i] == utils.alpha and self.pi[i+1] == utils.alpha:
            return utils.BLOCK
        elif self.pi[i] == utils.alpha or self.pi[i+1] == utils.alpha:
            return utils.BREAKPOINT

        if self.pi[i+1] - self.pi[i] != 1:
            return utils.BREAKPOINT
        if self.bpi[i+1] != self.biota[max(self.pi[i], self.pi[i+1])]:
            if self.bpi[i+1] > self.biota[max(self.pi[i], self.pi[i+1])]:
                return utils.OVERCHARGED
            else:
                return utils.UNDERCHARGED
        return utils.BLOCK

    def updateTBreakpoints(self):
        self.breakpoints = []
        self.overcharged = []
        self.undercharged = []
        self.intergenicBreakpoints = []
        for i in range(len(self.pi) - 1):
            adj = self.hasTBreakpoint(i)
            if adj == utils.BREAKPOINT:
                self.breakpoints.append(i)
                self.intergenicBreakpoints.append(i)
            elif adj == utils.OVERCHARGED:
                self.overcharged.append(i)
                self.intergenicBreakpoints.append(i)
            elif adj == utils.UNDERCHARGED:
                self.undercharged.append(i)
                self.intergenicBreakpoints.append(i)

        self.updateAlphabets()
        self.updateInverse()

    def score(self):
        return len(self.alphabetPhi) + len(self.intergenicBreakpoints)
