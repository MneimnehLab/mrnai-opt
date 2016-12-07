#!/usr/bin/python

import numpy as np
import itertools
import copy
import time
import datetime
import sys
import argparse
from collections import defaultdict


INFINITY = sys.maxint

def errorMsgs(*args):
    for msg in args:
        sys.stderr.write("ERROR: " + msg + "\n")


class SubpotEnum(object):
    def __init__(self, sequences, energies_file, loopAround, partMatrixFile, eps_frac, oneSidedGap, gapSize):
        self.sequences      = sequences
        self.energies_file  = energies_file
        self.loopAround     = loopAround
        self.ns             = [len(seq) for seq in self.sequences]
        self.partMatrixFile = partMatrixFile
        self.eps_frac       = eps_frac
        self.oneSidedGap    = oneSidedGap
        self.gapSize        = gapSize

        self.ns_plus_1      = map(lambda x: x+1, self.ns)
        self.levels         = len(self.ns)
        self.numRNAs        = len(self.ns)        
        self.RNAs           = list(range(len(self.ns)))

        if loopAround:
            self.levelsEnum = list(range(self.levels))
            self.RNAs.append(0)
        else:
            self.levelsEnum = list(range(self.levels-1))

        self.stack             = []
        self.acceptableConfigs = []

        # read all data files
        self.readEnergies()
        self.readPartial(self.partMatrixFile)
        
        # find subopts
        self.findSubOptimalConfigs()


    def readPartial(self, filename):
        self.partialOpt = {}
        mn = 0

        try:
            with open(filename) as matrx:
                for line in matrx:
                    data    = line.split()
                    indices = map(int, data[:-1])
                    score   = float(data[-1])
                    self.partialOpt[tuple(indices)] = score
                    mn      = min(mn, score)

            zeroIndex = tuple([0] * len(indices))
            self.partialOpt[zeroIndex] = 0.0
            self.opt = mn
            
        except IOError:
            errorMsgs("Cannot read partial weights file '" + filename + "'!",
                      "Please supply filename using -m")
            sys.exit(1)


    def getBoundary(self, config):
        smallest = {k:INFINITY for k,val in enumerate(self.ns_plus_1)}
        # smallest = {k:val for k,val in enumerate(ns_plus_1)}

        for (l,i,j,w1,w2) in config:
            smallest[l] = min(smallest[l], i-w1)
            smallest[(l+1)%len(self.ns_plus_1)] = min(smallest[(l+1)%len(self.ns_plus_1)], j-w2)

        return [smallest[l] for l in range(len(self.ns_plus_1))]


    def getFirstInLevel(self, config):
        firstWins = {k:None for k in range(len(self.ns_plus_1))}

        for win in config:
            (l,i,j,w1,w2) = win
            if firstWins[l]:
                firstWins[l] = min(firstWins[l], win)
            else:
                firstWins[l] = win

        return firstWins


    def preprocessWindows(self, weightsMatrix, allWindows, partialOpt, numRNAs):
        # for each window w, create a list of all windows w' that appear after it
        # and sort that list by weight(w') + partialOpt(w'.boundary)
        nextWins = defaultdict(list)
        weights = weightsMatrix
        for (l,i,j,w1,w2) in allWindows:
            thisList = []
            for (lA,iA,jA,w1A,w2A) in allWindows:
                if lA != l:
                    continue

                # if winA appears after win then add to list
                if iA < i-w1 and jA < j-w2:
                    actualBoundary = [0] * numRNAs
                    actualBoundary[l] = iA
                    actualBoundary[(l+1) % numRNAs] = jA

                    score = partialOpt[tuple(actualBoundary)] + weights[(lA,iA,jA,w1A,w2A)]

                    thisList.append( ((lA,iA,jA,w1A,w2A), score) )

            # thisList.sort(key = lambda x: x[1] , reverse = True )
            thisList.sort(key = lambda x: x[1])

            nextWins[(l,i,j,w1,w2)] = thisList
        
        return nextWins


    def findSubOptimalConfigs(self):
        gap = self.gapSize
        
        t1 = time.time()
        st = datetime.datetime.fromtimestamp(t1).strftime('%Y-%m-%d %H:%M:%S')
        print "Started at %s" % st
        
        self.eps = -1*self.eps_frac * float(self.opt)
        upperBound = self.opt + self.eps
        sys.stderr.write("Optimal energy = " + str(self.opt) + "\n")
        sys.stderr.write("Finding structures in range " + str(self.opt) + " to " + str(upperBound) + "\n")
        
        allWindows = self.windows
        
        # nextWins = preprocessWindows(weightsMatrix, allWindows, partialOpt, numRNAs)
        
        weights = self.weightsMatrix

        self.initStackWithWins()
        
        while len(self.stack) > 0:
            config = self.stack.pop()

            totalWeight = sum(weights[window_] for window_ in config)    
            if totalWeight <= self.opt + self.eps:
                self.acceptableConfigs.append((config, totalWeight))
                
            configBoundary = self.getBoundary(config)
            remainingWindows = self.getRemainingWindows(config, configBoundary)
            
            while len(remainingWindows) > 0:
                window = remainingWindows.pop()
            
                (l,i,j,w1,w2) = window

                # check if it's the last terminal
                # i.e., there is no window win' below this such that both left pegs of win' touch the boundary
                if self.isTerminal(window, config, configBoundary, l):
                    
                    # At this point we should check whether adding this edge is viable 
                    # i.e., solution is within suboptimal range
                    
                    # Before that, we need to construct the left "boundary" of this structure!
                    configWithWin = config + [window]
                    boundaryWithThisWin = self.getBoundary(configWithWin)
                    leftOverBoundary = [min(x-1, y)  for (x,y) in zip(boundaryWithThisWin, self.ns)]

                    bestPossible = totalWeight + weights[window] + \
                            self.partialOpt[tuple(leftOverBoundary)]

                    if bestPossible <= self.opt + self.eps:
                        # push this new config, terminals on the stack
                        updatedConfig = config + [window]
                        self.stack.append(updatedConfig)
                        
                    else:
                        # should we terminate search here?
                        pass
        
        self.acceptableConfigs.sort(key=lambda x:x[1], reverse=False)
        
        print "\n\n#Acceptable configs = "
        for config,weight in self.acceptableConfigs:
            print config, "\t",weight

        timeElapsed = time.time() - t1
        sys.stderr.write("Time elapsed = " + str(timeElapsed) + "\n")
    
    def initStackWithWins(self):
        for window in self.windows:
            actualBoundary = self.getBoundary([window])
            leftOverBoundary = [min(x-1, y)  for (x,y) in zip(actualBoundary, self.ns)]

            if self.partialOpt[tuple(leftOverBoundary)] + self.weightsMatrix[window] <= self.opt + self.eps:
                # stack.append( ([window], {}, {}) ) 
                self.stack.append( [window]) 
    

    def getRemainingWindows(self, config, configBoundary):
        remainingWindows = []
            
        finalWins = self.getFirstInLevel(config)
        
        # to find remaining windows, look at windows on this level to the left of the boundary of the config
        for level in self.levelsEnum:
            bi = configBoundary[self.RNAs[level]]
            bj = configBoundary[self.RNAs[level+1]]
            
            if finalWins[level] is not None:
                fl, fi, fj, fw1, fw2 = finalWins[level]

                if self.oneSidedGap == True:
                    candidateWindows = [(l,i,j,w1,w2) for (l,i,j,w1,w2) in self.windowsInLevel[level] \
                            if ((i < fi-self.gapSize-fw1 and j < fj-fw2) or (i < fi-fw1 and j < fj-self.gapSize-fw2))\
                            and i < bi and j < bj]
                else:
                    candidateWindows = [(l,i,j,w1,w2) for (l,i,j,w1,w2) in self.windowsInLevel[level] if i < fi-self.gapSize-fw1 and j < fj-self.gapSize-fw2 and i < bi and j < bj]
                
            else:
                candidateWindows = [(l,i,j,w1,w2) for (l,i,j,w1,w2) in self.windowsInLevel[level] if i < bi and j < bj]

            remainingWindows.extend(candidateWindows)

        return remainingWindows


    def isTerminal(self, window, config, configBoundary, l):
        for level in self.levelsEnum[l+2:]:
            for win in config:
                if win[0] == level:
                    bi = configBoundary[self.RNAs[level]]
                    bj = configBoundary[self.RNAs[level+1]]
                    
                    if (win[1]-win[3], win[2]-win[4]) == (bi, bj):
                        return False

        return True


    def readEnergies(self):
        ns = self.ns
        maxSize = max(ns)+1
        nums = len(ns) if self.loopAround else len(ns)-1
        self.weightsMatrix = np.zeros((nums, maxSize, maxSize, 26, 26), dtype=np.float32)
        self.windows = []
        self.windowsInLevel = defaultdict(list)
    

        sys.stderr.write("Reading energies... ")
        weightsFilename = self.energies_file
        try:
            with open(weightsFilename) as weightsFile:
                for line in weightsFile:
                    if len(line) == 0: break
                    
                    data = line.split()
                    data[:5] = map(int, data[:5])

                    level, i, j = data[0], data[2], data[4]
                    w1 = data[2] - data[1]
                    w2 = data[4] - data[3]
                    
                    if w1 != w2 or w1 == 0: continue

                    deltaEnergy = float(data[5])

                    if deltaEnergy >= 0:
                        continue

                    self.weightsMatrix[(level, i, j, w1, w2)] = deltaEnergy
                    self.windows.append((level, i, j, w1, w2))
                    self.windowsInLevel[level].append((level, i, j, w1, w2))
    
        
        except IOError:
            errorMsgs("Incorrect file path for energies: " + weightsFilename + "!",\
                        "Please supply correct path using -f")
            sys.exit(1)

        sys.stderr.write("Done.\n")


def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--partial-matrix', dest='matrix', default='partial_matrix')
    parser.add_argument('-f', '--energies-file', dest='energies_file', default='energies')
    parser.add_argument('-e', '--epsilon', dest='eps', default=0.25, type=float)
    parser.add_argument('-g', '--gap-sides', dest='oneSidedGap', default=1, type=int)
    parser.add_argument('-l', '--loop-around', dest='loopAround', default=1, type=int)
    args = parser.parse_args()

    inp         = raw_input("")
    sequences   = inp.strip().split("&")
    loopAround  = args.loopAround == 1   
    eps_frac    = args.eps
    oneSidedGap = args.oneSidedGap
    gapSize     = 1
    energies_file  = args.energies_file     
    partMatrixFile = args.matrix
    
    

    subopt = SubpotEnum(sequences, energies_file, loopAround, partMatrixFile, eps_frac, oneSidedGap, gapSize)

if __name__ == "__main__":
    main()
