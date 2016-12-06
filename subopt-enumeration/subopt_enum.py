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
gap = 1

          
# This returs a list of windows
# A window is represented by a tuple (l,i,j,w1,w2) where i and j are hte END points of hte window, 
# and l is the level, w is the size of the window ... therefore [(i, i+w1), (j, j+w2)]
# def getAllWindows(weights, ns):    
#     return [win for win in weights if weights[win] < 0]


def errorMsgs(*args):
    for msg in args:
        sys.stderr.write("ERROR: " + msg + "\n")


def readPartial(filename):
    partialOpt = {}
    mn = 0

    try:
        with open(filename) as matrx:
            for line in matrx:
                data    = line.split()
                indices = map(int, data[:-1])
                score   = float(data[-1])
                partialOpt[tuple(indices)] = score
                mn      = min(mn, score)

        zeroIndex = tuple([0] * len(indices))
        partialOpt[zeroIndex] = 0.0
        return partialOpt, mn

    except IOError:
        errorMsgs("Cannot read partial weights file '" + filename + "'!",
                  "Please supply filename using -m")
        sys.exit(1)


def getBoundary(config, ns_plus_1):
    smallest = {k:INFINITY for k,val in enumerate(ns_plus_1)}
    # smallest = {k:val for k,val in enumerate(ns_plus_1)}

    for (l,i,j,w1,w2) in config:
        smallest[l] = min(smallest[l], i-w1)
        smallest[(l+1)%len(ns_plus_1)] = min(smallest[(l+1)%len(ns_plus_1)], j-w2)

    return [smallest[l] for l in range(len(ns_plus_1))]


def getFirstInLevel(config, ns_plus_1):
    firstWins = {k:None for k in range(len(ns_plus_1))}

    for win in config:
        (l,i,j,w1,w2) = win
        if firstWins[l]:
            firstWins[l] = min(firstWins[l], win)
        else:
            firstWins[l] = win

    return firstWins


def findSubOptimalConfigs(weightsMatrix, windows, partMatrixFile, ns, eps_frac, loopAround=True, oneSidedGap=True):
    global gap
    
    t1 = time.time()
    st = datetime.datetime.fromtimestamp(t1).strftime('%Y-%m-%d %H:%M:%S')
    print "Started at %s" % st
    
    partialOpt, opt = readPartial(partMatrixFile)
    eps = -1*eps_frac * float(opt)
    upperBound = opt + eps
    sys.stderr.write("Optimal energy = " + str(opt) + "\n")
    sys.stderr.write("eps frac =" + str(eps_frac) + "\n")
    sys.stderr.write("eps =" + str(eps) + "\n")
    sys.stderr.write("Finding structures in range " + str(opt) + " to " + str(upperBound) + "\n")
    sys.stderr.write("Wrap around = " + str(loopAround) + "\n")
    
    
    setAllWindows = windows
    allWindows = setAllWindows
    
    
    ns_plus_1 = map(lambda x: x+1, ns)
    
    levels  = len(ns)
    numRNAs = len(ns)
    
    RNAs = list(range(len(ns)))

    if loopAround:
        levelsEnum = list(range(levels))
        RNAs.append(0)
    else:
        levelsEnum = list(range(levels-1))

    windowsInLevel = {}
    for level in levelsEnum:
        windowsInLevel[level] = []
    
    # for l in range(levels-1):
    for l in levelsEnum:
        # print "l =", l
        #for i in range(1,(ns[l]+1)) for j in range(1,(ns[l+1]+1)) for w1 in range(1,26) for w2 in range(1,26) if i > w1 and j > w2]:
        #windowsInLevel[l].append((l,i,j,w1,w2))
        windowsInLevel[l].extend(   [(l_,i,j,w1,w2) for (l_,i,j,w1,w2) in allWindows if l_ == l]   )
        if l == 3:
            print [(l_,i,j,w1,w2) for (l_,i,j,w1,w2) in allWindows if l_ == l]
    
    
    
    '''
    for l in range(3):
        for el in edgesInLevel[l]:
            print el
    '''
        
    '''
    With this insertion method, we are enforcing the following order on edges:
    (1,1) (1,2)... (1,n) (2,1), (3,1)... (n,1) (2,2)....(2,n) (3,2) (4,2)... (n,2) and so on
    This should reduce the number of lookups but it's not helping!
    
    for l in range(levels-1):
        for i in range(1,ns[l]+1):
            for j in range(i,ns[l+1]+1):
                
                edgesInLevel[l].append((l,i,j))
                
            for j in range(i+1,n+1):
                
                edgesInLevel[l].append((l,j,i))
    '''       


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
                # print (lA,iA,jA,w1A,w2A)
                actualBoundary = [0] * levels
                actualBoundary[l] = iA
                actualBoundary[(l+1) % numRNAs] = jA

                # print actualBoundary
                score = partialOpt[tuple(actualBoundary)] + weights[(lA,iA,jA,w1A,w2A)]

                thisList.append( ((lA,iA,jA,w1A,w2A), score) )

        # thisList.sort(key = lambda x: x[1] , reverse = True )
        thisList.sort(key = lambda x: x[1])

        nextWins[(l,i,j,w1,w2)] = thisList
        # print "nextWins for ", (l,i,j,w1,w2), " are: ", thisList
    
    # print "created list"
    
    
    

    stack = []

    for window in allWindows:
        actualBoundary = getBoundary([window], ns_plus_1)
        # print actualBoundary, [window]
        leftOverBoundary = [min(x-1, y)  for (x,y) in zip(actualBoundary, ns)]

        if partialOpt[tuple(leftOverBoundary)] + weights[window] <= opt + eps:
            stack.append( ([window], {}, {}) ) 

    stackLens = []
    added = []
    acceptableConfigs = []
    
    while len(stack) > 0:
        
        # print '\n\n\n'
        # print '------------------------------------------------'
        # print '\n\n\n'
        # print stack
        addedCount = 0
        stackLens.append(len(stack))
        
        config, terminals, finals = stack.pop()
        
        totalWeight = sum(weights[window_] for window_ in config)
        
        if totalWeight <= opt + eps:
            acceptableConfigs.append((config, totalWeight))
            
        remainingWindows = []
        
        configBoundary = getBoundary(config, ns_plus_1)
        finalWins = getFirstInLevel(config, ns_plus_1)
        
        # print 'curr config =', config
        # print "configBoundary =", configBoundary

        # to find remaining windows, look at windows on this level to the left of the boundary of the config
        for level in levelsEnum:
            bi = configBoundary[RNAs[level]]
            bj = configBoundary[RNAs[level+1]]
            
            if finalWins[level] is not None:
                fl, fi, fj, fw1, fw2 = finalWins[level]

                if oneSidedGap == True:
                    candidateWindows = [(l,i,j,w1,w2) for (l,i,j,w1,w2) in windowsInLevel[level] \
                            if ((i < fi-gap-fw1 and j < fj-fw2) or (i < fi-fw1 and j < fj-gap-fw2))\
                            and i < bi and j < bj]
                else:
                    candidateWindows = [(l,i,j,w1,w2) for (l,i,j,w1,w2) in windowsInLevel[level] if i < fi-gap-fw1 and j < fj-gap-fw2 and i < bi and j < bj]
                
            else:
                candidateWindows = [(l,i,j,w1,w2) for (l,i,j,w1,w2) in windowsInLevel[level] if i < bi and j < bj]

            remainingWindows.extend(candidateWindows)

       
        # print 'all rem wins =', remainingWindows
        while len(remainingWindows) > 0:
            window = remainingWindows.pop()
        
            (l,i,j,w1,w2) = window

            # print 'rem win =', window
            
            # check if it's the last terminal
            # i.e., there is no window win' below this such that both left pegs of win' touch the boundary
            
            isLastTerm = True
            # for level in range(l+2, levels-1):
            for level in levelsEnum[l+2:]:
                # print '**lev =', level
                for win in config:
                    if win[0] == level:
                        # bi, bj = configBoundary[level:level+2]
                        bi = configBoundary[RNAs[level]]
                        bj = configBoundary[RNAs[level+1]]
                        
                        if (win[1]-win[3], win[2]-win[4]) == (bi, bj):
                            isLastTerm = False
                            # print 'win =', win, ' touches boundary'
                            break
                    
            if isLastTerm:
                # print "for config", config, " ... ", window, "is last term"
                
                '''
                At this point we should check whether adding this edge is viable 
                i.e., solution is within suboptimal range
                
                Before that, we need to construct the left "boundary" of this structure!
                From above, we have the boundary of current configuration
                Add edge to it and compute actual boundary
                '''
                

                configWithWin = config + [window]
                boundaryWithThisWin = getBoundary(configWithWin, ns_plus_1)
                leftOverBoundary = [min(x-1, y)  for (x,y) in zip(boundaryWithThisWin, ns)]


                # print 'configWithWin =', configWithWin
                # print 'leftOverBoundary =', leftOverBoundary
                # print 'boundaryWithThisWin =', boundaryWithThisWin
                # print 'totalWeight =', totalWeight
                # print 'partialOpt[tuple(leftOverBoundary)] =', partialOpt[tuple(leftOverBoundary)]
                # print 'weights[window] =', weights[window]

                all = totalWeight + \
                        weights[window] + \
                        partialOpt[tuple(leftOverBoundary)]

                # print 'all =', all
                
                
                if all <= opt + eps:
                    # print "     ",
                    # print "fits opt condition",
                    # print "     ",
                    # print "ADDDDDD TERM"
                    
                    
                    #update config - cannot use ".add()" because it is by reference
                    config2 = list(config)
                    config2.append(window)
                    addedCount += 1
                    
                    #update terminals
                    terminals2 = terminals.copy()
                    terminals2[l] = window
                    
                    finals2 = finals.copy()
                    finals2[l] = window
                    
                    # if l > 0: terminals2[l-1] = ()
                    # if l < levels-1: terminals2[l+1] = ()
                    
                    #push this new config, terminals on the stack
                    stack.append( (config2, terminals2, finals2) )
                    # print "       ->    ", stack
                
                else:
                    # print "no more good wins, so breaking"
                    pass
            else:
                # print 'not last term'
                pass
                
            # print ''
            
            
        added.append(addedCount)
     
    
    acceptableConfigs.sort(key=lambda x:x[1], reverse=False)
    
    print "\n\n#Acceptable configs = "
    for config,weight in acceptableConfigs:
        print config, "\t",weight


    t2 = time.time()
    st2 = datetime.datetime.fromtimestamp(t2).strftime('%Y-%m-%d %H:%M:%S')
    
    # print "#End Acceptable"
    

    # print "Ended at   %s" % st2
    # print " "

    
    timeElapsed = t2 - t1
    sys.stderr.write("Time elapsed = ", str(timeElapsed),"\n")
    

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--partial-matrix', dest='matrix', default='partial_matrix')
    parser.add_argument('-f', '--energies-file', dest='energies_file', default='energies')
    parser.add_argument('-e', '--epsilon', dest='eps', default=0.25, type=float)
    parser.add_argument('-g', '--gap-sides', dest='oneSidedGap', default=1, type=int)
    parser.add_argument('-l', '--loop-around', dest='loopAround', default=1, type=int)
    args = parser.parse_args()


    # inp = sys.stdin.readlines()
    # sequences = inp[0].strip().split("&")

    inp = raw_input("")
    sequences = inp.strip().split("&")

    # print sequences

    ns = [len(seq) for seq in sequences]

    # read window weights file
    
    maxSize = max(ns)+1

    if args.loopAround == 1:
        weightsMatrix = np.zeros((len(ns), maxSize, maxSize, 26, 26), dtype=np.float32)
    else:
        weightsMatrix = np.zeros((len(ns)-1, maxSize, maxSize, 26, 26), dtype=np.float32)

    windows = []
    
    
    sys.stderr.write("Reading energies... ")
    weightsFilename = args.energies_file
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

                deltaEnergy = float(data[7])

                if deltaEnergy >= 0:
                    continue

                weightsMatrix[(level, i, j, w1, w2)] = deltaEnergy
                windows.append((level, i, j, w1, w2))
    
    except IOError:
        sys.stderr.write("\nIncorrect file path for energies! \nPlease supply correct path using -f\n")
        return

    sys.stderr.write("Done.\n")

    
    

    findSubOptimalConfigs(weightsMatrix, windows, args.matrix, ns, args.eps, \
            args.loopAround == 1 , (args.oneSidedGap == 1))


    
    

if __name__ == "__main__":
    main()
