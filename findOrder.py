#!/usr/bin/python

# TODO: Clean the directories in alignment that are not optimal
# capture control+c and stop the threads?

from scipy import cluster
import numpy as np
import os
from multiprocessing import Pool
import itertools

def getCount(path, value="nbAlis"):
    # return the number of alignments
    # can also return total number of aligned bases
    # or other values in the future (% (matching) bases?)
    if value== "nbAlis":
        return len(os.listdir(path))
    elif  value== "aliLengths":
        total = 0
        for ali in os.listdir(path):
            total += int(ali.split("_")[2])
        return total


def getNamesList(pairwiseDirsPath):
    namesSet = set()
    for pair in os.listdir(pairwiseDirsPath):
        namesSet.update(pair.split("_"))

    return list(namesSet)


def returnBestPutativePerm(pairwiseDirsPath, value = "nbAlis", measure= "score"):
    ''' 
    Returns the most likely permutation 
    '''

    namesList = getNamesList(pairwiseDirsPath)
    tempNamesHash = dict([(namesList[i],i) for i in range(len(namesList))])

    distMatrix =  np.matrix([[0 for i in range(len(namesList))] for j in range(len(namesList))])

    for i in range(len(namesList)-1):
        
        for j in range(i+1, len(namesList)):
            if os.path.exists(os.path.join(pairwiseDirsPath, namesList[i]+"_"+namesList[j])):
                distMatrix[i,j]=getCount(os.path.join(pairwiseDirsPath, namesList[i]+"_"+namesList[j]), value)
                distMatrix[j,i]= distMatrix[i,j]
            else:
                distMatrix[i,j]=getCount(os.path.join(pairwiseDirsPath, namesList[j]+"_"+namesList[i]), value)
                distMatrix[j,i]=distMatrix[i,j]
    #linkaegeMatrix = cluster.hierarchy.average(distMatrix)
    order=[] 
    print distMatrix

    maxValue = distMatrix[distMatrix > 0].max()
    order.append(namesList[np.where(distMatrix == maxValue)[0][0,0]])
    order.append(namesList[np.where(distMatrix == maxValue)[1][0,0]])


    tempNamesHash.pop(order[0])
    tempNamesHash.pop(order[1])

    while len(tempNamesHash) > 0 : 
        if measure  == "score":
            measureSoFar = 0
        # else:
        # measureSoFar =  float("inf")
        index = -1
        for i in tempNamesHash.keys():
            avgDist = 1.0 * sum([distMatrix[namesList.index(x),tempNamesHash[i]] for x in order])/len(order) 
        
            if avgDist > measureSoFar:
                index = tempNamesHash[i]
                measureSoFar = avgDist
        order.append(namesList[index]) 
        tempNamesHash.pop(namesList[index])
    return order



def testPermutation(pairwiseDirsPath, genomesDir, order,  measure= "score", minScore=None):
    # TODO: Clean and re-write this to use Seanome's find_csr instead os.system
    # if minScore not None, stop if measure is exceeded
    
    seanomeDir = "/home/mahdi/work/tobo_lab/new_pipline/Seanome"
    seanomeCommand = "python %s/Seanome.py --threads %s find_csr -a %s -g %s -f %s -o %s  -l %s" 


    # given that we have the pairwiseDirsPath, we don't need to run the seeding 
    # stored in order[0] and order[1]
    if os.path.exists(os.path.join( pairwiseDirsPath, order[0]+"_"+order[1])):
        outputPrefix =  order[0]+"_"+order[1]
        alisDir =  os.path.join(pairwiseDirsPath, outputPrefix)
    else:
        outputPrefix =  order[1]+"_"+order[0]
        alisDir =  os.path.join(pairwiseDirsPath, outputPrefix)

    for sp2 in order[2:]:
        outputPrefix = "%s_%s" % (outputPrefix, sp2) 
        genomeFile = "%s/%s_masked.fasta" % (genomesDir, sp2)
        fastaOut = "fasta/%s" % outputPrefix
        alisOut = "alignments/%s" % outputPrefix
        os.system(seanomeCommand %  (seanomeDir, 18, alisDir, genomeFile, fastaOut, alisOut, 90) )

        total =  getCount(alisOut, value="nbAlis")

        if minScore and total <= minScore:
            print "broke after %s" % alisOut
            break
        alisDir =  alisOut

    return total


def returnBestTree(pairwiseDirsPath,  value = "nbAlis"):
    namesSet = set() 
    for pair in os.listdir(pairwiseDirsPath):
        namesSet.update(pair.split("_"))
    namesList = list(namesSet)

    distMatrix =  np.matrix([[0 for i in range(len(namesList))] for j in range(len(namesList))])

    for i in range(len(namesList)-1):
        for j in range(i+1, len(namesList)):
            if os.path.exists(os.path.join(pairwiseDirsPath, namesList[i]+"_"+namesList[j])):
                distMatrix[i,j]=getCount(os.path.join(pairwiseDirsPath, namesList[i]+"_"+namesList[j]), value)
            else:
                distMatrix[i,j]=getCount(os.path.join(pairwiseDirsPath, namesList[j]+"_"+namesList[i]), value)

    linkaegeMatrix = cluster.hierarchy.average(distMatrix)
    
    print distMatrix
    print linkaegeMatrix
    print namesList
    upperNodes={}
    upperNodesCount = len(namesList) -1
    order = []
    for i in range(linkaegeMatrix.shape[0]):
        if upperNodes.has_key(int(linkaegeMatrix[i][0])):
            sp1 =  str(upperNodes[int(linkaegeMatrix[i][0])])
        else:
            sp1 = namesList[ int(linkaegeMatrix[i][0])]

        if upperNodes.has_key(int(linkaegeMatrix[i][1])):
            sp2 =  str(upperNodes[int(linkaegeMatrix[i][1])])
        else:
            sp2 = namesList[ int(linkaegeMatrix[i][1])]

        upperNodesCount+= 1
        upperNodes[upperNodesCount] = sp1+"_"+sp2
        order.append([sp1,sp2])
        print "iteration %s: %s vs. %s " % (i, sp1, sp2) 
    return order




def explorePermutations(pairwiseDirsPath, genomesDir, value ="nbAlis",  measure= "score", maxTries=20):
    
   order = returnBestPutativePerm(pairwiseDirsPath, value, measure)
   total = testPermutation(pairwiseDirsPath, genomesDir, order,  minScore=None)
   nbTries=1
   for newOrder in itertools.permutations(order):
       newTotal = testPermutation(pairwiseDirsPath, genomesDir, newOrder,  minScore=total)
       if newTotal > total:
           total = newTotal
           order = newOrder
       if nbTries == maxTries:
           break
       nbTries+=1

   print "The best value so far is %s" % total
   print "The best order so far is %s" % order
   





   



#value = "nbAlis"
#pairwiseDirsPath = "."
#print returnBestPutativePerm(pairwiseDirsPath,  value)
