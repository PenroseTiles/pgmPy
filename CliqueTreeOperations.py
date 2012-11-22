import numpy as np
from CliqueTree import *

def createCliqueTree( factorList):
    """ return a Clique Tree object given a list of factors """

    V=getUniqueVar(factorList)
    totalVars=len(V)
    cardinality=np.zeros(len(V)).tolist()
    for i in range(len(V)):
        for j in range(len(factorList)):
            try:
                indx= factorList[j].getVar().tolist().index( V[i] )
                cardinality[i]=factorList[j].getCard().tolist()[indx]
                break
            except:
                continue

    edges=np.zeros( (totalVars, totalVars))

    for f in factorList:
        variableList=f.getVar()
        for j in range(len(variableList) ):
            for k in range (len(variableList) ):
                edges[ variableList[j]-1, variableList[k]-1 ]=1

    (nrows,nedges)=np.shape(edges)            

    C=CliqueTree()
    C.setCard( cardinality )
    C.setEdges(np.zeros( (totalVars, totalVars)))
    C.setFactorList(factorList)
    print 'length of factorList: ', len(factorList)
    print C.toString()
    cliquesConsidered = 0

    while cliquesConsidered < len(V):
        bestClique = 0
        bestScore = sys.maxint
        for i in range(nrows):
            score=np.sum( edges[i,:] )
            if score > 0 and score < bestScore:
                bestScore = score
                bestClique = i+1
        cliquesConsidered+=1
    
        edges, factorList=C.eliminateVar(bestClique, edges, factorList)

    return C
# <codecell>


