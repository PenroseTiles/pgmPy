import numpy as np
from CliqueTree import *
from FactorOperations import *

def createCliqueTree( factorList):
    """ return a Clique Tree object given a list of factors
        it peforms VE and returns the clique tree the VE
        ordering defines. See Chapter 9 of Friedman and Koller
        Probabilistic Graphical Models"""

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
    #print 'length of factorList: ', len(factorList)
    #print C.toString()
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


def PruneTree ( C ):
    """ prune a clique tree by determing if neighboring cliques are 
        supersets of each other. E.g.: [A,B,E] -- [A,B] -- [A,D] 
        pruned: [A,B,E] -- [A,D] """

    ctree_edges=C.getEdges()
    (nrows,ncols)=np.shape( ctree_edges )
    totalNodes=nrows
    Cnodes=C.getNodeList()
    #print 'Cnodes: ', Cnodes
    toRemove=[]
    #print range( totalNodes )

    for i in range ( totalNodes ):
        if i in toRemove: continue
        #np.nonzero returns tuple, hence the [0]
        #we collect the neighbors of the ith clique
        neighborsI = np.nonzero ( ctree_edges[i,:] )[0].tolist()
        for c in range ( len(neighborsI) ):
            j= neighborsI[c]
            assert ( i != j), 'i cannot equal j: PruneTree'
            if j in toRemove: continue
            #here is where we look for superset neighboring nodes in the CTree
            if sum (  [ x in Cnodes[j] for x in Cnodes[i] ]  ) == len( Cnodes[i] ):
                for nk in neighborsI:
                    cnodes_i = set ( Cnodes[i] )
                    cnodes_nk= set ( Cnodes[nk] )
                    if len( list ( set.intersection( cnodes_i, cnodes_nk) ) ) == len (Cnodes[i]):
                        neighborsI_set=set( neighborsI )
                        nk_set=set( [nk] )
                        ctree_edges [ list( neighborsI_set - nk_set ), nk ] = 1
                        ctree_edges [  nk, list( neighborsI_set - nk_set )] = 1
                        break
                ctree_edges[i,:]=0
                ctree_edges[:,i]=0
                toRemove.append(i)
    toKeep =  list ( set ( range( totalNodes ) ) - set ( toRemove ) )
    for indx in toRemove:
        Cnodes[indx]=[]
    Cnodes=[ item for item in Cnodes if len(item) > 0 ]
    ctree_edges= ctree_edges[np.ix_(toKeep, toKeep)]

    C.setNodeList( Cnodes )
    C.setEdges( ctree_edges )

    #return the pruned tree with the updated nodes and edges
    return C

def CliqueTreeObserveEvidence ( C, E ):
    """ given a CliqueTree object C and list of values E, which represent evidence, update
        the factors of the cliqueTree C to reflect the observed evidence.
        Note that ObserveEvidence in FactorOperations assumes E is a Nx2 matrix,
        here we build the Nx2 matrix by assuing the jth index of E is the evidence
        for the variable j"""
    factorList= C.getFactorList()
    for j in range ( len (E)):
        if E[j] > 0:
            factorList=ObserveEvidence( factorList, np.array(np.matrix( [ j+1, E[j]] ) ) )
    C.setFactorList(factorList)
    #return the new CliqueTree object with the updated evidence
    return C


