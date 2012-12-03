import numpy as np
from CliqueTree import *
from FactorOperations import *
import pprint
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


def CliqueTreeInitialPotential( C ):
    """ given a tree, calculate the initial potentials for each of the cliques
        the factors in the updated clique list are FActor objects"""

    N= C.getNodeCount()
    totalFactorCount=C.getFactorCount()

    nodeList=C.getNodeList()
    factorList=C.getFactorList()

    cliqueList=[ Factor( [], [], [], str(i) )  for i in range(N)  ]
    #edges=np.zeros( (N,N) )

    """ First assign the factors to appropriate cliques
    based on the skeleton cliqueTree cTree"""

    factorsUsed=np.zeros( totalFactorCount, dtype=int).tolist()

    for i in range(N):
        cliqueList[i].setVar( nodeList[i] )
        F=[]
        
        for j in range( len(factorList) ):
            if len( factorList[j].getVar().tolist() ) == len ( list( set.intersection ( set(cliqueList[i].getVar().tolist() ), set( factorList[j].getVar().tolist() ) ) ) ):
        
                if factorsUsed[j] == 0:
                    F.append( factorList[j] )
                    factorsUsed[j] = 1
    #print F
        F= [ f.getFactor() for f in F ]
        cliqueList[i]=ComputeJointDistribution ( F )

    C.setNodeList(cliqueList)

    return C

def getNextClique(P, messages):

    """ we need to come up wih a proper message passing order. A clique is ready to pass
        messages upward once its recieved all downstream messages from its neighbor (and vice versa)
        its ready to transmit downstream once it recieves all its upstream messages

        the ith clique C_i is ready to transmit to its neighbor C_j when C_i recieves all its
        messages from neigbors except C_j. In cTree message passing, each message is passed
        once.  To get the process started we start with our initial potential cTree, P
        and an empty matrix of factors, representing messages passed between the nodes on the clique
        tree """
    i=j=-1
    edges=P.getEdges()
    #print edges
    (nrow, ncol) = np.shape(edges)

    for r in range(nrow):
        
        #we want to ignore nodes with only one neighbor
        #becuae they are ready to pass messages
        if np.sum(edges[r,:] ) == 1:
            continue 

        foundmatch=0

        for c in range(ncol):
            if  edges[r,c] == 1 and messages[r,c].getVarCount()  == 0:
                #list of indices indicating neighbors or r
                #print 'r,c: ', r, ' ', c
                #print 'edges[r,c]: ', edges[r,c]
                Nbs=np.nonzero(edges[:,r])[0]
                #print 'Nbs before:', Nbs
                Nbs=Nbs[np.nonzero(Nbs!= c)[0]]
                #print 'Nbs after: ', Nbs
                allnbmp=1 #neighbors messages passed?
                
                #find all of r's neighbors have sent messages *to* r
                for z in range( len(Nbs) ):
                    #print messages[Nbs[z],r].getVarCount()
                    if messages[ Nbs[z],r].getVarCount()  == 0:
                        allnbmp=0

                if allnbmp == 1:
                    foundmatch=1
                    break
        print
        if foundmatch==1:
            #sys.stderr.write("found match!\n")
            i=r
            j=c
            break
        print
    return (i,j)


def CliqueTreeCalibrate( P, isMax=False):
    """ this function performs sum-product or max-product algorithm for clique tree calibration.
        P is the CliqueTree object. isMax is a boolean flag that when set to True performs Max-Product
        instead of the default Sum-Product. The function returns a calibrated clique tree in which the
        values of the factors is set to final calibrated potentials. """

    if isMax == True:
        pass

    ctree_edges=P.getEdges()
    ctree_cliqueList=P.getNodeList()
    N=P.getNodeCount() #Ni is the total number of nodes (cliques) in cTree

    #dummyFactor=Factor( [], [], [], 'factor')
    #set up messsages to be passed
    #MESSAGES[i,j] represents the message going from clique i to clique j
    #MESSAGES will be a matrix of Factor objects
    MESSAGES=np.tile( Factor( [], [], [], 'factor'), (N,N))



    """While there are ready cliques to pass messages between, keep passing
    messages. Use GetNextCliques to find cliques to pass messages between.
    Once you have clique i that is ready to send message to clique
    j, compute the message and put it in MESSAGES(i,j).
    Remember that you only need an upward pass and a downward pass."""

    """ leaf nodes are ready to pass messages right away
        so we initialize MESSAGES with leaf message factors
        recall, a node is a leave if row sum is equal to 1"""
    for row in range(N):
        rowsum= np.sum( ctree_edges[row,:] )
        if rowsum ==1 :
            #Returns a tuple of arrays, one for each dimension, we want the first, hence the [0]
            leafnode=np.nonzero( ctree_edges[row,:] )[0].tolist()[0]
            #I discovered NumPy set operations http://docs.scipy.org/doc/numpy/reference/routines.set.html
            marginalize=np.setdiff1d( ctree_cliqueList[row].getVar(),  ctree_cliqueList[leafnode].getVar() ).tolist()
            sepset=np.intersect1d( ctree_cliqueList[row].getVar(), ctree_cliqueList[leafnode].getVar() ).tolist()

            """ if isMax, this is sumproduct, so we do factor marginalization """
            if isMax == 0:
                #MESSAGES(row,leafnode)=FactorMarginalization(P.cliqueList(row),marginalize);
                MESSAGES[row,leafnode]=FactorMarginalization(ctree_cliqueList[row], marginalize )
                if np.sum( MESSAGES[row,leafnode].getVal() ) != 1:
                    newVal=MESSAGES[row,leafnode].getVal() / np.sum( MESSAGES[row,leafnode].getVal() )
                    MESSAGES[row,leafnode].setVal(newVal)
            else:
                pass

    """ now that the leaf messages are initialized, we begin with the rest of the clique tree
    now we do a single pass to arrive at the calibrated clique tree. We depend on
    GetNextCliques to figure out which nodes i,j pass messages to each other"""

    while True:
        (i,j)=getNextClique(P,MESSAGES)
        if sum ( [ i, j] ) == 0:
            break

        """ similiar to above, we figure out the sepset and what variables to marginalize out
        between the two cliques"""
        marginalize=np.setdiff1d( ctree_cliqueList[i].getVar(),  ctree_cliqueList[j].getVar() ).tolist()
        sepset=np.intersect1d( ctree_cliqueList[i].getVar(), ctree_cliqueList[j].getVar() ).tolist()

        """ find all the incoming neighbors, except j """
        Nbs=np.nonzero( ctree_edges[:,i])
        Nbs_minusj=[ elem for elem in Nbs if elem !=  j ]
        #see numpy for matlab users http://www.scipy.org/NumPy_for_Matlab_Users
        # these are incoming messages to the ith clique
        Nbsfactors=MESSAGES[np.ix_(Nbs_minusj, i)]

        """ this is sum/product """
        if isMax == 0:
            if len(Nbsfactors) == 1
                Nbsproduct=FactorProduct( Nbsfactors[0], IdentityFactor(Nbsfactors[0]) )
            else:
                Nbsproduct=ComputeJointDistribution( Nbsfactors )

            #now mulitply wiht the clique factor
            CliqueNbsproduct=FactorProduct( Nbsproduct, ctree_cliqueList[i] )
            CliqueMarginal= FactorMarginalization ( CliqueNbsproduct, marginalize )
            #normalize the marginal
            newVal=CliqueMarginal.getVal() / np.sum( CliqueMarginal.getVal().getVal() )
            CliqueMarginal.setVal( newVal )
            MESSAGES[i,j] = CliqueMarginal
        else:
            pass