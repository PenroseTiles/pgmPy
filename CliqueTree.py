from Factor import *
from FactorOperations import *
import numpy as np
import networkx as nx
import pdb

class CliqueTree(object):
    'represent a Clique tree'

    def __init__(self, nodeList=[], edges=[], factorList=[],evidence=[]):
        self.nodeList=nodeList
        self.edges=edges
        self.factorList=factorList
        self.evidence=evidence
        self.card= []
        self.factorInds= []


    def toString(self):
        
        print 'nodes: ', self.nodeList
        print 'card: ', self.card
        print 'factorList: ', len( self.factorList)
        print 'factorInds:', self.factorInds
        print 'edges:\n',self.edges
        
        
    def setNodeList(self,nodeList):
        self.nodeList=nodeList

    def setEdges(self,edges):
        self.edges=edges

    def setFactorList(self,factorList):
        self.factorList=factorList
        #self.factorInds= len( factorList ) * [None]

    def setEvidence(self,evidence):
        self.evidence=evidence


    def setCard(self, cardinality):
        self.card = cardinality

    def getNodeList(self):
        return self.nodeList

    def getNodeCount(self):
        return len(self.nodeList)


    def getEdges(self):
        return self.edges

    def getFactorList(self):
        return self.factorList

    def getFactorCount(self):
        return len(self.factorList)
    

    def getEvidence(self):
        return self.evidence

    def incorporateEvidence(self):
        for j in range ( len(self.evidence)):
            k=j+1
            if self.evidence[j] > 0:
                self.factorList=ObserveEvidence(self.factorList, np.matrix([[k, self.evidence[j] ]] ) )
    

    def eliminateVar(self,Z,E,factorList):
        """ a variable elimination function
            based on https://github.com/indapa/PGM/blob/master/Prog4/EliminateVar.m

            Z is the variable to be eliminated. We base this code on the matlab file
            linked to above as well as the Sum-product VE pseudo code in Koller and Friedman
            page 298

            E is a numpy 2d matrix representing adjacency matrix of variables
            It represents the induced VE graph
            Once a variable is eliminated, its edges are removed from E

            """

        useFactors = []#the index of the factor that contains the variable Z
        scope = []

        #print 'Z: ', Z

        
        #get a list containining the index in self.factorLlist of factors
        #that contain the variable Z to be eliminated
        # get the scope of variables from the factors that contain variable Z
        for i in range (len(factorList)):
            
            if Z in factorList[i].getVar().tolist():
                useFactors.append(i)#the ith factor is being currently involved in elimination
                scope=list(set.union(set(scope), factorList[i].getVar().tolist() ))

        
        # update edge map
        """ These represent the induced edges for the VE graph.
         once the variable Z is eliminated, its edges are removed from the graph
         but in the process of elimination, we create a new factor. This
         introduces fill edges (see pg. 307 Koller and Friedman)
         Z is one based, but the indices in E are zero based, hence Z-1
         also the variable names in scope are 1 based, so we subtract 1 when
         updating the induced VE graph """

        for i in range ( len(scope)):
            for j in range ( len(scope)):
                if i != j:
                    E[  scope[i]-1, scope[j]-1   ]=1
                    E[ scope[j]-1,  scope[i]-1   ]=1
        E[Z-1,:]=0
        E[:,Z-1]=0

        #G=nx.from_numpy_matrix(E)
        #print 'induced graph edges:\n', (G.edges())
        #nx.draw_shell(G)
        #plt.show()

        
        #these are the indices of factorList which are not involved in VE
        unusedFactors= list( set.difference ( set(range(len(factorList))), set(useFactors)    )   )
        
        newF=None
        #check first if there are any unused factors left!
        if len(unusedFactors) > 0:
            newF=len(unusedFactors)*[None]
            newmap=np.zeros(max(unusedFactors)+1,dtype=int).tolist()
            
            #newF is a new factor list, we populate it first
            #with the unused factors
            #newmap is maps the new location of ith unusedFactor
            for i in range( len(unusedFactors)):
                newF[i]=factorList[ unusedFactors[i] ]
                newmap[ unusedFactors[i] ]= i
            
        #print 'newmap ', newmap,"\n"
        #print 'length of newmap: ', len(newmap), "\n"

        newFactor = Factor( [], [], [], 'newFactor')

        #we multiple in all the factors that contain the variable Z
        for i in range( len (useFactors)):
            newFactor = FactorProduct(newFactor,factorList[ useFactors[i] ])
        

        #then we marginalize Z out and obtain a new factor
        #then append it the end of newF, the new factor list
        newFactor = FactorMarginalization( newFactor,[Z] )
        #print 'newFactor: ',newFactor
        #newF(length(nonUseFactors)+1) = newFactor;
        if newFactor != None:
            newF.append ( newFactor )

        

        if newF != None:
            factorList=newF
        #return E

        ########################################################################
        """ the remaining code builds the edges of the clique tree """

        """ add new node with the factors that contain the variable Z
            adding a  new node represents new clique.
            The scope of every factor generated during the variable elimination process is a clique pg. 309 Koller & Friedman """

        self.nodeList.append ( scope )
        
        #newC is the total number of nodes in the clique tree
        newC=len( self.nodeList )
        #print 'newC: ', newC

        #factorInds are individual factors with one variable ... I think
        self.factorInds.append ( len(unusedFactors) + 1  )
        

        #print 'range( newC -1) ', range( newC-1  )
        #print 'factorInds: ', self.factorInds
        #print 'useFactors: ', useFactors
        #pdb.set_trace()
        """ we update the edges of the clique tree """
        for i in range( newC -1 ):
            
            #if self.factorInds [ i ] -1 in useFactors:
            #there was the off by onoe erorr - the values in factorInds
            #were one-based, need to subtract 1
            if self.factorInds [ i ] -1  in useFactors:
                
                self.edges[ i, newC-1 ] = 1
                self.edges [ newC-1, i ] = 1
                self.factorInds[ i ]  = 0
            else:
                if self.factorInds [i] != 0:
                    #print 'i: ', i
                    #print 'factorInds: ', self.factorInds
                    #print 'newmap: ', newmap
                    #print 'newmap [ self.factorInds[i] -1: ', newmap [ self.factorInds[i] -1 ]
                    #print 'self.factorInds[ i ]  = newmap [ self.factorInds[i] - 1  ] + 1 '
                    if len(unusedFactors) > 0:
                        #self.factorInds[ i ]  = newmap [ self.factorInds[i] -1  ] +1
                        self.factorInds[ i ]  = newmap [ self.factorInds[i] -1  ] +1
                        #self.factorInds[ i ]  = newmap [ self.factorInds[i]   ]
                    
        #print 'factorInds right before returning: ', self.factorInds
        return E, factorList

        
        
        
        
        
        


        
