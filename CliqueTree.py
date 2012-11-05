from Factor import *
from FactorOperations import *
import numpy as np 
class CliqueTree(object):
    'represent a Clique tree'

    def __init__(self, nodeList=[], edges=[], factorList=[],evidence=[]):
        self.nodeList=nodeList
        self.edges=edges
        self.factorList=factorList
        self.evidence=evidence
        self.card= []
        self.factorInds= len( factorList ) * [None]


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
        self.factorInds= len( factorList ) * [None]

    def setEvidence(self,evidence):
        self.evidence=evidence


    def setCard(self, cardinality):
        self.card = cardinality

    def getNodeList(self):
        return self.nodeList


    def getEdges(self):
        return self.edges

    def getFactorList(self):
        return self.factorList

    def getEvidence(self):
        return self.evidence

    def eliminateVar(self,Z):
        """ a variable elimination function
            based on https://github.com/indapa/PGM/blob/master/Prog4/EliminateVar.m

            Z is the variable to be eliminated. We base this code on the matlab file
            linked to above as well as the Sum-product VE pseudo code in Koller and Friedman
            page 298

            """

        useFactors = []#the index of the factor that contains the variable Z
        scope = []

        #get a list containining the index in self.factorLlist of factors
        #that contain the variable Z to be eliminated
        # get the scope of variables from the factors that contain variable Z
        for i in range (len(self.factorList)):
            if Z in self.factorList[i].getVar().tolist():
                useFactors.append(i)#the ith factor is being currently involved in elimination
                scope=list(set.union(set(scope), self.factorList[i].getVar().tolist() ))

        print 'scope: ', scope
        print 'useFactors: ', useFactors
        print 'Z: ', Z


        # update edge map
        # These represent the induced edges for the VE graph.
        #once the variable Z is eliminated, its edges are removed from the graph
        # but in the process of elimination, we create a new factor. This
        #inroduces fill edges pg. 307 Koller and Friedman
        for i in range ( len(scope)):
            for j in range ( len(scope)):
                if i!=j:
                    self.edges[i,j]=1
                    self.edges[j,i]=1
        self.edges[Z-1,:]=0
        self.edges[:,Z-1]=0

        

        #these are teh indices of factorList which are not involved in VE
        unusedFactors= list( set.difference ( set(range(len(self.factorList))), set(useFactors)    )   )
        print 'unusedFactors:', unusedFactors
        print 'useFactors ', useFactors
        print 'length of unused factors: ', len(unusedFactors)

        newF=len(unusedFactors)*[None]
        newmap=np.zeros(max(unusedFactors)+1,dtype=int).tolist()
        print 'newmap initially: ', newmap
        for i in range( len(unusedFactors)):
            newF[i]=self.factorList[ unusedFactors[i] ]
            #print unusedFactors[i]
            newmap[ unusedFactors[i] ]= i
            
        print 'newmap ', newmap,"\n"
        print 'length of newmap: ', len(newmap), "\n"

        newFactor = Factor( [], [], [], 'newFactor')

        for i in range( len (useFactors)):
            newFactor = FactorProduct(newFactor,self.factorList[ useFactors[i] ])
        print 'newFactor:\n ', newFactor
        #for z in range ( len(newF) ):
        #    print z, newF[z].getVar()


        newFactor = FactorMarginalization( newFactor,[Z] )
        #newF(length(nonUseFactors)+1) = newFactor;
        newF.append ( newFactor )
        #print 'length of newF after FactorMarginaliztion: ', len(newF)
        #for z in range ( len(newF) ):
        #    print z, newF[z].getVar()
        #update the factorList, after eliminating the variable and generating the newFactor
        self.factorList=newF

        self.nodeList.append ( scope ) #add to the node list of the clique tree
                                       #the factors  that contained the variable to be eliminated
        print 'self.nodeList ', self.nodeList
        #newC is the total number of nodes in the clique tree
        newC=len( self.nodeList )
        #print 'newC: ', newC

        #print self.edges

        #print 'factorInds: ',self.factorInds
        #self.factorInds[ newC -1] = len(unusedFactors) + 1
        self.factorInds.append ( len(unusedFactors) + 1  )


        print "last for loop ...\n"
        for i in range( newC -1 ):
            print 'i: ', i
            if self.factorInds [ i-1 ] in useFactors:
                self.edges[ i-1, newC-1 ] = 1
                self.edges [ newC-1, i-1 ] = 1
                self.factorInds[ i-1 ] = 0
            else:
                if self.factorInds [i] != None:
                    print 'factorInds[i]: ', self.factorInds[i]
                    self.factorInds[ i ]  = newmap [ self.factorInds[i] -1 ]
                    


        
        
        
        #print scope
        
        #self.nodeList.append( scope )
        #print self.nodeList
        #for i in range ( len( unusedFactors) ):
        #    print i
        #    newmap [ unusedFactors[i]-1 ] = i+1
        
        


        
