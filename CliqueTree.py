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

    def setNodeList(self,nodeList):
        self.nodeList=nodeList

    def setEdges(self,edges):
        self.edges=edges

    def setFactorList(self,factorList):
        self.factorList=factorList

    def setEvidence(self,evidence):
        self.evidence=evidence

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
            based on https://github.com/indapa/PGM/blob/master/Prog4/EliminateVar.m """

        useFactors = []#the index of the factor that contains the variable Z
        scope = []

        for i in range (len(self.factorList)):
            if Z in self.factorList[i].getVar().tolist():
                useFactors.append(i)#the ith factor is being currently involved in elimination
                scope=list(set.union(set(scope), self.factorList[i].getVar().tolist() ))

        print 'scope: ', scope
        print 'useFactors: ', useFactors
        print 'Z: ', Z


        # update edge map
        # These represent the induced edges for the VE graph.

        for i in range ( len(scope)):
            for j in range ( len(scope)):
                if i!=j:
                    self.edges[i,j]=1
                    self.edges[j,i]=1
        self.edges[Z-1,:]=0
        self.edges[:,Z-1]=0

        #nonUseFactors = setdiff(1:length(F),[useFactors]);
        #list( set.difference( set( range( len(self.factorList)) , set(useFactors) ) )

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
            
        print 'newmap: ', newmap

        newFactor = Factor( [], [], [], 'newFactor')

        for i in range( len (useFactors)):
            newFactor = FactorProduct(newFactor,self.factorList[ useFactors[i] ])
        print newFactor

        newFactor = FactorMarginalization( newFactor,[Z] )
        #newF(length(nonUseFactors)+1) = newFactor;
        newF.append ( newFactor )
        

        #for i in range ( len( unusedFactors) ):
        #    print i
        #    newmap [ unusedFactors[i]-1 ] = i+1
        
        


        
