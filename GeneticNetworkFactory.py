from Factor import *
from FactorOperations import *
from PedigreeFactors import *
import itertools
import numpy as np
"""" Still not sure how this is going to work
     This class is a factory for generating a genetic network
     If we consider each location in the genome independent
     we generate a new network for each position along an interval.
     A network is a collection of factors, so we return a python
     list of Pedigree factors (either GenotypeAlleleFreqFactor for founders
     or GenotypeGivenParentsFactor for non-founders. Each phenotype is conditionally
     independent given its genotoype, so each member of the pedigree has a
     PhenotypeGivenGenotypeFactor"""

class GeneticNetworkFactory(object):
    

    def __init__(self, pedfile, alphaList, allelefreq, chrom, position):
        #parse pedfile
        self.alphaList=alphaList
        self.allelefreq=allelefreq
        self.totalAlleles=len(allelefreq)
        self.chrom=chrom
        self.pos=position

        self.pedigree=Pedfile(pedfile)
        self.pedigree.parsePedfile()
        self.pedlist=self.pedigree.getPedList()


        #list of factors that will comprise the Genetic network
        self.totalFactors=self.pedigree.getTotalSize() * 2
        self.factorList=self.totalFactors*[None]

        
    def constructNetwork(self):
        totalPeople=self.pedigree.getTotalSize()
        for i in range( totalPeople ):
    
            if self.pedlist[i].isFounder():
                #print self.pedlist[i].getid()
                self.factorList[i]=GenotypeAlleleFreqFactor(self.allelefreq,self.pedlist[i].getid(),self.pedlist[i].getid())
                #factorList(i)=genotypeGivenAlleleFreqsFactor(alleleFreqs,i);
            else:
                #3print self.pedlist[i].getParents(), self.pedlist[i].getid()
                #GenotypeGivenParentsFactor(2,"bart","homer","marge","""Bart | Homer, Marge """)
                self.factorList[i]=GenotypeGivenParentsFactor(self.totalAlleles, self.pedlist[i].getid(), self.pedlist[i].getParents()[0], self.pedlist[i].getParents()[1], "child|Father,Child")
        
            self.factorList[i+totalPeople]=PhenotypeGivenGenotypeFactor(self.alphaList, i, i+totalPeople, "|".join([ str(i), str(i+totalPeople)]) )

    def getFactorList(self):
        return self.factorList


