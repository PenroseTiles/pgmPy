import sys
import numpy as np
from Factor import *
from PGMcommon import *
from FactorOperations import *
import itertools

class PhenotypeFactor (object):
    """ represents a factor that encodes Pr(phenotype|genotype)
        for purposes here the variable to the left of the conditioning
        bar is the first variable in the PhenotypeFactor var list. """

    def __init__(self, isDominant, genotypeVar, phenotypeVar, name):

        #instantiate a Factor object
        phenotype = Factor( [phenotypeVar, genotypeVar], [2, 3], [], name )
    
        phenotype.setVal( np.zeros ( np.prod(phenotype.getCard())).tolist() )
        #this enumerates the values the factor can take
        # since there are 2x3 cardinality, 6 possible assignments
        assignments=IndexToAssignment( np.arange(np.prod(phenotype.getCard())), phenotype.getCard() )
        val=val = np.zeros(np.prod(phenotype.getCard() ))
        (nrows,ncols)=np.shape(assignments)

        for i in range(np.prod([2,3])):
         #if its dominant, if you have at least one copy, you have the phenotype
            (pheno,geno)=assignments[i]
            if isDominant==1:
                if pheno ==1: #affected
                    if geno ==1 or geno ==2:
                        val[i]=1
                    else:
                        val[i]=0
                else:#uneffected
                    if geno == 3:
                        val[i]=1


            if isDominant == 0:
                if pheno == 1:
                    if geno==3:
                        val[i]=1
                else:
                    if geno ==1 or geno == 2:
                        val[i]=1


        phenotype.setVal( val.tolist() )

        self.phenotype=phenotype

    def __str__(self):
        return self.phenotype.__str__()


class GenotypeAlleleFreqFactor (object):
    """ construct a factor that has the probability of each genotype
        given allele frequencies Pr(genotype|allele_freq)"""

    def __init__(self, allelefreqs, genotypeVar, name):
        self.allelefreq=allelefreqs
        #number of alleles == number of allele frequencies passed in
        numAlleles=len(allelefreqs)
        self.allelesToGenotypes=None
        self.genotypesToAlleles=None
        self.genotypeFactor=None

        #map alleles to genotypes and genotyeps to alleles
        (self.allelesToGenotypes, self.genotypesToAlleles)=generateAlleleGenotypeMappers(numAlleles)
        (ngenos,ploidy)=np.shape(self.genotypesToAlleles)


        self.genotypeFactor = Factor( [genotypeVar], [], [], name)
        #the cardinality of the factor is the number of genotypes
        self.genotypeFactor.setCard( [ngenos] )

        #set the values to zero initially
        values=np.zeros( (np.prod(self.genotypeFactor.getCard()))).tolist()
        
        for i in range (ngenos):
            alleles=self.genotypesToAlleles[i,:].tolist()
            

            if alleles[0] == alleles[1]:
                values[i]= np.prod( [ allelefreqs[j] for j in alleles ])
                
            else:
               values[i]= np.prod( [ allelefreqs[j] for j in alleles ]) * 2
        
        self.genotypeFactor.setVal( values )

    def __str__(self):
        return self.genotypeFactor.__str__()