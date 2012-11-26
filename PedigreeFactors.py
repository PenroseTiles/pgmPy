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



    def getVar(self):
        return self.phenotype.getVar()
    def getCard(self):
        return self.phenotype.getCard()
    def getVal(self):
        return self.phenotype.getVal()

    def getFactor(self):
        return self.phenotype

    def __str__(self):
        return self.phenotype.__str__()



class PhenotypeGivenGenotypeFactor(object):
    """ construct factor of phenotype|genotype
        #prob of being effected, given the ith genotype
        #alphaList[i] is the prob of being effected given the ith genotype """
    def __init__(self,alphaList, phenotypeVar, genotypeVar , name):
        self.phenotypeFactor=Factor( [ phenotypeVar, genotypeVar], [], [], name)
        self.alpha=np.array ( alphaList)

        ngenotypes=len(alphaList)
        self.phenotypeFactor.setCard( [2, ngenotypes])

        values=[x for x in range( np.prod(self.phenotypeFactor.getCard()))]

        for i in range( len(alphaList )):

            values[i]=alphaList[i]
            values[i+1]=1-alphaList[i]
        ctr=0
        alphas=2*len(alphaList)*[None]
        for i in range(len(alphaList)):
            alphas[ctr]=alphaList[i];
            ctr=ctr+1
            alphas[ctr]=1-alphaList[i];
            ctr=ctr+1

        values=alphas
        self.phenotypeFactor.setVal( values)

    def getVar(self):
        return self.phenotypeFactor.getVar()
    def getCard(self):
        return self.phenotypeFactor.getCard()
    def getVal(self):
        return self.phenotypeFactor.getVal()
    def setVal(self,val):
        self.phenotypeFactor.setVal(val)
    def getFactor(self):
        return self.phenotypeFactor


    def __str__(self):
        return self.phenotypeFactor.__str__()


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


    def getVar(self):
        return self.genotypeFactor.getVar()
    def getCard(self):
        return self.genotypeFactor.getCard()
    def getVal(self):
        return self.genotypeFactor.getVal()
    def setVal(self,val):
        self.genotypeFactor.setVal(val)
    def getFactor(self):
        return self.genotypeFactor


    def __str__(self):
        return self.genotypeFactor.__str__()

class GenotypeGivenParentsFactor (object):
    """ construct factor that has prob of genotype of child given both parents
        Pr(g_child| g_mother, g_father """

    def __init__(self,numAlleles, genotypeVarChild, genotypeVarParentOne, genotypeVarParentTwo, name):
        self.genotypeFactor =  Factor( [genotypeVarChild, genotypeVarParentOne, genotypeVarParentTwo ], [ ], [ ], name)

        #map alleles to genotypes and genotyeps to alleles
        (self.allelesToGenotypes, self.genotypesToAlleles)=generateAlleleGenotypeMappers(numAlleles)

        (ngenos,ploidy)=np.shape(self.genotypesToAlleles)
        

        
        self.genotypeFactor.setCard([ ngenos,ngenos,ngenos ] )
        #set the values to zero initially
        values=np.zeros( (np.prod(self.genotypeFactor.getCard()))).tolist()

        #iterate thru variable assignments to random variables
        #assign probablities based on Punnet square crosses
        assignments=IndexToAssignment( np.arange(np.prod(self.genotypeFactor.getCard())), self.genotypeFactor.getCard() )-1
        for z in range( np.prod(self.genotypeFactor.getCard() ) ):
            curr_assign= assignments[z]
            childAssignment=int(curr_assign[0])

            parent1gametes= self.genotypesToAlleles[curr_assign[1],:]
            parent2gametes= self.genotypesToAlleles[curr_assign[2],:]
            #print 'parental gametes: ', parent1gametes, parent2gametes
            #print 'child assignment: ', childAssignment
            #list of tuples containing list of zygote(genotype) tuples
            zygote_list=list(itertools.product(parent1gametes,parent2gametes))
            punnet_freq=[  self.allelesToGenotypes[zygote[0],zygote[1]] for zygote in zygote_list ]
            histc={}
            hist=[]
            for g in range( ngenos):
                histc[g]=0.
            for x in punnet_freq:
                histc[x]+=1.
            #print histc.values()
            for g in range (ngenos):
                hist.append ( histc[g] )
            #print punnet_freq
            hist=(np.array ( hist)) /4
            #print 'hist:', hist
            #print zygote_list
            values[z]=hist[childAssignment]

        self.genotypeFactor.setVal( values )

    def getVar(self):
        return self.genotypeFactor.getVar()
    def getCard(self):
        return self.genotypeFactor.getCard()
    def getVal(self):
        return self.genotypeFactor.getVal()
    def setVal(self, val):
        self.genotypeFactor.setVal(val)

    def getFactor(self):
        return self.genotypeFactor
    def genotypeSlice(self):
        pass
        #see this http://stackoverflow.com/q/4257394/1735942

    def __str__(self):
        return self.genotypeFactor.__str__()


class ChildCopyGivenParentalsFactor(object):
    """ this represents a de-coupled factor
        given a parents two haplotypes, returns
        factor whose values are the probablity
        of inheriting (grand)paternal or (grand)maternal
        haplotype. This allows for some more flexibility
        in modeling inheritance, rather than clumping
        a single parent's haplotype into a genotype
        i.e. GenotypeGivenParentsFactor """

    def __init__(self, numAlleles, geneCopyVarChild, geneCopyHapOne, geneCopyHapTwo):
        self.numalleles=numAlleles
        self.hapone=geneCopyVarChild
        self.haptwo=geneCopyHapTwo

        #geneCopyFactor = struct('var', [], 'card', [], 'val', []);
        self.geneCopyFactor=Factor( [geneCopyVarChild, geneCopyHapOne, geneCopyHapTwo ], [], [], 'child|hap1,hap2')
        self.geneCopyFactor.setCard( [self.numalleles,self.numalleles,self.numalleles ])
        values=np.zeros( np.prod([ self.numalleles,self.numalleles,self.numalleles])).tolist()
        #this keeps track of what posiiton you are in the values list
        index=0
        #the number of iterations thru the nested for loops should be equal to numallels^3

        for i in range(numAlleles):
        #iterate through alleles from
        #grand(paternal) haplotype
            for j in range(numAlleles):
            #iterate through alleles from
            #grand(maternal) haplotype
                for k in range(numAlleles):
                #iterate thru child alleles
                    print i, j, k
                    if j==k:#child has grandmotherhap
                        if i==k:#grandfatherhap is the same
                            values[index]=1
                        else:
                            values[index]=.5
                    elif i==k:#child has grandfather hap
                        values[index]=.5
                    else:
                        pass
                    index+=1
        #print values
        self.geneCopyFactor.setVal( values )

    def getVar(self):
        return self.geneCopyFactor.getVar()
    def getCard(self):
        return self.geneCopyFactor.getCard()
    def getVal(self):
        return self.geneCopyFactor.getVal()
    def getFactor(self):
        return self.geneCopyFactor
    def __str__(self):
        return self.geneCopyFactor.__str__()
  
class ChildCopyGivenFreqFactor(object):
    """ for a founder, its particular haplotype is proprortional to the
    given allelel freq of the locus. This factor is part of the decoupled
    Bayesian Genetic network , along with ChildCopyGivenParentalsFactor"""
    
    def __init__(self, alleleFreqs, geneCopyVar):
        numAlleles = len(alleleFreqs)
        self.geneCopyFactor=Factor( [geneCopyVar], [], [], 'founderHap')
        self.geneCopyFactor.setCard ( [numAlleles])
        self.geneCopyFactor.setVal( alleleFreqs )
        #geneCopyFactor = struct('var', [], 'card', [], 'val', [])
        #geneCopyFactor.var(1) = geneCopyVar;
        #geneCopyFactor.card(1) = numAlleles;
        #geneCopyFactor.val = alleleFreqs';


    def getVar(self):
        return self.geneCopyFactor.getVar()
    def getCard(self):
        return self.geneCopyFactor.getCard()
    def getVal(self):
        return self.geneCopyFactor.getVal()
    def getFactor(self):
        return self.genCopyFactor
    def __str__(self):
        return self.geneCopyFactor.__str__()

class phenotypeGivenHaplotypesFactor(object):
    """ factor represents Pr(phenotype| paternal haplotype, maternal haplotype)
    very similiar to PhenotypeGivenGenotypeFactor, but we are de-coupling into
    paternal and maternal alleles rather than genotype"""

    def __init__(self, alphaList, numAlleles, geneCopyVarOne, geneCopyVarTwo, phenotypeVar):
        
        self.numalleles=numAlleles
        self.alphaList=alphaList
        self.phenotypeFactor=Factor([phenotypeVar,geneCopyVarOne, geneCopyVarTwo], [], [], 'phenotype| geneCopy1, geneCopy2')

        ngenos=len(alphaList)
        self.phenotypeFactor.setCard( [ 2, numAlleles, numAlleles])
        #phenotypeFactor.val = zeros(1, prod(phenotypeFactor.card));
        values=np.zeros( (1, np.prod(self.phenotypeFactor.getCard()))).flatten().tolist()

        affectedAlphas=alphaList
        unaffectedAlphas=[ 1- alpha for alpha in alphaList]


        (allelesToGenotypes, genotypesToAlleles) = generateAlleleGenotypeMappers(numAlleles)
        assignments=IndexToAssignment( np.arange(np.prod(self.phenotypeFactor.getCard())), self.phenotypeFactor.getCard() )-1
        for z in range( np.prod(self.phenotypeFactor.getCard() ) ):
            curr_assign= assignments[z]
            curr_assign=assignments[z]
            genotype_num=allelesToGenotypes[curr_assign[1], curr_assign[2]]
            if curr_assign[0] == 0:
                values[z] = affectedAlphas[genotype_num]
            else:
                values[z] = unaffectedAlphas[genotype_num]
        self.phenotypeFactor.setVal( values )


            #genotype_num=allelesToGenotypes(assignment(2), assignment(3));


    def getVar(self):
        return self.phenotypeFactor.getVar()
    def getCard(self):
        return self.phenotypeFactor.getCard()
    def getVal(self):
        return self.phenotypeFactor.getVal()
    def getFactor(self):
        return self.phenotypeFactor
    def __str__(self):
        return self.phenotypeFactor.__str__()

    def __str__(self):
        return self.phenotypeFactor.__str__()


##########################

class Ped(object):
    """ represents a pedigree record in a Ped file http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped """

    def __init__(self,famid='.', indv='.', paternal='.', maternal='.', sex='.', phenotype='.'):
        """ attributes of a Ped object """
        self.famid=famid
        self.individ=indv
        self.pid=paternal
        self.mid=maternal
        self.sex=sex
        self.pheno=phenotype


    def setfamid(self,famid):
        self.famid=famid

    def setindvid(self,indv):
        self.individ=indv

    def setmid(self,mid):
         self.mid=mid

    def setpid(self,pid):
         self.pid=pid

    def setsex(self,sex):
         self.sex=sex

    def setpheno(self,pheno):
         self.pheno=pheno



    def getfamid(self):
        return self.famid

    def getid(self):
        return self.individ

    def getmid(self):
        return self.mid

    def getpid(self):
        return self.pid

    def getsex(self):
        return self.sex

    def getpheno(self):
        return self.pheno


    def isFounder(self):
        return self.pid == '0' and self.mid == '0'

    def getParents(self):

        return ( self.pid, self.mid)


    def __str__(self):
        return "\t".join( [ self.famid, self.individ, self.pid, self.mid, self.sex, self.pheno] )


class Pedfile(object):
    """ a Pedfile object has a list of Ped  objects """

    def __init__(self,filename):
        self.filename=filename
        self.fh=open(self.filename, 'r')
        self.pedlist=[]

    def parsePedfile(self):
        """ given a filehandle to a *.ped file read its contents and populate the list pedlist with Ped objects """
        for line in self.fh:
            fields=line.strip().split('\t')
            (famid, indv, pid, mid, sex, pheno)=fields[0:6]
            self.pedlist.append( Ped(famid, indv, pid, mid, sex, pheno) )

    
    def returnFounders(self):
        """ return the founders in a ped file (those with unknown paternal and maternids """
        founders=[]

        for pedobj in self.pedlist:
            if pedobj.getpid() == "0":
                founders.append(pedobj)

        return founders

    def returnFounderIds(self):
        """ return the indiv ids of the founders in the ped file"""
        founderids=[]
        for pedobj in self.pedlist:
            if pedobj.getpid() == "0":
                founderids.append( pedobj.getid() )
        return founderids

    def returnNonFounderIds(self):
        """ return the indiv ids of the founders in the ped file"""
        nonfounderids=[]
        for pedobj in self.pedlist:
            if pedobj.getpid() != "0":
                nonfounderids.append( pedobj.getid() )
        return nonfounderids

    def returnNonFounders(self):
        """ return the founders in a ped file (those with unknown paternal and maternids """
        nonfounders=[]

        for pedobj in self.pedlist:
            if pedobj.getpid() != "0":
                nonfounders.append(pedobj)

        return nonfounders

    def returnIndivids(self):
        """ return a list of indvi ids from the ped file """
        #samplelist=[]
        return [ ped.getid() for ped in self.pedlist]

    def getPedList(self):
        return self.pedlist
    def getTotalSize(self):
        return len(self.pedlist)

    def yieldMembers(self):
        for pedobj in self.pedlist:
            yield pedobj

    def __str__(self):
        return "\n".join( [ x.__str__() for x in self.pedlist ] )
