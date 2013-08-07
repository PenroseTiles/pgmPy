#!/usr/bin/env python
from Factor import *
import numpy as np
from PGMcommon import *
import sys
import itertools
import common
import pdb

def IndexToAssignment( I, D):

    """ given and index I (a row vector representing the indices of values a factor object's val field
        and D, an array representing the cadinality of variables in a factor object, this function produces
        a matrix of assignments, one assignment per row. See https://github.com/indapa/PGM/blob/master/Prog1/IndexToAssignment.m """

    a=np.reshape ( np.arange(np.prod(D)).repeat(len(D)), (np.prod(D),len(D)))
    

    b=tmp=list( D[:-1] )
    tmp.insert(0,1)
    tmp =np.cumprod ( np.array (tmp) )
    b=np.tile( np.cumprod(b), (len(I), 1))
    #print b

    #print np.floor ( a /b )
    c=np.tile ( D, ( len(I), 1) )

    assignment = np.mod ( np.floor( a/b), c)  +1
    return assignment


def AssignmentToIndex ( A, D):
    """ I = AssignmentToIndex(A, D) converts an assignment, A, over variables
        with cardinality D to an index into the .val vector for a factor.
        If A is a matrix then the function converts each row of A to an index.
        See https://github.com/indapa/PGM/blob/master/Prog1/AssignmentToIndex.m """
        
    D=D.flatten(0) #turn array into vector (note that this forces a copy), see http://www.scipy.org/NumPy_for_Matlab_Users#head-fd74115e6798fbf3a628094a55d1cb2b2b5cdd3c
    I=np.array( [] )
    (nrowA,ncolA)=np.shape(A)

    if nrowA== 1 or ncolA ==1: #if assginments are 1 row or 1 col
        #sys.stderr.write("if block ...\n")
        b=tmp=list( D[:-1] )
        tmp.insert(0,1)
        
        tmp =np.cumprod ( np.array (tmp) )
        tmp=(np.array(np.matrix(tmp)))
        #print "tmp: ", tmp
        
        a_flat=np.array ( np.matrix( A.flatten(0) ).transpose() )
        #print "a flat: ", a_flat
        I= ( tmp * (a_flat-1) ) + 1
        return I
        

    else:
        #sys.stderr.write("else block ...\n")
        b=tmp=list( D[:-1] )
        tmp.insert(0,1)
        tmp =np.cumprod ( np.array (tmp) )
        tmp = np.tile( tmp, (nrowA,1) )
        #print tmp
        #print (A-1)
        I= np.sum( np.multiply(tmp, (A-1)), 1) + 1

       
    return np.array( np.matrix( I ).transpose()  )


def SetValueOfAssignment( F, A, v, Vorder=None):
    """ % SetValueOfAssignment Sets the value of a variable assignment in a factor.
%
%   F = SetValueOfAssignment(F, A, v) sets the value of a variable assignment,
%   A, in factor F to v. The order of the variables in A are assumed to be the
%   same as the order in F.var.
%
%   F = SetValueOfAssignment(F, A, v, VO) sets the value of a variable
%   assignment, A, in factor F to v. The order of the variables in A are given
%   by the vector VO. See https://github.com/indapa/PGM/blob/master/Prog1/SetValueOfAssignment.m  """

    

    if Vorder == None:
        indx=AssignmentToIndex( A, F.getCard() )
    else:
        sys.stderr.write("assumes the order of variables in A are the sayme as in F.var ...\n")
        pass

    #http://stackoverflow.com/a/5183720, How to make List from Numpy Matrix in Python
    #http://stackoverflow.com/a/8373103, numpy function to set elements of array to a value given a list of indices
    indices=np.array(indx-1).flatten().tolist()
    zeros=np.zeros(len(A))
    zeros[indices]=v
    F.setVal( zeros.tolist() )

def GetValueOfAssignment( F, A, Vorder = None ):
    """ % GetValueOfAssignment Gets the value of a variable assignment in a factor.
%
%   v = GetValueOfAssignment(F, A) returns the value of a variable assignment,
%   A, in factor F. The order of the variables in A are assumed to be the
%   same as the order in F.var.
%
%   v = GetValueOfAssignment(F, A, VO) gets the value of a variable assignment,
%   A, in factor F. The order of the variables in A are given by the vector VO. See https://github.com/indapa/PGM/blob/master/Prog1/GetValueOfAssignment.m """

    if Vorder  == None:
        indx= AssignmentToIndex ( A, F.getCard() )
    else:
        sys.stderr.write("The order of the variables in A are assumed to be the same as the order in F var\n")
        pass

    indices=np.array(indx-1).flatten().tolist()
    return np.array ( np.matrix ( F.getVal()[indices] ))

def FactorProduct ( A, B):
    """ FactorProduct Computes the product of two factors.
%       C = FactorProduct(A,B) computes the product between two factors, A and B,
%       where each factor is defined over a set of variables with given dimension.
%       The factor data structure has the following fields:
%       .var    Vector of variables in the factor, e.g. [1 2 3]
%       .card   Vector of cardinalities corresponding to .var, e.g. [2 2 2]
%       .val    Value table of size prod(.card)
%
%       See also FactorMarginalization  IndexToAssignment,
%       AssignmentToIndex, and https://github.com/indapa/PGM/blob/master/Prog1/FactorProduct.m """

    #print "A: ", A
    #print "===="
    #print "B: ", B
    C=Factor()

   #check for empty factors
    if len( A.getVar() ) == 0 :
        sys.stderr.write("A factor is empty!\n")
        return B
    if len( B.getVar() ) == 0:
        sys.stderr.write("B factor is empty!\n")
        return A


    #check of  variables that in both A and B have the same cardinality
    #print 'A.getVar():  ', A.getVar()
    #print 'B.getVar(): ',B.getVar()
    #setA= set( A.getVar() )
    #setB= set( B.getVar() )
    #intersect=np.array( list( setA.intersection(setB)))
    intersect=np.intersect1d( A.getVar(), B.getVar() ).tolist()
    #print "Intersection of variables in FactorProduct ", intersect
    #print "A var: ",  A.getVar()
    #print "B var: ",  B.getVar()

    #if the intersection of variables in the two factors
    #is non-zero, then make sure they have the same cardinality
    if len(intersect) > 0:
        #iA=np.nonzero(intersect - A.getVar()==0)[0].tolist() # see this http://stackoverflow.com/a/432146, return the index of something in an array?
        iA=getIndex( A.getVar(), intersect )
        #print "iA: ", iA
        #iB=np.nonzero(intersect - B.getVar()==0)[0].tolist()
        iB = getIndex (  B.getVar(), intersect )
        #print "iB: ", iB

        # check to see if any of the comparisons in the  array resulting from  of a.getCard()[iA] == b.getCard()[iB] 
        # are all False. If so print an error and exit
        if len( np.where( A.getCard()[iA].all() == B.getCard()[iB].all() ==False)[0].tolist() ) > 0:
            sys.stderr.write("dimensionality mismatch in factors!\n")
            sys.exit(1)

    #now set the variables of C to the union of variables in factors A and B
    #print 'setA ' ,setA
    #print 'setB ', setB
    #print list( setA.union(setB) )
    C.setVar( np.union1d ( A.getVar(), B.getVar() ).tolist()  )
    #C.setVar ( list( setA.union(setB) ) )
    mapA=isMember(A.getVar(), C.getVar() )
    mapB=isMember(B.getVar(), C.getVar() )

    

    #Set the cardinality of variables in C
    C.setCard( np.zeros( len(C.getVar())).tolist() )
    C.getCard()[mapA]=A.getCard()
    C.getCard()[mapB]=B.getCard()

    #intitialize the values of the factor C to be zero
    C.setVal( np.zeros(np.prod(C.getCard())).tolist() )

    #some helper indices to tell what indices of A and B values to multiply
    assignments=IndexToAssignment( np.arange(np.prod(C.getCard())), C.getCard() ) #get the assignment of values of C
    indxA=AssignmentToIndex(  assignments[:,mapA], A.getCard())-1 # re-arrange the assignment of C, to what it would be in factor  A
    indxB=AssignmentToIndex(  assignments[:,mapB], B.getCard())-1 # re-arange the assignment of C to what it would be in  factorB

    

    c_val=A.getVal()[indxA.flatten().tolist()] * B.getVal()[indxB.flatten().tolist()] #now that we have the index into A.val and B.val vector, multiply them to factor product
    C.setVal ( c_val.tolist() )

    return C

def FactorMarginalization(A,V):
    """   FactorMarginalization Sums given variables out of a factor.
          B = FactorMarginalization(A,V) computes the factor with the variables
          in V summed out. The factor data structure has the following fields:
          .var    Vector of variables in the factor, e.g. [1 2 3]
          .card   Vector of cardinalities corresponding to .var, e.g. [2 2 2]
          .val    Value table of size prod(.card)

          The resultant factor should have at least one variable remaining or this
          function will throw an error.   See also FactorProduct, IndexToAssignment , and AssignmentToIndex
          Based on matlab code found here: https://github.com/indapa/PGM/blob/master/Prog1/FactorMarginalization.m """

    #the resulting factor after marginalizing out variables in python list V that are in 
    #the factor A
    B=Factor()

    #check for empy factor or variable list
    if len( A.getVar() ) == 0 or len(V) == 0:
        return A

    #construct the variables of the marginalized factor by 
    #computing the set difference between A.var and V
    #These variables in the difference set will be the scope of the new factor
    setA=set( A.getVar() )
    setV=set(V)
    Bvar=np.array( list( setA.difference(setV)))
    mapB=isMember(Bvar, A.getVar()) #indices of the variables of the new factor in the original factor A
    #print mapB,  Bvar

    #check to see if the new factor has empty scope
    if len(Bvar) == 0:
        sys.stderr.write("FactorMarginalization:Error, resultant factor has empty scope...\n")
        return None
    #set the marginalized factor's variable scope and cardinality
    B.setVar( Bvar.tolist() )
    B.setCard( A.getCard()[mapB] )
    B.setVal( np.zeros(np.prod(B.getCard())).tolist() )

    #compute some helper indices
    assignments=IndexToAssignment ( np.arange(np.prod(A.getCard()) ), A.getCard() )
    #indxB tells which values in A to sum together when marginalizing out the variable(s) in B
    indxB=AssignmentToIndex( assignments[:,mapB], B.getCard())-1

    #accum is a numpy implementation of matlab accumarray
    #accumarray sums data in each group
    #here the group(s) are defined in indxB
    #indxB is a map to tell which value in A.val to map the sum to
    #see http://blogs.mathworks.com/loren/2008/02/20/under-appreciated-accumarray/
    marginal_vals=accum(indxB, A.getVal() )
    
    #set the marginal values to the new factor with teh variable(s) in V summed(marginalized) out
    B.setVal( marginal_vals.tolist() )
    return B



def ObserveEvidence (INPUTS, EVIDENCE):

    """   ObserveEvidence Modify a vector of factors given some evidence.
          F = ObserveEvidence(INPUTS, EVIDENCE) sets all entries in the vector of factors,INPUTS,
          that are not consistent with the evidence, E, to zero. F is a vector of
          factors, each a data structure with the following fields:
          .var    Vector of variables in the factor, e.g. [1 2 3]
          .card   Vector of cardinalities corresponding to .var, e.g. [2 2 2]
          .val    Value table of size prod(.card)
          EVIDENCE is an N-by-2 matrix, where each row consists of a variable/value pair.
          Variables are in the first column and values are in the second column.   """
    (nrows, ncols)=np.shape(EVIDENCE)
    #total_factors=len(INPUTS)
    #iterate through evidence
    for i in range(nrows):
        variable=EVIDENCE[i,0]
        value=EVIDENCE[i,1]
        #print 'var: ', variable, 'value: ', value
        if int(value) == 0:
            print "Evidence is not set for variable: ', variable, ' in evidence matrix.\n"
            continue

        for factor in INPUTS:
            #the following returns a list
            indx=np.where( factor.getVar() == variable )[0].tolist()
            if indx: #if the indx is not empty, it contains the index value  of the evidence variable in factor.val array
                indx=indx[0]
                
            if value > factor.getCard()[indx] or value < 0:
                sys.stderr.write("invalid evidene for variable X_'" + str(variable) + " = " + str(value) + "\n")
                sys.exit(1)

            #get the assignments of variables for the factor
            assignments=IndexToAssignment( np.arange(np.prod( factor.getCard() )), factor.getCard() )
            # now get the indices in  the assignments that don't agree with the observed  value (evidence)
            mask=np.where( assignments[:,indx] != value )[0].tolist()
            #we are going to update the val array for the current factor
            newvals=factor.getVal()
            #set the mask indices to zero and reset the val array of the factor
            newvals[mask]=0
            factor.setVal( newvals.tolist() )

            #now check to see the validity of the updated values of the factor
            #given the observed evidence. We cannot have all zero values for the factor!
            zeroIndices=np.where ( factor.getVal() == 0)[0].tolist()
            if len(zeroIndices) == len (factor.getVal() ):
                sys.stderr.write("All variable values are zero, which is not possible.\n")
    return INPUTS

def ComputeJointDistribution(INPUTS):
    """ ComputeJointDistribution Computes the joint distribution defined by a set of given factors

    Joint = ComputeJointDistribution(INPUTS) computes the joint distribution
    defined by a set of given factors

    Joint is a factor that encapsulates the joint distribution given by INPUTS
    INPUTS is a vector of Factor objects containing the factors defining the distribution

    """

    totalFactors = len(INPUTS)
    #check for empty list of INPUTS

    if totalFactors== 0:
        sys.stderr.write("Empty factor list given as input\n")
        return Factor( [], [], [] )
 
    else:
        # see http://docs.python.org/library/functions.html#reduce for description of Python reduce function
        return reduce(lambda x, y: FactorProduct(x,y), INPUTS)


def ComputeMarginal(V, F, E):
    """
        ComputeMarginal Computes the marginal over a set of given variables
        M = ComputeMarginal(V, F, E) computes the marginal over variables V
        in the distribution induced by the set of factors F, given evidence E

        M is a factor containing the marginal over variables V

        V is a vector containing the variables in the marginal e.g. [1 2 3] for X_1, X_2 and X_3.
        i.e. a result of FactorMarginalization

        F is a vector of factors (struct array) containing the factors
        defining the distribution

        E is an N-by-2 matrix, each row being a variable/value pair.
        Variables are in the first column and values are in the second column.
        If there is no evidence, pass in the empty matrix [] for E.

    """
    totalFactors=len(F)
    #reshape a 1d array to 1 x ncol array
    #since ObserveEvidence requires Nx2 array, we reshape to a 2 column array
    #see http://stackoverflow.com/a/12576163 for reshaping 1d array to 2d array
    EVIDENCE= np.reshape( np.array ( E ), (-1,2) )
    #print np.shape(EVIDENCE)
    
    if totalFactors == 0:
        sys.stderr.write("empty factor list given as input.\n")
        return Factor( [], [], [])
    # union of all variables in list of factors F
    variableList=[] # a list of of lists, where each element is a list containing the variables of the factor in F
    for factor in F:
        var=factor.getVar().tolist()
        variableList.append( var )

    #get the union of variables across all the factor in F
    #see this http://stackoverflow.com/a/2151553, Pythonic Way to Create Union of All Values Contained in Multiple Lists
    union_variables = set().union(*variableList)
    #print union_variables
    #v contains the variables not in the list of variables in the marginal
    v=list( union_variables.difference(V) )
   
    # compute the joint distribution, but then reduce it, given the evidence
    # ComputeJointDistribution returns a factor, but ObserveEvidence expects a list
    # of factors as the first argument, so hence the need for brackets [ ]
    # ObserveEvidence returns a list, but we want the first element so thats why the [0]
    jointE= ObserveEvidence ( [ComputeJointDistribution ( F )], EVIDENCE )[0]

    #now we need to re-normaize the joint, since observe evidence doesn't do it for us
    jointE_normalizedVal = jointE.getVal()/np.sum( jointE.getVal() )
    jointE.setVal( jointE_normalizedVal.tolist() )

    return FactorMarginalization ( jointE, v)
    

def IdentityFactor( F ):
    return Factor ( F.getVar().tolist(), F.getCard().tolist(), np.ones( np.prod( F.getCard() ) ), F.getName()+ '_identity' )



def SumProductEliminateVar(z, factorList):

    """ this is a non-graph based  sum-product variable elimination function
    z is a variable to eliminate
       pass in a list of factors

       1. figure out which factor contain z in their variable scope
       2. figure out which factors don't contain z in their scope
       3. mulitply in all factors that have z
       4. sum out z (marginalize) and return new factor with variable z eliminated"""

    useFactors = []# list  of factors that contains the variable Z
    unusedFactors=[] #list of factors that don't contain variable Z
    scope = []

    """get a list containining the index in self.factorLlist of factors
       that contain the variable Z to be eliminated
       get the scope of variables from the factors that contain variable Z """
    for fi in factorList:
        if z in fi.getVar().tolist():
            useFactors.append(fi)#the ith factor is being currently involved in elimination
            scope=list(set.union(set(scope), fi.getVar().tolist() ))
        else:
            unusedFactors.append( fi )

    #for f in useFactors:
    #    print 'useFactor: ', f
    #print '==='
    #for f in  unusedFactors:
    #    print 'unusedFactor: ', f

    psiFactor= ComputeJointDistribution ( useFactors )
    tauFactor=FactorMarginalization( psiFactor,[z] )

    #print 'psiFactor: ', psiFactor
    #print 'tauFactor: ', tauFactor
    return unusedFactors + [ tauFactor ]

def SumProductVE ( Z, F ):

    """ A wrapper function for SumProductEliminateVar
        sum-product variable elimination based on pseudocode algorithm 9.1 in Koller and Friedman
        We are giving a list of variables to eliminate in Z (in the order we want to
        elimiinate them) and a list of factors F
        eliminate each one getting getting the marginal distribution of the last variable in the list
        Z. """

    for z in Z:
        F=SumProductEliminateVar(z, F)
    return reduce(lambda x, y: FactorProduct(x,y), F)


def FactorMaxMarginalization( A, V ):
    """ computes the factor with the variables in V *maxed* out.
        The resulting factor will have all the variables in A minus
        those variables in V. This is quite similiar to FactorMarginalization, but rather then summing out variables in V
        we take the max. In the code, this translates passing np.max as the function to accum
        See section  13.2 in Koller and Friedman  for more information"""

    B=Factor()
    #check for empy factor or variable list
    if len( A.getVar() ) == 0 or len(V) == 0:
        return A
    Bvar=np.setdiff1d( A.getVar(), V)
    mapB=isMember(Bvar, A.getVar())

    if len(Bvar) == 0:
        sys.stderr.write("FactorMaxMarginalization: Error, resultant factor has empty scope...\n")
        return np.max (A.getVal() )
    #set the marginalized factor's variable scope and cardinality
    B.setVar( Bvar.tolist() )
    B.setCard( A.getCard()[mapB] )
    B.setVal( np.zeros(np.prod(B.getCard())).tolist() )

    #compute some helper indices
    assignments=IndexToAssignment ( np.arange(np.prod(A.getCard()) ), A.getCard() )
    #indxB tells which values in A to sum together when marginalizing out the variable(s) in B
    indxB=AssignmentToIndex( assignments[:,mapB], B.getCard())-1

    #here we pass in the function np.max
    #NumPy and Python are awesome
    max_vals=accum(indxB, A.getVal(), np.max )
    B.setVal( max_vals.tolist() )

    return B


def MaxProductEliminateVar(z, factorList):

    """ this is a non-graph based  MAX-product variable elimination function
        z is a variable to eliminate
        pass in a list of factors

       1. figure out which factor contain z in their variable scope
       2. figure out which factors don't contain z in their scope
       3. mulitply in all factors that have z
       4. max marginalize out z  and return new factor with variable z eliminated"""

    useFactors = []# list  of factors that contains the variable Z
    unusedFactors=[] #list of factors that don't contain variable Z
    scope = []

    """get a list containining the index in self.factorLlist of factors
       that contain the variable Z to be eliminated
       get the scope of variables from the factors that contain variable Z """
    for fi in factorList:
        if z in fi.getVar().tolist():
            useFactors.append(fi)#the ith factor is being currently involved in elimination
            scope=list(set.union(set(scope), fi.getVar().tolist() ))
        else:
            unusedFactors.append( fi )

   
    """ psiFactor is an intermediate factor, prior to max-marginalization """
    psiFactor= ComputeJointDistribution ( useFactors )
    tauFactor=FactorMaxMarginalization( psiFactor,[z] )

    
    """ we return tuple consisting of
    1. a list factors that are unused, plus the result of max-marginal
    such that the variable z is not eliminated from the list of factors remaining.
    2. For traceback, we return the intermediate factor generated in the process of eliminating
    variable z. """
    return unusedFactors + [ tauFactor ], psiFactor

def TracebackMAP(FI, Z):
    """We are back-tracing to the most probable assingments here ...
    See psuedo-code in Koller and Friedman page 557

    In order to return the most probable assignment from MaxProductVE
    we take in a list of intermediate factors, FI, that were generated in the process
    of MaxProductVE. Z is the same elimination ordering used as MaxProductVE.
    We traceback our steps by iterating in reverse order the elimination ordering Z.

    Following arguments  section 13.2.2 page 558 in Koller and Friedman, as one eliminates
    variables you cannot determine their maximizing value. But you can compute their 'conditional'
    maximizing value - their max value given the other variables not eliminate yet. Once
    the last variable is eliminated, we can traceback to get the maximzing value of the remaining
    variables. Hence the reason for iterating thru the elimination ordering Z in reverse order
    
    Returns a python dictionary with key: variable value: variable assignment in the MAP"""

    z_i=Z[-1]
    f_i=FI[-1]
    Z=Z[:-1]
    FI=FI[:-1]

    #print 'z_i:', z_i
    #print 'f_i:', f_i

    values=f_i.getVal().tolist()
    fidx= IndexToAssignment( np.arange( np.prod( f_i.getCard()  ) ), f_i.getCard() )
    maxidx=values.index(max(values))
    

    maxed_vars={} #key variable, value: max assignment value

    #print 'variable: ', z_i
    #print 'max value: ', fidx.flatten()[maxidx]
    maxed_vars[z_i]=int(fidx.flatten()[maxidx])
    #print maxed_vars
    #print
    for (z, f) in itertools.izip( reversed(Z), reversed(FI) ):
        #print z
        #print f
        #print 'setdiff: ', np.setdiff1d( f_i.getVar(), [z]).tolist()
        
        variables=np.setdiff1d( f_i.getVar(), [z]).tolist()
        evidence=[ [v, maxed_vars[v] ] for v in variables]
        #print 'Evidence: ',evidence
        f=ObserveEvidence( [f], np.matrix(evidence) )[0]
        #print f
        values=f.getVal().tolist()
        fidx= IndexToAssignment( np.arange( np.prod( f.getCard()  ) ), f.getCard() )
        #print fidx
        #print max(values)
        maxidx=values.index(max(values))

        maxed_vars[z]=int(fidx.flatten()[maxidx])
        #print 'variable: ', z
        #print 'max value: ', fidx.flatten()[maxidx]
        
        #print
    #print maxed_vars
    return maxed_vars



def MaxDecoding ( F ):
    """ F is a list of max marginal factors passed in. The factors have a variable scope over a single variable only
        So no backtracing is inovlved, we just get the index of the highest number in the value array.
        The code here is based on https://github.com/indapa/PGM/blob/master/Prog4/MaxDecoding.m """
    ASSIGNMENTS=[]
    for f in F:
        values=f.getVal().tolist()
        ASSIGNMENTS.append ( values.index( max(values) ) )
    return ASSIGNMENTS
    
def MaxDecodingNonUniq ( F ):
    """ F is a list of max marginal factors passed in. We don't assume that there is a unique
        max value. So we get the indices of the non-uniq max value as a tuple and add it to """
    ASSIGNMENTS=[]
    for f in F:
        values=f.getVal().tolist()
        maxvalue=max(values)
        """ if the maxvalue is duplicated, we get the indices of where it resides in the value array """
        if common.isMaxDuplicated(values):
            dup_indices_list=[dup for dup in sorted(common.list_duplicates(values)) ]
            dup_values= [ x for (x, y) in dup_indices_list ]
            dup_indices= [ y for (x, y) in dup_indices_list ]
            non_uniq_max_indices=tuple(dup_indices [ dup_values.index(maxvalue) ])
            ASSIGNMENTS.append ( non_uniq_max_indices )
        else:
            ASSIGNMENTS.append( values.index(maxvalue))
    return ASSIGNMENTS


def posterior_genotypes_values(factorList, ALPHABET,samplenames,bedstring,fh):
    """ given the factorlist of posterior marginals, the genotype alphabet, samplenames,bedstring position,
        and prettybase file handle, print to file the posterior marginals for all 10 possibel genotypes
        for each sample. """
    genotype_factors=factorList[0:len(samplenames)]
    sample_factorObj_zip=zip(samplenames, genotype_factors)
    #print bedstring
    for sample, f in sample_factorObj_zip:
        #print sample, ": "
        #values=f.getVal().tolist()
        
        prob_val_normalized=( lognormalize( f.getVal() ) )
        #print prob_val_normalized.tolist()
        #genotype_probZip=zip(ALPHABET,values)
        posteriors=[]
        #print prob_val_normalized.tolist()
        for posterior_val in prob_val_normalized.tolist():
        #for posterior_val in values:
            posteriors.append(str(posterior_val))
            #posteriors.append(str(round(posterior_val,5) ))
        #print posteriors
        gstring="\t".join(posteriors)
        #print gstring
        outstring="\t".join([bedstring, sample,gstring])
        
        fh.write(outstring + "\n")
        

    
def MaxProductVE ( Z, F ):

    """ A wrapper function for MaxProductEliminateVar
        sum-product variable elimination based on pseudocode algorithm 9.1 in Koller and Friedman
        We are giving a list of variables to eliminate in Z (in the order we want to
        elimiinate them) and a list of factors F
        eliminate each one getting getting the marginal distribution of the last variable in the list
        Z.

        Returns the probabliity of the MAP configuration as well as the variable assignments of the MAP configuration"""
    intermediateMaxFactors=[]
    for z in Z:
        (F, intermediateFactor)=MaxProductEliminateVar(z, F)
        intermediateMaxFactors.append ( intermediateFactor )
       
        #intermediateMaxFactors.append ( intermediateFactor )
   
    #MaxDecodingBT( intermediateMaxFactors, Z )
    bt_results=TracebackMAP( intermediateMaxFactors, Z )
    return (reduce(lambda x, y: FactorProduct(x,y), F), bt_results)


def FactorSum ( A, B):
    """ FactorSum Computes the sum of two factors.
%       Similiar to FactorProduct
        We would use this in log space where multiplication becomes addition
%       Based on the code here https://github.com/indapa/PGM/blob/master/Prog4/FactorSum.m """


    C=Factor()

   #check for empty factors
    if len( A.getVar() ) == 0 :
        sys.stderr.write("A factor is empty!\n")
        return B
    if len( B.getVar() ) == 0:
        sys.stderr.write("B factor is empty!\n")
        return A


    #check of  variables that in both A and B have the same cardinality
    #print 'A.getVar():  ', A.getVar()
    #print 'B.getVar(): ',B.getVar()
    #setA= set( A.getVar() )
    #setB= set( B.getVar() )
    #intersect=np.array( list( setA.intersection(setB)))
    intersect=np.intersect1d( A.getVar(), B.getVar() ).tolist()
    #print "Intersection of variables in FactorProduct ", intersect
    #print "A var: ",  A.getVar()
    #print "B var: ",  B.getVar()

    #if the intersection of variables in the two factors
    #is non-zero, then make sure they have the same cardinality
    if len(intersect) > 0:
        #iA=np.nonzero(intersect - A.getVar()==0)[0].tolist() # see this http://stackoverflow.com/a/432146, return the index of something in an array?
        iA=getIndex( A.getVar(), intersect )
        #print "iA: ", iA
        #iB=np.nonzero(intersect - B.getVar()==0)[0].tolist()
        iB = getIndex (  B.getVar(), intersect )
        #print "iB: ", iB

        # check to see if any of the comparisons in the  array resulting from  of a.getCard()[iA] == b.getCard()[iB]
        # are all False. If so print an error and exit
        if len( np.where( A.getCard()[iA].all() == B.getCard()[iB].all() ==False)[0].tolist() ) > 0:
            sys.stderr.write("dimensionality mismatch in factors!\n")
            sys.exit(1)

    #now set the variables of C to the union of variables in factors A and B
    #print 'setA ' ,setA
    #print 'setB ', setB
    #print list( setA.union(setB) )
    C.setVar( np.union1d ( A.getVar(), B.getVar() ).tolist()  )
    #C.setVar ( list( setA.union(setB) ) )
    mapA=isMember(A.getVar(), C.getVar() )
    mapB=isMember(B.getVar(), C.getVar() )



    #Set the cardinality of variables in C
    C.setCard( np.zeros( len(C.getVar())).tolist() )
    C.getCard()[mapA]=A.getCard()
    C.getCard()[mapB]=B.getCard()

    #intitialize the values of the factor C to be zero
    C.setVal( np.zeros(np.prod(C.getCard())).tolist() )

    #some helper indices to tell what indices of A and B values to multiply
    assignments=IndexToAssignment( np.arange(np.prod(C.getCard())), C.getCard() ) #get the assignment of values of C
    indxA=AssignmentToIndex(  assignments[:,mapA], A.getCard())-1 # re-arrange the assignment of C, to what it would be in factor  A
    indxB=AssignmentToIndex(  assignments[:,mapB], B.getCard())-1 # re-arange the assignment of C to what it would be in  factorB
    #print 'indxA ', indxA
    #print 'indxB ', indxB


    c_val=A.getVal()[indxA.flatten().tolist()] + B.getVal()[indxB.flatten().tolist()] #now that we have the index into A.val and B.val vector, multiply them to factor product
    C.setVal ( c_val.tolist() )

    return C


def LogFactor( F ):
    """ return a factor whose values are the  natural log of the orginal factor F  """
    
    return Factor ( F.getVar().tolist(), F.getCard().tolist(), np.log ( F.getVal() ).tolist(), F.getName() )


def ExpFactorNormalize ( logF ):
        """ exponentiate a factor to probablity space from log space
        Since moving from log space to probablity space incures a decrease in dynamic range,
        factors should be normalized before applying the transform. One trick we use here is
        to shift every entry by the maximum entry. For example: phi[i] = exp{logPhi[i] -c}
        The value of c is max(logPhi). This type of transformation ensures the resulting factor
        has a maximum entry of 1 and prevents overflow. See page 360 of Koller and Friedman text"""
        logPhi=logF.getVal()
        phi=np.exp(logPhi-np.max(logPhi) )
        logF.setVal( phi )
        return logF

def ExpFactor( logF ):
    """ similiar to above, but don't normalize  """
    logPhi=logF.getVal()
    phi=np.exp( logPhi )
    logF.setVal( phi )
    return logF




def FactorDiv ( A, B):
    """ FactorProduct Computes the dividend of two factors.
%       Similiar to Factor Product, but if we divide 0/0, return 0
    see page 365 in Koller and Friedman for definition of FactorDivision """

    #print "A: ", A
    #print "===="
    #print "B: ", B
    C=Factor()

   #check for empty factors
    if len( A.getVar() ) == 0 :
        sys.stderr.write("A factor is empty!\n")
        return B
    if len( B.getVar() ) == 0:
        sys.stderr.write("B factor is empty!\n")
        return A


    #check of  variables that in both A and B have the same cardinality
    #print 'A.getVar():  ', A.getVar()
    #print 'B.getVar(): ',B.getVar()
    #setA= set( A.getVar() )
    #setB= set( B.getVar() )
    #intersect=np.array( list( setA.intersection(setB)))
    intersect=np.intersect1d( A.getVar(), B.getVar() ).tolist()
    #print "Intersection of variables in FactorProduct ", intersect
    #print "A var: ",  A.getVar()
    #print "B var: ",  B.getVar()

    #if the intersection of variables in the two factors
    #is non-zero, then make sure they have the same cardinality
    if len(intersect) > 0:
        #iA=np.nonzero(intersect - A.getVar()==0)[0].tolist() # see this http://stackoverflow.com/a/432146, return the index of something in an array?
        iA=getIndex( A.getVar(), intersect )
        #print "iA: ", iA
        #iB=np.nonzero(intersect - B.getVar()==0)[0].tolist()
        iB = getIndex (  B.getVar(), intersect )
        #print "iB: ", iB

        # check to see if any of the comparisons in the  array resulting from  of a.getCard()[iA] == b.getCard()[iB] 
        # are all False. If so print an error and exit
        if len( np.where( A.getCard()[iA].all() == B.getCard()[iB].all() ==False)[0].tolist() ) > 0:
            sys.stderr.write("dimensionality mismatch in factors!\n")
            sys.exit(1)

    #now set the variables of C to the union of variables in factors A and B
    #print 'setA ' ,setA
    #print 'setB ', setB
    #print list( setA.union(setB) )
    C.setVar( np.union1d ( A.getVar(), B.getVar() ).tolist()  )
    #C.setVar ( list( setA.union(setB) ) )
    mapA=isMember(A.getVar(), C.getVar() )
    mapB=isMember(B.getVar(), C.getVar() )

    

    #Set the cardinality of variables in C
    C.setCard( np.zeros( len(C.getVar())).tolist() )
    C.getCard()[mapA]=A.getCard()
    C.getCard()[mapB]=B.getCard()

    #intitialize the values of the factor C to be zero
    C.setVal( np.zeros(np.prod(C.getCard())).tolist() )

    #some helper indices to tell what indices of A and B values to multiply
    assignments=IndexToAssignment( np.arange(np.prod(C.getCard())), C.getCard() ) #get the assignment of values of C
    indxA=AssignmentToIndex(  assignments[:,mapA], A.getCard())-1 # re-arrange the assignment of C, to what it would be in factor  A
    indxB=AssignmentToIndex(  assignments[:,mapB], B.getCard())-1 # re-arange the assignment of C to what it would be in  factorB
    
    numerator=A.getVal()[indxA.flatten().tolist()]
    denominator=B.getVal()[indxB.flatten().tolist()]
    
    #print numerator
    #print denominator
    #print zip(numerator, denominator)
    val= map( lambda x: common.zerodiv_tuple(x), zip(numerator,denominator)  )
    #print val
    C.setVal ( val )
    
    return C