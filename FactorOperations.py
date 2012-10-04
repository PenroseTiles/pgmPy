#!/usr/bin/env python
from Factor import *
import numpy as np
from PGMcommon import *
import sys



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

    C=Factor()

   #check for empty factors
    if len( A.getVar() ) == 0 :
        sys.stderr.write("A factor is empty!\n")
        return B
    if len( B.getVar() ) == 0:
        sys.stderr.write("B factor is empty!\n")
        return A


    #check of  variables that in both A and B have the same cardinality
    setA= set( A.getVar() )
    setB= set( B.getVar() )
    intersect=np.array( list( setA.intersection(setB)))

    #if the intersection of variables in the two factors
    #is non-zero, then make sure they have the same cardinality
    if len(intersect) > 0:
        iA=np.nonzero(intersect - A.getVar()==0)[0].tolist() # see this http://stackoverflow.com/a/432146, return the index of something in an array?
        iB=np.nonzero(intersect - B.getVar()==0)[0].tolist()

        # check to see if any of the comparisons in the  array resulting from  of a.getCard()[iA] == b.getCard()[iB] 
        # are all False. If so print an error and exit
        if len( np.where( A.getCard()[iA] == B.getCard()[iB] ==False)[0].tolist() ) > 0:
            sys.stderr.write("dimensionality mismatch in factors!\n")
            sys.exit(1)

    #now set the variables of C to the union of variables in factors A and B
    C.setVar ( list( setA.union(setB) ) )
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
    print mapB,  Bvar

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
  #if totalFactors ==2:
  #    return FactorProduct ( INPUTS[0], INPUTS[1] )
    else:
        # see http://docs.python.org/library/functions.html#reduce for description of Python reduce function
        return reduce(lambda x, y: FactorProduct(x,y), INPUTS)





