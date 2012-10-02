#!/usr/bin/env python
from Factor import *
import numpy as np
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