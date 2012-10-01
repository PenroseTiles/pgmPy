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
        sys.stderr.write("if block ...\n")
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
