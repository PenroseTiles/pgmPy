#!/usr/bin/env python
from Factor import *
import numpy as np

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

    
