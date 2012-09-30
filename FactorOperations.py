#!/usr/bin/env python
from Factor import *
import numpy as np

def IndexToAssignment( I, D):
    a=np.reshape ( np.arange(np.prod(D)).repeat(len(D)), (np.prod(D),len(D)))
    #print np.shape(a)
    #print a

    b=tmp=list( D[:-1] )
    tmp.insert(0,1)
    tmp =np.cumprod ( np.array (tmp) )
    b=np.tile( np.cumprod(b), (len(I), 1))
    #print b

    #print np.floor ( a /b )
    c=np.tile ( D, ( len(I), 1) )

    assignment = np.mod ( np.floor( a/b), c)  +1
    return assignment

    #repmat(cumprod([1, D(1:end - 1)]), length(I), 1)
